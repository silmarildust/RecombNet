from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path
from typing import Sequence

import numpy as np

from .core import (
    ARGResult,
    ExtendedBarcodeResult,
    analyze_group,
    build_arg_graph,
    compute_extended_barcode,
    load_fasta,
    load_manifest,
    render_arg_graph,
    segregating_sites,
)


def _read_positions(path: str | None) -> np.ndarray | None:
    if not path:
        return None
    values: list[int] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                values.append(int(line))
    return np.array(values, dtype=int)


def _write_bars(prefix: str, result: ExtendedBarcodeResult) -> None:
    out = Path(prefix + ".bars.tsv")
    with out.open("w", encoding="utf-8") as f:
        f.write("birth\tdeath\tlength\tbirth_edge\tdeath_triangle\tgenerator_edges\n")
        for interval in result.intervals:
            birth_edge = f"{interval.birth_edge[0]}-{interval.birth_edge[1]}"
            death_triangle = "" if interval.death_triangle is None else "-".join(map(str, interval.death_triangle))
            generator = ",".join(f"{a}-{b}" for a, b in interval.generator_edges)
            death = "inf" if interval.death is None else f"{interval.death:.6f}"
            length = "inf" if np.isinf(interval.length) else f"{interval.length:.6f}"
            f.write(f"{interval.birth:.6f}\t{death}\t{length}\t{birth_edge}\t{death_triangle}\t{generator}\n")


def _write_b1(prefix: str, result: ExtendedBarcodeResult) -> None:
    out = Path(prefix + ".b1.tsv")
    mat = np.asarray(result.count_matrix, dtype=int)
    with out.open("w", encoding="utf-8") as f:
        f.write("start\tend\tb1\n")
        for i in range(mat.shape[0]):
            for j in range(i + 1, mat.shape[1]):
                if mat[i, j] > 0:
                    f.write(f"{i}\t{j}\t{mat[i, j]}\n")


def _plot_group_summary(prefix: str, result: ExtendedBarcodeResult) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    births = [iv.birth for iv in result.intervals if np.isfinite(iv.birth)]
    axes[0].hist(births, bins=min(25, max(1, len(births))))
    axes[0].set_title("H1 birth-time histogram")
    axes[0].set_xlabel("birth")
    axes[0].set_ylabel("count")

    mat = np.asarray(result.count_matrix, dtype=int)
    axes[1].imshow(mat, aspect="auto", interpolation="nearest")
    axes[1].set_title("H1 support-count matrix")
    axes[1].set_xlabel("site index")
    axes[1].set_ylabel("site index")

    fig.suptitle("RecombNet summary")
    fig.savefig(prefix + ".png", dpi=200)
    plt.close(fig)


def _scan_single(args: argparse.Namespace) -> None:
    labels, seqs = load_fasta(args.fasta)
    positions = _read_positions(args.positions)
    if positions is not None and len(positions) != len(seqs[0]):
        raise SystemExit(
            f"positions file has {len(positions)} entries, but the FASTA alignment has {len(seqs[0])} columns"
        )

    poss, compressed = segregating_sites(seqs, exclude_compatible=args.exclude_compatible)
    if poss.size == 0:
        raise SystemExit("No segregating sites found after filtering")

    if positions is not None:
        positions = positions[poss]

    print(f"Read {len(seqs)} sequences with {len(poss)} segregating sites")

    result = compute_extended_barcode(
        compressed,
        max_sites=args.sites,
        window=args.window,
        min_sites=args.truncate,
        positions=positions,
        labels=labels,
        model=args.model,
        exclude_compatible=False,
        legacy=False,
    )
    assert isinstance(result, ExtendedBarcodeResult)

    _write_bars(args.output, result)
    _write_b1(args.output, result)
    with open(args.output + ".pkl", "wb") as f:
        pickle.dump(result, f)

    manifest = {
        "mode": "scan",
        "input": args.fasta,
        "output_prefix": args.output,
        "labels": labels,
        "num_sequences": len(seqs),
        "num_segregating_sites": int(len(poss)),
        "window": args.window,
        "sites": args.sites,
        "truncate": args.truncate,
        "model": args.model,
        "exclude_compatible": args.exclude_compatible,
    }
    with open(args.output + ".json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    if not args.no_graphics:
        _plot_group_summary(args.output, result)

    print(f"Wrote {args.output}.bars.tsv, {args.output}.b1.tsv, {args.output}.pkl, {args.output}.json")


def _scan_arg(args: argparse.Namespace) -> None:
    manifest_rows = load_manifest(args.manifest)
    groups = []
    for row in manifest_rows:
        groups.append(
            analyze_group(
                row["path"],
                group_id=row["group_id"],
                location=row["location"],
                lineages=row["lineages"],
                positions=_read_positions(args.positions) if args.positions else None,
                max_sites=args.sites,
                window=args.window,
                min_sites=args.truncate,
                model=args.model,
                exclude_compatible=args.exclude_compatible,
            )
        )
        print(f"Analysed {row['group_id']} ({row['location']})")

    arg = render_arg_graph(
        groups,
        output_prefix=args.output,
        top_events=args.top_events,
        min_interval_length=args.min_interval_length,
        merge_shared_lineages=True,
    )

    summary = {
        "mode": "arg",
        "manifest": args.manifest,
        "output_prefix": args.output,
        "groups": [
            {
                "group_id": g.group_id,
                "location": g.location,
                "lineages": g.lineages,
                "path": g.path,
                "num_sequences": len(g.labels),
                "num_segregating_sites": int(len(g.segregating_positions)),
                "num_intervals": len(g.extended.intervals),
            }
            for g in groups
        ],
        "graphml": arg.graphml_path,
        "png": arg.png_path,
        "top_events": args.top_events,
        "min_interval_length": args.min_interval_length,
        "sites": args.sites,
        "truncate": args.truncate,
        "window": args.window,
        "model": args.model,
        "exclude_compatible": args.exclude_compatible,
    }
    with open(args.output + ".arg.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote {arg.png_path}, {arg.graphml_path}, and {args.output}.arg.json")


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="recombnet",
        description="RecombNet: recombination screening and visual ARG generation for aligned SARS-CoV-2 sequences",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    scan = sub.add_parser("scan", help="analyse a single aligned FASTA file")
    scan.add_argument("fasta", help="aligned FASTA file (Nextclade-aligned is fine)")
    scan.add_argument("-o", "--output", default="out", help="output prefix (default: out)")
    scan.add_argument("-p", "--positions", default="", help="file with genomic positions, one integer per alignment column")
    scan.add_argument("-s", "--sites", type=int, default=5, help="maximum subset size (default: 5)")
    scan.add_argument("-t", "--truncate", type=int, default=2, help="minimum subset size (default: 2)")
    scan.add_argument("-w", "--window", type=int, default=13, help="sliding window width over segregating sites (default: 13)")
    scan.add_argument("-e", "--exclude-compatible", action="store_true", help="exclude fully compatible segregating sites")
    scan.add_argument("-n", "--no-graphics", action="store_true", help="do not write a PNG summary")
    scan.add_argument("--model", choices=["hamming", "jc69", "k80", "t92"], default="hamming", help="distance model")
    scan.set_defaults(func=_scan_single)

    arg = sub.add_parser("arg", help="build a visual ARG from a manifest of location files")
    arg.add_argument("manifest", help="TSV/CSV manifest with file/path, location, and lineages columns")
    arg.add_argument("-o", "--output", default="out", help="output prefix (default: out)")
    arg.add_argument("-p", "--positions", default="", help="optional genomic positions file applied to every alignment")
    arg.add_argument("-s", "--sites", type=int, default=5, help="maximum subset size (default: 5)")
    arg.add_argument("-t", "--truncate", type=int, default=2, help="minimum subset size (default: 2)")
    arg.add_argument("-w", "--window", type=int, default=13, help="sliding window width over segregating sites (default: 13)")
    arg.add_argument("-e", "--exclude-compatible", action="store_true", help="exclude fully compatible segregating sites")
    arg.add_argument("--model", choices=["hamming", "jc69", "k80", "t92"], default="hamming", help="distance model")
    arg.add_argument("--top-events", type=int, default=12, help="maximum H1 events to draw per location file")
    arg.add_argument("--min-interval-length", type=float, default=0.0, help="ignore intervals shorter than this length")
    arg.set_defaults(func=_scan_arg)

    parsed = parser.parse_args(argv)
    parsed.func(parsed)


if __name__ == "__main__":
    main()
