from __future__ import annotations

from dataclasses import dataclass, asdict, field
from itertools import combinations
from pathlib import Path
from typing import Sequence, Iterator
import csv
import json
import math
import os
import pickle
import re

import numpy as np


# -----------------------------
# Data structures
# -----------------------------


@dataclass(frozen=True)
class H1Interval:
    """A single H1 persistence interval."""

    birth: float
    death: float | None
    birth_edge: tuple[int, int]
    death_triangle: tuple[int, int, int] | None = None
    generator_edges: tuple[tuple[int, int], ...] = ()

    @property
    def length(self) -> float:
        if self.death is None or math.isinf(self.death):
            return math.inf
        return float(self.death - self.birth)


@dataclass
class BarcodeResult:
    """Barcode summary for one metric space / one subsample."""

    dimensions: list[int]
    births: list[float]
    lengths: list[float]
    generators: list[list[tuple[tuple[int, int], float]]]
    intervals: list[H1Interval] = field(default_factory=list)
    simplices: list[dict] = field(default_factory=list)

    def as_legacy(self) -> list:
        return [self.dimensions, self.births, self.lengths, self.generators]


@dataclass
class ExtendedBarcodeResult:
    """Aggregated result across many subsamples."""

    intervals: list[H1Interval]
    count_matrix: np.ndarray
    labels: list[str] | None = None
    positions: np.ndarray | None = None
    windows: list[tuple[int, ...]] = field(default_factory=list)
    barcodes: list[BarcodeResult] = field(default_factory=list)

    def as_legacy(self) -> list:
        return [
            [bc.as_legacy() for bc in self.barcodes],
            self.labels,
            self.positions,
            [
                [(int(i), int(j)), int(self.count_matrix[i, j])]
                for i in range(self.count_matrix.shape[0])
                for j in range(i + 1, self.count_matrix.shape[1])
                if self.count_matrix[i, j] > 0
            ],
            self.count_matrix,
        ]


@dataclass
class GroupAnalysis:
    """All outputs for one location file."""

    group_id: str
    location: str
    lineages: list[str]
    path: str
    labels: list[str]
    seqs: list[str]
    segregating_positions: np.ndarray
    compressed: list[str]
    extended: ExtendedBarcodeResult


@dataclass
class ARGNode:
    node_id: str
    node_type: str
    label: str
    location: str | None = None
    lineage: str | None = None
    support: int = 0
    birth: float | None = None
    death: float | None = None
    span: str | None = None


@dataclass
class ARGResult:
    """Graph-ready ARG summary."""

    nodes: list[ARGNode]
    edges: list[tuple[str, str, dict]]
    groups: list[GroupAnalysis]
    graphml_path: str | None = None
    png_path: str | None = None


class DistanceMatrix:
    """A lightweight wrapper matching the old callable distance object."""

    def __init__(self, matrix: Sequence[Sequence[float]], max_index: int | None = None):
        self.m = np.asarray(matrix, dtype=float)
        self.n = self.m.shape[0]
        self.maxi = int(max_index if max_index is not None else self.n - 1)

    def __len__(self) -> int:
        return self.n

    def __call__(self, x: int, y: int) -> float:
        return float(self.m[x][y])


# -----------------------------
# FASTA / manifest handling
# -----------------------------


def load_fasta(path: str | os.PathLike[str]) -> tuple[list[str], list[str]]:
    """Read a minimal FASTA file; all sequences must be aligned."""
    labels: list[str] = []
    seqs: list[str] = []
    current: list[str] = []
    have_header = False

    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if have_header:
                    seqs.append("".join(current).upper())
                    current = []
                labels.append(line[1:].strip())
                have_header = True
            else:
                if not have_header:
                    raise ValueError("FASTA record missing header before sequence data")
                current.append(line.replace(" ", "").upper())
    if have_header:
        seqs.append("".join(current).upper())

    if not seqs:
        raise ValueError("No sequences found in FASTA file")

    length = len(seqs[0])
    if any(len(s) != length for s in seqs):
        raise ValueError("All sequences must be aligned and have the same length")
    return labels, seqs


def _split_list(value: str) -> list[str]:
    parts = [p.strip() for p in re.split(r"[;,|]+", value or "") if p.strip()]
    return parts


def load_manifest(path: str | os.PathLike[str]) -> list[dict[str, object]]:
    """Load a TSV/CSV manifest describing per-location sequence files.

    Expected columns: file/path, location, lineages. Optional: group_id.
    """
    path = str(path)
    with open(path, "r", encoding="utf-8") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        except csv.Error:
            dialect = csv.get_dialect("excel-tab")
        reader = csv.DictReader(f, dialect=dialect)
        rows = list(reader)

    out: list[dict[str, object]] = []
    for idx, row in enumerate(rows):
        norm = {k.strip().lower(): (v.strip() if isinstance(v, str) else v) for k, v in row.items() if k is not None}
        file_key = norm.get("file") or norm.get("path") or norm.get("fasta") or norm.get("alignment")
        location = norm.get("location") or norm.get("country") or norm.get("site")
        lineages = norm.get("lineages") or norm.get("lineage") or norm.get("pango") or ""
        if not file_key or not location:
            raise ValueError(f"Manifest row {idx+2} must contain at least file/path and location columns")
        group_id = norm.get("group_id") or norm.get("id") or Path(str(file_key)).stem
        out.append(
            {
                "group_id": str(group_id),
                "path": str(file_key),
                "location": str(location),
                "lineages": _split_list(str(lineages)),
            }
        )
    return out


# -----------------------------
# Site filtering / distances
# -----------------------------


def segregating_sites(seqs: Sequence[str], exclude_compatible: bool = False) -> tuple[np.ndarray, list[str]]:
    """Return segregating-site indices and compressed haplotypes."""
    arr = np.array([list(s) for s in seqs], dtype="U1")
    positions: list[int] = []
    compressed = ["" for _ in seqs]

    for idx, col in enumerate(arr.T):
        if len(set(col)) <= 1:
            continue
        if "-" in col:
            continue

        compatible = False
        if exclude_compatible:
            for other in arr.T:
                haplotypes = {a + b for a, b in zip(col, other)}
                if len(haplotypes) > 3:
                    compatible = True
                    break
        else:
            compatible = True

        if compatible:
            positions.append(idx)
            for i, base in enumerate(col):
                compressed[i] += base

    return np.array(positions, dtype=int), compressed


def transverse(str1: str, str2: str) -> tuple[int, int]:
    transitions = {"AG", "GA", "CT", "TC"}
    transversions = {"AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"}
    ctransitions = 0
    ctransversions = 0
    for m1, m2 in zip(str1, str2):
        p = (m1 + m2).upper()
        if p in transitions:
            ctransitions += 1
        elif p in transversions:
            ctransversions += 1
    return ctransitions, ctransversions


def _distance(str1: str, str2: str, model: str = "hamming", gc_fraction: float | None = None) -> float:
    n = len(str1)
    diffs = sum(a != b for a, b in zip(str1, str2))
    p = diffs / float(max(n, 1))
    if model == "hamming":
        return float(diffs)
    if model == "jc69":
        val = 1.0 - (4.0 * p / 3.0)
        return float("inf") if val <= 0 else float(-0.75 * math.log(val))
    if model == "k80":
        ts, tv = transverse(str1, str2)
        p_ts = ts / float(max(n, 1))
        p_tv = tv / float(max(n, 1))
        a = 1.0 - (2.0 * p_ts) - p_tv
        b = 1.0 - (2.0 * p_tv)
        if a <= 0 or b <= 0:
            return float("inf")
        return float(-0.5 * math.log(a) - 0.25 * math.log(b))
    if model == "t92":
        if gc_fraction is None:
            raise ValueError("gc_fraction is required for the T92 model")
        ts, tv = transverse(str1, str2)
        p_ts = ts / float(max(n, 1))
        p_tv = tv / float(max(n, 1))
        h = 2.0 * gc_fraction * (1.0 - gc_fraction)
        a = 1.0 - (p_ts / (h if h else 1.0)) - p_tv
        b = 1.0 - (2.0 * p_tv)
        if a <= 0 or b <= 0 or h == 0:
            return float("inf")
        return float(-h * math.log(a) - 0.5 * (1.0 - h) * math.log(b))
    raise ValueError(f"Unknown distance model: {model}")


def dist_list(seqs: Sequence[str], positions: tuple[int, int] | None = None, model: str = "hamming") -> DistanceMatrix:
    n = len(seqs)
    mat = np.zeros((n, n), dtype=float)
    gc_fraction = None
    if model == "t92":
        x = "".join(seqs)
        gc = sum(c in "GgCc" for c in x)
        at = sum(c in "AaTt" for c in x)
        gc_fraction = 0.5 if gc + at == 0 else gc / float(gc + at)
    for i in range(n):
        for j in range(i, n):
            d = _distance(seqs[i], seqs[j], model=model, gc_fraction=gc_fraction)
            mat[i, j] = mat[j, i] = d
    return DistanceMatrix(mat, max_index=n - 1 if positions is None else positions[1] - positions[0])


# -----------------------------
# Persistent homology core (up to H1)
# -----------------------------


def _edge_key(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


@dataclass(frozen=True)
class _Simplex:
    vertices: tuple[int, ...]
    dim: int
    filtration: float


@dataclass
class _Column:
    simplex_index: int
    dim: int
    filtration: float
    value: int


def _build_simplices(dist: np.ndarray, threshold: float) -> tuple[list[_Simplex], dict[tuple[int, ...], int]]:
    n = dist.shape[0]
    simplices: list[_Simplex] = []
    for i in range(n):
        simplices.append(_Simplex((i,), 0, 0.0))

    edges: list[tuple[tuple[int, int], float]] = []
    for i, j in combinations(range(n), 2):
        d = float(dist[i, j])
        if d <= threshold:
            edges.append((_edge_key(i, j), d))
    edges.sort(key=lambda x: (x[1], x[0]))
    for e, d in edges:
        simplices.append(_Simplex(e, 1, d))

    edge_set = {s.vertices for s in simplices if s.dim == 1}
    tri: list[tuple[tuple[int, int, int], float]] = []
    for i, j, k in combinations(range(n), 3):
        e1 = _edge_key(i, j)
        e2 = _edge_key(i, k)
        e3 = _edge_key(j, k)
        if e1 in edge_set and e2 in edge_set and e3 in edge_set:
            d = max(dist[i, j], dist[i, k], dist[j, k])
            if d <= threshold:
                tri.append(((i, j, k), float(d)))
    tri.sort(key=lambda x: (x[1], x[0]))
    for t, d in tri:
        simplices.append(_Simplex(t, 2, d))

    simplices.sort(key=lambda s: (s.filtration, s.dim, s.vertices))
    index = {s.vertices: i for i, s in enumerate(simplices)}
    return simplices, index


def _boundary_bitset(simplex: _Simplex, index: dict[tuple[int, ...], int]) -> int:
    if simplex.dim == 0:
        return 0
    bits = 0
    verts = simplex.vertices
    for face in combinations(verts, simplex.dim):
        bits ^= 1 << index[tuple(face)]
    return bits


def compute_barcode(
    distances2: DistanceMatrix | np.ndarray | Sequence[Sequence[float]],
    skeleton: int = 2,
    threshold: float | None = None,
    legacy: bool = True,
) -> BarcodeResult | list:
    if skeleton < 2:
        raise ValueError("RecombNet requires skeleton >= 2")

    dist = np.asarray(getattr(distances2, "m", distances2), dtype=float)
    if threshold is None:
        finite = dist[np.isfinite(dist)]
        threshold = float(np.max(finite)) if finite.size else 0.0

    simplices, index = _build_simplices(dist, threshold)
    columns: list[_Column] = []
    for idx, s in enumerate(simplices):
        columns.append(_Column(idx, s.dim, s.filtration, _boundary_bitset(s, index)))

    reduced: dict[int, int] = {}
    low_to_col: dict[int, int] = {}
    pair: dict[int, int] = {}

    for col in columns:
        value = col.value
        while value:
            low = value.bit_length() - 1
            if low in low_to_col:
                value ^= reduced[low_to_col[low]]
            else:
                break
        col.value = value
        reduced[col.simplex_index] = value
        if value:
            low = value.bit_length() - 1
            low_to_col[low] = col.simplex_index
            pair[low] = col.simplex_index

    edge_births = {
        idx: simplices[idx].filtration
        for idx, s in enumerate(simplices)
        if s.dim == 1 and reduced[idx] == 0
    }

    intervals: list[H1Interval] = []
    dimensions: list[int] = []
    births: list[float] = []
    lengths: list[float] = []
    generators: list[list[tuple[tuple[int, int], float]]] = []

    edge_by_index = {idx: simplices[idx].vertices for idx, s in enumerate(simplices) if s.dim == 1}

    for edge_idx, birth_time in sorted(edge_births.items(), key=lambda x: (x[1], x[0])):
        death_idx = pair.get(edge_idx)
        death_time = simplices[death_idx].filtration if death_idx is not None else None
        gen_edges: tuple[tuple[int, int], ...] = ()
        if death_idx is not None:
            support: list[tuple[int, int]] = []
            bitset = reduced[death_idx]
            for idx, edge in edge_by_index.items():
                if bitset & (1 << idx):
                    support.append(edge)
            gen_edges = tuple(sorted(support))
        interval = H1Interval(
            birth=float(birth_time),
            death=float(death_time) if death_time is not None else None,
            birth_edge=edge_by_index[edge_idx],
            death_triangle=simplices[death_idx].vertices if death_idx is not None and simplices[death_idx].dim == 2 else None,
            generator_edges=gen_edges,
        )
        intervals.append(interval)
        dimensions.append(1)
        births.append(float(birth_time))
        lengths.append(interval.length)
        generators.append([(e, 1.0) for e in gen_edges])

    result = BarcodeResult(
        dimensions=dimensions,
        births=births,
        lengths=lengths,
        generators=generators,
        intervals=intervals,
        simplices=[{"vertices": s.vertices, "dim": s.dim, "filtration": s.filtration} for s in simplices],
    )
    return result.as_legacy() if legacy else result


# -----------------------------
# Subsampling / aggregation
# -----------------------------


def subsample(
    seqs: Sequence[str],
    max_sites: int,
    window: int,
    min_sites: int,
) -> Iterator[tuple[tuple[int, ...], list[str]]]:
    if not seqs:
        return
    n_sites = len(seqs[0])
    window = max(1, min(window, n_sites))
    max_sites = max(1, max_sites)
    min_sites = max(1, min_sites)
    positions = list(range(n_sites))

    for start in range(0, n_sites - window + 1):
        stop = start + window
        window_positions = positions[start:stop]
        upper = min(max_sites, len(window_positions))
        for k in range(min_sites, upper + 1):
            for combo in combinations(window_positions, k):
                compressed = ["" for _ in seqs]
                for q in combo:
                    for i, s in enumerate(seqs):
                        compressed[i] += s[q]
                yield tuple(combo), compressed


def compute_extended_barcode(
    seqs: Sequence[str],
    max_sites: int,
    window: int,
    min_sites: int,
    cores: int = 1,
    positions: np.ndarray | None = None,
    labels: list[str] | None = None,
    model: str = "hamming",
    exclude_compatible: bool = False,
    legacy: bool = False,
) -> ExtendedBarcodeResult | list:
    if exclude_compatible:
        raise NotImplementedError(
            "exclude_compatible is handled at the site-filtering stage; pass filtered sequences instead"
        )
    if not seqs:
        raise ValueError("No sequences provided")

    n_sites = len(seqs[0])
    count_matrix = np.zeros((n_sites, n_sites), dtype=int)
    barcodes: list[BarcodeResult] = []
    intervals: list[H1Interval] = []
    windows: list[tuple[int, ...]] = []

    for combo, compressed in subsample(seqs, max_sites=max_sites, window=window, min_sites=min_sites):
        windows.append(combo)
        dist = dist_list(compressed, model=model)
        bc = compute_barcode(dist, legacy=False)
        assert isinstance(bc, BarcodeResult)
        barcodes.append(bc)

        if combo and bc.intervals:
            i, j = combo[0], combo[-1]
            count_matrix[min(i, j), max(i, j)] += len(bc.intervals)
        intervals.extend(bc.intervals)

    result = ExtendedBarcodeResult(
        intervals=intervals,
        count_matrix=count_matrix,
        labels=labels,
        positions=positions,
        windows=windows,
        barcodes=barcodes,
    )
    return result.as_legacy() if legacy else result


# -----------------------------
# Utility summaries
# -----------------------------


def b1_matrix(rr: ExtendedBarcodeResult | np.ndarray | list) -> list[list[int]]:
    if isinstance(rr, ExtendedBarcodeResult):
        return rr.count_matrix.tolist()
    if isinstance(rr, np.ndarray):
        return rr.astype(int).tolist()
    if isinstance(rr, list) and rr and isinstance(rr[-1], np.ndarray):
        return rr[-1].astype(int).tolist()
    raise TypeError("Unsupported result type for b1_matrix")


def breaking_points(rr: ExtendedBarcodeResult | np.ndarray | list) -> list[int]:
    mat = np.asarray(b1_matrix(rr), dtype=int)
    return [i for i in range(max(0, mat.shape[0] - 1)) if mat[i, i + 1] == 1]


# -----------------------------
# Group analysis and ARG graph construction
# -----------------------------


def analyze_group(
    fasta_path: str | os.PathLike[str],
    group_id: str,
    location: str,
    lineages: Sequence[str],
    *,
    positions: np.ndarray | None = None,
    max_sites: int = 5,
    window: int = 13,
    min_sites: int = 2,
    model: str = "hamming",
    exclude_compatible: bool = False,
) -> GroupAnalysis:
    labels, seqs = load_fasta(fasta_path)
    poss, compressed = segregating_sites(seqs, exclude_compatible=exclude_compatible)
    if poss.size == 0:
        raise ValueError(f"No segregating sites found in {fasta_path}")
    if positions is not None:
        if len(positions) != len(seqs[0]):
            raise ValueError("positions file length does not match the alignment length")
        positions = positions[poss]
    ext = compute_extended_barcode(
        compressed,
        max_sites=max_sites,
        window=window,
        min_sites=min_sites,
        positions=positions,
        labels=labels,
        model=model,
        exclude_compatible=False,
        legacy=False,
    )
    assert isinstance(ext, ExtendedBarcodeResult)
    return GroupAnalysis(
        group_id=group_id,
        location=location,
        lineages=[x for x in lineages if x],
        path=str(fasta_path),
        labels=labels,
        seqs=seqs,
        segregating_positions=poss,
        compressed=compressed,
        extended=ext,
    )


def _group_key(group: GroupAnalysis) -> str:
    return f"{group.location}::{group.group_id}"


def _event_span(interval: H1Interval, group: GroupAnalysis) -> str:
    if group.extended.positions is not None and len(group.extended.positions) > 0:
        b = min(interval.birth_edge)
        e = max(interval.birth_edge)
        return f"{int(group.extended.positions[b])}-{int(group.extended.positions[e])}"
    return f"{min(interval.birth_edge)}-{max(interval.birth_edge)}"


def build_arg_graph(
    groups: Sequence[GroupAnalysis],
    *,
    top_events: int = 12,
    min_interval_length: float = 0.0,
    merge_shared_lineages: bool = True,
) -> tuple[list[ARGNode], list[tuple[str, str, dict]]]:
    """Build a graph that visualizes recombination signals across locations and lineages.

    The graph is a practical, visual ARG: location nodes connect to lineage nodes,
    and each location-specific recombination event becomes an event node linked to the
    location and its lineages.
    """
    nodes: list[ARGNode] = []
    edges: list[tuple[str, str, dict]] = []
    seen: set[str] = set()

    def add_node(node: ARGNode) -> None:
        if node.node_id not in seen:
            seen.add(node.node_id)
            nodes.append(node)

    lineage_nodes: dict[str, str] = {}
    for group in groups:
        gid = _group_key(group)
        add_node(ARGNode(node_id=f"group:{gid}", node_type="group", label=group.location, location=group.location, support=len(group.labels)))
        for lineage in group.lineages:
            lid = f"lineage:{lineage}"
            if lineage not in lineage_nodes:
                lineage_nodes[lineage] = lid
                add_node(ARGNode(node_id=lid, node_type="lineage", label=lineage, lineage=lineage))
            if merge_shared_lineages:
                edges.append((lid, f"group:{gid}", {"kind": "present_in"}))

        intervals = [iv for iv in group.extended.intervals if iv.length >= min_interval_length]
        intervals.sort(key=lambda iv: (-(iv.length if math.isfinite(iv.length) else 1e18), iv.birth))
        for idx, interval in enumerate(intervals[:top_events]):
            span = _event_span(interval, group)
            eid = f"event:{gid}:{idx}"
            add_node(
                ARGNode(
                    node_id=eid,
                    node_type="event",
                    label=f"event {idx+1}",
                    location=group.location,
                    support=len(interval.generator_edges),
                    birth=interval.birth,
                    death=interval.death,
                    span=span,
                )
            )
            edges.append((f"group:{gid}", eid, {"kind": "recombination_signal", "birth": interval.birth, "length": interval.length, "span": span}))
            for lineage in group.lineages:
                edges.append((f"lineage:{lineage}", eid, {"kind": "associated_lineage"}))

    return nodes, edges


def render_arg_graph(
    groups: Sequence[GroupAnalysis],
    *,
    output_prefix: str,
    top_events: int = 12,
    min_interval_length: float = 0.0,
    merge_shared_lineages: bool = True,
) -> ARGResult:
    import networkx as nx
    import matplotlib.pyplot as plt

    nodes, edges = build_arg_graph(
        groups,
        top_events=top_events,
        min_interval_length=min_interval_length,
        merge_shared_lineages=merge_shared_lineages,
    )

    G = nx.Graph()
    for node in nodes:
        G.add_node(node.node_id, **asdict(node))
    for u, v, attrs in edges:
        G.add_edge(u, v, **attrs)

    # Manual layered layout: groups on the left, lineages in the middle, events on the right.
    pos: dict[str, tuple[float, float]] = {}
    groups_nodes = [n for n in nodes if n.node_type == "group"]
    lineage_nodes = [n for n in nodes if n.node_type == "lineage"]
    event_nodes = [n for n in nodes if n.node_type == "event"]

    def spaced_y(count: int) -> list[float]:
        if count <= 1:
            return [0.0]
        return list(np.linspace(0.95, 0.05, count))

    for node, y in zip(groups_nodes, spaced_y(len(groups_nodes))):
        pos[node.node_id] = (0.05, y)
    for node, y in zip(lineage_nodes, spaced_y(len(lineage_nodes))):
        pos[node.node_id] = (0.5, y)
    for node, y in zip(event_nodes, spaced_y(len(event_nodes))):
        pos[node.node_id] = (0.95, y)

    fig, ax = plt.subplots(figsize=(16, max(7, 0.45 * max(len(nodes), 8))))
    node_colors = []
    node_shapes = {"group": "s", "lineage": "o", "event": "D"}
    for node in nodes:
        if node.node_type == "group":
            node_colors.append("#7da0fa")
        elif node.node_type == "lineage":
            node_colors.append("#7bc67e")
        else:
            node_colors.append("#f7b267")

    # Draw edges by type.
    group_edges = [(u, v) for u, v, a in edges if a.get("kind") == "present_in"]
    signal_edges = [(u, v) for u, v, a in edges if a.get("kind") == "recombination_signal"]
    lineage_edges = [(u, v) for u, v, a in edges if a.get("kind") == "associated_lineage"]

    if group_edges:
        nx.draw_networkx_edges(G, pos, edgelist=group_edges, width=1.4, alpha=0.55, ax=ax)
    if lineage_edges:
        nx.draw_networkx_edges(G, pos, edgelist=lineage_edges, width=0.9, alpha=0.35, style="dashed", ax=ax)
    if signal_edges:
        nx.draw_networkx_edges(G, pos, edgelist=signal_edges, width=2.2, alpha=0.8, edge_color="#d62828", ax=ax)

    # Draw nodes one shape at a time.
    for node_type, shape in node_shapes.items():
        subset = [n.node_id for n in nodes if n.node_type == node_type]
        if not subset:
            continue
        subset_colors = [
            "#7da0fa" if G.nodes[n]["node_type"] == "group" else "#7bc67e" if G.nodes[n]["node_type"] == "lineage" else "#f7b267"
            for n in subset
        ]
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=subset,
            node_shape=shape,
            node_color=subset_colors,
            node_size=1400 if node_type == "group" else 1100 if node_type == "lineage" else 1250,
            edgecolors="black",
            linewidths=1.0,
            ax=ax,
        )

    labels = {n.node_id: n.label for n in nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9, ax=ax)

    # Add event annotations.
    for node in event_nodes:
        txt = []
        if node.location:
            txt.append(node.location)
        if node.span:
            txt.append(f"span {node.span}")
        if node.birth is not None:
            txt.append(f"birth {node.birth:.3f}")
        if node.death is not None:
            txt.append(f"death {node.death:.3f}")
        if txt and node.node_id in pos:
            ax.text(pos[node.node_id][0] + 0.015, pos[node.node_id][1] - 0.03, "\n".join(txt), fontsize=7, ha="left", va="top")

    ax.set_axis_off()
    ax.set_title("RecombNet visual ARG", fontsize=18)
    fig.tight_layout()
    png_path = output_prefix + ".arg.png"
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    graphml_path = output_prefix + ".arg.graphml"
    nx.write_graphml(G, graphml_path)

    return ARGResult(nodes=nodes, edges=edges, groups=list(groups), graphml_path=graphml_path, png_path=png_path)


# -----------------------------
# Serialization helpers
# -----------------------------


def save_result(obj, path: str | os.PathLike[str]) -> None:
    path = str(path)
    if path.endswith(".json"):
        with open(path, "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=2, default=lambda o: asdict(o) if hasattr(o, "__dataclass_fields__") else str(o))
    else:
        with open(path, "wb") as f:
            pickle.dump(obj, f)


def load_result(path: str | os.PathLike[str]):
    with open(path, "rb") as f:
        return pickle.load(f)
