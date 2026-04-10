# RecombNet

RecombNet is a Python 3 package for recombination screening with persistent homology and **visual ancestral recombination graph (ARG)** generation.

It is a modernized replacement for the old TARGet-style workflow and is designed to work on **Nextclade-aligned SARS-CoV-2 sequences**.

## What RecombNet does

RecombNet supports two related workflows:

1. **Single-file scan**
   - Read one aligned FASTA file.
   - Find segregating sites.
   - Scan windows of sites.
   - Compute H1 persistent homology.
   - Write barcode summaries and a compact plot.

2. **Visual ARG generation from multiple locations and lineages**
   - Read a manifest that lists one aligned FASTA file per location.
   - Analyse each location file.
   - Combine the recombination signals into a visual ARG.
   - Export a PNG figure and a GraphML network file.

The ARG here is a **visual, data-driven reconstruction** based on recombination signals from persistent homology. It is meant for exploratory analysis and reporting, not for clinical or forensic claims.

## Installation

From the project directory:

```bash
pip install -e .
```

This installs the command-line tool:

```bash
recombnet
```

## Input requirements

### 1) Aligned FASTA files from Nextclade

Each input FASTA must already be aligned. Nextclade output is a good fit.

All sequences in a file must have the same length.

Example:

```fasta
>seq_001
ATGCC...-
>seq_002
ATGTC...-
```

### 2) One file per location

For the ARG workflow, prepare **one aligned FASTA file per location**.

Example layout:

```text
data/
  china.fasta
  philippines.fasta
  singapore.fasta
  south_korea.fasta
  united_states.fasta
  united_kingdom.fasta
  austria.fasta
  denmark.fasta
  germany.fasta
```

### 3) Manifest file for ARG mode

RecombNet uses a TSV or CSV manifest with at least these columns:

- `file` or `path` — path to the aligned FASTA file
- `location` — location name
- `lineages` — comma-, semicolon-, or pipe-separated lineage labels

Optional columns:

- `group_id` — a custom label for the file/group

Example TSV:

```tsv
group_id	file	location	lineages
china_panel	data/china.fasta	China	XBC.1;BA.2;B.1.617.2
philippines_panel	data/philippines.fasta	Philippines	XBC.1;BA.2;B.1.617.2
singapore_panel	data/singapore.fasta	Singapore	XBC.1;BA.2;B.1.617.2
south_korea_panel	data/south_korea.fasta	South Korea	XBC.1;BA.2;B.1.617.2
usa_asia_panel	data/united_states.fasta	United States of America	XBC.1;BA.2;B.1.617.2
uk_panel	data/united_kingdom.fasta	United Kingdom	XBE;BA.5.2;BE.4.1
usa_europe_panel	data/united_states_uk.fasta	United States of America	XBE;BA.5.2;BE.4.1
austria_panel	data/austria.fasta	Austria	XBZ;BA.5.2.1;EF.1.3
denmark_panel	data/denmark.fasta	Denmark	XBZ;BA.5.2.1;EF.1.3
germany_panel	data/germany.fasta	Germany	XBZ;BA.5.2.1;EF.1.3
```

Important note: the lineages listed above are the **group labels for the ARG workflow**. If you have one large location file, you can either keep the whole file in one manifest row or split the data into filtered per-panel FASTA files before running RecombNet.

## Installation check

After installation, confirm that the command exists:

```bash
recombnet --help
```

You should see two subcommands:

- `scan` — single-file analysis
- `arg` — visual ARG generation from a manifest

---

# 1. Single-file scan mode

Use this when you want persistent-homology summaries for one aligned FASTA file.

```bash
recombnet scan data/china.fasta -o china_scan
```

Useful options:

```bash
recombnet scan data/china.fasta \
  -o china_scan \
  -s 5 \
  -t 2 \
  -w 13 \
  --model hamming
```

### Single-file outputs

RecombNet writes:

- `china_scan.bars.tsv` — H1 intervals
- `china_scan.b1.tsv` — H1 support-count matrix
- `china_scan.pkl` — saved Python object
- `china_scan.json` — run metadata
- `china_scan.png` — summary figure

---

# 2. Visual ARG mode

This is the main workflow for your SARS-CoV-2 use case.

## Recommended workflow for the SARS-CoV-2 data

You said your data are aligned by Nextclade and organized as **one file per location**.

Use a manifest that groups those files into the lineage/location panels you want to compare:

### Panel A
- **Lineages**: `XBC.1`, `BA.2`, `B.1.617.2`
- **Locations**: `China`, `Philippines`, `Singapore`, `South Korea`, `United States of America`

### Panel B
- **Lineages**: `XBE`, `BA.5.2`, `BE.4.1`
- **Locations**: `United Kingdom`, `United States of America`

### Panel C
- **Lineages**: `XBZ`, `BA.5.2.1`, `EF.1.3`
- **Locations**: `Austria`, `Denmark`, `Germany`

If you are analysing all three panels, make **three manifest files** or one combined manifest with clear `group_id` labels.

## Example manifest files

### `manifest_panel_a.tsv`

```tsv
group_id	file	location	lineages
china	data/china.fasta	China	XBC.1;BA.2;B.1.617.2
philippines	data/philippines.fasta	Philippines	XBC.1;BA.2;B.1.617.2
singapore	data/singapore.fasta	Singapore	XBC.1;BA.2;B.1.617.2
south_korea	data/south_korea.fasta	South Korea	XBC.1;BA.2;B.1.617.2
united_states_a	data/united_states_a.fasta	United States of America	XBC.1;BA.2;B.1.617.2
```

### `manifest_panel_b.tsv`

```tsv
group_id	file	location	lineages
united_kingdom	data/united_kingdom.fasta	United Kingdom	XBE;BA.5.2;BE.4.1
united_states_b	data/united_states_b.fasta	United States of America	XBE;BA.5.2;BE.4.1
```

### `manifest_panel_c.tsv`

```tsv
group_id	file	location	lineages
austria	data/austria.fasta	Austria	XBZ;BA.5.2.1;EF.1.3
denmark	data/denmark.fasta	Denmark	XBZ;BA.5.2.1;EF.1.3
germany	data/germany.fasta	Germany	XBZ;BA.5.2.1;EF.1.3
```

## Run the ARG builder

```bash
recombnet arg manifest_panel_a.tsv -o panel_a_arg
recombnet arg manifest_panel_b.tsv -o panel_b_arg
recombnet arg manifest_panel_c.tsv -o panel_c_arg
```

You can also run a combined manifest if you prefer one merged network.

## ARG outputs

For each run, RecombNet writes:

- `PREFIX.arg.png` — the visual ARG figure
- `PREFIX.arg.graphml` — the network in GraphML format
- `PREFIX.arg.json` — a summary of the groups and settings

## How to read the ARG figure

The figure uses three layers:

- **Location nodes** — each aligned FASTA file / geographic location
- **Lineage nodes** — the lineages listed for that file in the manifest
- **Event nodes** — persistent-homology recombination signals detected in that location file

Edges indicate:

- which lineages were included in a location file,
- which recombination events were detected in that file,
- and how the ARG-like visual summary connects them.

A denser event cluster usually means more persistent H1 signal in the corresponding location file.

---

# 3. Using genomic coordinates from Nextclade

If you have a file with genomic coordinates for the aligned columns, pass it with `-p`:

```bash
recombnet scan data/china.fasta -o china_scan -p positions.txt
```

For ARG mode, you can reuse the same positions file if every location file uses the same alignment coordinate system.

The positions file must contain **one integer per alignment column**.

Example:

```text
1
2
3
4
5
...
```

If no positions file is provided, RecombNet uses alignment-column indices.

---

# 4. Recommended parameter choices

For SARS-CoV-2 Nextclade-aligned data, a good starting point is:

```bash
-s 5 -t 2 -w 13 --model hamming
```

You can make the scan more conservative by increasing `-t` or decreasing `-w`.

If you want to reduce clutter in the ARG visualization, raise `--min-interval-length`:

```bash
recombnet arg manifest_panel_a.tsv -o panel_a_arg --min-interval-length 0.01
```

You can also reduce the number of displayed events per group with `--top-events`:

```bash
recombnet arg manifest_panel_a.tsv -o panel_a_arg --top-events 8
```

---

# 5. Python API

You can use RecombNet directly from Python:

```python
from recombnet import load_fasta, segregating_sites, compute_extended_barcode, analyze_group, render_arg_graph

labels, seqs = load_fasta("data/china.fasta")
poss, compressed = segregating_sites(seqs)
result = compute_extended_barcode(compressed, max_sites=5, window=13, min_sites=2, legacy=False)
print(result.count_matrix)
```

ARG workflow from Python:

```python
from recombnet import analyze_group, render_arg_graph

g1 = analyze_group(
    "data/china.fasta",
    group_id="china",
    location="China",
    lineages=["XBC.1", "BA.2", "B.1.617.2"],
    max_sites=5,
    window=13,
    min_sites=2,
)

g2 = analyze_group(
    "data/philippines.fasta",
    group_id="philippines",
    location="Philippines",
    lineages=["XBC.1", "BA.2", "B.1.617.2"],
    max_sites=5,
    window=13,
    min_sites=2,
)

render_arg_graph([g1, g2], output_prefix="demo_arg")
```

---

# 6. Troubleshooting

### "All sequences must be aligned and have the same length"
Your FASTA file contains sequences with different lengths. Re-run Nextclade alignment or trim them to the same reference coordinate system.

### "No segregating sites found after filtering"
The file may be too homogeneous, too short, or overly filtered. Try another location file or relax the site filtering.

### The ARG looks too busy
Use one or more of these:

- lower `--top-events`
- raise `--min-interval-length`
- split the analysis into separate manifests
- run one geographic panel at a time

### The input file already has Nextclade output
That is fine. RecombNet only needs the aligned FASTA file. If you also have Nextclade metadata TSVs, keep them alongside the FASTA files so you can map sequence labels back to lineage annotations.

---

# 7. Command summary

Single-file scan:

```bash
recombnet scan INPUT.fasta -o PREFIX
```

ARG mode:

```bash
recombnet arg MANIFEST.tsv -o PREFIX
```

## License

RecombNet is released under GPL-3.0-or-later for compatibility with the original project.
