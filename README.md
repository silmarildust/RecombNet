# RecombNet

RecombNet is a Python 3 toolkit for reconstructing and visualizing ancestral recombination graphs (ARGs) from aligned SARS-CoV-2 sequence data using persistent homology.

It is designed for **Nextclade-aligned CSV files** named in the format:

```text
LINEAGE_LOCATION.csv
```

Examples:

```text
XBC.1_China.csv
BA.2_Philippines.csv
B.1.617.2_United_States_of_America.csv
```

## What RecombNet does

RecombNet supports two workflows:

1. **Single-file scan**
   - Read one aligned FASTA or CSV file.
   - Find segregating sites.
   - Compute persistent homology.
   - Write barcode summaries and a small summary plot.

2. **Visual ARG generation**
   - Read a manifest of CSV files, one file per lineage-location pair.
   - Infer metadata from the filename or from manifest columns.
   - Analyse each file separately.
   - Combine the recombination signals into one visual ARG.
   - Export a PNG figure, GraphML network, and JSON summary.

The ARG output is a **visual, exploratory reconstruction** based on topological recombination signals. It is useful for comparing clusters and lineage/location patterns, but it is not a clinical or phylogenetic ground truth.

## Installation

From the project directory:

```bash
pip install -e .
```

This installs the command-line tool:

```bash
recombnet
```

## Input format

### 1) CSV files

RecombNet now accepts Nextclade-style `.csv` files directly.

The package looks for an explicit aligned-sequence column first. If your CSV has one, RecombNet will use it automatically. If your CSV does **not** contain an explicit sequence column, pass a reference FASTA with `--reference` so RecombNet can reconstruct a best-effort aligned sequence from the Nextclade mutation columns.

The filename should follow:

```text
LINEAGE_LOCATION.csv
```

For your project, that means files such as:

```text
XBC.1_China.csv
XBC.1_Philippines.csv
XBC.1_Singapore.csv
XBC.1_South_Korea.csv
XBC.1_United_States_of_America.csv
XBE_United_Kingdom.csv
XBE_United_States_of_America.csv
XBZ_Austria.csv
XBZ_Denmark.csv
XBZ_Germany.csv
```

### 2) Required columns in the Nextclade CSV

Your file headers already include the key metadata columns RecombNet uses:

- `seqName` for sequence labels
- `Nextclade_pango` for lineage annotation when present
- `substitutions`, `deletions`, `insertions`, `missing`, and related fields for reconstruction when no explicit sequence column exists

If the CSV contains a direct sequence column, RecombNet will detect it automatically. If not, use:

```bash
--reference data/reference.fasta
```

with the same SARS-CoV-2 reference used for the Nextclade alignment.

### 3) Manifest file for ARG mode

You can provide a minimal manifest with just a `path` column, because RecombNet infers lineage and location from the filename.

Example TSV:

```tsv
path
XBC.1_China.csv
XBC.1_Philippines.csv
XBC.1_Singapore.csv
XBC.1_South_Korea.csv
XBC.1_United_States_of_America.csv
```

You may also include optional columns such as `group_id`, `location`, and `lineages` if you want to override the filename-based inference.

## Single-file scan mode

Use this when you want persistent-homology summaries for one aligned file.

```bash
recombnet scan data/XBC.1_China.csv -o china_scan --reference data/reference.fasta
```

If your CSV already contains an explicit aligned-sequence column, you can omit `--reference`.

Useful options:

```bash
recombnet scan data/XBC.1_China.csv \
  -o china_scan \
  --reference data/reference.fasta \
  -s 5 \
  -t 2 \
  -w 13
```

### Single-file outputs

RecombNet writes:

- `china_scan.bars.tsv` — H1 intervals
- `china_scan.b1.tsv` — H1 support-count matrix
- `china_scan.pkl` — saved Python object
- `china_scan.json` — run metadata
- `china_scan.png` — summary figure

## Visual ARG mode

This is the main workflow for your SARS-CoV-2 analysis.

### Your study design

Use three separate ARG runs, one for each lineage/location panel:

**Panel A**
- Lineages: `XBC.1`, `BA.2`, `B.1.617.2`
- Locations: `China`, `Philippines`, `Singapore`, `South Korea`, `United States of America`

**Panel B**
- Lineages: `XBE`, `BA.5.2`, `BE.4.1`
- Locations: `United Kingdom`, `United States of America`

**Panel C**
- Lineages: `XBZ`, `BA.5.2.1`, `EF.1.3`
- Locations: `Austria`, `Denmark`, `Germany`

### Example manifest files

#### `manifest_panel_a.tsv`

```tsv
path
XBC.1_China.csv
XBC.1_Philippines.csv
XBC.1_Singapore.csv
XBC.1_South_Korea.csv
XBC.1_United_States_of_America.csv
```

#### `manifest_panel_b.tsv`

```tsv
path
XBE_United_Kingdom.csv
XBE_United_States_of_America.csv
```

#### `manifest_panel_c.tsv`

```tsv
path
XBZ_Austria.csv
XBZ_Denmark.csv
XBZ_Germany.csv
```

### Run the ARG builder

```bash
recombnet arg manifest_panel_a.tsv -o panel_a_arg --reference data/reference.fasta
recombnet arg manifest_panel_b.tsv -o panel_b_arg --reference data/reference.fasta
recombnet arg manifest_panel_c.tsv -o panel_c_arg --reference data/reference.fasta
```

### ARG outputs

For each run, RecombNet writes:

- `PREFIX.arg.png` — the visual ARG figure
- `PREFIX.arg.graphml` — the network in GraphML format
- `PREFIX.arg.json` — a summary of the groups and settings

## How the ARG figure is organized

The visual ARG uses layered nodes:

- **Group nodes** represent each input file / lineage-location pair.
- **Lineage nodes** represent the lineages associated with the file.
- **Event nodes** represent the strongest persistent-H1 recombination signals detected in that file.

Edges indicate:

- which lineage appears in which group,
- which recombination event was detected in that group,
- and how the event is associated with the group’s lineage labels.

## Using `--sequence-column`

If your CSV has a direct aligned-sequence column with a non-standard name, point RecombNet to it explicitly:

```bash
recombnet scan XBC.1_China.csv -o china_scan --sequence-column alignedSequence
```

The same option is available in `arg` mode.

## Troubleshooting

### "No usable aligned sequences found"

This usually means the CSV does not include an explicit sequence column and you did not provide `--reference`.

### "Sequences are not aligned to the same length"

All input sequences must be aligned and have equal length.

### "positions file length does not match the alignment length"

The genomic positions file must contain one integer per alignment column.

## Development notes

RecombNet is a Python 3 rewrite of a legacy TARGet-style workflow. The package was modernized to:

- remove Python 2-only dependencies,
- support CSV-based Nextclade workflows,
- infer lineage/location from file naming conventions,
- preserve the persistent-homology-driven barcode analysis,
- and generate visual ARGs suitable for comparative SARS-CoV-2 analysis.

## License

GPL-3.0-or-later
