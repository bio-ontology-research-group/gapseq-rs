# gapseq-rs

Rust reimplementation of [gapseq](https://github.com/jotech/gapseq).
For a detailed comparison with the original R/bash implementation, see
[COMPARISON.md](COMPARISON.md).

## What it does

gapseq-rs reconstructs genome-scale metabolic models from bacterial
proteomes. Given a protein FASTA, it predicts metabolic pathways, detects
transporters, assembles a draft stoichiometric model, infers a growth
medium, and gap-fills the model so it can simulate growth.

```
gapseq doall genome.faa.gz -f output/
```

This produces a gap-filled SBML model that loads directly in COBRApy,
COBRAToolbox, or any SBML-compatible tool.

## Install

### Prerequisites

An external sequence aligner (pick one):

| Tool | Install |
|------|---------|
| BLAST+ | `apt install ncbi-blast+` or `conda install -c bioconda blast` |
| DIAMOND | `apt install diamond` or `conda install -c bioconda diamond` |
| MMseqs2 | `apt install mmseqs2` or `conda install -c bioconda mmseqs2` |

### Reference data

```bash
git clone --depth 1 https://github.com/jotech/gapseq.git
```

Then download the reference sequence database:

```bash
gapseq update-sequences -t Bacteria
```

### Build from source

```bash
git clone https://github.com/bio-ontology-research-group/gapseq-rs.git
cd gapseq-rs
cargo build --release
# Binary: target/release/gapseq
```

## Quick start

```bash
# Full reconstruction pipeline (find → transport → draft → medium → fill)
gapseq --data-dir gapseq/dat doall genome.faa.gz -f output/ -A diamond

# Step by step
gapseq --data-dir gapseq/dat find -p all -A diamond -o output/ genome.faa
gapseq --data-dir gapseq/dat find-transport -A diamond -o output/ genome.faa
gapseq --data-dir gapseq/dat draft -r output/*-Reactions.tbl -t output/*-Transporter.tbl -o output/
gapseq --data-dir gapseq/dat medium -m output/*-draft.gmod.cbor -p output/*-Pathways.tbl
gapseq --data-dir gapseq/dat fill output/*-draft.gmod.cbor -n output/*-medium.csv -r output/*-Reactions.tbl -o output/
```

### Output files

| File | Contents |
|------|----------|
| `*-all-Reactions.tbl` | Per-reaction homology hits + pathway context |
| `*-all-Pathways.tbl` | Pathway completeness predictions |
| `*-Transporter.tbl` | Detected transporters |
| `*-draft.gmod.cbor` | Draft model (native format) |
| `*-draft.xml` | Draft model (SBML L3V1 + FBC2 + groups) |
| `*-medium.csv` | Predicted growth medium |
| `*-filled.gmod.cbor` | Gap-filled model (native format) |
| `*-filled.xml` | Gap-filled model (SBML) |
| `*-filled-added.tsv` | Reactions added during gap-filling |

## Subcommands

| Command | Description |
|---------|-------------|
| `doall` | Full pipeline: find → transport → draft → medium → fill |
| `find` | Pathway and reaction detection |
| `find-transport` | Transporter detection |
| `draft` | Build a draft metabolic model |
| `medium` | Rule-based growth medium inference |
| `fill` | Iterative gap-filling (pFBA + KO essentiality) |
| `fba` | FBA / pFBA on an existing model |
| `adapt` | Add/remove reactions or force growth on compounds |
| `pan` | Build a pan-draft model from multiple drafts |
| `update-sequences` | Sync reference sequence database from Zenodo |
| `convert` | Convert between CBOR and JSON model formats |
| `export-sbml` | Export a model as SBML |

Run any command with `-h` for full option documentation.

## Documentation

| Document | Contents |
|----------|----------|
| [User guide](docs/user-guide.md) | Install, quick-start, per-subcommand recipes, troubleshooting |
| [CLI reference](docs/cli-reference.md) | Every flag of every subcommand |
| [Architecture](docs/architecture.md) | Crate dependency graph, data flow, LP plumbing |
| [Feature matrix](docs/feature-matrix.md) | R source → Rust module mapping, status per feature |
| [Porting notes](docs/porting-notes.md) | Intentional deviations from upstream gapseq |
| [Comparison](COMPARISON.md) | Performance benchmarks and feature comparison with upstream |

## License

GPL-3.0-or-later — same as [gapseq](https://github.com/jotech/gapseq).

## Citation

If you use gapseq-rs, please cite the original gapseq paper:

> Zimmermann J, Kaleta C, Özbek Ö, et al. gapseq: informed prediction of
> bacterial metabolic pathways and reconstruction of accurate metabolic
> models. *Genome Biology* 22, 81 (2021).
> https://doi.org/10.1186/s13059-021-02295-1
