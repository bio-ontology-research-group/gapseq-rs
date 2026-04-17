# gapsmith

Rust reimplementation of [gapseq](https://github.com/jotech/gapseq).
For a detailed comparison with the original R/bash implementation, see
[COMPARISON.md](COMPARISON.md).

## What it does

gapsmith reconstructs genome-scale metabolic models from bacterial
proteomes. Given a protein FASTA, it predicts metabolic pathways, detects
transporters, assembles a draft stoichiometric model, infers a growth
medium, and gap-fills the model so it can simulate growth.

```
gapsmith doall genome.faa.gz -f output/
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

Plus a C++ toolchain (`cmake`, `gcc`/`clang`) for the bundled HiGHS LP
solver — any Linux / macOS system with dev-tools installed will do.

### Option 1: pre-built binary (fastest)

```bash
# Pick the right target for your OS/arch; see Releases page.
TARGET=x86_64-unknown-linux-gnu
VER=$(curl -s https://api.github.com/repos/bio-ontology-research-group/gapsmith/releases/latest \
  | grep '"tag_name"' | head -1 | sed 's/.*"\(v[^"]*\)".*/\1/')
curl -L https://github.com/bio-ontology-research-group/gapsmith/releases/download/$VER/gapsmith-$VER-$TARGET.tar.gz | tar xz
cd gapsmith-$VER-$TARGET
./gapsmith --version
```

Each release tarball bundles the binary + curated data tables. See the
[releases page](https://github.com/bio-ontology-research-group/gapsmith/releases).

### Option 2: cargo install (directly from git)

```bash
cargo install --git https://github.com/bio-ontology-research-group/gapsmith.git gapsmith-cli
```

Installs `gapsmith` into `~/.cargo/bin/`. You still need the `data/`
curation tables — clone the repo or grab them from a release tarball.

### Option 3: build from source

```bash
git clone https://github.com/bio-ontology-research-group/gapsmith.git
cd gapsmith
cargo build --release
# Binary: target/release/gapsmith, curated data in ./data/
```

### Reference data

Three parts, fetched independently:

1. **Curation tables** (subex, medium rules, biomass templates, …) —
   vendored in this repo under `data/`. ~1 MB. Auto-used when running
   from a checkout; bundled inside release tarballs.
2. **Large public reference tables** (SEED reactions + metabolites,
   MNXref cross-refs, ~65 MB) — fetched on demand from upstream gapseq's
   GitHub mirror:

    ```bash
    gapsmith update-data -o path/to/dat
    ```
3. **Sequence database** (per-reaction FASTAs, ~2 GB) — downloaded
   from Zenodo on demand:

    ```bash
    gapsmith update-sequences -D path/to/dat/seq -t Bacteria
    ```

After that you have a complete data directory and no longer need any
upstream gapseq checkout. Point all subsequent invocations at it with
`--data-dir path/to/dat`.

License-restricted data (MetaCyc pathways, KEGG, BiGG, BRENDA, VMH) is
left opt-in; a forthcoming `--accept-license` flag will gate loading
those.

## Quick start

```bash
# Full reconstruction pipeline (find → transport → draft → medium → fill)
gapsmith --data-dir path/to/dat doall genome.faa.gz -f output/ -A diamond

# Step by step
gapsmith --data-dir path/to/dat find -p all -A diamond -o output/ genome.faa
gapsmith --data-dir path/to/dat find-transport -A diamond -o output/ genome.faa
gapsmith --data-dir path/to/dat draft -r output/*-Reactions.tbl -t output/*-Transporter.tbl -o output/
gapsmith --data-dir path/to/dat medium -m output/*-draft.gmod.cbor -p output/*-Pathways.tbl
gapsmith --data-dir path/to/dat fill output/*-draft.gmod.cbor -n output/*-medium.csv -r output/*-Reactions.tbl -o output/
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
| `update-data` | Fetch the large public reference tables (SEED, MNXref) |
| `convert` | Convert between CBOR and JSON model formats |
| `export-sbml` | Export a model as SBML |

Run any command with `-h` for full option documentation.

## Documentation

Full documentation is published at
**https://bio-ontology-research-group.github.io/gapsmith/**.

Local copies:

| Document | Contents |
|----------|----------|
| [User guide](docs/user-guide.md) | Install, quick-start, per-subcommand recipes, troubleshooting |
| [CLI reference](docs/cli-reference.md) | Every flag of every subcommand |
| [Architecture](docs/architecture.md) | Crate dependency graph, data flow, LP plumbing |
| [Feature matrix](docs/feature-matrix.md) | R source → Rust module mapping, status per feature |
| [Porting notes](docs/porting-notes.md) | Intentional deviations from upstream gapseq |
| [Performance](docs/performance.md) | Shipped optimisations, benchmarks, semantic-parity results |
| [Comparison](COMPARISON.md) | Performance benchmarks and feature comparison with upstream |

## License

GPL-3.0-or-later — same as [gapseq](https://github.com/jotech/gapseq).

## Citation

If you use gapsmith, please cite the original gapseq paper:

> Zimmermann J, Kaleta C, Özbek Ö, et al. gapseq: informed prediction of
> bacterial metabolic pathways and reconstruction of accurate metabolic
> models. *Genome Biology* 22, 81 (2021).
> https://doi.org/10.1186/s13059-021-02295-1
