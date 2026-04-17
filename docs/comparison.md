# Comparison with upstream gapseq

gapsmith is a Rust reimplementation of
[gapseq](https://github.com/jotech/gapseq) (R/bash, ~9k LOC). This
document covers performance benchmarks and feature differences.

## Benchmarks

Wall-clock comparison on four real bacterial proteomes. Hardware:
56-core Xeon, 128 GB RAM, Debian 13, NVMe SSD. Same aligner binary
(NCBI `blastp` 2.17.0 via bioconda). Same reference sequence database
(Zenodo v1.4, record 16908828).

### Test genomes

| Organism | Accession | Proteins |
|----------|-----------|-------:|
| *Candidatus Blochmannia floridanus* | GCF_000007725.1 | 517 |
| *Bacillus subtilis* 168 | GCF_000009045.1 | 4,237 |
| *Escherichia coli* K-12 MG1655 | GCF_000005845.2 | 4,300 |
| *Salmonella enterica* Typhimurium LT2 | GCF_000006945.2 | 4,554 |

### Results

| Genome | Stage | gapseq (R) | gapsmith | Speedup |
|--------|-------|----------:|----------:|--------:|
| *B. floridanus* (517) | `find -p all` | 117 s | 34 s | **3.5×** |
| *B. floridanus* (517) | `find-transport` | 9 s | 8 s | **1.2×** |
| *B. subtilis* (4,237) | `find -p all` | 205 s | 73 s | **2.8×** |
| *B. subtilis* (4,237) | `find-transport` | 25 s | 14 s | **1.8×** |
| *E. coli* (4,300) | `find -p all` | 218 s | 76 s | **2.9×** |
| *E. coli* (4,300) | `find-transport` | 27 s | 14 s | **1.9×** |
| *S. Typhimurium* (4,554) | `find -p all` | 211 s | 76 s | **2.8×** |
| *S. Typhimurium* (4,554) | `find-transport` | 27 s | 14 s | **1.9×** |

gapsmith also uses 35–40% less peak memory (e.g. 498 MB vs 786 MB for
*E. coli* `find`).

### Stages without R baseline

The `draft`, `medium`, and `fill` stages could not be benchmarked
against upstream R gapseq because the `cobrar` R package (which replaced
the archived `sybil`) failed to install on the benchmark host due to a
libsbml/libxml2 ABI conflict in the conda R 4.5 environment. This is a
packaging issue on the test rig, not a limitation of either tool.

gapsmith timings for these stages (*E. coli* K-12):

| Stage | Wall-time | Peak RSS |
|-------|----------:|---------:|
| `draft` | 0.6 s | 152 MB |
| `medium` | 0.1 s | 47 MB |
| `fill` (Steps 1 + 2 + 2b) | 52 s | 103 MB |

These stages are fast in gapsmith because the LP solver (HiGHS) runs
in-process via `good_lp`, rather than shelling out through R's `cobrar`
wrapper to an external GLPK/CPLEX binary.

### Reproducing the benchmarks

```bash
# Requires: blastp, Rscript (with data.table, stringr, Biostrings),
# and gapsmith binary on PATH.
bash tools/bench/run_bench.sh <genomes_dir> <results_dir>
python3 tools/bench/aggregate.py <results_dir>
```

## Feature comparison

| Feature | gapseq (R) | gapsmith |
|---------|:----------:|:---------:|
| `find` (pathway detection) | ✅ | ✅ byte-identical output |
| `find-transport` | ✅ | ✅ row-count parity |
| `draft` (model assembly) | ✅ | ✅ 0 libSBML errors |
| `medium` (rule-based inference) | ✅ | ✅ byte-identical output |
| `fill` (4-phase gap-filling) | ✅ | ✅ Steps 1–4 |
| `adapt` (add/remove reactions) | ✅ | ✅ (EC/KEGG resolution deferred) |
| `pan` (pan-draft union) | ✅ | ✅ |
| `doall` (end-to-end pipeline) | ✅ | ✅ |
| `update-sequences` (Zenodo sync) | ✅ | ✅ |
| EC/TC conflict resolution | ✅ | ❌ (affects <1% of reactions) |
| MIRIAM cross-ref annotations | ✅ | ❌ (SBML loads without them) |
| HMM-based taxonomy prediction | ✅ | ❌ (use `--taxonomy` flag) |
| Precomputed alignment input | ❌ | ✅ `--aligner precomputed` |
| Batch-cluster alignment | ❌ | ✅ `batch-align` |
| In-process LP solver | ❌ (external GLPK) | ✅ (HiGHS, bundled) |
| FBA/pFBA subcommand | ❌ | ✅ `fba` |
| Native CBOR model format | ❌ (RDS) | ✅ (replaces RDS) |
| Single static binary | ❌ | ✅ |

## Key differences

### Solver

gapseq uses R's `cobrar` package (or the older `sybil`) which wraps
GLPK or CPLEX. gapsmith uses [HiGHS](https://highs.dev/) via
[`good_lp`](https://crates.io/crates/good_lp), statically linked at
build time. No runtime solver dependency.

### Model format

gapseq stores models as R `.RDS` files. gapsmith uses CBOR (`.gmod.cbor`)
as its native format — compact, fast to load, language-agnostic. Both
tools emit SBML for interchange.

### Alignment

Both tools shell out to the same external aligner binaries (blastp,
diamond, mmseqs2). gapsmith additionally supports `--aligner precomputed`
(skip the aligner, read a user-supplied TSV) and `batch-align` (cluster
N genomes with mmseqs2, align once, expand per-genome).

### Known deferred items

See [docs/porting-notes.md](docs/porting-notes.md) for the full list.
The main gaps:

- **EC/TC conflict resolution** — affects <1% of multi-EC-annotated genes.
- **MIRIAM cross-refs** — SBML emits ModelSEED id only; COBRApy round-trip still works.
- **HMM-based taxonomy** — CLI requires `--taxonomy Bacteria|Archaea` instead of auto-detecting.
- **`adapt` EC/KEGG resolution** — use direct SEED reaction ids or MetaCyc pathway ids.
