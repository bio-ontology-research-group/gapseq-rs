# gapsmith

gapsmith is a Rust reimplementation of
[gapseq](https://github.com/jotech/gapseq). It reconstructs
genome-scale metabolic models from bacterial proteomes: pathway
detection, transporter detection, draft model assembly, growth-medium
inference, and iterative gap-filling — all as a single static binary.

## What it does

Given a protein FASTA, gapsmith:

1. **Finds pathways** — aligns proteins against per-reaction reference
   FASTAs and scores pathway completeness (`gapsmith find`).
2. **Finds transporters** — matches against a TCDB-derived substrate
   table (`gapsmith find-transport`).
3. **Assembles a draft model** — combines the above into a
   stoichiometrically consistent metabolic network (`gapsmith draft`).
4. **Infers a growth medium** — applies boolean rules over detected
   reactions to predict minimal-medium compounds (`gapsmith medium`).
5. **Gap-fills** — an iterative pFBA + knock-out loop adds the minimal
   set of extra reactions needed for biomass growth (`gapsmith fill`).

Each step is independently callable, or chain them end-to-end with
`gapsmith doall`.

## One-liner

```bash
gapsmith doall genome.faa.gz -f output/
```

This produces a gap-filled SBML model that loads directly in COBRApy,
COBRAToolbox, or any SBML-compatible tool.

## Why a rewrite?

Upstream gapseq is R + bash (~9k LOC). gapsmith replaces the R stack —
including the `cobrar` wrapper and external GLPK/CPLEX — with:

- In-process [HiGHS](https://highs.dev/) via `good_lp` for FBA/pFBA
- Native CBOR model format (replaces R's RDS)
- Single static binary, no R install
- `--aligner precomputed` mode for re-using one alignment across many
  genomes
- `batch-align` for clustering N genomes with mmseqs2 and amortising
  the per-reaction alignment

Benchmarks show **~3× faster** pathway detection on real bacterial
proteomes and 35–40% lower peak memory. See
[Comparison with upstream gapseq](comparison.md) for details.

## Scaling out: metagenomes and multi-genome runs

gapsmith also ships **optional, additive** features for community and
multi-genome workloads:

- **`--gspa-run`** — consume protein clusters + alignments produced
  upstream by [gspa](https://github.com/bio-ontology-research-group/gspa)
  so thousands of genomes share one alignment.
- **`doall-batch`** — rayon-parallel driver with a `--shard i/N`
  selector for SLURM array jobs (scales from one machine to 1 M genomes).
- **`community per-mag`** — per-MAG FBA under a shared (union) medium.
  Linear in N MAGs.
- **`community cfba`** — compose N drafts into one community model with
  a shared `_e0` pool and a weighted-sum biomass objective; optional
  balanced-growth constraint for classical cFBA semantics.

See [Multi-genome & metagenome workflows](multi-genome.md) for the full
recipe.

## Where to start

- **Install & first run**: [User guide](user-guide.md)
- **Every CLI flag**: [CLI reference](cli-reference.md)
- **Metagenomes / multi-genome runs**: [Multi-genome workflows](multi-genome.md)
- **How the crates fit together**: [Architecture](architecture.md)
- **R → Rust module mapping**: [Feature matrix](feature-matrix.md)
- **Intentional deviations from upstream**: [Porting notes](porting-notes.md)
