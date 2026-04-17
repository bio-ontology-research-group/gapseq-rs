# Architecture

A bird's-eye view of how gapsmith is organised, how data flows through
it, and what each crate is responsible for.

---

## 1. Crate dependency graph

```
             gapsmith-core            ← Types only, no I/O / DB / solver
              ↑  ↑  ↑
              │  │  │
    ┌─────────┘  │  └────────────────┐
    │            │                   │
  gapsmith-io   gapsmith-db           gapsmith-sbml
    │            │                   ↑
    └─┬──────────┘                   │
      │                              │
      ▼                              │
  gapsmith-align ──┐                   │
      │          │                   │
      ▼          │                   │
  gapsmith-find    │                   │
      │          │                   │
      ├─▶ gapsmith-transport           │
      │                              │
      ▼                              │
  gapsmith-draft ──────────────────────┤
      │                              │
      ▼                              │
  gapsmith-fill ◀─── gapsmith-medium    │
      │                              │
      └──────────────▶ gapsmith-cli ◀──┘
```

Arrows point from "depended on" toward "depends on". `gapsmith-core` is
the leaf: every other crate builds on its `Model` / `Reaction` /
`Metabolite` / `StoichMatrix` types.

---

## 2. Data flow for `gapsmith doall`

```text
  genome.faa.gz
        │
        ▼
  ┌──────────────┐                      ┌──────────────────────┐
  │ gapsmith find  │ ──── aligner (shell) │  dat/seq/<tax>/      │
  │ -p all       │          ─ BLAST/    │  rev/ unrev/ rxn/    │
  └──────┬───────┘            diamond/  └──────────────────────┘
         │                    mmseqs2
         │   *-Reactions.tbl  + *-Pathways.tbl
         ▼
  ┌──────────────────────┐
  │ gapsmith find-transport│ ─── same aligner, subex.tbl + tcdb.fasta
  └──────┬───────────────┘
         │   *-Transporter.tbl
         ▼
  ┌────────────────────┐
  │ gapsmith draft       │ ─── seed_reactions_corrected.tsv, biomass JSON
  └──────┬─────────────┘
         │   *-draft.gmod.cbor  (CBOR) + *-draft.xml  (SBML)
         ▼
  ┌──────────────────┐
  │ gapsmith medium    │ ─── medium_prediction_rules.tsv + *-Pathways.tbl
  │ (auto)           │
  └──────┬───────────┘
         │   *-medium.csv
         ▼
  ┌────────────────────┐
  │ gapsmith fill        │ ─── HiGHS LP (in-process via good_lp)
  │ 4-phase suite      │
  └──────┬─────────────┘
         │   *-filled.gmod.cbor  +  *-filled.xml  +  *-filled-added.tsv
         ▼
  (COBRApy, COBRAToolbox, cobrar, …)
```

Every arrow between stages is plain filesystem I/O. No daemon, no IPC.
Each subcommand is individually re-runnable; if you change the medium
you only re-run `fill`.

---

## 3. Model representation

[`gapsmith_core::Model`] is the single source of truth:

```rust
pub struct Model {
    pub annot: ModelAnnot,           // provenance metadata
    pub compartments: Vec<Compartment>,
    pub mets: Vec<Metabolite>,
    pub rxns: Vec<Reaction>,
    pub genes: Vec<GeneId>,
    pub s: StoichMatrix,              // sprs CsMat in CSC
}
```

Serialised via serde as:

- **CBOR** (`.gmod.cbor`) — native, compact, fast load. The gapsmith
  internal format; replaces upstream's `.RDS`.
- **JSON** (`.json`) — human-readable; `gapsmith convert` converts both
  ways. Same semantic content as CBOR.
- **SBML** (`.xml`) — Level 3 Version 1 + FBC 2 + groups 1. Written by
  [`gapsmith_sbml::write_sbml`]; loads cleanly in COBRApy / COBRAToolbox
  / cobrar.

Metabolite ids use the `cpd00001_c0` convention (compartment baked into
the id) for SBML SId compliance. Reaction ids likewise — `rxn00001_c0`,
`EX_cpd00001_e0`, `bio1`, etc.

---

## 4. LP plumbing

The gap-filler does a lot of LPs — dozens per Step 1 / 2 / 2b / 3 / 4
run, sometimes thousands in a `--full-suite`. Key design decisions:

### Split-flux encoding

Every reaction `r` becomes two variables `vp_r, vn_r ≥ 0`. The net
flux is `v_r = vp_r − vn_r`. This makes `|v_r| = vp_r + vn_r` linear —
essential for pFBA.

### Bound translation

| model `[lb, ub]`           | `vp_r` upper | `vn_r` upper |
|----------------------------|---|---|
| `lb ≥ 0, ub ≥ 0`           | `ub`           | `0`          |
| `lb ≤ 0, ub ≥ 0`           | `ub`           | `−lb`        |
| `lb ≤ 0, ub ≤ 0`           | `0`            | `−lb`        |

Strict-direction reactions additionally get `vp_r ≥ lb` (forward-forced)
or `vn_r ≥ −ub` (backward-forced).

### Mass balance

Walk the CSC matrix once column-by-column (`O(nnz)` total). For each
non-zero `(row, col, coef)`, accumulate `coef · (vp[col] − vn[col])`
onto row `row`'s expression. The naive per-row scan would be
`O(m · n)` — ~45 B iterations on a 3k×5k model. The `O(nnz)` walk
solves it in under 100 ms on ecoli.

### Solver

`good_lp` 1.15 with HiGHS backend. HiGHS is statically linked via
CMake at build time — no runtime HiGHS dependency. Optional CBC
fallback behind `--features cbc`: `pfba_heuristic` retries once with
CBC after HiGHS exhausts the tolerance / coefficient ladder.

---

## 5. Alignment layer

`gapsmith-align` exposes a single [`Aligner`] trait:

```rust
pub trait Aligner {
    fn run(&self, query: &Path, targets: &[&Path], opts: &AlignOpts)
        -> Result<Vec<Hit>, AlignError>;
}
```

Implementors — one struct each — shell out to the matching binary and
normalise its output to a [`Hit`] with the gapseq-canonical column
set (qseqid / pident / evalue / bitscore / qcovs / stitle / sstart /
send).

- [`BlastpAligner`] builds an indexed BLAST DB in a temp dir per call.
- [`DiamondAligner`] mirrors gapseq's `--more-sensitive` default.
- [`Mmseqs2Aligner`] uses the explicit `createdb → search →
  convertalis` 4-command pipeline (not `easy-search` — see
  [porting-notes](porting-notes.md) for why).
- [`PrecomputedTsvAligner`] skips all shelling out and reads a
  user-supplied TSV.
- [`BatchClusterAligner`] is the Rust-only feature: mmseqs2-clusters
  N input genomes, runs one alignment, expands cluster membership
  back to per-genome hits.
- [`GspaRunAligner`] is the multi-genome / metagenome entry point:
  given a [gspa](multi-genome.md#1-consuming-a-gspa-manifest---gspa-run)
  manifest (`clusters.tsv` + `alignment/*.tsv` + `genomes.tsv`) and a
  target genome id, it fans every rep-level hit onto that genome's
  cluster members. No alignment is re-run — gspa clustered across N
  genomes once, upstream.

---

## 6. Candidate-pool / gap-fill

The 4-phase gap-fill suite (`gapsmith-fill::run_suite`) is the heaviest
computation in the pipeline. Sketch of one fill call:

```
  draft model  +  medium CSV  +  Reactions.tbl (bitscore weights)
                         │
                         ▼
          ┌──────────────────────────────┐
          │ apply_medium(draft)          │
          │ attach EX_<target>_c0 obj    │
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ build_full_model             │  ← every approved SEED rxn
          │   draft + candidates         │    not already present,
          │   deduped by stoich hash     │    sorted by core status
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ (optional) detect_futile_    │  ← parallel LP probe,
          │   cycles  +  drop            │    `--prune-futile`
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ pfba_heuristic               │  ← tolerance ladder,
          │   min_growth ≥ k             │    `1e-6 → 1e-9`
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ extract utilized candidates  │
          │ add to draft                 │
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ KO essentiality loop         │  ← zero bounds, recheck FBA,
          │   serial, core-first         │    drop non-essential
          └──────────────┬───────────────┘
                         ▼
          ┌──────────────────────────────┐
          │ Step 2: biomass components   │  ← per substrate, MM_glu
          ├──────────────────────────────┤
          │ Step 2b: aerobic/anaerobic   │
          ├──────────────────────────────┤
          │ Step 3: ESP1-5 energy screen │  (`--full-suite`)
          ├──────────────────────────────┤
          │ Step 4: fermentation screen  │  (`--full-suite`)
          └──────────────────────────────┘
                         │
                         ▼
                    filled model
```

---

## 6b. Multi-genome and community paths

The single-genome flow above is unchanged. Three additional, opt-in
paths sit alongside it.

### 6b.1 Alignment reuse via gspa

```
   N proteomes ───▶  [gspa, upstream]  ──▶  gspa-run directory
                     (mmseqs2 linclust,      clusters.tsv
                      diamond vs. ref)       alignment/*.tsv
                                             genomes.tsv
                           │
     ┌─────────────────────┼─────────────────────┐
     ▼                     ▼                     ▼
  gapsmith find         find-transport         doall-batch
   --gspa-run …          --gspa-run …           --gspa-run …
     │                     │                     │
     └─ reads clusters.tsv + alignment/*.tsv, expands rep
        hits onto each target genome's cluster members
        (gapsmith_align::GspaRunAligner).
```

No alignment is run inside gapsmith when `--gspa-run` is set. Per-genome
`find` still applies its bitscore / coverage / completeness gates on the
fanned-out hits, so the downstream pipeline is unchanged.

### 6b.2 `doall-batch` parallel driver

```
   genomes.tsv  ──▶ sort + shard (`--shard i/N`) ──▶ rayon pool
                                                       │
                                                       ├─▶ doall genome_1
                                                       ├─▶ doall genome_2
                                                       ├─▶ ...
                                                       └─▶ doall genome_k
   out-dir/
     <genome_1>/
     <genome_2>/
     ...
     doall-batch-errors.tsv  ← present only with --continue-on-error
```

Each worker calls `doall::run_cli` in-process with its `--threads`
pinned to 1 so rayon owns the parallelism. SLURM array tasks get a
disjoint `--shard i/N` — no inter-task coordination, no shared state.

### 6b.3 Community modes

```
   N drafts (gmod.cbor) ──┐
                          │
              ┌───────────┴───────────┐
              ▼                       ▼
   community per-mag           community cfba
   (gapsmith-fill::             (gapsmith-fill::
    community::                  community::
     union_medium,                compose_models,
     per_mag_weights)             add_community_biomass)
              │                       │
              ▼                       ▼
   per-mag-growth.tsv +        community.gmod.cbor +
   community-medium.csv        community-fba.tsv
   (N independent FBAs)        (1 big LP across the
                                block-diagonal model)
```

`compose_models` suffixes every private metabolite and reaction with
`__<organism>`, leaves extracellular (`_e0`) ids **un**suffixed so they
collapse into one shared pool, and takes the most permissive
`EX_*_e0` bounds across organisms. `add_community_biomass` appends a
`community_biomass` pseudo-metabolite plus a draining `bio_community`
reaction (`Σ w_k · bio_k`), optionally adding one equality per organism
(`v(bio_k) == v(bio_community)`) for balanced-growth cFBA. The regular
[`fba::fba`] solves the resulting LP with no solver-side changes.

---

## 7. Testing surface

Every algorithmic component has unit tests (~170 total, `cargo test
--workspace`). On top, three integration layers:

1. **Shell-parity tests** (`crates/gapsmith-align/tests/*_parity.rs`)
   — run the real BLAST / diamond / mmseqs2 binaries and diff against
   our wrappers' output.
2. **R-parity tests** (`crates/gapsmith-find/tests/complex_parity.rs`,
   `pipeline_parity.rs`, `crates/gapsmith-transport/tests/parity.rs`)
   — run the actual R gapseq via `Rscript` and diff column-by-column.
3. **SBML validator** (`tools/validate_sbml.py`) — uv-managed venv
   with `python-libsbml` + `cobra`; runs libSBML consistency checks +
   COBRApy round-trip on emitted SBML. 0 errors on every model.

---

## 8. File-system layout

```
gapsmith/
  crates/
    gapsmith-core/               # types
    gapsmith-io/                 # CBOR/JSON + path resolver
    gapsmith-db/                 # dat/*.tsv loaders
    gapsmith-sbml/               # SBML writer
    gapsmith-align/              # aligner trait + 6 backends (incl. gspa-run)
    gapsmith-find/               # pathway + reaction detection
    gapsmith-transport/          # transporter detection
    gapsmith-draft/              # draft model builder
    gapsmith-medium/             # rule-based medium inference
    gapsmith-fill/               # FBA / pFBA / gap-filler / suite / community
    gapsmith-cli/                # clap dispatch + every subcommand
  tools/
    bench/                     # R-vs-rs benchmark runner
    validate_sbml.py           # libSBML + COBRApy validator
    r_complex_detection.R      # R parity harness
    .sbml-validate/            # uv venv (python-libsbml + cobra)
  docs/
    user-guide.md              # end-user walk-through
    cli-reference.md           # exhaustive flag list
    multi-genome.md            # metagenome, gspa integration, cFBA
    feature-matrix.md          # R source → Rust module pointers
    porting-notes.md           # intentional deviations from upstream
    architecture.md            # this file
  Cargo.toml                   # workspace manifest
  README.md                    # status + benchmarks
```
