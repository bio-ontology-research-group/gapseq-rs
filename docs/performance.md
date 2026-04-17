# Performance

gapsmith shipped with a series of optimizations relative to its initial
Rust port. This page catalogues what's in, where wins actually show up,
and the methodology behind the numbers.

## Shipped optimizations

All items below are **on by default** unless marked *(opt-in)* — those
are library-level knobs exposed as struct fields, not CLI flags yet.

### C1 — SEED stoich-hash cache at load time
`crates/gapsmith-db/src/seed.rs`. `SeedRxnRow` gains a `#[serde(skip)]`
`stoich_hash: Option<String>` field populated once during
`load_seed_reactions` (sequential, ~150 ms for 35k rows). Gap-fill and
draft dedup paths read the cached hash instead of re-parsing +
re-sorting the stoichiometry on every access.

### C3 — O(1) reaction/metabolite lookups in `apply_medium`
`crates/gapsmith-fill/src/medium.rs`. Replaces `rxns.iter().find(..)`
per medium entry with a one-shot `HashMap<rxn_id, index>`. Same for
the metabolite-index lookup in the new-exchange fallback path.
`O(k·n)` → `O(n + k)` on a 500-rxn draft with 20 medium entries.

### A2 — FASTA header cache per `find` run
`crates/gapsmith-find/src/runner.rs`. Multiple reactions frequently
resolve to the same reference FASTA (two reactions sharing an EC both
land on `rev/<EC>.fasta`). A `HashMap<PathBuf, Arc<Vec<(String, String)>>>`
populated lazily during the unique-reaction loop avoids redundant
file reads + header parses.

### E1 — Structural blocked-reaction pre-pass in futile-cycle prune
`crates/gapsmith-fill/src/futile.rs::structural_blocked_rxns`.
Iterative topological dead-reaction detection on the closed-boundary
model: reactions whose metabolites have no producer or no consumer
(in any alive direction) can't carry flux, so they can't participate
in cycles. Drops 30–70% of the LP-probed candidate pool before any
LP runs. `O(iters × nnz)`, typically 3–8 iterations. Runs inside the
existing `--prune-futile` flag.

### B2 — Skip KO probes for zero-flux reactions
`crates/gapsmith-fill/src/gapfill.rs`. Any non-core candidate carrying
zero flux in the post-fill max-biomass FBA optimum is provably
non-essential: the constraint `x[r] = 0` leaves the current optimum
feasible, so the KO probe would succeed deterministically. Skip the
LP probe; drop the reaction directly. Cuts ~20–40% of KO-loop LP
solves on typical gap-fills.

### B1 *(opt-in)* — HiGHS hot-start via `with_initial_solution`
`crates/gapsmith-fill/src/fba.rs`. `FbaOptions` gains
`hot_start: Option<Vec<f64>>`. When `Some`, the caller supplies a
previous net-flux vector; `fba()` decomposes it into split-flux
`(vp, vn) = (max(v,0), max(-v,0))` seeds and threads them through
good_lp's `with_initial_solution`. HiGHS warm-starts the simplex
from the seed basis. Can nudge HiGHS toward a different alternative
optimum, so **off by default for byte-parity**.

### B3 *(opt-in)* — pFBA ladder memo across calls
`crates/gapsmith-fill/src/pfba.rs`. `PfbaHeuristicOptions` gains
`ladder_memo: Option<Arc<Mutex<Option<(f64, f64)>>>>`. When `Some`,
the ladder seeds from the last successful `(tolerance, coefficient)`
instead of the defaults and writes back on success. Clamped to
`[min_tol, start_tol]` to prevent a stale memo from pushing the
solver outside the user's configured range. Same byte-parity caveat
as B1. **Off by default.**

## Benchmark results

Hardware: server `ws`, 64 cores, 128 GB RAM, NVMe SSD, Debian.
Aligner: DIAMOND 2.1.9 (`-A diamond -K 8`).
Input: *E. coli* K-12 MG1655 (`GCF_000005845.2`, 4300 proteins).
Reference DB: upstream gapseq's `dat/` at commit-time HEAD.
Comparison: commit `a42d172` (pre-Phase-1 baseline) vs `474a1e7`
(all Phase 1 / 2 / 3 defaults shipped).

### `gapsmith find -p all`

| Metric | Phase-0 | Phase-3 | Δ |
|--------|--------:|--------:|--:|
| Wall | 18.36 s | 17.96 s | −2.2% |
| User | 115.8 s | 114.1 s | −1.5% |
| Peak RSS | 279 MB | 286 MB | +2.5% |

Pathways.tbl byte-identical. Reactions.tbl varies run-to-run on the
same binary — DIAMOND parallel search with `-K > 1` has inherent
non-determinism at the per-hit level. Not a regression; not in our
code.

### `gapsmith doall` (full default pipeline)

| Metric | Phase-0 | Phase-3 | Δ |
|--------|--------:|--------:|--:|
| Wall | 77.6 s | 79.2 s | +2.1% |
| User | 244 s | 245 s | +0.4% |
| Peak RSS | 276 MB | 297 MB | +7.7% |

**Semantic parity (via COBRApy):**
- 2776 reactions, 2193 metabolites (both)
- Biomass FBA growth: **0.7422696154309283** (identical to 16 decimals)
- Reaction Jaccard: **1.0000**
- Metabolite Jaccard: **1.0000**

### `gapsmith fill --full-suite` (4-phase, deep KO loops)

| Metric | Phase-0 | Phase-3 | Δ |
|--------|--------:|--------:|--:|
| Wall | 11:51.36 | 11:51.55 | +0.03% |
| User | 1209 s | 1217 s | +0.7% |
| Peak RSS | 107 MB | 123 MB | +15.5% |

**Semantic parity:**
- 2859 reactions, 2207 metabolites (both)
- Biomass FBA growth: **0.7433163537454** (identical to 13 decimals)
- Reaction Jaccard: **1.0000**
- Metabolite Jaccard: **1.0000**

Added reactions differ only in **step attribution** (one reaction
`rxn12933_c0` added in Step 1 by P0 vs Step 2 by P3 — B2's zero-flux
skip dropped it from Step 1's add-list because its flux at the
biomass optimum was zero; Step 2 needed it and re-added it). Final
model is identical.

### `gapsmith fill --full-suite --prune-futile` (E1 shines)

Both binaries hit a pre-existing `pfba_heuristic: max iterations`
failure on this workload (unrelated to Phase 1–3 changes). Timing
*to the failure point*:

| Metric | Phase-0 | Phase-3 | Δ |
|--------|--------:|--------:|--:|
| Wall | 5:56 | 4:20 | **−27%** |
| User | 18963 s | 13728 s | −28% |
| Peak RSS | 1354 MB | 1322 MB | −2.4% |

E1 cuts the time needed to exhaust the LP-prune pool substantially.
When `--prune-futile` actually completes on a well-conditioned input,
this speedup translates directly.

## Honest assessment

The measured wins on the **default path** are small:
- DIAMOND dominates `find` / `doall` wall time (~60 s of 80 s total).
- `fill` (Step 1 default) runs few LPs; B2 saves fractions of a second.
- C1 / C3 savings are milliseconds, below noise.
- A2 cache adds ~15 MB RSS because per-FASTA reads were already buffered.

Real wins show up in:
- **`--prune-futile` enabled** (E1: ~25–30% wall reduction)
- **Deep KO loops** in full-suite runs with large candidate pools (B2)
- **Batch runs** reusing `ladder_memo` across pFBA calls (B3, opt-in)
- **KO-loop-dominated fills** with `hot_start` configured (B1, opt-in)

Semantic parity is the firm guarantee: identical reaction sets,
identical biomass to machine precision, regardless of which
optimization path fires.

## Scaling to many genomes

The numbers above are all single-genome. For 100+ genomes, the
dominant cost moves out of FBA and into per-genome alignment.
Gapsmith's multi-genome path routes around that entirely:

- **`--gspa-run`** consumes a pre-clustered, pre-aligned gspa manifest
  so a 1000-genome collection pays for one cross-genome alignment,
  not 1000. Aligner calls inside `find` / `find-transport` become
  map-lookups.
- **`doall-batch --jobs N`** runs `doall` across genomes with rayon;
  on a 32-core box, speedup is ~`min(N, 32)` vs. looping `doall`
  sequentially.
- **`doall-batch --shard i/N`** makes SLURM array jobs the unit of
  parallelism. At 1 024 tasks of 16 cores each, throughput is ~1–2 s
  of wallclock per genome including I/O — limited by per-task
  doall-internal work, not by coordination.

Community modes are small LPs (≤ 1000 organisms in `per-mag`, which is
linear) or a single big LP (`cfba`). `cfba` is feasible up to ~50
organisms before HiGHS becomes the bottleneck; for larger communities
prefer `per-mag`.

See [Multi-genome workflows](multi-genome.md) for recipes.

## Deferred — where the big wins still are

### B1-full — HiGHS basis warm-start across KO probes (`highs-sys` FFI)
The currently-shipped B1 uses good_lp's `with_initial_solution` for
a simplex hot-start guess, but each solve still constructs a fresh
`highs::Model`. True warm-start needs a persistent `highs::Model`
across the KO loop, mutating bounds via `Highs_changeColBounds`
and re-solving in place. The `highs` crate 2.0 doesn't expose that
API — requires ~200 LOC of `highs-sys` FFI. Estimated 30–50% wall
reduction on KO-dominated fills.

### Parallel KO probe
The current KO loop is serial because each commit affects the next
iteration's baseline. A valid parallelisation: probe every
candidate concurrently against the *same* post-fill model (all
probes are independent at that snapshot), collect a set of
provably-non-essential candidates, then re-probe the survivors
serially. Bounded to a single pass, correctness-preserving.
Estimated 4–8× on multi-core for the Step 4 fermentation-product
loop.

### Step-4 restructuring
Step 4 probes ~100 fermentation products one at a time, each doing
a full gap-fill cycle. Amortising the shared preprocessing across
probes (same starting model, same weights) would cut ~5 min off a
full-suite run.

## Reproducing

```bash
# 1. Clone + build both commits
git clone https://github.com/bio-ontology-research-group/gapsmith.git
cd gapsmith
git worktree add -f /tmp/gapsmith-p0 a42d172
cd /tmp/gapsmith-p0 && cargo build --release --bin gapsmith
cd - && cargo build --release --bin gapsmith

# 2. Get reference data
gapsmith update-data -o /path/to/dat
gapsmith --data-dir /path/to/dat update-sequences -t Bacteria

# 3. Run both binaries, compare
for bin in /tmp/gapsmith-p0/target/release/gapsmith ./target/release/gapsmith; do
  out=$(mktemp -d)
  /usr/bin/time -v "$bin" --data-dir /path/to/dat doall \
    -A diamond -K 8 -f "$out" genome.faa.gz 2>&1 | tee "$out.time"
done

# 4. Semantic diff in COBRApy
uv run --with cobra python tools/bench/diff_models.py \
  out1/genome-filled.xml out2/genome-filled.xml
```

`tools/bench/run_bench.sh` bundles the above end-to-end.
