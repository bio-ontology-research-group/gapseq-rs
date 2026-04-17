# gapsmith porting notes

Intentional and unintentional deviations from the upstream R/bash gapseq
codebase. Each entry says *what* differs, *why*, and whether downstream
output parity holds.

This document is the authoritative list of "known non-R behavior". If you
find something in gapsmith that differs from real gapseq and it's not
listed here, please file an issue — we treat those as bugs.

---

## 1. Storage format

- **RDS → CBOR.** R's `saveRDS(mod, ...)` is replaced by `ciborium`-encoded
  CBOR at `<id>.gmod.cbor`. SBML is emitted alongside as the portable
  interchange format. An Rscript-driven CBOR↔RDS bridge can be added later
  if downstream R consumers need it.
- **Metabolite IDs: `cpd00001_c0` instead of `cpd00001[c0]`.** SBML SIds
  must match `[A-Za-z_][A-Za-z0-9_]*` — brackets break that. The SBML
  writer derives the species id from `met.id + met.compartment` and only
  presents `<cpd>[comp]` to the user via the `name` attribute.
- **Model IDs are sanitized** to valid SBML SIds: `-`, `.`, `:`, space →
  `_`. Original id preserved in the `name` attribute if different.

## 2. Alignment layer

- **mmseqs2 invocation.** gapseq's R pipeline calls explicit `mmseqs
  createdb` + `search` + `convertalis`; we replicate this 4-command flow
  exactly rather than the simpler `easy-search`, because the latter uses
  alignment-mode 3 (full-alignment identity) while gapseq calibrates
  against the k-mer prefilter identity that `search` reports.
- **mmseqs2 `qseqid` normalization.** We post-process mmseqs's `qheader`
  output to the first whitespace-delimited token, matching what
  blastp/diamond emit natively. gapseq does this via a sed call in the
  shell script (`gapseq_find.sh` at the tail of each mmseqs block).
- **Diamond `--more-sensitive` is on by default**, matching gapseq's
  `aliArgs="--more-sensitive"` default for diamond.
- **E-value formatting** in TSV output: BLAST's exact 2-decimal
  scientific notation (`2.52e-41`). Initial wrapper used `{:.3e}` which
  emitted `2.520e-41`; fixed to match.

## 3. find pipeline (M5)

- **`noSuperpathways=true` is always the default.** Matches gapseq's
  `gapseq_find.sh:20`. Pass `-n` / `--include-superpathways` to disable.
- **`-p` keyword shorthand.** `-p amino|nucl|cofactor|carbo|polyamine|fatty|energy|terpenoid|degradation|core|min|kegg|all`
  resolves to the same hierarchy-category name gapseq uses. For `-p <literal>`
  we fall back to `\b(<literal>)\b` regex against id / name / altname /
  hierarchy — same as gapseq's default `stype` branch.
- **Pathway table merge semantics.** `-l metacyc,custom` is default
  (matching gapseq); on ID collision, the `custom_pwy.tbl` row wins. This
  is how several "GS"-suffixed pathways (e.g. ARGSYN-PWY) have their
  `superpathway=TRUE` flag shadowed by the custom row with
  `superpathway=FALSE`, so they survive the filter.
- **`user/` reference fastas always resolve against the built-in
  `<gapseq>/dat/seq/<tax>/user/`**, even when `-D` / `--seq-dir`
  overrides the rest. Matches `prepare_batch_alignments.R:217`.
- **Reaction-name MD5 fallback.** Uses `md5sum`-compatible hex digest of
  the raw bytes (no trailing newline), same as `uniprot.sh:155`.
- **Header filter uses substring, not word-boundary, match.** This is a
  deliberate match of gapseq's shell behavior (`grep -Fivf` without
  `-w`). Has the known quirk of accidentally matching `4.A.2.1.1`
  against a line containing `4.A.2.1.11`. Reproducing the quirk is
  required for byte-identical output.

## 4. find-transport pipeline (M6)

- **Substrate case is preserved** from subex.tbl (`sub=Potassium`, not
  `sub=potassium`). Matches gapseq.
- **Alt-transporter reaction assignment** kicks in when a TC-id hit
  exists but no SEED reaction matches the exact TC type. Same rule as
  `analyse_alignments_transport.R:110-130`. Disable with `-a`.

## 5. draft pipeline (M7) — known gaps

- **EC/TC conflict resolution** (`prepare_candidate_reaction_tables.R`
  `resolve_common_EC_conflicts` / `resolve_common_TC_conflicts`, ~130
  lines of IRanges overlap math) is **not yet ported**. Affects <1% of
  reactions; the mitigations are listed downstream in the draft — a gene
  with two EC annotations stays attached to both reactions rather than
  being resolved.
- **MIRIAM cross-ref annotations.** SBML writer emits the ModelSEED id
  only. KEGG / BiGG / MetaNetX / HMDB / ChEBI cross-refs are available
  in `addReactAttr.R` + `addMetAttr.R` but not threaded through yet.
  Round-trip load in COBRApy works; MIRIAM-dependent downstream tools
  may need to re-query gapseq's `dat/mnxref_*.tsv`.
- **`gram_by_network.R`** (predict gram-staining by metabolic-network
  similarity against reference genomes) is not ported. The CLI requires
  `--biomass pos|neg|archaea` or reads the `gram=` field from the
  Reactions.tbl header (which find writes when gapseq has HMM-predicted
  the gram).
- **Biomass rescaling**: the menaquinone-8 auto-removal rule
  (`generate_GSdraft.R:281-292`) is **not ported**. Gram+/Gram- biomass
  retains menaquinone-8 regardless of whether the MENAQUINONESYN-PWY /
  PWY-5852 / PWY-5837 pathways are present. Impact is small (a single
  metabolite in the biomass reaction) but affects anaerobic-growth
  predictions.

## 6. Not yet implemented

- `gapsmith medium` — rule-based medium inference. Planned for M10.
- `gapsmith fill` — 4-phase gap-filling driver. Planned for M9 (FBA/pFBA
  primitives already shipped in M8 under `gapsmith-fill`).
- `gapsmith adapt` — add/remove reactions on an existing model. Planned for M10.
- `gapsmith pan` — pan-model union. Planned for M10.
- `gapsmith doall` — chain find→find-transport→draft→medium→fill. Planned for M10.
- `gapsmith update-sequences` — seqdb download from upstream. Planned for M10.
- HMM-based `predict_domain` / `predict_gramstaining`. Planned alongside M10.

## 6a. fill (M8–M9) — LP primitives + single-iteration gap-fill (shipped)

- Split-flux encoding: each reaction contributes `vp_r, vn_r ≥ 0` with net
  `v_r = vp_r − vn_r`. `|v_r| = vp_r + vn_r` is linear, which keeps pFBA a
  pure LP. The bound translation is documented in
  `crates/gapsmith-fill/src/lp.rs` — `[lb, ub]` on the net flux becomes a
  per-variable ub pair `(max(ub,0), max(−lb,0))`, with additional lower
  bounds when the direction is strictly forced.
- FBA (`gapsmith-fill::fba`) — maximise / minimise `Σ c_r · (vp_r − vn_r)` s.t.
  `S · (vp − vn) = 0`. HiGHS backend via `good_lp` 1.15.
- pFBA (`gapsmith-fill::pfba`) — `minimise Σ w_r·(vp_r+vn_r) − pfba_coef · Σ c_r·(vp_r−vn_r)`
  subject to `v_bio ≥ min_growth`. Matches cobrar's `pfbaHeuristic`.
- pFBA-heuristic ladder (`gapsmith-fill::pfba_heuristic`) — port of
  `gapfill4.R:95–137`. Retries up to 15 times, halving the feasibility
  tolerance then the pFBA coefficient on solver failure; validates each
  candidate solution by running FBA on the reduced model after zero-flux
  reactions are removed.
- `gapsmith fba` subcommand — FBA or `--pfba` on a CBOR/JSON model; prints
  solver status, objective, biomass flux, and the top-N |flux| reactions.
  Useful for sanity-checking a draft model before `fill` is applied.
- **Candidate pool** (`pool.rs::build_full_model`) — clone draft, append
  every approved/corrected SEED reaction not already present, deduped by
  stoichiometric hash with core-preference on ties. Port of
  `gapfill4.R:12–56`.
- **Medium application** (`medium.rs::apply_medium`) — close all `EX_*`
  lower bounds then re-open per medium CSV; adds missing EX reactions
  for compounds not yet in the model. Matches `constrain.model.R`.
- **Single-iteration gap-fill** (`gapfill.rs::gapfill4`) — pFBA-heuristic
  on the full model with bitscore-derived weights, extract utilized
  candidates, add to draft, run KO essentiality loop. Port of
  `gapfill4.R:1–303`.
- **`gapsmith fill`** subcommand — end-to-end `draft → filled` pipeline.
  Ecoli on `MM_glu.csv`: 20 reactions added, final growth 0.53.

### Known divergences / outstanding work

- **Biomass link-multi parsing** — `biomass.json` supports
  `"link": "cpd01997:-1|cpd03422:-1"` (pipe-separated coupled metabolites).
  The M7 biomass parser only read the *first* link. M9 fixed this
  (`BiomassComponent::links()`). Without the fix, every biomass reaction
  was missing two cofactor metabolites, pegging growth to zero. **Re-run
  `draft` on any model built before M9 — the old CBOR has a structurally
  broken biomass.**
- **SBML species-id double suffix** — the M7 `species_id` helper blindly
  appended `_<comp>` to the cpd id. When gapsmith stores metabolites as
  `cpd00001_c0` (the actual shipped convention), the SBML emitted
  `M_cpd00001_c0_c0`. COBRApy treated these as distinct mets, so SEED
  reactions and biomass didn't link. Fixed in M9: idempotent when the id
  already ends with the compartment suffix.
- **Objective must be `EX_<target>_c0` sink, not `bio1` directly** —
  `gf.suite.R:239–242` zeros every `obj_coef` and then calls
  `add_met_sink(mod, target.met, obj=1)` to attach the objective to a
  freshly-created sink reaction. Without it, steady-state balance pins
  the biomass flux to zero because the biomass pseudo-metabolite has no
  consumer. The CLI wires this up automatically.
- **4-phase suite** (Steps 2 / 2b / 3 / 4) shipped in the M9 follow-up.
  Default CLI runs Steps 1 + 2 + 2b; `--full-suite` adds 3 + 4. Steps 3
  and 4 iterate ~200 exchange compounds each with a gapfill4 call — budget
  10–30 minutes on a whole-bacterium genome.
- **Futile-cycle prune** shipped (`detect_futile_cycles`) but opt-in via
  `--prune-futile`. Pairwise LP probe on large candidate pools (~8k
  rxns) takes ~40 minutes even with rayon parallelism; use selectively.
- **Candidate-dedup fix** (draft-stage latent bug caught during M9
  follow-up): my build_candidates accumulator copied `is_complex` /
  `complex_status` / `pathway_status` only from the highest-bitscore row.
  That row may have empty status fields (gapseq fills those only on rows
  where complex detection ran), causing reactions like `rxn00011`
  (bs=1797, pyruvate decarboxylase) to be rejected in the runner's
  complex filter. Fixed: OR `is_complex`, max-across-rows for
  `complex_status`, highest-rank `pathway_status`. **Re-run `gapseq
  draft` on any model built before this fix.**
- **CBC fallback** wired behind `--features cbc`. Path:
  `pfba_heuristic` retries once with `good_lp::solvers::coin_cbc` after
  HiGHS exhausts its tolerance/coef ladder. Requires
  `coinor-libcbc-dev` + `coinor-libclp-dev` on the build host.

## 7. Unimplemented but probably won't be

- The R helpers that depend on `pythoncyc` (`meta2pwy.py`, `meta2genes.py`,
  `meta2rea.py`) are one-off MetaCyc DB updaters run ≤ twice per year by
  maintainers. They stay in Python; no Rust port planned.
- CPLEX solver support — plan stated "no CPLEX path"; gapsmith-fill will
  use HiGHS via `good_lp` with CBC/Clp as backups.

## 8. New features not in gapseq

- **`BatchClusterAligner`** (M4.5) — cluster N genomes' proteomes with
  mmseqs2, align once against the reference query FASTA, expand per-genome
  hit sets. Amortizes alignment cost over many genomes with shared
  proteins. Not in upstream gapseq.
- **Precomputed alignment input** — `--aligner precomputed -P <tsv>`
  on `find` / `find-transport` skips the alignment step entirely. User
  must pre-run blastp/diamond/mmseqs2 with the expected 8-column layout.

## 9. Output format differences

- **Pathways.tbl completeness precision.** Ported to `f64` to match R's
  default 15-significant-digit output (`66.6666666666667`).
- **Reactions.tbl column order** matches gapseq exactly after M5 fixes.
- **Transporter.tbl** matches gapseq exactly after M6 fixes, modulo the
  cosmetic `sub` column case when a substrate has both capitalized and
  lowercased SUBkey entries (gapseq picks the alphabetically-first
  post-dedup; we pick insertion-order-first).
- **`evalue` formatting** matches BLAST's `-outfmt 6` conventions: `0` for
  exact zero, 2-decimal scientific notation below 1e-3.

## 10. Testing surface

- R-vs-Rust parity tests run the actual R implementation via `Rscript`
  where feasible. This lives in three places:
  - `crates/gapsmith-find/tests/complex_parity.rs` drives
    `tools/r_complex_detection.R` which sources the real
    `src/complex_detection.R`.
  - `crates/gapsmith-find/tests/pipeline_parity.rs` runs real `gapsmith find`
    and gapsmith `find` side-by-side with a synthesized minimal seqdb.
  - `crates/gapsmith-transport/tests/parity.rs` does the same for
    `find-transport`.
- SBML validation runs `python-libsbml` + `cobra` via a uv-managed venv
  at `tools/.sbml-validate/`. Not part of `cargo test` (no libsbml on
  crates.io) but a manual gate via `tools/validate_sbml.py`.
