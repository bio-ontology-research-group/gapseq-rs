# Feature matrix

Exhaustive list of everything gapsmith implements, one row per feature,
with pointers to both the upstream R source and the Rust module.

Status legend:

- ✅ — implemented and tested against real gapseq where feasible
- 🆕 — Rust-only feature (no upstream equivalent)
- ⚠️ — shipped but intentionally deviates from upstream (see [porting-notes.md](porting-notes.md))
- ❌ — deferred / intentionally not ported

---

## 1. Subcommands

| Subcommand | R source | Rust module | Status |
|---|---|---|---|
| `gapsmith test` | — | `gapsmith-cli/src/commands/test.rs` | ✅ |
| `gapsmith find` | `src/gapseq_find.sh` + `src/*.R` | `gapsmith-find/` + `gapsmith-cli/src/commands/find.rs` | ✅ byte-identical on `PWY-6587` & `amino` |
| `gapsmith find-transport` | `src/transporter.sh` + `src/analyse_alignments_transport.R` | `gapsmith-transport/` + `gapsmith-cli/src/commands/find_transport.rs` | ✅ TC-set + row-count identical |
| `gapsmith draft` | `src/generate_GSdraft.R` + `src/prepare_candidate_reaction_tables.R` | `gapsmith-draft/` + `gapsmith-cli/src/commands/draft.rs` | ✅ SBML validates (0 libSBML errors) |
| `gapsmith medium` | `src/predict_medium.R` | `gapsmith-medium/` + `gapsmith-cli/src/commands/medium.rs` | ✅ byte-identical on ecoli |
| `gapsmith fill` | `src/gf.suite.R` + `src/gapfill4.R` | `gapsmith-fill/` + `gapsmith-cli/src/commands/fill.rs` | ✅ 4-phase suite + KO loop |
| `gapsmith adapt` | `src/adapt.R` + `src/gf.adapt.R` | `gapsmith-cli/src/commands/adapt.rs` | ⚠️ EC / KEGG / name resolution deferred |
| `gapsmith pan` | `src/pan-draft.R` | `gapsmith-cli/src/commands/pan.rs` | ✅ union + binary table |
| `gapsmith doall` | `src/doall.sh` | `gapsmith-cli/src/commands/doall.rs` | ✅ end-to-end on ecore in 2m47s |
| `gapsmith update-sequences` | `src/update_sequences.sh` | `gapsmith-cli/src/commands/update_sequences.rs` | ✅ Zenodo sync + md5 diff |
| `gapsmith convert` | — | `gapsmith-cli/src/commands/convert.rs` | 🆕 CBOR ↔ JSON round-trip |
| `gapsmith example-model` | — | `gapsmith-cli/src/commands/example_model.rs` | 🆕 toy model fixture |
| `gapsmith db inspect` | — | `gapsmith-cli/src/commands/db.rs` | 🆕 reference-data row-count dump |
| `gapsmith export-sbml` | `cobrar::writeSBMLmod` | `gapsmith-cli/src/commands/export_sbml.rs` | 🆕 CBOR → SBML |
| `gapsmith align` | — | `gapsmith-cli/src/commands/align.rs` | 🆕 debug-wrap for a single aligner |
| `gapsmith batch-align` | — | `gapsmith-cli/src/commands/batch_align.rs` | 🆕 cluster N genomes + single alignment |
| `gapsmith doall-batch` | — | `gapsmith-cli/src/commands/doall_batch.rs` | 🆕 rayon + SLURM-shard parallel `doall` across N genomes |
| `gapsmith community per-mag` | — | `gapsmith-cli/src/commands/community.rs` | 🆕 per-MAG FBA under a shared (union) medium |
| `gapsmith community cfba` | — | `gapsmith-cli/src/commands/community.rs` | 🆕 compose N drafts, weighted-sum biomass objective |
| `gapsmith fba` | — | `gapsmith-cli/src/commands/fba.rs` | 🆕 FBA / pFBA standalone |

---

## 2. Core algorithms

### 2.1 Alignment layer (`gapsmith-align`)

| Feature | R source | Rust module |
|---|---|---|
| BLASTp wrapper | `gapseq_find.sh` blastp block | `blast.rs::BlastpAligner` |
| tBLASTn wrapper | same, for `-n nucl` | `blast.rs::TblastnAligner` |
| DIAMOND wrapper | `gapseq_find.sh` diamond block | `diamond.rs::DiamondAligner` |
| mmseqs2 wrapper (full pipeline) | `gapseq_find.sh` mmseqs block | `mmseqs2.rs::Mmseqs2Aligner` |
| Precomputed TSV input | — | `precomputed.rs::PrecomputedTsvAligner` 🆕 |
| Batch-cluster (N genomes → 1 alignment) | — | `batch.rs::BatchClusterAligner` 🆕 |
| gspa-run manifest reader (cluster-aware hit expansion) | — | `gspa.rs::{GspaManifest, GspaRunAligner}` 🆕 |
| 2-decimal scientific e-value format | BLAST `-outfmt 6` native | `tsv.rs` |

### 2.2 `find` pipeline (`gapsmith-find`)

| Feature | R source | Rust module |
|---|---|---|
| Pathway table loader (meta / kegg / seed / custom) | `gapseq_find.sh:520-532` | `gapsmith-db::PathwayTable` |
| metacyc + custom merge (custom-wins-on-id) | same | same |
| Keyword-shorthand resolution (`amino`, `carbo`, ...) | `gapseq_find.sh:40-60` | `pathways.rs::MatchMode::Hierarchy` |
| Reference FASTA resolver (`user/` → `rxn/` → `rev/EC` → `unrev/EC` → md5) | `prepare_batch_alignments.R:150-234` | `seqfile.rs` |
| Complex-subunit detection | `complex_detection.R` | `complex.rs` (R-parity on 9 cases) |
| Hit classification with exception table | `analyse_alignments.R:108-189` | `classify.rs` |
| Pathway completeness scoring (`f64` precision) | `filter_pathways.R:10-34` | `pathways.rs::score` |
| `dbhit` lookup (EC + altEC + MetaCyc id + enzyme name) | `getDBhit.R:60-130` | `dbhit.rs` |
| `noSuperpathways=true` default | `gapseq_find.sh:20` | `find::FindOptions` |
| Word-boundary-less header filter (matches shell `grep -Fivf`) | `gapseq_find.sh` | `seqfile.rs` ⚠️ intentional |
| Output writers (`Reactions.tbl`, `Pathways.tbl`) | same | `output.rs` |

### 2.3 `find-transport` pipeline (`gapsmith-transport`)

| Feature | R source | Rust module |
|---|---|---|
| `subex.tbl` substrate filter | `transporter.sh:140-280` | `filter.rs` |
| TC-id parsing + type canonicalisation | `analyse_alignments_transport.R:1-188` | `tc.rs` |
| Substrate resolution (tcdb_all + FASTA header fallback) | same | `runner.rs` |
| Alt-transporter reaction assignment (gated by `--nouse-alternatives`) | `analyse_alignments_transport.R:110-130` | `runner.rs` |
| Substrate-case preservation (gapseq emits `sub=Potassium`) | shell behaviour | `data.rs` |

### 2.4 Draft model builder (`gapsmith-draft`)

| Feature | R source | Rust module |
|---|---|---|
| Candidate selection (bitscore ≥ cutoff OR pathway support) | `prepare_candidate_reaction_tables.R` + `generate_GSdraft.R:55-100` | `candidate.rs` |
| Stoichiometric hash dedup | `generate_rxn_stoich_hash.R` | `stoich_hash.rs` |
| Best-status-across-rows (OR `is_complex`, max `complex_status`, highest-rank `pathway_status`) | implicit in R's data.table merges | `candidate.rs::build_candidates` ⚠️ explicit |
| Biomass JSON parser (single + pipe-separated multi-link) | `parse_BMjson.R:1-107` | `biomass.rs` + `gapsmith-db::BiomassComponent::links` |
| Biomass cofactor mass-rescaling | `generate_GSdraft.R:281-292` | `biomass.rs` ⚠️ menaquinone-8 auto-removal deferred |
| GPR composition (and / or tree, "subunit undefined" edge cases) | `get_gene_logic_string.R` | `gpr.rs` |
| Diffusion + exchange expansion | `add_missing_exRxns.R:1-156` | `exchanges.rs` |
| Conditional transporter additions (butyrate, IPA, PPA, phloretate) | `generate_GSdraft.R` | `runner.rs::add_conditional_transporters` |
| SBML ID sanitiser (`-/./:/space → _`) | — | `builder.rs` 🆕 |
| Cytosolic met-id format (`cpd00001_c0` not `cpd00001[c0]`) | — | `builder.rs` ⚠️ for SBML SId compliance |

### 2.5 FBA / pFBA solver (`gapsmith-fill`)

| Feature | R source | Rust module |
|---|---|---|
| Split-flux LP encoding (`vp, vn ≥ 0`) | implicit in `cobrar::pfbaHeuristic` | `lp.rs::SplitFluxLp` |
| FBA | `cobrar::fba` | `fba.rs::fba` |
| pFBA (single call) | `cobrar::pfbaHeuristic` | `pfba.rs::pfba` |
| pFBA-heuristic tolerance ladder (15 iters, `1e-6 → 1e-9`, pFBA-coef relaxation) | `gapfill4.R:95-137` | `pfba.rs::pfba_heuristic` |
| HiGHS solver | — (R uses glpk/cplex) | `good_lp` 1.15 + `highs-sys` |
| CBC fallback | — | `pfba.rs::pfba_cbc` (feature-gated `cbc`) 🆕 |
| Row-expression builder (`O(nnz)`) | implicit | `fba.rs::build_row_exprs` 🆕 performance |

### 2.6 Gap-filling (`gapsmith-fill`)

| Feature | R source | Rust module |
|---|---|---|
| `gapfill4` single-iteration driver | `gapfill4.R:1-303` | `gapfill.rs::gapfill4` |
| Candidate pool (`draft + all approved SEED, stoich-hash deduped`) | `construct_full_model.R` + `gapfill4.R:12-56` | `pool.rs::build_full_model` |
| `rxnWeights` derivation from bitscores | `prepare_candidate_reaction_tables.R:222-228` | `pool.rs::rxn_weight` |
| KO essentiality loop (serial, core-first, highest-weight-first) | `gapfill4.R:247-280` | `gapfill.rs::gapfill4` |
| Medium application (close all EX, open per-medium, add missing EX) | `constrain.model.R` | `medium.rs::apply_medium` |
| Environment overrides (`env_highH2.tsv`) | `adjust_model_env.R` | `medium.rs::apply_environment_file` |
| Step 1 (user medium + biomass target) | `gf.suite.R:244-258` | `suite.rs::run_suite` |
| Step 2 (per-biomass-component on MM_glu + carbon sources) | `gf.suite.R:285-372` | `suite.rs::step2` |
| Step 2b (aerobic / anaerobic variant) | `gf.suite.R:377-464` | `suite.rs::run_suite` |
| Step 3 (energy-source screen with ESP1-5) | `gf.suite.R:480-581` | `suite.rs::step3` |
| Step 4 (fermentation-product screen) | `gf.suite.R:585-683` | `suite.rs::step4` |
| Target-met sink as objective | `add_met_sink` in `add_missing_exRxns.R:56-72` | `suite.rs::add_target_sink_obj` |
| Futile-cycle detector (parallel pairwise LP probe) | recent upstream `cccbb6f0` | `futile.rs::detect_futile_cycles` (opt-in `--prune-futile`) |
| Community model composition (block-diagonal, shared `_e0`) | — | `community.rs::compose_models` 🆕 |
| Weighted-sum community biomass + optional balanced-growth | — | `community.rs::add_community_biomass` 🆕 |
| Union-medium + per-MAG weights (community `per-mag` mode) | — | `community.rs::{union_medium, per_mag_weights, weighted_growth}` 🆕 |

### 2.7 Medium inference (`gapsmith-medium`)

| Feature | R source | Rust module |
|---|---|---|
| Rules-table loader | `predict_medium.R:46` | `rules.rs::load_rules` |
| Boolean-expression evaluator (`\| & ! < > == <= >=`) | `eval(parse(text=))` | `boolexpr.rs::eval` |
| Counting-rule support (`a + b + c < 3`) | same (R int arithmetic) | `boolexpr.rs::parse_sum` |
| Cross-rule dedup + mean flux | `predict_medium.R:84-86` | `predict.rs::predict_medium` |
| Saccharides / Organic acids category dedup | `predict_medium.R:88-92` | `predict.rs::predict_medium` |
| Manual flux overrides | `predict_medium.R:94-114` | `predict.rs::parse_manual_flux` |
| Proton balancer | `predict_medium.R:121-132` | `predict.rs::predict_medium` |

### 2.8 Serialisation (`gapsmith-sbml`, `gapsmith-io`)

| Feature | R source | Rust module |
|---|---|---|
| CBOR round-trip | — | `gapsmith-io::{read,write}_model_cbor` 🆕 |
| JSON round-trip | — | `gapsmith-io::{read,write}_model_json` 🆕 |
| SBML L3V1 + FBC2 + groups writer | `cobrar::writeSBMLmod` | `gapsmith-sbml::write_sbml` |
| SBML SId idempotent on mets with compartment suffix | — | `writer.rs::species_id` 🆕 bugfix |
| Streaming via `quick-xml` | — | `writer.rs` 🆕 no libSBML dep |
| SBML consistency validation | libSBML native | `tools/validate_sbml.py` (libSBML + COBRApy) |

### 2.9 Reference-data loaders (`gapsmith-db`)

| Feature | R source | Rust module |
|---|---|---|
| `seed_reactions_corrected.tsv` | `data.table::fread` | `seed.rs::load_seed_reactions` |
| `seed_metabolites_edited.tsv` | same | `seed.rs::load_seed_metabolites` |
| MNXref cross-refs (`mnxref_*.tsv`) | same | `mnxref.rs` |
| `meta_pwy.tbl` / `kegg_pwy.tbl` / `seed_pwy.tbl` / `custom_pwy.tbl` | same | `pathway.rs` |
| `subex.tbl` | same | `subex.rs` |
| `tcdb.tsv` | same | `tcdb.rs` |
| `exception.tbl` | same | `exception.rs` |
| `medium_prediction_rules.tsv` | same | `gapsmith-medium::rules` |
| `complex_subunit_dict.tsv` | same | `complex.rs` |
| Biomass JSON (Gram+, Gram-, archaea, user custom) | `parse_BMjson.R` | `biomass.rs` |
| SEED stoichiometry parser (`-1:cpd00001:0:0:"H2O";...`) | `parse_BMjson.R:21-29` | `stoich_parse.rs` |

---

## 3. New Rust-only features

| Feature | Rust module | Motivation |
|---|---|---|
| Precomputed alignment input (`--aligner precomputed -P <tsv>`) | `gapsmith-align::PrecomputedTsvAligner` | Skip per-genome BLAST when the user pre-runs diamond / mmseqs2 at batch scale |
| BatchClusterAligner (`gapsmith batch-align`) | `gapsmith-align::BatchClusterAligner` | Amortise alignment cost over N genomes via one mmseqs2 cluster + single alignment |
| gspa-run manifest reader (`--gspa-run <dir>`) | `gapsmith-align::GspaRunAligner` | Consume precomputed cluster-rep hits from the upstream [gspa](https://github.com/bio-ontology-research-group/gspa) pipeline; fans rep hits onto per-genome members |
| `gapsmith doall-batch` | `gapsmith-cli::commands::doall_batch` | Rayon + SLURM-array-friendly driver for reconstructing 100 → 1 M genomes in one batch |
| `gapsmith community per-mag` | `gapsmith-fill::community` + CLI | Shared-medium per-MAG FBA for metagenomes with 50+ MAGs |
| `gapsmith community cfba` | `gapsmith-fill::community` + CLI | Full community LP (block-diagonal compose, weighted-sum biomass, optional balanced-growth) |
| In-process LP (HiGHS via good_lp) | `gapsmith-fill` | Replaces R cobrar's shelled-out glpk / cplex; faster warm-starts |
| Optional CBC fallback backend | `gapsmith-fill::pfba_cbc` (`--features cbc`) | When HiGHS exhausts the tolerance ladder on pathological LPs |
| CBOR native format | `gapsmith-io` | Fast, compact, stdlib-free; replaces R's RDS |
| `gapsmith fba` subcommand | `gapsmith-cli` | Standalone FBA / pFBA without shelling into R |
| `gapsmith convert` subcommand | `gapsmith-cli` | CBOR ↔ JSON round-trip for inspection |
| `gapsmith db inspect` subcommand | `gapsmith-cli` | Smoke-test the reference data directory |
| `gapsmith export-sbml` subcommand | `gapsmith-cli` | Write an arbitrary CBOR model as SBML |

---

## 4. Known gaps (deferred)

| Gap | Upstream location | Workaround / plan |
|---|---|---|
| EC / TC conflict resolution (IRanges overlap math) | `prepare_candidate_reaction_tables.R::resolve_common_{EC,TC}_conflicts` | Affects <1 % of multi-EC annotations. Plan: port when a user case needs it. |
| MIRIAM cross-ref annotations (KEGG / BiGG / MetaNetX / HMDB / ChEBI) | `addReactAttr.R` + `addMetAttr.R` | SBML emits ModelSEED id only; round-trip in COBRApy still works. |
| HMM-based taxonomy / gram prediction | `predict_domain.R`, `predict_gramstaining.R` | CLI requires explicit `--taxonomy Bacteria|Archaea`; gram defaults to `neg`. |
| Gene-name MD5 fallback in seqfile resolver | `uniprot.sh:179` | Common-case MD5 fallback is ported; gene-name branch rarely fires. |
| Menaquinone-8 auto-removal (gated on MENAQUINONESYN-PWY / PWY-5852 / PWY-5837) | `generate_GSdraft.R:281-292` | Bio1 includes cpd15500 regardless; affects anaerobic predictions marginally. |
| `gram_by_network.R` (predict gram by metabolic-network similarity) | same | Requires explicit `-b pos|neg|archaea` instead. |
| `adapt` EC / KEGG / enzyme-name resolution | `adapt.R::ids2seed` strategies 3–7 | Direct SEED + pathway id resolution works; user can pre-resolve via `gapsmith find`. |
| `pan` weight medianing (`custom_median`) | `pan-draft_functions.R` | Pan-model emits without merged `rxnWeights` metadata; `gapsmith fill` on a pan-draft needs the source Reactions.tbl. |
| CPLEX solver support | — | Plan explicitly calls for HiGHS + optional CBC; no CPLEX path. |
| MetaCyc DB updaters (`meta2pwy.py`, `meta2genes.py`, `meta2rea.py`) | upstream Python helpers | Run once per year by maintainers; kept in Python. |

---

## 5. Testing surface

| Test suite | Tests | What it asserts |
|---|---|---|
| `gapsmith-core` unit | 18 | Type invariants, serde round-trips, stoichiometric matrix construction. |
| `gapsmith-io` unit | 5 | CBOR / JSON round-trip, data-dir auto-detect. |
| `gapsmith-db` unit | 18 | Every reference-data parser on realistic inputs. |
| `gapsmith-sbml` unit + integration | 2 + 1 | SBML writer emits every FBC2 / groups element; libSBML validates cleanly. |
| `gapsmith-align` unit + smoke + parity | 18 + 4 + 3 | Aligner trait, precomputed TSV, gspa-run manifest + fan-out, BLAST / diamond / mmseqs2 shell parity. |
| `gapsmith-find` unit + smoke + parity | 36 + 1 + 2 | Pathway scoring, complex detection (R-parity on 9 cases), `find -p PWY-6587` and `-p amino` byte-identical against real gapseq. |
| `gapsmith-transport` unit + parity | 7 + 1 | TC parsing, substrate resolution, end-to-end row+TC-id parity against real gapseq. |
| `gapsmith-draft` unit + smoke | 10 + 1 | Biomass rescaling, GPR composition, stoich dedup, conditional transporters. |
| `gapsmith-fill` unit + textbook + smoke | 20 + 5 + 1 | FBA / pFBA / pFBA-heuristic on toys; community compose + cFBA on 2-organism toy; `gapfill4` end-to-end on ecoli draft. |
| `gapsmith-medium` unit | 14 | Boolean-expression evaluator (incl. counting rules), rule loader, cross-rule dedup, proton balance. |
| `gapsmith-cli` integration | 6 + 1 | CBOR↔JSON round-trip end-to-end via the binary; SLURM-shard parser. |
| **Total** | **~170** | |

Run the full suite:

```bash
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
```

---

## 6. LOC breakdown

| Crate | `src/` LOC | `tests/` LOC |
|---|---:|---:|
| `gapsmith-core` | ~1 000 | — |
| `gapsmith-io` | ~330 | — |
| `gapsmith-db` | ~1 600 | — |
| `gapsmith-sbml` | ~870 | ~220 |
| `gapsmith-align` | ~1 250 | ~560 |
| `gapsmith-find` | ~2 700 | ~380 |
| `gapsmith-transport` | ~1 040 | ~115 |
| `gapsmith-draft` | ~1 670 | ~80 |
| `gapsmith-fill` | ~2 150 | ~400 |
| `gapsmith-medium` | ~550 | — |
| `gapsmith-cli` | ~2 400 | ~60 |
| **Total** | **~17 000** | **~2 100** |
