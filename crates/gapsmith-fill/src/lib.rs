//! LP and gap-filling for gapsmith.
//!
//! This crate hosts every linear-program the reconstruction pipeline
//! invokes, from the simplest FBA on a small toy network up to the full
//! 4-phase gap-fill driver over 10 000+-reaction candidate pools.
//!
//! # Pipeline at a glance
//!
//! 1. **FBA / pFBA primitives** ([`fba`], [`mod@pfba`]) — split-flux LP with
//!    `vp, vn ≥ 0` so `|v|` stays linear. Backend: `good_lp` + HiGHS.
//! 2. **pFBA-heuristic tolerance ladder** ([`pfba_heuristic`]) — retries
//!    up to 15 times, halving the feasibility tolerance then the pFBA
//!    coefficient on each failure. Port of `gapfill4.R:95–137`.
//! 3. **Medium application** ([`apply_medium`]) — closes every `EX_*`
//!    lower bound then reopens per-medium entries. Also adds missing
//!    exchange reactions for compounds not yet in the model.
//! 4. **Candidate pool** ([`build_full_model`]) — clones the draft and
//!    appends every SEED-approved reaction not already present, deduped
//!    by stoichiometric hash. This is the "full model" gap-filling
//!    searches against.
//! 5. **Single-iteration gap-fill** ([`gapfill4`]) — run pFBA-heuristic
//!    on the full model, extract candidate reactions that carried flux,
//!    add them to the draft, KO-essentiality-prune. Port of
//!    `gapfill4.R:1–303`.
//! 6. **4-phase suite** ([`run_suite`]) — orchestrates Steps 1 / 2 / 2b
//!    / 3 / 4 from `gf.suite.R`. Step 1 is the user-medium biomass fill,
//!    Steps 2 / 2b iterate per biomass substrate, Steps 3 / 4 do the
//!    energy-source and fermentation-product screens.
//! 7. **Futile-cycle prune** ([`detect_futile_cycles`]) — rayon-parallel
//!    pairwise LP probe that drops candidate reactions saturating at
//!    `0.99 · max_flux` with all boundaries closed. Opt-in because it's
//!    expensive on large pools (~40 min on ~8 k candidates).
//!
//! # Solver choice
//!
//! HiGHS is the default via `good_lp`'s `highs` feature; it's statically
//! linked via CMake at build time, so the binary has no runtime solver
//! dependency. Build with `--features cbc` to enable a CBC fallback
//! triggered when `pfba_heuristic`'s tolerance ladder exhausts itself.

pub mod error;
pub mod fba;
pub mod futile;
pub mod gapfill;
pub mod lp;
pub mod medium;
pub mod pfba;
pub mod pool;
pub mod suite;

pub use error::{FillError, SolveStatus};
pub use fba::{fba, FbaOptions, FbaSolution};
pub use gapfill::{drop_reactions, gapfill4, GapfillOptions, GapfillReport};
pub use lp::SplitFluxLp;
pub use medium::{apply_environment_file, apply_medium, read_medium, MediumEntry, MediumError};
pub use pfba::{pfba, pfba_heuristic, PfbaHeuristicOptions, PfbaOptions, PfbaSolution};
pub use pool::{
    apply_medium_to_full, build_full_model, pfba_weights, read_weights_from_reactions_tbl,
    rxn_weight, strip_compartment, HitRecord, RxnWeights, DEFAULT_BCORE, DEFAULT_DUMMY_WEIGHT,
    DEFAULT_HIGH_EVI,
};
pub use futile::{detect_futile_cycles, FutileOptions};
pub use suite::{run_suite, SuiteOptions, SuiteReport};
