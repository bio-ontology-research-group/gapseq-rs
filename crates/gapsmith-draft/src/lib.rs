//! Draft model construction.
//!
//! Port of `src/generate_GSdraft.R`. Produces a [`gapsmith_core::Model`]
//! from the `*-Reactions.tbl` and `*-Transporter.tbl` emitted by the
//! find / find-transport pipeline. The core algorithm:
//!
//! 1. Load the find + transport hit tables.
//! 2. Select every SEED reaction whose associated homology hit clears
//!    `high.evi.rxn.BS` (default 200) OR whose parent pathway is
//!    predicted present.
//! 3. Filter the reactions database to `gapseq.status ∈ {approved, corrected}`
//!    and intersect with the selected SEED reaction set.
//! 4. Dedupe via stoichiometric hash (`generate_rxn_stoich_hash.R`).
//! 5. Build the [`gapsmith_core::Model`]: add compartments, metabolites,
//!    reactions, biomass reaction, exchanges + diffusion.
//! 6. Emit CBOR (and optionally SBML) via `gapsmith_io` / `gapsmith_sbml`.
//!
//! Scope notes (simpler than the R original):
//!
//! - EC / TC conflict resolution (`prepare_candidate_reaction_tables.R`
//!   `resolve_common_EC_conflicts` / `resolve_common_TC_conflicts`) is
//!   deferred — covers <1% of reactions and requires IRanges-style
//!   overlap math. Document as a known gap.
//! - The conditional transporter inclusions (e.g. butyrate, IPA, PPA)
//!   are implemented as a simple rule table.
//! - The per-metabolite / per-reaction CVTerms / annotation block is
//!   limited to the ModelSEED identifier; full MIRIAM expansion is
//!   deferred.

pub mod biomass;
pub mod builder;
pub mod candidate;
pub mod exchanges;
pub mod gpr;
pub mod reactions_tbl;
pub mod runner;
pub mod stoich_hash;

pub use biomass::{parse_biomass_json, BiomassError, BiomassSpec};
pub use builder::{build_model, BuilderOptions};
pub use candidate::{CandidateReaction, CandidateTable};
pub use exchanges::{add_missing_diffusion, add_missing_exchanges};
pub use gpr::build_gpr_string;
pub use reactions_tbl::{
    read_reactions_tbl, read_transporter_tbl, ReactionRow, TransporterRow,
};
pub use runner::{run, DraftError, DraftOptions, DraftReport};
pub use stoich_hash::rxn_stoich_hash;
