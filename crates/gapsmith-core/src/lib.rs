//! Core types for gapsmith.
//!
//! Every other crate depends on this one. Kept deliberately small: only the
//! in-memory representation of a metabolic model plus the supporting ID,
//! stoichiometry, and GPR types. No I/O, no solvers, no databases.
//!
//! # Overview
//!
//! [`Model`] is the central type — an analogue of cobrar's `ModelOrg` S4
//! class. It bundles:
//!
//! - Compartments (cytosol `c0`, extracellular `e0`, periplasm `p0`).
//! - [`Metabolite`]s (compound id + compartment + chemical formula).
//! - [`Reaction`]s (id, name, bounds, objective coefficient, GPR, EC list,
//!   gapseq-specific `gs.origin` provenance tag, SEED curation status).
//! - A sparse [`StoichMatrix`] in CSC layout (`sprs` under the hood),
//!   one column per reaction, one row per metabolite.
//! - A [`ModelAnnot`] bundle of provenance metadata that travels with the
//!   model through CBOR / SBML round-trips.
//!
//! # ID newtypes
//!
//! [`CpdId`], [`RxnId`], [`GeneId`] wrap `Arc<str>`. Cloning is cheap and
//! the types are mutually incompatible so a metabolite id can't
//! accidentally be used where a reaction id is expected.
//!
//! # Serialisation
//!
//! Every public type derives `Serialize` + `Deserialize`. CBOR (via
//! `ciborium`) is the primary storage format; JSON works too. See the
//! `gapsmith-io` crate for the round-trip helpers.

pub mod compartment;
pub mod gpr;
pub mod id;
pub mod metabolite;
pub mod model;
pub mod reaction;
pub mod stoich;

pub use compartment::{Compartment, CompartmentId};
pub use gpr::{Gpr, GprParseError};
pub use id::{CpdId, GeneId, RxnId};
pub use metabolite::Metabolite;
pub use model::{Model, ModelAnnot, ModelError};
pub use reaction::{Reaction, Reversibility, SeedStatus};
pub use stoich::StoichMatrix;
