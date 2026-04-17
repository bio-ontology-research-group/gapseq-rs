//! Compartment identifiers.
//!
//! gapseq models use three compartments: cytosol (`c0`), extracellular (`e0`),
//! and periplasm (`p0`). See `src/construct_full_model.R:8`.

use serde::{Deserialize, Serialize};

/// Compact compartment index. Stored as `u8` so per-metabolite compartment
/// membership is a cheap vector of bytes.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct CompartmentId(pub u8);

impl CompartmentId {
    pub const CYTOSOL: CompartmentId = CompartmentId(0);
    pub const EXTRACELLULAR: CompartmentId = CompartmentId(1);
    pub const PERIPLASM: CompartmentId = CompartmentId(2);
}

/// Compartment metadata carried on the model.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Compartment {
    pub id: String,   // e.g. "c0"
    pub name: String, // e.g. "Cytosol"
}

impl Compartment {
    pub fn cytosol() -> Self {
        Self { id: "c0".into(), name: "Cytosol".into() }
    }
    pub fn extracellular() -> Self {
        Self { id: "e0".into(), name: "Extracellular".into() }
    }
    pub fn periplasm() -> Self {
        Self { id: "p0".into(), name: "Periplasm".into() }
    }
    pub fn default_three() -> Vec<Compartment> {
        vec![Self::cytosol(), Self::extracellular(), Self::periplasm()]
    }
}
