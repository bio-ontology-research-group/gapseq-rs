//! Metabolite metadata.

use crate::{CompartmentId, CpdId};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Metabolite {
    pub id: CpdId,
    pub name: String,
    pub formula: Option<String>,
    #[serde(default)]
    pub charge: i32,
    pub compartment: CompartmentId,
    /// MetaNetX cross-reference (from `dat/mnxref_seed.tsv`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mnx: Option<String>,
}

impl Metabolite {
    pub fn new(id: impl Into<CpdId>, name: impl Into<String>, compartment: CompartmentId) -> Self {
        Self {
            id: id.into(),
            name: name.into(),
            formula: None,
            charge: 0,
            compartment,
            mnx: None,
        }
    }
}
