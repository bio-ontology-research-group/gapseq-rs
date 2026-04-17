//! Reaction metadata.
//!
//! The stoichiometry itself lives in the model-level `StoichMatrix`; this
//! struct carries only per-reaction attributes that the solver and output
//! writers consume.

use crate::RxnId;
use serde::{Deserialize, Serialize};

/// SEED reactions carry a curation status used by gapseq to filter the
/// candidate pool. See `src/correct_seed_rxnDB.R`.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SeedStatus {
    Approved,
    Corrected,
    NotAssessed,
    Removed,
    /// Status not recorded (synthetic reactions: exchange, demand, biomass).
    #[default]
    None,
}

impl SeedStatus {
    pub fn is_usable(self) -> bool {
        matches!(self, SeedStatus::Approved | SeedStatus::Corrected)
    }
}

/// Reversibility marker.
///
/// Matches the `>` / `<` / `=` codes used throughout `dat/seed_reactions_corrected.tsv`.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Reversibility {
    Forward,
    Backward,
    Reversible,
}

impl Reversibility {
    pub fn from_code(c: char) -> Option<Self> {
        match c {
            '>' => Some(Self::Forward),
            '<' => Some(Self::Backward),
            '=' => Some(Self::Reversible),
            _ => None,
        }
    }

    pub fn code(self) -> char {
        match self {
            Self::Forward => '>',
            Self::Backward => '<',
            Self::Reversible => '=',
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Reaction {
    pub id: RxnId,
    pub name: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ec: Vec<String>,
    pub lb: f64,
    pub ub: f64,
    #[serde(default)]
    pub obj_coef: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gpr_raw: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subsystem: Option<String>,
    #[serde(default)]
    pub seed_status: SeedStatus,
    #[serde(default)]
    pub is_exchange: bool,
    #[serde(default)]
    pub is_biomass: bool,
    /// gapseq gs.origin code (see `generate_GSdraft.R`): 0 high-evidence SEED,
    /// 5 conditional transporter, 6 biomass, 7 exchange, 8 diffusion, 9 pathway-only.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gs_origin: Option<i8>,
    /// Bitscore from homology search (NaN / None if absent).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bitscore: Option<f32>,
    /// pFBA weight (populated by the gap-filler).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub weight: Option<f32>,
}

impl Reaction {
    pub fn new(id: impl Into<RxnId>, name: impl Into<String>, lb: f64, ub: f64) -> Self {
        Self {
            id: id.into(),
            name: name.into(),
            ec: Vec::new(),
            lb,
            ub,
            obj_coef: 0.0,
            gpr_raw: None,
            subsystem: None,
            seed_status: SeedStatus::None,
            is_exchange: false,
            is_biomass: false,
            gs_origin: None,
            bitscore: None,
            weight: None,
        }
    }

    pub fn reversibility(&self) -> Reversibility {
        match (self.lb < 0.0, self.ub > 0.0) {
            (true, true) => Reversibility::Reversible,
            (false, true) => Reversibility::Forward,
            (true, false) => Reversibility::Backward,
            // lb == 0 and ub == 0: treat as forward-blocked.
            (false, false) => Reversibility::Forward,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reversibility_roundtrip() {
        for c in ['>', '<', '='] {
            let r = Reversibility::from_code(c).unwrap();
            assert_eq!(r.code(), c);
        }
        assert!(Reversibility::from_code('?').is_none());
    }

    #[test]
    fn reversibility_from_bounds() {
        assert_eq!(Reaction::new("r", "r", -1000.0, 1000.0).reversibility(), Reversibility::Reversible);
        assert_eq!(Reaction::new("r", "r", 0.0, 1000.0).reversibility(), Reversibility::Forward);
        assert_eq!(Reaction::new("r", "r", -1000.0, 0.0).reversibility(), Reversibility::Backward);
    }
}
