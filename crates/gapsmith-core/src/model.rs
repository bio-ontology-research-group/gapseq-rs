//! In-memory metabolic model.
//!
//! Analogue of cobrar's `ModelOrg` S4 class. Replaces RDS serialization with
//! serde (CBOR/JSON — see `gapsmith-io`).

use crate::{Compartment, CpdId, GeneId, Metabolite, Reaction, RxnId, StoichMatrix};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Provenance and version metadata that travels with the model.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ModelAnnot {
    /// Genome / organism identifier (usually the FASTA basename).
    pub id: String,
    /// Human-readable name.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// gapseq version that produced the model.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gapsmith_version: Option<String>,
    /// Reference sequence database version.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seqdb_version: Option<String>,
    /// Taxonomic domain (`Bacteria`, `Archaea`, ...).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tax_domain: Option<String>,
    /// Gram staining (`pos`, `neg`, `na`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gram: Option<String>,
    /// Free-form notes (merged into SBML `<notes>` on export).
    #[serde(default)]
    pub notes: Vec<String>,
}

#[derive(Debug, thiserror::Error)]
pub enum ModelError {
    #[error("unknown reaction id `{0}`")]
    UnknownRxn(String),
    #[error("unknown metabolite id `{0}`")]
    UnknownMet(String),
    #[error("duplicate reaction id `{0}`")]
    DuplicateRxn(String),
    #[error("duplicate metabolite id `{0}`")]
    DuplicateMet(String),
    #[error(
        "stoichiometric matrix shape {got:?} does not match model dimensions {expected:?}"
    )]
    ShapeMismatch { got: (usize, usize), expected: (usize, usize) },
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Model {
    pub annot: ModelAnnot,
    pub compartments: Vec<Compartment>,
    pub mets: Vec<Metabolite>,
    pub rxns: Vec<Reaction>,
    /// Genes referenced by any reaction's GPR.
    #[serde(default)]
    pub genes: Vec<GeneId>,
    /// Sparse stoichiometric matrix in CSC (one column per reaction).
    pub s: StoichMatrix,
}

impl Model {
    /// Empty model with the three standard compartments (`c0`, `e0`, `p0`).
    pub fn new(id: impl Into<String>) -> Self {
        Self {
            annot: ModelAnnot { id: id.into(), ..Default::default() },
            compartments: Compartment::default_three(),
            mets: Vec::new(),
            rxns: Vec::new(),
            genes: Vec::new(),
            s: StoichMatrix::zeros(0, 0),
        }
    }

    pub fn rxn_count(&self) -> usize {
        self.rxns.len()
    }
    pub fn met_count(&self) -> usize {
        self.mets.len()
    }

    /// Lookup map from reaction id → index. Built on demand (O(n)); callers
    /// that need many lookups should cache the result.
    pub fn rxn_index(&self) -> HashMap<RxnId, usize> {
        self.rxns
            .iter()
            .enumerate()
            .map(|(i, r)| (r.id.clone(), i))
            .collect()
    }

    pub fn met_index(&self) -> HashMap<CpdId, usize> {
        self.mets
            .iter()
            .enumerate()
            .map(|(i, m)| (m.id.clone(), i))
            .collect()
    }

    /// Sanity-check that the stoichiometric matrix shape matches the model.
    pub fn check_shape(&self) -> Result<(), ModelError> {
        let expected = (self.mets.len(), self.rxns.len());
        let got = (self.s.rows(), self.s.cols());
        if got == expected {
            Ok(())
        } else {
            Err(ModelError::ShapeMismatch { got, expected })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CompartmentId, Reversibility};

    fn toy_model() -> Model {
        // 3 mets × 2 rxns: A -> B, B -> C.
        let mut m = Model::new("toy");
        m.mets.push(Metabolite::new("cpdA", "A", CompartmentId::CYTOSOL));
        m.mets.push(Metabolite::new("cpdB", "B", CompartmentId::CYTOSOL));
        m.mets.push(Metabolite::new("cpdC", "C", CompartmentId::CYTOSOL));

        let mut r1 = Reaction::new("r1", "A -> B", 0.0, 1000.0);
        r1.obj_coef = 0.0;
        m.rxns.push(r1);
        m.rxns.push(Reaction::new("r2", "B -> C", 0.0, 1000.0));

        m.s = StoichMatrix::from_triplets(
            3,
            2,
            vec![(0, 0, -1.0), (1, 0, 1.0), (1, 1, -1.0), (2, 1, 1.0)],
        );
        m
    }

    #[test]
    fn shape_check() {
        let m = toy_model();
        m.check_shape().unwrap();
    }

    #[test]
    fn bad_shape_detected() {
        let mut m = toy_model();
        m.rxns.push(Reaction::new("r3", "C -> sink", 0.0, 1000.0));
        assert!(m.check_shape().is_err());
    }

    #[test]
    fn reaction_indexing() {
        let m = toy_model();
        let idx = m.rxn_index();
        assert_eq!(idx[&RxnId::new("r1")], 0);
        assert_eq!(idx[&RxnId::new("r2")], 1);
    }

    #[test]
    fn reversibility_from_bounds() {
        let m = toy_model();
        assert_eq!(m.rxns[0].reversibility(), Reversibility::Forward);
    }
}
