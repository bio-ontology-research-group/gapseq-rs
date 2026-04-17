//! `dat/complex_subunit_dict.tsv` loader.
//!
//! Columns: `rxn, subunit_synonym, subunit`.
//!
//! Used during alignment post-processing to map subunit synonyms (case-insensitive,
//! substring search) to canonical subunit names, so that a complex can be
//! considered "present" when enough of its subunits have good blast hits
//! (`src/complex_detection.R:10–37`).

use crate::common::{csv_err, io_err, DbError};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplexSubunitEntry {
    pub rxn: String,
    pub subunit_synonym: String,
    pub subunit: String,
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct ComplexSubunitTable {
    pub rows: Vec<ComplexSubunitEntry>,
    /// Inverse lookup: `rxn → [(synonym, canonical)]` built at load time.
    pub by_rxn: HashMap<String, Vec<(String, String)>>,
}

impl ComplexSubunitTable {
    pub fn load(path: impl AsRef<Path>) -> Result<Self, DbError> {
        let path = path.as_ref();
        let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .quoting(false)
            .flexible(true)
            .from_reader(f);
        let mut rows = Vec::new();
        for rec in rdr.deserialize::<ComplexSubunitEntry>() {
            rows.push(rec.map_err(|e| csv_err(path, e))?);
        }
        let mut by_rxn: HashMap<String, Vec<(String, String)>> = HashMap::new();
        for r in &rows {
            by_rxn
                .entry(r.rxn.clone())
                .or_default()
                .push((r.subunit_synonym.clone(), r.subunit.clone()));
        }
        tracing::info!(path = %path.display(), rows = rows.len(), "loaded complex subunit dict");
        Ok(Self { rows, by_rxn })
    }

    pub fn for_rxn(&self, rxn: &str) -> Option<&[(String, String)]> {
        self.by_rxn.get(rxn).map(|v| v.as_slice())
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
}
