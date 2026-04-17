//! `dat/subex.tbl` loader.
//!
//! Columns: `name, metacyc, vmh, seed, group`. Maps substrate names to
//! exchange-reaction IDs across three namespaces. Used by the transporter
//! pipeline (`src/transporter.sh`, `src/seed_transporter.R`).

use crate::common::{csv_err, io_err, DbError};
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubexRow {
    pub name: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub metacyc: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub vmh: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub seed: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub group: String,
}

pub fn load(path: impl AsRef<Path>) -> Result<Vec<SubexRow>, DbError> {
    let path = path.as_ref();
    let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .quoting(false)
        .flexible(true)
        .from_reader(f);
    let mut out = Vec::new();
    for rec in rdr.deserialize::<SubexRow>() {
        out.push(rec.map_err(|e| csv_err(path, e))?);
    }
    tracing::info!(path = %path.display(), rows = out.len(), "loaded subex");
    Ok(out)
}
