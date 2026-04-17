//! `dat/tcdb_substrates.tbl` loader.
//!
//! Format: two columns (no header), TAB-separated.
//!
//! ```text
//! 2.A.1.28.4<TAB>CHEBI:5651;ferroheme b
//! 1.A.11.4.1<TAB>CHEBI:7435;ammonium|CHEBI:7434;ammonia|...
//! ```
//!
//! The substrate column is a `|`-separated list of `CHEBI:<id>;<name>` pairs.

use crate::common::{csv_err, io_err, DbError};
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TcdbSubstrateRow {
    pub tc_id: String,
    pub substrates: Vec<TcdbSubstrate>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TcdbSubstrate {
    pub chebi: String,
    pub name: String,
}

pub fn load_substrates(path: impl AsRef<Path>) -> Result<Vec<TcdbSubstrateRow>, DbError> {
    let path = path.as_ref();
    let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .quoting(false)
        .flexible(true)
        .from_reader(f);

    let mut out = Vec::new();
    for rec in rdr.records() {
        let rec = rec.map_err(|e| csv_err(path, e))?;
        if rec.len() < 2 {
            continue;
        }
        let tc_id = rec.get(0).unwrap_or("").to_string();
        let raw = rec.get(1).unwrap_or("");
        let subs: Vec<TcdbSubstrate> = raw
            .split('|')
            .filter_map(|pair| {
                let (chebi, name) = pair.split_once(';')?;
                Some(TcdbSubstrate { chebi: chebi.to_string(), name: name.to_string() })
            })
            .collect();
        out.push(TcdbSubstrateRow { tc_id, substrates: subs });
    }
    tracing::info!(path = %path.display(), rows = out.len(), "loaded tcdb substrates");
    Ok(out)
}
