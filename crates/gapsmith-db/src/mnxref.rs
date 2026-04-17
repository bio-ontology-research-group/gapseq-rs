//! MNXref cross-reference loaders.
//!
//! gapseq consults two companion files:
//!
//! - `dat/mnxref_seed.tsv` — 5 cols: `MNX_ID`, `Balance`, `EC`, `Source` (SEED
//!   id), `kegg`.
//! - `dat/mnxref_seed-other.tsv` — 3 cols: `MNX_ID`, `seed`, `other`
//!   (`<prefix>:<id>`).
//!
//! The `dat/mnxref_reac_xref.tsv` raw MNXref dump is large (~240k lines) and
//! starts with a `###`-banner and `#`-prefixed license headers. We don't load
//! it by default — the seed-specific tables above already express the joins
//! gapseq actually uses.

use crate::common::{csv_err, io_err, DbError};
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MnxrefSeed {
    #[serde(rename = "MNX_ID")]
    pub mnx_id: String,
    #[serde(rename = "Balance", default, skip_serializing_if = "String::is_empty")]
    pub balance: String,
    #[serde(rename = "EC", default, skip_serializing_if = "String::is_empty")]
    pub ec: String,
    #[serde(rename = "Source", default, skip_serializing_if = "String::is_empty")]
    pub source: String,
    #[serde(rename = "kegg", default, skip_serializing_if = "String::is_empty")]
    pub kegg: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MnxrefSeedOther {
    #[serde(rename = "MNX_ID")]
    pub mnx_id: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub seed: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub other: String,
}

impl MnxrefSeedOther {
    /// Split `<prefix>:<id>` out of the `other` column. Returns `None` if
    /// the column is empty or has no `:`.
    pub fn split_other(&self) -> Option<(&str, &str)> {
        self.other.split_once(':')
    }
}

pub fn load_mnxref_seed(path: impl AsRef<Path>) -> Result<Vec<MnxrefSeed>, DbError> {
    let path = path.as_ref();
    let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .quoting(false)
        .flexible(true)
        .from_reader(f);
    let mut out = Vec::with_capacity(4_500);
    for rec in rdr.deserialize::<MnxrefSeed>() {
        match rec {
            Ok(r) => out.push(r),
            Err(e) => return Err(csv_err(path, e)),
        }
    }
    tracing::info!(path = %path.display(), rows = out.len(), "loaded mnxref_seed");
    Ok(out)
}

pub fn load_mnxref_seed_other(path: impl AsRef<Path>) -> Result<Vec<MnxrefSeedOther>, DbError> {
    let path = path.as_ref();
    let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .quoting(false)
        .flexible(true)
        .from_reader(f);
    let mut out = Vec::with_capacity(12_000);
    for rec in rdr.deserialize::<MnxrefSeedOther>() {
        match rec {
            Ok(r) => out.push(r),
            Err(e) => return Err(csv_err(path, e)),
        }
    }
    tracing::info!(path = %path.display(), rows = out.len(), "loaded mnxref_seed_other");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_mnxref_seed() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("s.tsv");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "MNX_ID\tBalance\tEC\tSource\tkegg").unwrap();
        writeln!(f, "MNXR94683\ttrue\t\trxn07919\t").unwrap();
        let r = load_mnxref_seed(&p).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].source, "rxn07919");
    }

    #[test]
    fn split_other_works() {
        let row = MnxrefSeedOther {
            mnx_id: "MNXR100003".into(),
            seed: "rxn11552".into(),
            other: "RXN-12078".into(), // no colon
        };
        assert!(row.split_other().is_none());
        let row = MnxrefSeedOther {
            mnx_id: "x".into(),
            seed: "y".into(),
            other: "bigg:R_EX_glc".into(),
        };
        assert_eq!(row.split_other(), Some(("bigg", "R_EX_glc")));
    }
}
