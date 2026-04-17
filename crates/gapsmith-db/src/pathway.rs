//! Pathway table loader — handles `meta_pwy.tbl`, `kegg_pwy.tbl`,
//! `seed_pwy.tbl`, `custom_pwy.tbl`.
//!
//! The four files are structurally identical apart from column count: some
//! have a trailing `spont` column, some don't. We parse positionally so the
//! same loader handles all variants.
//!
//! Columns (max schema, 14):
//! `id, name, altname, hierarchy, taxrange, reaId, reaEc, keyRea, reaName,
//!  reaNr, ecNr, superpathway, status, spont`.

use crate::common::{csv_err, io_err, DbError};
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum PwySource {
    MetaCyc,
    Kegg,
    Seed,
    Custom,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathwayRow {
    pub id: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub altname: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub hierarchy: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub taxrange: String,
    /// Comma-separated reaction ids (MetaCyc RXN ids, KEGG R-numbers, or SEED rxn ids).
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub rea_id: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub rea_ec: String,
    /// Comma-separated "key reaction" ids (presence ⇒ pathway considered present).
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub key_rea: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub rea_name: String,
    #[serde(default)]
    pub rea_nr: u32,
    #[serde(default)]
    pub ec_nr: u32,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub superpathway: String,
    /// Free-form status string (usually `TRUE` / `FALSE` / empty).
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub status: String,
    /// Comma-separated ids of reactions considered spontaneous. Empty when
    /// the source file lacks the `spont` column.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub spont: String,
    pub source: PwySource,
}

impl PathwayRow {
    pub fn rea_ids(&self) -> Vec<&str> {
        self.rea_id.split(',').map(str::trim).filter(|s| !s.is_empty()).collect()
    }
    pub fn ec_list(&self) -> Vec<&str> {
        self.rea_ec.split(',').map(str::trim).filter(|s| !s.is_empty()).collect()
    }
    pub fn key_rea_list(&self) -> Vec<&str> {
        self.key_rea.split(',').map(str::trim).filter(|s| !s.is_empty()).collect()
    }
    pub fn spont_list(&self) -> Vec<&str> {
        self.spont.split(',').map(str::trim).filter(|s| !s.is_empty()).collect()
    }
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct PathwayTable {
    pub source: Option<PwySource>,
    pub rows: Vec<PathwayRow>,
}

impl PathwayTable {
    pub fn load(path: impl AsRef<Path>, source: PwySource) -> Result<Self, DbError> {
        let path = path.as_ref();
        let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .quoting(false)
            .flexible(true)
            .from_reader(f);
        let headers = rdr.headers().map_err(|e| csv_err(path, e))?.clone();
        let col = |name: &str| headers.iter().position(|h| h.trim() == name);
        let c = Cols {
            id: col("id").unwrap_or(0),
            name: col("name").unwrap_or(1),
            altname: col("altname"),
            hierarchy: col("hierarchy"),
            taxrange: col("taxrange"),
            rea_id: col("reaId"),
            rea_ec: col("reaEc"),
            key_rea: col("keyRea"),
            rea_name: col("reaName"),
            rea_nr: col("reaNr"),
            ec_nr: col("ecNr"),
            superpathway: col("superpathway"),
            status: col("status"),
            spont: col("spont"),
        };
        let mut rows = Vec::new();
        for rec in rdr.records() {
            let rec = rec.map_err(|e| csv_err(path, e))?;
            rows.push(PathwayRow {
                id: rec.get(c.id).unwrap_or("").to_string(),
                name: rec.get(c.name).unwrap_or("").to_string(),
                altname: c.altname.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                hierarchy: c.hierarchy.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                taxrange: c.taxrange.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                rea_id: c.rea_id.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                rea_ec: c.rea_ec.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                key_rea: c.key_rea.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                rea_name: c.rea_name.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                rea_nr: c.rea_nr.and_then(|i| rec.get(i).and_then(|s| s.trim().parse().ok())).unwrap_or(0),
                ec_nr: c.ec_nr.and_then(|i| rec.get(i).and_then(|s| s.trim().parse().ok())).unwrap_or(0),
                superpathway: c.superpathway.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                status: c.status.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                spont: c.spont.and_then(|i| rec.get(i)).unwrap_or("").to_string(),
                source,
            });
        }
        tracing::info!(path = %path.display(), rows = rows.len(), ?source, "loaded pathway table");
        Ok(Self { source: Some(source), rows })
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
}

struct Cols {
    id: usize,
    name: usize,
    altname: Option<usize>,
    hierarchy: Option<usize>,
    taxrange: Option<usize>,
    rea_id: Option<usize>,
    rea_ec: Option<usize>,
    key_rea: Option<usize>,
    rea_name: Option<usize>,
    rea_nr: Option<usize>,
    ec_nr: Option<usize>,
    superpathway: Option<usize>,
    status: Option<usize>,
    spont: Option<usize>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_meta_pwy_schema() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("m.tsv");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(
            f,
            "id\tname\taltname\thierarchy\ttaxrange\treaId\treaEc\tkeyRea\treaName\treaNr\tecNr\tsuperpathway\tstatus\tspont"
        )
        .unwrap();
        writeln!(
            f,
            "PWY-1\tExample\talt\th\ttax\trxn1,rxn2\t1.1.1.1\trxn1\tex\t2\t1\tFALSE\tTRUE\trxn2"
        )
        .unwrap();
        let t = PathwayTable::load(&p, PwySource::MetaCyc).unwrap();
        assert_eq!(t.rows.len(), 1);
        let r = &t.rows[0];
        assert_eq!(r.id, "PWY-1");
        assert_eq!(r.rea_ids(), vec!["rxn1", "rxn2"]);
        assert_eq!(r.key_rea_list(), vec!["rxn1"]);
        assert_eq!(r.spont_list(), vec!["rxn2"]);
        assert_eq!(r.rea_nr, 2);
    }

    #[test]
    fn parses_kegg_pwy_without_spont() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("k.tsv");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(
            f,
            "id\tname\taltname\thierarchy\ttaxrange\treaId\treaEc\tkeyRea\treaName\treaNr\tecNr\tsuperpathway\tstatus"
        )
        .unwrap();
        writeln!(
            f,
            "map00010\tGlycolysis\t\tkegg;Metabolism\t\tR01061\t1.2.1.12\t\t\t1\t1\tFALSE\tTRUE"
        )
        .unwrap();
        let t = PathwayTable::load(&p, PwySource::Kegg).unwrap();
        assert_eq!(t.rows.len(), 1);
        assert!(t.rows[0].spont.is_empty());
    }
}
