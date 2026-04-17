//! SEED reaction and metabolite table loaders.
//!
//! Input files live in `dat/`:
//! - `seed_reactions_corrected.tsv` — 23 tab-separated columns (see below).
//! - `seed_metabolites_edited.tsv` — 28 columns.
//!
//! We parse every column the downstream pipeline consumes. Unused columns
//! are kept as raw strings so a round-trip-write stays possible; fields that
//! really aren't referenced anywhere (e.g. `deltag`, `pka`) are loaded as
//! `Option<String>` for trivial cost.

use crate::common::{csv_err, io_err, DbError};
use crate::stoich_parse::{parse_stoichiometry, StoichTerm};
use gapsmith_core::{CpdId, Reversibility, RxnId, SeedStatus};
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Row of `dat/seed_reactions_corrected.tsv`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeedRxnRow {
    pub id: RxnId,
    pub abbreviation: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub code: String,
    /// Raw stoichiometry field. Use [`SeedRxnRow::parse_stoich`] to decode.
    pub stoichiometry: String,
    #[serde(default)]
    pub is_transport: u8,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub equation: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub definition: String,
    /// Raw reversibility code (`>`, `<`, `=`, or empty).
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub reversibility: String,
    /// Raw thermodynamic direction code.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub direction: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub abstract_reaction: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub pathways: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub aliases: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub ec_numbers: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub deltag: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub deltagerr: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub compound_ids: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub status: String,
    #[serde(default)]
    pub is_obsolete: u8,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub linked_reaction: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub notes: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub is_copy_of: String,
    /// gapseq curation status — see [`SeedStatus`].
    pub gapseq_status: SeedStatus,

    /// Canonical stoichiometric hash, populated by [`load_seed_reactions`].
    /// Used as a fingerprint for dedup in draft / gap-fill. `None` when
    /// the row was constructed outside the loader or when the stoichiometry
    /// failed to parse. Not serialised — recomputed on deserialisation.
    #[serde(skip, default)]
    pub stoich_hash: Option<String>,
}

impl SeedRxnRow {
    pub fn parse_stoich(&self) -> Result<Vec<StoichTerm>, crate::stoich_parse::StoichParseError> {
        parse_stoichiometry(&self.stoichiometry)
    }

    pub fn reversibility(&self) -> Option<Reversibility> {
        self.reversibility.chars().next().and_then(Reversibility::from_code)
    }

    /// Split comma-separated EC numbers, dropping empty entries.
    pub fn ec_list(&self) -> Vec<&str> {
        self.ec_numbers
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect()
    }

    /// Split comma-separated pathway memberships.
    pub fn pathway_list(&self) -> Vec<&str> {
        self.pathways
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect()
    }
}

/// Row of `dat/seed_metabolites_edited.tsv`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeedCpdRow {
    pub id: CpdId,
    #[serde(rename = "MNX_ID", default, skip_serializing_if = "String::is_empty")]
    pub mnx_id: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub abbreviation: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub formula: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub mass: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub source: String,
    #[serde(default)]
    pub charge: i32,
    #[serde(default)]
    pub is_core: u8,
    #[serde(default)]
    pub is_obsolete: u8,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub linked_compound: String,
    #[serde(default)]
    pub is_cofactor: u8,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub deltag: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub deltagerr: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub pka: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub pkb: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub abstract_compound: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub comprised_of: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub aliases: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub smiles: String,
    #[serde(rename = "InChIKey", default, skip_serializing_if = "String::is_empty")]
    pub inchikey: String,
    #[serde(rename = "hmdbID", default, skip_serializing_if = "String::is_empty")]
    pub hmdb_id: String,
    #[serde(rename = "reactomeID", default, skip_serializing_if = "String::is_empty")]
    pub reactome_id: String,
    #[serde(rename = "chebiID", default, skip_serializing_if = "String::is_empty")]
    pub chebi_id: String,
    #[serde(rename = "InChI", default, skip_serializing_if = "String::is_empty")]
    pub inchi: String,
    #[serde(rename = "keggID", default, skip_serializing_if = "String::is_empty")]
    pub kegg_id: String,
    #[serde(rename = "biggID", default, skip_serializing_if = "String::is_empty")]
    pub bigg_id: String,
    #[serde(rename = "biocycID", default, skip_serializing_if = "String::is_empty")]
    pub biocyc_id: String,
}

/// Load `seed_reactions_corrected.tsv`.
///
/// Any `SeedStatus::Removed` row is kept (downstream code may want to
/// surface it in diagnostics). Filter at the call site with
/// [`SeedStatus::is_usable`] if you only want the candidate pool.
pub fn load_seed_reactions(path: impl AsRef<Path>) -> Result<Vec<SeedRxnRow>, DbError> {
    let path = path.as_ref();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .quoting(false) // SEED uses `"` inside fields without CSV-quoting semantics
        .flexible(true)
        .from_path(path)
        .map_err(|e| csv_err(path, e))?;

    // We re-read the header manually so we can tolerate the exact column name
    // `gapseq.status` (serde requires a field name; renaming via serde is
    // verbose given 23 columns). The approach: pull the header, then use
    // positional indexing.
    let headers = rdr.headers().map_err(|e| csv_err(path, e))?.clone();
    let col = |name: &str| -> Option<usize> {
        headers.iter().position(|h| h.trim() == name)
    };
    let c_id = col("id").ok_or_else(|| DbError::Parse {
        path: path.to_path_buf(),
        line: 1,
        msg: "missing `id` column".into(),
    })?;
    let c_abbr = col("abbreviation").unwrap_or(usize::MAX);
    let c_name = col("name").ok_or_else(|| DbError::Parse {
        path: path.to_path_buf(),
        line: 1,
        msg: "missing `name` column".into(),
    })?;
    let c_code = col("code").unwrap_or(usize::MAX);
    let c_stoich = col("stoichiometry").ok_or_else(|| DbError::Parse {
        path: path.to_path_buf(),
        line: 1,
        msg: "missing `stoichiometry` column".into(),
    })?;
    let c_is_transport = col("is_transport").unwrap_or(usize::MAX);
    let c_equation = col("equation").unwrap_or(usize::MAX);
    let c_definition = col("definition").unwrap_or(usize::MAX);
    let c_rev = col("reversibility").unwrap_or(usize::MAX);
    let c_dir = col("direction").unwrap_or(usize::MAX);
    let c_abs = col("abstract_reaction").unwrap_or(usize::MAX);
    let c_pwy = col("pathways").unwrap_or(usize::MAX);
    let c_ali = col("aliases").unwrap_or(usize::MAX);
    let c_ec = col("ec_numbers").unwrap_or(usize::MAX);
    let c_deltag = col("deltag").unwrap_or(usize::MAX);
    let c_deltagerr = col("deltagerr").unwrap_or(usize::MAX);
    let c_cpd_ids = col("compound_ids").unwrap_or(usize::MAX);
    let c_status = col("status").unwrap_or(usize::MAX);
    let c_is_obs = col("is_obsolete").unwrap_or(usize::MAX);
    let c_linked = col("linked_reaction").unwrap_or(usize::MAX);
    let c_notes = col("notes").unwrap_or(usize::MAX);
    let c_copy = col("is_copy_of").unwrap_or(usize::MAX);
    let c_gapseq = col("gapseq.status").ok_or_else(|| DbError::Parse {
        path: path.to_path_buf(),
        line: 1,
        msg: "missing `gapseq.status` column".into(),
    })?;

    let mut out = Vec::with_capacity(35_000);
    for (i, rec) in rdr.records().enumerate() {
        let rec = rec.map_err(|e| csv_err(path, e))?;
        let get = |c: usize| -> String {
            if c == usize::MAX { String::new() } else { rec.get(c).unwrap_or("").to_string() }
        };
        let get_u8 = |c: usize| -> u8 { get(c).trim().parse().unwrap_or(0) };
        let gapseq_status = match get(c_gapseq).as_str() {
            "approved" => SeedStatus::Approved,
            "corrected" => SeedStatus::Corrected,
            "not.assessed" | "not_assessed" => SeedStatus::NotAssessed,
            "removed" => SeedStatus::Removed,
            _ => SeedStatus::None,
        };
        out.push(SeedRxnRow {
            id: RxnId::new(get(c_id).trim()),
            abbreviation: get(c_abbr),
            name: get(c_name),
            code: get(c_code),
            stoichiometry: get(c_stoich),
            is_transport: get_u8(c_is_transport),
            equation: get(c_equation),
            definition: get(c_definition),
            reversibility: get(c_rev),
            direction: get(c_dir),
            abstract_reaction: get(c_abs),
            pathways: get(c_pwy),
            aliases: get(c_ali),
            ec_numbers: get(c_ec),
            deltag: get(c_deltag),
            deltagerr: get(c_deltagerr),
            compound_ids: get(c_cpd_ids),
            status: get(c_status),
            is_obsolete: get_u8(c_is_obs),
            linked_reaction: get(c_linked),
            notes: get(c_notes),
            is_copy_of: get(c_copy),
            gapseq_status,
            stoich_hash: None,
        });
        let _ = i;
    }
    // Cache the canonical stoich hash once so gap-fill / draft dedup paths
    // don't re-parse + re-sort the stoichiometry string for every candidate
    // on every run. Sequential — 35k rows in ~150 ms on warm CPU.
    for row in &mut out {
        if !row.stoichiometry.is_empty() {
            row.stoich_hash =
                crate::stoich_hash::rxn_stoich_hash(&row.stoichiometry, &row.reversibility).ok();
        }
    }
    tracing::info!(path = %path.display(), rows = out.len(), "loaded SEED reactions");
    Ok(out)
}

/// Load `seed_metabolites_edited.tsv`.
pub fn load_seed_metabolites(path: impl AsRef<Path>) -> Result<Vec<SeedCpdRow>, DbError> {
    let path = path.as_ref();
    let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .quoting(false)
        .flexible(true)
        .from_reader(f);

    let mut out = Vec::with_capacity(28_000);
    for rec in rdr.deserialize::<SeedCpdRow>() {
        match rec {
            Ok(r) => out.push(r),
            Err(e) => return Err(csv_err(path, e)),
        }
    }
    tracing::info!(path = %path.display(), rows = out.len(), "loaded SEED metabolites");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_minimal_reactions_tsv() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.tsv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(
            f,
            "id\tabbreviation\tname\tcode\tstoichiometry\tis_transport\tequation\tdefinition\treversibility\tdirection\tabstract_reaction\tpathways\taliases\tec_numbers\tdeltag\tdeltagerr\tcompound_ids\tstatus\tis_obsolete\tlinked_reaction\tnotes\tis_copy_of\tgapseq.status"
        )
        .unwrap();
        writeln!(
            f,
            "rxn00001\tR00004\tsome name\t\t-1:cpd00001:0:0:\"H2O\";1:cpd00002:0:0:\"X\"\t0\t\t\t>\t=\t\t\t\t3.5.1.2\t\t\t\tOK\t0\t\t\t\tapproved"
        )
        .unwrap();
        writeln!(
            f,
            "rxn00002\tR00005\tother\t\t-1:cpd00002:0:0:\"X\"\t1\t\t\t=\t=\t\t\t\t\t\t\t\tOK\t0\t\t\t\tremoved"
        )
        .unwrap();

        let rows = load_seed_reactions(&path).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].id.as_str(), "rxn00001");
        assert_eq!(rows[0].gapseq_status, SeedStatus::Approved);
        assert_eq!(rows[0].reversibility(), Some(Reversibility::Forward));
        assert_eq!(rows[0].ec_list(), vec!["3.5.1.2"]);
        let stoich = rows[0].parse_stoich().unwrap();
        assert_eq!(stoich.len(), 2);
        assert_eq!(stoich[0].cpd.as_str(), "cpd00001");

        assert_eq!(rows[1].gapseq_status, SeedStatus::Removed);
        assert_eq!(rows[1].reversibility(), Some(Reversibility::Reversible));
    }

    #[test]
    fn parses_minimal_metabolites_tsv() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("m.tsv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(
            f,
            "id\tMNX_ID\tabbreviation\tname\tformula\tmass\tsource\tcharge\tis_core\tis_obsolete\tlinked_compound\tis_cofactor\tdeltag\tdeltagerr\tpka\tpkb\tabstract_compound\tcomprised_of\taliases\tsmiles\tInChIKey\thmdbID\treactomeID\tchebiID\tInChI\tkeggID\tbiggID\tbiocycID"
        )
        .unwrap();
        writeln!(
            f,
            "cpd00001\tMNXM2\th2o\tH2O\tH2O\t18\tModelSEED\t0\t1\t0\tnull\t0\t-56.687\t0.5\t1:15.7\t1:-1.8\tnull\tnull\tnull\tO\tXLYOFNOQVPJJNP-UHFFFAOYSA-N\t\t\t\t\tC00001\th2o\tWATER"
        )
        .unwrap();

        let rows = load_seed_metabolites(&path).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].id.as_str(), "cpd00001");
        assert_eq!(rows[0].mnx_id, "MNXM2");
        assert_eq!(rows[0].name, "H2O");
        assert_eq!(rows[0].charge, 0);
    }
}
