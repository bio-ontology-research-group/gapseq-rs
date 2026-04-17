//! `dat/medium_prediction_rules.tsv` loader.
//!
//! One row per nutrient rule. Columns:
//! - `Nutrient` — human-readable name
//! - `cpd.id` — SEED compound id (cpdNNNNN)
//! - `rule` — boolean expression over pathways / reactions / compounds
//! - `maxFlux` — flux to assign when the rule evaluates true
//! - `proton.balance` — whether this entry contributes to proton balance
//! - `Comment`, `Category` — metadata

use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum RulesError {
    #[error("i/o error on `{path}`: {source}")]
    Io {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("medium_prediction_rules.tsv is missing column `{0}`")]
    MissingColumn(&'static str),
    #[error("malformed row at line {line}: {message}")]
    BadRow { line: usize, message: String },
}

#[derive(Debug, Clone)]
pub struct MediumRule {
    pub nutrient: String,
    pub cpd_id: String,
    pub rule: String,
    pub max_flux: f64,
    pub proton_balance: bool,
    pub category: String,
}

pub fn load_rules(path: &Path) -> Result<Vec<MediumRule>, RulesError> {
    // The real `dat/medium_prediction_rules.tsv` contains Latin-1
    // accented characters in a Comment column (e.g. "β-D-xylose"). Treat
    // invalid UTF-8 as lossy rather than fail hard.
    let bytes =
        std::fs::read(path).map_err(|e| RulesError::Io { path: path.to_path_buf(), source: e })?;
    let text = String::from_utf8_lossy(&bytes).into_owned();
    let mut lines = text.lines();
    let header = lines.next().ok_or(RulesError::MissingColumn("Nutrient"))?;
    let cols: Vec<&str> = header.split('\t').collect();
    let ix = |name: &str| -> Result<usize, RulesError> {
        cols.iter()
            .position(|c| *c == name)
            .ok_or(RulesError::MissingColumn(match name {
                "Nutrient" => "Nutrient",
                "cpd.id" => "cpd.id",
                "rule" => "rule",
                "maxFlux" => "maxFlux",
                "proton.balance" => "proton.balance",
                "Category" => "Category",
                _ => "",
            }))
    };
    let i_nut = ix("Nutrient")?;
    let i_cpd = ix("cpd.id")?;
    let i_rule = ix("rule")?;
    let i_flux = ix("maxFlux")?;
    let i_proton = ix("proton.balance")?;
    let i_cat = ix("Category")?;

    let mut out = Vec::new();
    for (n, line) in lines.enumerate() {
        let line_no = n + 2;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        let need = i_nut.max(i_cpd).max(i_rule).max(i_flux).max(i_proton).max(i_cat) + 1;
        if parts.len() < need {
            // Partial trailing rows exist in the real file (trailing empty
            // fields) — gapseq tolerates these.
            continue;
        }
        let cpd = parts[i_cpd].trim();
        if !cpd.starts_with("cpd") {
            continue;
        }
        let max_flux: f64 = match parts[i_flux].trim() {
            "" | "NA" => continue,
            s => s.parse().map_err(|e| RulesError::BadRow {
                line: line_no,
                message: format!("maxFlux `{s}`: {e}"),
            })?,
        };
        let proton_balance = matches!(
            parts[i_proton].trim(),
            "TRUE" | "True" | "true" | "T" | "1"
        );
        out.push(MediumRule {
            nutrient: parts[i_nut].trim().to_string(),
            cpd_id: cpd.to_string(),
            rule: parts[i_rule].trim().to_string(),
            max_flux,
            proton_balance,
            category: parts[i_cat].trim().to_string(),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn roundtrip_tiny_table() {
        let td = tempfile::tempdir().unwrap();
        let path = td.path().join("rules.tsv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(
            f,
            "Nutrient\tcpd.id\trule\tmaxFlux\tproton.balance\tComment\tCategory"
        )
        .unwrap();
        writeln!(f, "Water\tcpd00001\tTRUE\t100\tFALSE\tCore\tInorganics").unwrap();
        writeln!(f, "Glu\tcpd00027\trxn05226\t5\tTRUE\tsugar\tSaccharides").unwrap();
        drop(f);
        let rules = load_rules(&path).unwrap();
        assert_eq!(rules.len(), 2);
        assert_eq!(rules[0].cpd_id, "cpd00001");
        assert!((rules[0].max_flux - 100.0).abs() < 1e-9);
        assert!(!rules[0].proton_balance);
        assert!(rules[1].proton_balance);
    }
}
