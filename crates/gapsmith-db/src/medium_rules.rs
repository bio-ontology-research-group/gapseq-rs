//! `dat/medium_prediction_rules.tsv` loader.
//!
//! Columns: `Nutrient, cpd.id, rule, maxFlux, proton.balance, Comment, Category`.
//! The `rule` column is a Boolean expression compiled at medium-prediction
//! time (see `src/predict_medium.R:46–86`). For now we keep it verbatim —
//! parsing happens later in the `gapsmith-medium` crate.

use crate::common::{csv_err, io_err, DbError};
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MediumRule {
    pub nutrient: String,
    pub cpd_id: String,
    pub rule: String,
    /// `None` for rows that ship `NA` (e.g. "do not supply O2 under anaerobic").
    pub max_flux: Option<f64>,
    pub proton_balance: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub comment: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub category: String,
}

pub fn load(path: impl AsRef<Path>) -> Result<Vec<MediumRule>, DbError> {
    let path = path.as_ref();
    let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .quoting(false)
        .flexible(true)
        .from_reader(f);

    // `medium_prediction_rules.tsv` in the reference data is Latin-1, not
    // UTF-8 (old source material with `“ ”` smart quotes). We therefore read
    // `byte_records()` and convert each field via `String::from_utf8_lossy`.
    let headers_raw = rdr.byte_headers().map_err(|e| csv_err(path, e))?.clone();
    let headers: Vec<String> = headers_raw
        .iter()
        .map(|h| String::from_utf8_lossy(h).trim().to_string())
        .collect();
    let col = |name: &str| -> Option<usize> {
        headers.iter().position(|h| h == name)
    };
    let c_nut = col("Nutrient").ok_or_else(|| DbError::Parse {
        path: path.to_path_buf(),
        line: 1,
        msg: "missing `Nutrient` column".into(),
    })?;
    let c_cpd = col("cpd.id").ok_or_else(|| DbError::Parse {
        path: path.to_path_buf(),
        line: 1,
        msg: "missing `cpd.id` column".into(),
    })?;
    let c_rule = col("rule").ok_or_else(|| DbError::Parse {
        path: path.to_path_buf(),
        line: 1,
        msg: "missing `rule` column".into(),
    })?;
    let c_flux = col("maxFlux").ok_or_else(|| DbError::Parse {
        path: path.to_path_buf(),
        line: 1,
        msg: "missing `maxFlux` column".into(),
    })?;
    let c_proton = col("proton.balance").unwrap_or(usize::MAX);
    let c_comment = col("Comment").unwrap_or(usize::MAX);
    let c_cat = col("Category").unwrap_or(usize::MAX);

    let mut out = Vec::new();
    for rec in rdr.byte_records() {
        let rec = rec.map_err(|e| csv_err(path, e))?;
        let get = |c: usize| -> String {
            if c == usize::MAX {
                String::new()
            } else {
                rec.get(c)
                    .map(|b| String::from_utf8_lossy(b).trim().to_string())
                    .unwrap_or_default()
            }
        };
        let nutrient = get(c_nut);
        let cpd_id = get(c_cpd);
        let rule = get(c_rule);
        // Blank separator rows have every field empty — skip them silently.
        if nutrient.is_empty() && cpd_id.is_empty() && rule.is_empty() {
            continue;
        }
        let flux_raw = get(c_flux);
        let max_flux = match flux_raw.as_str() {
            "" | "NA" | "na" | "N/A" => None,
            other => Some(other.parse::<f64>().map_err(|_| DbError::Parse {
                path: path.to_path_buf(),
                line: rec.position().map(|p| p.line()).unwrap_or(0),
                msg: format!("maxFlux `{other}` is not a number"),
            })?),
        };
        out.push(MediumRule {
            nutrient,
            cpd_id,
            rule,
            max_flux,
            proton_balance: get(c_proton),
            comment: get(c_comment),
            category: get(c_cat),
        });
    }
    tracing::info!(path = %path.display(), rows = out.len(), "loaded medium rules");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_rules() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("r.tsv");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(
            f,
            "Nutrient\tcpd.id\trule\tmaxFlux\tproton.balance\tComment\tCategory"
        )
        .unwrap();
        writeln!(f, "Water\tcpd00001\tTRUE\t100\tFALSE\tCore medium compound\tInorganics").unwrap();
        writeln!(f, "O2\tcpd00007\tpwy1\tNA\tFALSE\tno O2\tInorganics").unwrap();
        writeln!(f, "\t\t\t\t\t\t").unwrap(); // blank separator
        writeln!(f, "Glc\tcpd00027\trxn1\t5\tTRUE\t\tSaccharides").unwrap();
        let rows = load(&p).unwrap();
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].cpd_id, "cpd00001");
        assert_eq!(rows[0].max_flux, Some(100.0));
        assert_eq!(rows[1].max_flux, None);
        assert_eq!(rows[2].cpd_id, "cpd00027");
    }
}
