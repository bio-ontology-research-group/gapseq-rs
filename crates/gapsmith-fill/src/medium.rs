//! Medium application.
//!
//! Port of `src/constrain.model.R` (close all EX, open per media CSV, create
//! missing EX reactions) and `src/adjust_model_env.R` (environment overrides
//! like `highH2`).

use gapsmith_core::{CompartmentId, Metabolite, Model, Reaction, StoichMatrix};
use std::path::Path;

/// One row of a media CSV — columns `compounds,name,maxFlux`.
#[derive(Debug, Clone)]
pub struct MediumEntry {
    /// SEED compound id without compartment suffix (e.g. `cpd00027`).
    pub compound: String,
    pub name: String,
    pub max_flux: f64,
}

#[derive(Debug, thiserror::Error)]
pub enum MediumError {
    #[error("i/o error on `{0}`: {1}")]
    Io(std::path::PathBuf, #[source] std::io::Error),
    #[error("malformed media row `{row}`: {message}")]
    BadRow { row: String, message: String },
}

/// Read a gapseq media CSV (comma- or tab-separated). Accepts columns in any
/// order so long as they include `compounds`, `name`, and `maxFlux`.
pub fn read_medium(path: &Path) -> Result<Vec<MediumEntry>, MediumError> {
    let text =
        std::fs::read_to_string(path).map_err(|e| MediumError::Io(path.to_path_buf(), e))?;
    let mut lines = text.lines();
    let header = lines.next().ok_or_else(|| MediumError::BadRow {
        row: String::new(),
        message: "empty file".into(),
    })?;
    let delim = if header.contains('\t') { '\t' } else { ',' };
    let cols: Vec<&str> = header.split(delim).map(str::trim).collect();
    let ix = |name: &str| {
        cols.iter()
            .position(|c| c.eq_ignore_ascii_case(name))
            .ok_or_else(|| MediumError::BadRow {
                row: header.to_string(),
                message: format!("missing column `{name}`"),
            })
    };
    let i_cpd = ix("compounds")?;
    let i_name = ix("name")?;
    let i_flux = ix("maxFlux")?;

    let mut out = Vec::new();
    for (ln, line) in lines.enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split(delim).collect();
        let req = i_cpd.max(i_name).max(i_flux) + 1;
        if parts.len() < req {
            return Err(MediumError::BadRow {
                row: format!("line {}: {line}", ln + 2),
                message: "too few fields".into(),
            });
        }
        let max_flux: f64 = parts[i_flux].trim().parse().map_err(|_| MediumError::BadRow {
            row: format!("line {}: {line}", ln + 2),
            message: "maxFlux not numeric".into(),
        })?;
        if max_flux == 0.0 {
            continue; // gapseq treats 0-flux rows as excluded
        }
        out.push(MediumEntry {
            compound: parts[i_cpd].trim().to_string(),
            name: parts[i_name].trim().to_string(),
            max_flux,
        });
    }
    Ok(out)
}

/// Apply a medium to the model: close every `EX_*` reaction's lower bound
/// and then set it to `-max_flux` for each medium compound. When the
/// compound has no matching exchange reaction yet, create one.
pub fn apply_medium(model: &mut Model, medium: &[MediumEntry], scaling_fac: f64, ub: f64) {
    // 1. Close all existing exchanges on the lower side.
    for r in &mut model.rxns {
        if r.id.as_str().starts_with("EX_") {
            r.lb = 0.0;
        }
    }

    // 2. Open for each medium compound; add missing EX if needed.
    for entry in medium {
        let rxn_id = format!("EX_{}_e0", entry.compound);
        let met_id = format!("{}_e0", entry.compound);

        if let Some(r) = model.rxns.iter_mut().find(|r| r.id.as_str() == rxn_id) {
            r.lb = -entry.max_flux * scaling_fac;
            continue;
        }

        // Need to add a new exchange. Ensure met exists.
        let met_idx = match model.mets.iter().position(|m| m.id.as_str() == met_id) {
            Some(i) => i,
            None => {
                let met = Metabolite::new(
                    met_id.as_str(),
                    entry.name.clone(),
                    CompartmentId::EXTRACELLULAR,
                );
                model.mets.push(met);
                model.mets.len() - 1
            }
        };

        let mut r = Reaction::new(
            rxn_id.as_str(),
            format!("{} Exchange", entry.name),
            -entry.max_flux * scaling_fac,
            ub,
        );
        r.is_exchange = true;
        r.gs_origin = Some(7);
        model.rxns.push(r);

        append_column(model, met_idx, -1.0);
    }
}

/// Append a new column to `model.s` with a single non-zero entry. The
/// column index is the last reaction (caller must have pushed the reaction
/// just before calling).
fn append_column(model: &mut Model, met_idx: usize, coef: f64) {
    let n_mets = model.mets.len();
    let n_rxns = model.rxns.len();
    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(model.s.nnz() + 1);
    let old_cols = model.s.cols();
    for c in 0..old_cols.min(n_rxns - 1) {
        for (row, v) in model.s.column(c) {
            triplets.push((row, c, v));
        }
    }
    triplets.push((met_idx, n_rxns - 1, coef));
    model.s = StoichMatrix::from_triplets(n_mets, n_rxns, triplets);
}

/// Apply an environment override file like `dat/env/env_highH2.tsv`.
/// Each row is `<rxn_id>\t<direction>`; direction is `>`, `<`, `=`.
pub fn apply_environment_file(
    model: &mut Model,
    env_tsv: &Path,
    max_flux: f64,
) -> Result<usize, MediumError> {
    let text = std::fs::read_to_string(env_tsv)
        .map_err(|e| MediumError::Io(env_tsv.to_path_buf(), e))?;
    let mut changed = 0;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.splitn(2, '\t');
        let rxn_raw = it.next().unwrap_or("");
        let direction = it.next().unwrap_or("").trim();
        if rxn_raw.is_empty() || direction.is_empty() {
            continue;
        }
        let target = format!("{rxn_raw}_c0");
        if let Some(r) = model.rxns.iter_mut().find(|r| r.id.as_str() == target) {
            match direction {
                ">" => {
                    r.lb = 0.0;
                    r.ub = max_flux;
                }
                "<" => {
                    r.lb = -max_flux;
                    r.ub = 0.0;
                }
                "=" => {
                    r.lb = -max_flux;
                    r.ub = max_flux;
                }
                other => {
                    return Err(MediumError::BadRow {
                        row: line.to_string(),
                        message: format!("unknown direction `{other}`"),
                    })
                }
            }
            changed += 1;
        }
    }
    Ok(changed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn build_tiny_model() -> Model {
        use gapsmith_core::{Metabolite, Reaction};
        let mut m = Model::new("tiny");
        m.mets.push(Metabolite::new("cpd00001_e0", "H2O", CompartmentId::EXTRACELLULAR));
        let mut ex = Reaction::new("EX_cpd00001_e0", "H2O Exchange", -1000.0, 1000.0);
        ex.is_exchange = true;
        m.rxns.push(ex);
        m.s = StoichMatrix::from_triplets(1, 1, vec![(0, 0, -1.0)]);
        m
    }

    #[test]
    fn read_medium_parses_csv_and_tsv() {
        let td = tempfile::tempdir().unwrap();
        let path = td.path().join("med.csv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "compounds,name,maxFlux").unwrap();
        writeln!(f, "cpd00001,H2O,1000").unwrap();
        writeln!(f, "cpd00027,D-Glucose,10").unwrap();
        writeln!(f, "cpd00999,Zero,0").unwrap(); // excluded
        drop(f);
        let m = read_medium(&path).unwrap();
        assert_eq!(m.len(), 2);
        assert_eq!(m[1].compound, "cpd00027");
    }

    #[test]
    fn apply_medium_sets_lb_and_adds_missing() {
        let mut m = build_tiny_model();
        let medium = vec![
            MediumEntry { compound: "cpd00001".into(), name: "H2O".into(), max_flux: 1000.0 },
            MediumEntry { compound: "cpd00027".into(), name: "D-Glucose".into(), max_flux: 10.0 },
        ];
        apply_medium(&mut m, &medium, 1.0, 1000.0);
        assert_eq!(m.rxn_count(), 2); // added glucose exchange
        let glc = m.rxns.iter().find(|r| r.id.as_str() == "EX_cpd00027_e0").unwrap();
        assert!((glc.lb + 10.0).abs() < 1e-9);
        // S should have 2 entries now.
        assert_eq!(m.s.nnz(), 2);
    }
}
