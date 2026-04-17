//! Rule-based medium prediction. Port of `src/predict_medium.R:42-135`.

use crate::boolexpr::{eval, BoolExprError};
use crate::rules::MediumRule;
use gapsmith_core::Model;
use gapsmith_db::SeedCpdRow;
use std::collections::{HashMap, HashSet};

#[derive(Debug, thiserror::Error)]
pub enum MediumPredictError {
    #[error("bool expr error on rule `{rule}`: {source}")]
    Expr { rule: String, #[source] source: BoolExprError },
    #[error("bad manual flux `{0}`")]
    BadManual(String),
}

#[derive(Debug, Clone)]
pub struct PredictedMedium {
    pub compounds: Vec<MediumLine>,
}

#[derive(Debug, Clone)]
pub struct MediumLine {
    pub cpd_id: String,
    pub name: String,
    pub max_flux: f64,
    pub category: String,
    pub charge: Option<f64>,
}

impl PredictedMedium {
    pub fn write_csv<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(w, "compounds,name,maxFlux")?;
        for c in &self.compounds {
            writeln!(w, "{},{},{}", c.cpd_id, c.name, fmt_flux(c.max_flux))?;
        }
        Ok(())
    }
}

fn fmt_flux(v: f64) -> String {
    if (v - v.round()).abs() < 1e-9 {
        format!("{}", v as i64)
    } else {
        format!("{v}")
    }
}

/// Manual flux overrides supplied as `cpd00027:10;cpd00007:18.5`.
pub fn parse_manual_flux(s: &str) -> Result<Vec<(String, f64)>, MediumPredictError> {
    if s.is_empty() {
        return Ok(Vec::new());
    }
    let mut out = Vec::new();
    for term in s.split(';').filter(|t| !t.is_empty()) {
        let (cpd, flux) = term
            .split_once(':')
            .ok_or_else(|| MediumPredictError::BadManual(term.to_string()))?;
        let flux: f64 = flux
            .trim()
            .parse()
            .map_err(|_| MediumPredictError::BadManual(term.to_string()))?;
        out.push((cpd.trim().to_string(), flux));
    }
    Ok(out)
}

/// Predict a growth medium for `model` given a set of predicted pathways
/// and the rules table. Mirrors `predict_medium.R`:
///
/// 1. For each rule, eval the boolean expression where an atom is TRUE iff
///    it matches a predicted pathway / a reaction in the model / a
///    compound in the model.
/// 2. Dedup by cpd id, averaging maxFlux across matching rules.
/// 3. In Saccharides / Organic acids categories, keep only the first rule
///    for each category (matches gapseq's preference for glucose over
///    other sugars).
/// 4. Apply manual flux overrides.
/// 5. Proton balance — if the net charge of `proton.balance=TRUE` rows is
///    negative, add enough H+ (cpd00067) to compensate.
pub fn predict_medium(
    model: &Model,
    predicted_pathways: &HashSet<String>,
    rules: &[MediumRule],
    manual_flux: &[(String, f64)],
    seed_cpds: &[SeedCpdRow],
) -> Result<PredictedMedium, MediumPredictError> {
    // -- Atoms present in the model/pathways --
    let rxn_ids: HashSet<String> = model
        .rxns
        .iter()
        .filter_map(|r| {
            let s = r.id.as_str();
            let trimmed = s.trim_end_matches("_c0").trim_end_matches("_e0").trim_end_matches("_p0");
            if trimmed.starts_with("rxn") { Some(trimmed.to_string()) } else { None }
        })
        .collect();
    let cpd_ids: HashSet<String> = model
        .mets
        .iter()
        .filter_map(|m| {
            let s = m.id.as_str();
            let trimmed = s.trim_end_matches("_c0").trim_end_matches("_e0").trim_end_matches("_p0");
            if trimmed.starts_with("cpd") { Some(trimmed.to_string()) } else { None }
        })
        .collect();

    // -- Evaluate each rule --
    struct Hit<'a> {
        rule: &'a MediumRule,
    }
    let mut hits: Vec<Hit> = Vec::new();
    for rule in rules {
        let r = &rule.rule;
        let ok = eval(r, |tok| {
            predicted_pathways.contains(tok)
                || rxn_ids.contains(tok)
                || cpd_ids.contains(tok)
        })
        .map_err(|e| MediumPredictError::Expr { rule: r.clone(), source: e })?;
        if ok {
            hits.push(Hit { rule });
        }
    }

    // -- Dedup by cpd id, mean of max_flux across matches --
    let mut by_cpd: HashMap<String, Vec<&Hit>> = HashMap::new();
    for h in &hits {
        by_cpd.entry(h.rule.cpd_id.clone()).or_default().push(h);
    }
    let mut merged: Vec<MediumEntryInternal> = Vec::new();
    for (cpd, group) in by_cpd {
        let mean_flux: f64 = group.iter().map(|h| h.rule.max_flux).sum::<f64>() / group.len() as f64;
        let first = group[0].rule;
        merged.push(MediumEntryInternal {
            cpd_id: cpd,
            nutrient: first.nutrient.clone(),
            category: first.category.clone(),
            max_flux: mean_flux,
            proton_balance: first.proton_balance,
        });
    }

    // Deterministic ordering — the R code keeps the original row order via
    // data.table merge. We emulate by sorting by the first-match position
    // in `rules`.
    let rule_pos: HashMap<&str, usize> = rules
        .iter()
        .enumerate()
        .map(|(i, r)| (r.cpd_id.as_str(), i))
        .collect();
    merged.sort_by_key(|m| *rule_pos.get(m.cpd_id.as_str()).unwrap_or(&usize::MAX));

    // -- Drop duplicates in Saccharides / Organic acids: keep the first --
    let mut seen: HashMap<String, usize> = HashMap::new();
    merged.retain(|m| {
        let count = seen.entry(m.category.clone()).or_insert(0);
        *count += 1;
        if matches!(m.category.as_str(), "Saccharides" | "Organic acids") {
            *count == 1
        } else {
            true
        }
    });

    // -- Apply manual flux overrides --
    for (cpd, flux) in manual_flux {
        if let Some(existing) = merged.iter_mut().find(|m| &m.cpd_id == cpd) {
            existing.max_flux = *flux;
        } else {
            merged.push(MediumEntryInternal {
                cpd_id: cpd.clone(),
                nutrient: cpd.clone(),
                category: String::new(),
                max_flux: *flux,
                proton_balance: false,
            });
        }
    }

    // -- Look up name + charge from SEED compound table --
    let meta: HashMap<&str, &SeedCpdRow> =
        seed_cpds.iter().map(|r| (r.id.as_str(), r)).collect();
    let mut out_lines = Vec::<MediumLine>::new();
    for m in &merged {
        let (name, charge) = match meta.get(m.cpd_id.as_str()) {
            Some(c) => (
                if c.name.is_empty() { m.nutrient.clone() } else { c.name.clone() },
                Some(f64::from(c.charge)),
            ),
            None => (m.nutrient.clone(), None),
        };
        out_lines.push(MediumLine {
            cpd_id: m.cpd_id.clone(),
            name,
            max_flux: m.max_flux,
            category: m.category.clone(),
            charge,
        });
    }

    // -- Proton balance --
    let mut balance_pairs: Vec<(f64, f64)> = Vec::new();
    for (i, m) in merged.iter().enumerate() {
        if m.proton_balance {
            if let Some(charge) = out_lines[i].charge {
                balance_pairs.push((m.max_flux, charge));
            }
        }
    }
    let net_charge: f64 = balance_pairs.iter().map(|(f, c)| f * c).sum();
    if net_charge < -1e-9 && !out_lines.iter().any(|l| l.cpd_id == "cpd00067") {
        let h_charge: f64 = meta
            .get("cpd00067")
            .map(|r| f64::from(r.charge))
            .unwrap_or(1.0);
        let denom = h_charge.abs().max(1.0);
        out_lines.push(MediumLine {
            cpd_id: "cpd00067".into(),
            name: "H+".into(),
            max_flux: -net_charge / denom,
            category: "Inorganics".into(),
            charge: Some(h_charge),
        });
    }

    Ok(PredictedMedium { compounds: out_lines })
}

#[derive(Debug, Clone)]
struct MediumEntryInternal {
    cpd_id: String,
    nutrient: String,
    category: String,
    max_flux: f64,
    proton_balance: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blank_model() -> Model {
        Model::new("m")
    }

    fn make_rules(rows: &[(&str, &str, &str, f64, bool, &str)]) -> Vec<MediumRule> {
        rows.iter()
            .map(|(n, c, r, f, pb, cat)| MediumRule {
                nutrient: n.to_string(),
                cpd_id: c.to_string(),
                rule: r.to_string(),
                max_flux: *f,
                proton_balance: *pb,
                category: cat.to_string(),
            })
            .collect()
    }

    #[test]
    fn always_included_rule() {
        let rules = make_rules(&[("Water", "cpd00001", "TRUE", 100.0, false, "Inorganics")]);
        let predicted = HashSet::new();
        let pm = predict_medium(&blank_model(), &predicted, &rules, &[], &[]).unwrap();
        assert_eq!(pm.compounds.len(), 1);
        assert_eq!(pm.compounds[0].cpd_id, "cpd00001");
    }

    #[test]
    fn pathway_gated_rule() {
        let rules = make_rules(&[
            ("Water", "cpd00001", "TRUE", 100.0, false, "Inorganics"),
            ("O2", "cpd00007", r#""PWY-7279""#, 12.5, false, "Inorganics"),
        ]);
        let with_pwy: HashSet<String> = ["PWY-7279".into()].into_iter().collect();
        let without: HashSet<String> = HashSet::new();
        let pm = predict_medium(&blank_model(), &with_pwy, &rules, &[], &[]).unwrap();
        assert!(pm.compounds.iter().any(|c| c.cpd_id == "cpd00007"));
        let pm = predict_medium(&blank_model(), &without, &rules, &[], &[]).unwrap();
        assert!(!pm.compounds.iter().any(|c| c.cpd_id == "cpd00007"));
    }

    #[test]
    fn saccharides_dedup_keeps_first() {
        let rules = make_rules(&[
            ("D-Glucose", "cpd00027", "TRUE", 5.0, true, "Saccharides"),
            ("D-Mannose", "cpd00138", "TRUE", 5.0, true, "Saccharides"),
        ]);
        let pm = predict_medium(&blank_model(), &HashSet::new(), &rules, &[], &[]).unwrap();
        assert_eq!(pm.compounds.len(), 1);
        assert_eq!(pm.compounds[0].cpd_id, "cpd00027");
    }

    #[test]
    fn manual_flux_overrides() {
        let rules = make_rules(&[("O2", "cpd00007", "TRUE", 12.5, false, "Inorganics")]);
        let manual = vec![("cpd00007".to_string(), 0.0)];
        let pm = predict_medium(&blank_model(), &HashSet::new(), &rules, &manual, &[]).unwrap();
        let o2 = pm.compounds.iter().find(|c| c.cpd_id == "cpd00007").unwrap();
        assert_eq!(o2.max_flux, 0.0);
    }

    #[test]
    fn manual_flux_adds_new_compound() {
        let rules = make_rules(&[("Water", "cpd00001", "TRUE", 100.0, false, "Inorganics")]);
        let manual = vec![("cpd99999".to_string(), 1.0)];
        let pm = predict_medium(&blank_model(), &HashSet::new(), &rules, &manual, &[]).unwrap();
        assert!(pm.compounds.iter().any(|c| c.cpd_id == "cpd99999"));
    }
}
