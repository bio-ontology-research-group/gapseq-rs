//! Biomass-template parser — port of `src/parse_BMjson.R`.
//!
//! Reads `dat/biomass/biomass_<domain>.json`, parses groups of
//! metabolites (DNA, RNA, Protein, Lipid, ...), rescales component
//! coefficients so each group's mass fraction adds up to the group's
//! target g/gDW, adds GAM terms, and returns a flat list of
//! `(cpd_id[comp], Scoef)` pairs ready for a `bio1` reaction.

use gapsmith_db::{biomass as db_bm, BiomassError as DbBiomassError, SeedCpdRow};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum BiomassError {
    #[error("biomass JSON load error: {0}")]
    Db(#[from] DbBiomassError),
    #[error("seed metabolite `{0}` referenced by biomass template not found")]
    MissingMetabolite(String),
    #[error("group `{group}` components sum to {got}, expected ~1")]
    BadGroupFractions { group: String, got: f64 },
    #[error("biomass template is missing required group `{0}`")]
    MissingGroup(&'static str),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiomassEntry {
    /// `cpd00115[c0]`-style id with compartment suffix.
    pub id: String,
    pub cpd: String,
    pub compartment: char,
    pub name: String,
    pub formula: String,
    pub scoef: f64,
    pub met_group: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiomassSpec {
    pub id: String,
    pub name: String,
    pub domain: String,
    pub entries: Vec<BiomassEntry>,
}

/// Load and expand a gapseq biomass JSON template.
pub fn parse_biomass_json(
    path: &Path,
    seed_cpds: &[SeedCpdRow],
) -> Result<BiomassSpec, BiomassError> {
    let tmpl = db_bm::BiomassTemplate::load(path)?;

    // Require DNA, RNA, Protein (matches `parse_BMjson.R:33-35`).
    let group_names: Vec<&str> = tmpl.met_groups.iter().map(|g| g.group_name.as_str()).collect();
    for required in ["DNA", "RNA", "Protein"] {
        if !group_names.contains(&required) {
            return Err(BiomassError::MissingGroup(match required {
                "DNA" => "DNA",
                "RNA" => "RNA",
                _ => "Protein",
            }));
        }
    }

    // Index seed metabolites by id.
    let met_by_id: HashMap<&str, &SeedCpdRow> =
        seed_cpds.iter().map(|r| (r.id.as_str(), r)).collect();

    // Accumulate entries by (cpd, compartment, met_group) — matching R's
    // merge + sum behaviour when multiple components contribute the same
    // id (e.g., coupled reactants via `link`).
    #[derive(Default)]
    struct Acc {
        coef: f64,
        mass: f64,
        formula: String,
        name: String,
        compartment: char,
    }
    let mut by_cpd: HashMap<(String, char, String), Acc> = HashMap::new();
    for g in &tmpl.met_groups {
        // Fraction renormalization: sum of per-group coefs should be 1
        // ("MOLFRACTION"). When not, rescale.
        let group_coef_sum: f64 = g.components.iter().map(|c| c.coef).sum();
        let fraction_norm =
            if g.unit_components.eq_ignore_ascii_case("MOLFRACTION") && group_coef_sum > 0.0 {
                group_coef_sum
            } else {
                1.0
            };

        // Collect expanded rows for this group.
        let mut grp_rows: Vec<(String, char, f64)> = Vec::new();
        for c in &g.components {
            let comp_ch = c.comp.chars().next().unwrap_or('c');
            grp_rows.push((c.id.clone(), comp_ch, c.coef / fraction_norm));
            for (linked_cpd, linked_coef) in c.links() {
                grp_rows.push((
                    linked_cpd.to_string(),
                    comp_ch,
                    linked_coef * c.coef / fraction_norm,
                ));
            }
        }

        // Merge with seed metabolites to pick up formula/mass; compute
        // denom = Σ (coef × mass), then scale Scoef = group.mass × coef /
        // denom × 1000.
        let mut denom = 0.0f64;
        let mut resolved: Vec<(String, char, f64, f64, String, String)> = Vec::new();
        for (cpd, comp, coef) in &grp_rows {
            let met = met_by_id
                .get(cpd.as_str())
                .ok_or_else(|| BiomassError::MissingMetabolite(cpd.clone()))?;
            let formula = clean_formula(&met.formula);
            let mass = formula_mass(&formula);
            denom += mass * *coef;
            resolved.push((
                cpd.clone(),
                *comp,
                *coef,
                mass,
                formula.clone(),
                met.name.clone(),
            ));
        }
        if denom == 0.0 {
            continue;
        }

        for (cpd, comp, coef, mass, formula, name) in resolved {
            let scoef = g.mass * coef / denom * 1000.0;
            let entry = by_cpd
                .entry((cpd.clone(), comp, g.group_name.clone()))
                .or_insert_with(|| Acc {
                    coef: 0.0,
                    mass,
                    formula: formula.clone(),
                    name: name.clone(),
                    compartment: comp,
                });
            entry.coef += scoef;
            entry.mass = mass;
            entry.formula = formula.clone();
            entry.name = name.clone();
            entry.compartment = comp;
        }
    }

    // Attach GAM terms.
    let gam = tmpl.energy_gam;
    let gam_terms: &[(&str, f64)] = &[
        ("cpd00001", gam),  // H2O
        ("cpd00002", gam),  // ATP
        ("cpd00008", -gam), // ADP
        ("cpd00009", -gam), // Phosphate
        ("cpd00067", -gam), // H+
    ];
    for (cpd, scoef) in gam_terms {
        let met = met_by_id
            .get(cpd)
            .ok_or_else(|| BiomassError::MissingMetabolite((*cpd).to_string()))?;
        let formula = clean_formula(&met.formula);
        let entry = by_cpd
            .entry(((*cpd).to_string(), 'c', "GAM".to_string()))
            .or_insert_with(|| Acc {
                coef: 0.0,
                mass: formula_mass(&formula),
                formula: formula.clone(),
                name: met.name.clone(),
                compartment: 'c',
            });
        entry.coef += *scoef;
    }

    // Now sum across met_group per (cpd, comp) and negate — R code does
    // `grp_tab[, Scoef := -Scoef]` because biomass draws metabolites as
    // consumption (negative) except the biomass pseudo-metabolite.
    let mut by_cpd_comp: HashMap<(String, char), (f64, String, String, String)> =
        HashMap::new();
    for ((cpd, comp, grp), acc) in &by_cpd {
        let entry = by_cpd_comp
            .entry((cpd.clone(), *comp))
            .or_insert_with(|| (0.0, String::new(), String::new(), String::new()));
        entry.0 += acc.coef;
        if entry.1.is_empty() {
            entry.1.clone_from(&acc.name);
            entry.2.clone_from(&acc.formula);
        }
        if !entry.3.is_empty() && !entry.3.contains(grp.as_str()) {
            entry.3.push_str(", ");
        }
        if !entry.3.contains(grp.as_str()) {
            entry.3.push_str(grp);
        }
    }

    let mut entries: Vec<BiomassEntry> = by_cpd_comp
        .into_iter()
        .map(|((cpd, comp), (coef, name, formula, grp))| BiomassEntry {
            id: format!("{cpd}[{comp}0]"),
            cpd,
            compartment: comp,
            name,
            formula,
            scoef: -coef,
            met_group: grp,
        })
        .collect();

    // Add the biomass pseudo-metabolite.
    entries.push(BiomassEntry {
        id: "cpd11416[c0]".into(),
        cpd: "cpd11416".into(),
        compartment: 'c',
        name: "Biomass".into(),
        formula: String::new(),
        scoef: 1.0,
        met_group: "Biomass".into(),
    });

    // Deterministic order by (met_group, id).
    entries.sort_by(|a, b| a.met_group.cmp(&b.met_group).then_with(|| a.id.cmp(&b.id)));

    Ok(BiomassSpec {
        id: tmpl.id,
        name: tmpl.name,
        domain: tmpl.domain,
        entries,
    })
}

/// Strip pseudo-residue markers `R`, `R1`, `R2`, ... from the formula —
/// matches `c_list[, formula := gsub("R|R[0-9]+","",formula)]`.
fn clean_formula(f: &str) -> String {
    let mut out = String::with_capacity(f.len());
    let bytes = f.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'R' {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            i = j;
            continue;
        }
        out.push(bytes[i] as char);
        i += 1;
    }
    out
}

/// Parse a chemical formula and return the monoisotopic mass in g/mol.
/// Mirrors `cobrar::mass()` — enough precision for biomass rescaling.
fn formula_mass(formula: &str) -> f64 {
    if formula.is_empty() {
        return 0.0;
    }
    let bytes = formula.as_bytes();
    let mut i = 0;
    let mut total = 0.0;
    while i < bytes.len() {
        let c = bytes[i] as char;
        if !c.is_ascii_alphabetic() {
            i += 1;
            continue;
        }
        // Element symbol: one uppercase, optional lowercase.
        let start = i;
        i += 1;
        while i < bytes.len() && (bytes[i] as char).is_ascii_lowercase() {
            i += 1;
        }
        let sym = &formula[start..i];
        // Count digits.
        let cnt_start = i;
        while i < bytes.len() && (bytes[i] as char).is_ascii_digit() {
            i += 1;
        }
        let n: u32 = if cnt_start == i {
            1
        } else {
            formula[cnt_start..i].parse().unwrap_or(1)
        };
        total += element_mass(sym) * f64::from(n);
    }
    total
}

/// Monoisotopic atomic masses for elements that occur in SEED
/// metabolites. Values from NIST. Missing elements contribute 0 — good
/// enough for biomass rescaling (denominator dominated by C/H/O/N).
fn element_mass(sym: &str) -> f64 {
    match sym {
        "H" => 1.007825,
        "C" => 12.0,
        "N" => 14.003074,
        "O" => 15.994915,
        "P" => 30.973762,
        "S" => 31.972071,
        "Na" => 22.989770,
        "Mg" => 23.985042,
        "K" => 38.963708,
        "Ca" => 39.962591,
        "Fe" => 55.934938,
        "Cl" => 34.968853,
        "Br" => 78.918338,
        "Mn" => 54.938050,
        "Co" => 58.933195,
        "Ni" => 57.935343,
        "Cu" => 62.929598,
        "Zn" => 63.929142,
        "Mo" => 97.905408,
        "I" => 126.904473,
        "F" => 18.998403,
        "B" => 11.009305,
        "Se" => 79.916521,
        "As" => 74.921596,
        _ => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clean_formula_strips_pseudo_residues() {
        assert_eq!(clean_formula("C10H12N5O7PR"), "C10H12N5O7P");
        assert_eq!(clean_formula("C10H12NR3"), "C10H12N");
        assert_eq!(clean_formula("C10H12"), "C10H12");
    }

    #[test]
    fn formula_mass_water() {
        let m = formula_mass("H2O");
        assert!((m - 18.010565).abs() < 0.001, "got {m}");
    }
}
