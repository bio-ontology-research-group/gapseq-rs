//! Stoichiometric hash — port of `src/generate_rxn_stoich_hash.R`.
//!
//! Produces a canonical `String` fingerprint of a SEED-format
//! stoichiometry + reversibility code. Used during the draft-build step
//! to drop duplicate reactions (same chemistry, different SEED id).
//!
//! Canonicalization rules (faithful to the R code):
//!
//! 1. Sort `(met_id[comp], coef)` pairs by composite key `met[comp]`.
//! 2. If `reversibility == "="` and the first coef is negative, invert
//!    every coef sign.
//! 3. If `reversibility == "<"`, invert every coef sign (treated as
//!    `">"` after inversion).
//! 4. Concatenate `coef:met[comp]` terms with `;`, append `#<rev>`.

use gapsmith_db::parse_stoichiometry;

pub fn rxn_stoich_hash(stoich: &str, rev: &str) -> Result<String, String> {
    let terms = parse_stoichiometry(stoich).map_err(|e| e.to_string())?;
    if terms.is_empty() {
        return Err("empty stoichiometry".into());
    }

    let mut rows: Vec<(String, i64)> = terms
        .iter()
        .map(|t| {
            (
                format!("{}[{}]", t.cpd.as_str(), t.compartment),
                t.coef as i64,
            )
        })
        .collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));

    let mut rev_out = rev.to_string();
    match rev {
        "=" => {
            if let Some((_, c)) = rows.first() {
                if *c < 0 {
                    for (_, c) in &mut rows {
                        *c = -*c;
                    }
                }
            }
        }
        "<" => {
            for (_, c) in &mut rows {
                *c = -*c;
            }
            rev_out = ">".into();
        }
        _ => {}
    }

    let body: Vec<String> = rows.iter().map(|(m, c)| format!("{c}:{m}")).collect();
    Ok(format!("{}#{rev_out}", body.join(";")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reversible_with_neg_first_is_flipped() {
        // Same reaction, different orientations (`=`) should hash the same.
        let a = rxn_stoich_hash(
            r#"-1:cpd00001:0:0:"H2O";1:cpd00002:0:0:"ATP""#,
            "=",
        )
        .unwrap();
        let b = rxn_stoich_hash(
            r#"1:cpd00001:0:0:"H2O";-1:cpd00002:0:0:"ATP""#,
            "=",
        )
        .unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn backward_converts_to_forward() {
        let a = rxn_stoich_hash(
            r#"-1:cpd00001:0:0:"H2O";1:cpd00002:0:0:"ATP""#,
            ">",
        )
        .unwrap();
        let b = rxn_stoich_hash(
            r#"1:cpd00001:0:0:"H2O";-1:cpd00002:0:0:"ATP""#,
            "<",
        )
        .unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn distinct_reactions_hash_differently() {
        let a = rxn_stoich_hash(r#"-1:cpd00001:0:0:"H2O";1:cpd00002:0:0:"ATP""#, ">").unwrap();
        let b = rxn_stoich_hash(r#"-1:cpd00001:0:0:"H2O";1:cpd00003:0:0:"ADP""#, ">").unwrap();
        assert_ne!(a, b);
    }
}
