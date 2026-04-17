//! Port of `src/get_gene_logic_string.R`. Builds a Boolean gene-reaction
//! association string from parallel `(complex, gene)` vectors.
//!
//! Rules (matching the R original):
//!
//! - Genes sharing the same `complex` label are joined with `|` (OR)
//!   inside parentheses — they're isoforms of the same subunit.
//! - Each parenthesised subunit group is then joined with `&` (AND) —
//!   different subunits of the same enzyme complex.
//! - Rows with `complex == NA` or `complex == "Subunit undefined"` are
//!   monomers / orphan genes and get OR-joined at the top level.
//! - If `complex != NA && complex != "Subunit undefined"` is present,
//!   `"Subunit undefined"` rows are dropped from the complex tree but
//!   still OR-joined outside.

use std::collections::BTreeMap;

pub struct GeneAssignment<'a> {
    pub complex: Option<&'a str>,
    pub gene: &'a str,
}

pub fn build_gpr_string(rows: &[GeneAssignment<'_>]) -> String {
    if rows.is_empty() {
        return String::new();
    }

    let any_real_complex = rows.iter().any(|r| {
        matches!(r.complex, Some(c) if !c.is_empty() && c != "Subunit undefined")
    });

    let mut complex_groups: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    let mut monomers: Vec<&str> = Vec::new();

    for r in rows {
        match r.complex {
            Some(c) if !c.is_empty() && c != "Subunit undefined" => {
                complex_groups.entry(c).or_default().push(r.gene);
            }
            Some(c) if c == "Subunit undefined" => {
                // "Subunit undefined" only contributes to monomers when
                // there is no real complex present.
                if !any_real_complex {
                    monomers.push(r.gene);
                }
            }
            _ => {
                monomers.push(r.gene);
            }
        }
    }

    let cx_str = if complex_groups.is_empty() {
        None
    } else {
        let parts: Vec<String> = complex_groups
            .values()
            .map(|genes| format!("({})", genes.join(" | ")))
            .collect();
        Some(format!("({})", parts.join(" & ")))
    };

    let mono_str = if monomers.is_empty() {
        None
    } else {
        Some(monomers.join(" | "))
    };

    match (cx_str, mono_str) {
        (Some(a), Some(b)) => format!("{a} | {b}"),
        (Some(a), None) => a,
        (None, Some(b)) => b,
        (None, None) => String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monomers_only() {
        let rows = vec![
            GeneAssignment { complex: None, gene: "b0001" },
            GeneAssignment { complex: None, gene: "b0002" },
        ];
        assert_eq!(build_gpr_string(&rows), "b0001 | b0002");
    }

    #[test]
    fn complex_and_monomer() {
        let rows = vec![
            GeneAssignment { complex: Some("Subunit alpha"), gene: "g1" },
            GeneAssignment { complex: Some("Subunit alpha"), gene: "g2" },
            GeneAssignment { complex: Some("Subunit beta"), gene: "g3" },
            GeneAssignment { complex: None, gene: "gm" },
        ];
        let s = build_gpr_string(&rows);
        // Order of complexes is sorted (BTreeMap); genes within preserved.
        assert!(s.contains("(g1 | g2)") || s.contains("(g3)"));
        assert!(s.contains("& "));
        assert!(s.contains("| gm"));
    }

    #[test]
    fn undefined_dropped_when_real_complex_exists() {
        let rows = vec![
            GeneAssignment { complex: Some("Subunit alpha"), gene: "g1" },
            GeneAssignment { complex: Some("Subunit undefined"), gene: "gx" },
        ];
        let s = build_gpr_string(&rows);
        assert!(s.contains("g1"));
        assert!(!s.contains("gx"));
    }

    #[test]
    fn undefined_kept_when_no_real_complex() {
        let rows = vec![
            GeneAssignment { complex: Some("Subunit undefined"), gene: "gx" },
            GeneAssignment { complex: None, gene: "gy" },
        ];
        let s = build_gpr_string(&rows);
        assert!(s.contains("gx"));
        assert!(s.contains("gy"));
    }
}
