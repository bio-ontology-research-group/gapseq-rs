//! Parser for SEED stoichiometry strings.
//!
//! Format (one term per reaction participant, `;`-separated):
//!
//! ```text
//! <coef>:<cpd_id>:<compartment>:<reserved>:"<name>"
//! ```
//!
//! Example (from `dat/seed_reactions_corrected.tsv`):
//!
//! ```text
//! -1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+"
//! ```
//!
//! Cross-ref: `src/construct_full_model.R:28–35` in the R reference
//! implementation.
//!
//! Design notes:
//!
//! - We split on the top-level `;` and `:` only — a `;` or `:` inside a
//!   double-quoted `name` field must be preserved. In practice the current
//!   gapseq data does not contain such cases, but the parser handles them
//!   anyway so future DB updates don't break us silently.
//! - The name field is optional (the whole trailing `:"..."` may be absent).
//! - Leading/trailing whitespace is trimmed.

use gapsmith_core::CpdId;

#[derive(Clone, Debug, PartialEq)]
pub struct StoichTerm {
    pub coef: f64,
    pub cpd: CpdId,
    pub compartment: u8,
    /// Reserved field seen as `0` in every real row. Preserved for fidelity.
    pub reserved: u8,
    /// Human-readable metabolite name as given in the stoichiometry field.
    /// May be empty if absent from the source.
    pub name: String,
}

#[derive(Debug, thiserror::Error)]
pub enum StoichParseError {
    #[error("term {index} (`{fragment}`): {msg}")]
    BadTerm {
        index: usize,
        fragment: String,
        msg: String,
    },
    #[error("empty stoichiometry string")]
    Empty,
}

/// Parse a full stoichiometry field.
pub fn parse_stoichiometry(s: &str) -> Result<Vec<StoichTerm>, StoichParseError> {
    let s = s.trim();
    if s.is_empty() {
        return Err(StoichParseError::Empty);
    }
    let mut out = Vec::new();
    for (i, raw_term) in split_top_level(s, ';').iter().enumerate() {
        let t = raw_term.trim();
        if t.is_empty() {
            continue;
        }
        out.push(parse_term(i, t)?);
    }
    if out.is_empty() {
        return Err(StoichParseError::Empty);
    }
    Ok(out)
}

fn parse_term(index: usize, t: &str) -> Result<StoichTerm, StoichParseError> {
    let parts = split_top_level(t, ':');
    if parts.len() < 3 {
        return Err(StoichParseError::BadTerm {
            index,
            fragment: t.to_string(),
            msg: format!("expected at least 3 `:`-separated fields, got {}", parts.len()),
        });
    }

    let coef: f64 = parts[0].trim().parse().map_err(|_| StoichParseError::BadTerm {
        index,
        fragment: t.to_string(),
        msg: format!("coefficient `{}` is not a number", parts[0]),
    })?;
    let cpd = CpdId::new(parts[1].trim());
    let compartment: u8 = parts[2].trim().parse().map_err(|_| StoichParseError::BadTerm {
        index,
        fragment: t.to_string(),
        msg: format!("compartment `{}` is not u8", parts[2]),
    })?;
    let reserved: u8 = if parts.len() >= 4 {
        parts[3].trim().parse().unwrap_or(0)
    } else {
        0
    };
    let name = if parts.len() >= 5 {
        // Remaining fields may contain `:` inside quoted names — rejoin.
        let joined = parts[4..].join(":");
        strip_quotes(joined.trim())
    } else {
        String::new()
    };
    Ok(StoichTerm { coef, cpd, compartment, reserved, name })
}

/// Split on `delim` at the top level — i.e. outside any `"..."` double-quoted
/// section. Quote state is reset by a non-escaped `"`.
fn split_top_level(s: &str, delim: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut buf = String::new();
    let mut in_quotes = false;
    let mut prev_backslash = false;
    for c in s.chars() {
        if prev_backslash {
            buf.push(c);
            prev_backslash = false;
            continue;
        }
        match c {
            '\\' => {
                buf.push(c);
                prev_backslash = true;
            }
            '"' => {
                in_quotes = !in_quotes;
                buf.push(c);
            }
            c if c == delim && !in_quotes => {
                out.push(std::mem::take(&mut buf));
            }
            _ => buf.push(c),
        }
    }
    out.push(buf);
    out
}

fn strip_quotes(s: &str) -> String {
    let bytes = s.as_bytes();
    if bytes.len() >= 2 && bytes.first() == Some(&b'"') && bytes.last() == Some(&b'"') {
        s[1..s.len() - 1].to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_row() {
        let s = r#"-1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+""#;
        let terms = parse_stoichiometry(s).unwrap();
        assert_eq!(terms.len(), 4);
        assert_eq!(terms[0].coef, -1.0);
        assert_eq!(terms[0].cpd.as_str(), "cpd00001");
        assert_eq!(terms[0].compartment, 0);
        assert_eq!(terms[0].name, "H2O");
        assert_eq!(terms[2].coef, 2.0);
        assert_eq!(terms[3].name, "H+");
    }

    #[test]
    fn five_term_real_row() {
        // Actual row observed in `dat/seed_reactions_corrected.tsv`.
        let s = r#"-1:cpd00001:0:0:"H2O";-3:cpd00067:0:0:"H+";-1:cpd00742:0:0:"Allophanate";2:cpd00011:0:0:"CO2";2:cpd00013:0:0:"NH3""#;
        let terms = parse_stoichiometry(s).unwrap();
        assert_eq!(terms.len(), 5);
        assert_eq!(terms.iter().map(|t| t.coef).sum::<f64>(), -1.0 - 3.0 - 1.0 + 2.0 + 2.0);
    }

    #[test]
    fn colon_inside_quoted_name_preserved() {
        let s = r#"-1:cpd00001:0:0:"X: Y: Z";1:cpd00002:0:0:"W""#;
        let terms = parse_stoichiometry(s).unwrap();
        assert_eq!(terms.len(), 2);
        assert_eq!(terms[0].name, "X: Y: Z");
    }

    #[test]
    fn semicolon_inside_quoted_name_preserved() {
        let s = r#"-1:cpd00001:0:0:"A; B";1:cpd00002:0:0:"C""#;
        let terms = parse_stoichiometry(s).unwrap();
        assert_eq!(terms.len(), 2);
        assert_eq!(terms[0].name, "A; B");
        assert_eq!(terms[1].cpd.as_str(), "cpd00002");
    }

    #[test]
    fn missing_name_is_ok() {
        let terms = parse_stoichiometry("-1:cpd00001:0:0").unwrap();
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].name, "");
    }

    #[test]
    fn fractional_coef() {
        let terms = parse_stoichiometry("-0.5:cpd00001:0:0:\"A\"").unwrap();
        assert_eq!(terms[0].coef, -0.5);
    }

    #[test]
    fn empty_is_error() {
        assert!(matches!(parse_stoichiometry("").unwrap_err(), StoichParseError::Empty));
        assert!(matches!(parse_stoichiometry("   ").unwrap_err(), StoichParseError::Empty));
        assert!(matches!(parse_stoichiometry(";;;").unwrap_err(), StoichParseError::Empty));
    }

    #[test]
    fn bad_term_reports_index() {
        let err = parse_stoichiometry("-1:cpd:0:0:X;bad").unwrap_err();
        assert!(matches!(err, StoichParseError::BadTerm { index: 1, .. }));
    }
}
