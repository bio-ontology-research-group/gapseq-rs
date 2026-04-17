//! Tiny boolean-expression parser used by medium-prediction rules.
//!
//! The grammar (matches the forms gapseq writes):
//!
//! ```text
//!   expr    := compare
//!   compare := sum ( ('<' | '>' | '<=' | '>=' | '==') sum )?
//!   sum     := or ('+' or)*
//!   or      := and ('|' and)*
//!   and     := not ('&' not)*
//!   not     := '!' not | atom
//!   atom    := '(' expr ')'
//!            | IDENT
//!            | '"' QUOTED '"'
//!            | NUMBER
//!            | 'TRUE' | 'FALSE'
//! ```
//!
//! Values are numeric (TRUE = 1.0, FALSE = 0.0); any non-zero result at
//! the top level is treated as `true`. The `+` operator lets gapseq count
//! the number of TRUE sub-expressions (e.g. `rxn02213 + rxn01255 < 3`).
//!
//! Identifiers are bare words (typically `cpdNNNNN` or `rxnNNNNN`).
//! Double-quoted strings are pathway ids (MetaCyc-style, e.g.
//! `"N2FIX-PWY"`).
//!
//! Evaluation is driven by a closure that returns `true` / `false` for each
//! atom. The caller decides whether a pathway / reaction / compound is
//! "present".

use std::iter::Peekable;
use std::str::CharIndices;

#[derive(Debug, thiserror::Error)]
pub enum BoolExprError {
    #[error("unexpected token at position {pos}: {message}")]
    UnexpectedToken { pos: usize, message: String },
    #[error("trailing input after expression at position {pos}")]
    TrailingInput { pos: usize },
    #[error("unterminated quoted string starting at {pos}")]
    UnterminatedQuote { pos: usize },
}

/// Evaluate a boolean expression.
///
/// `resolve(tok)` receives the atom text (for quoted strings the quotes are
/// stripped; for bare words it's the verbatim identifier) and returns its
/// truth value.
pub fn eval<F>(src: &str, mut resolve: F) -> Result<bool, BoolExprError>
where
    F: FnMut(&str) -> bool,
{
    let mut p = Parser::new(src);
    let v = p.parse_compare(&mut resolve)?;
    p.skip_ws();
    if let Some((pos, _)) = p.iter.peek().copied() {
        return Err(BoolExprError::TrailingInput { pos });
    }
    Ok(v.abs() > 0.5)
}

struct Parser<'a> {
    src: &'a str,
    iter: Peekable<CharIndices<'a>>,
}

impl<'a> Parser<'a> {
    fn new(src: &'a str) -> Self {
        Self { src, iter: src.char_indices().peekable() }
    }

    fn skip_ws(&mut self) {
        while let Some(&(_, c)) = self.iter.peek() {
            if c.is_whitespace() {
                self.iter.next();
            } else {
                break;
            }
        }
    }

    /// Top-level comparison: `sum (op sum)?`. When no operator is
    /// present the sum itself is the result (non-zero → true).
    fn parse_compare<F>(&mut self, r: &mut F) -> Result<f64, BoolExprError>
    where
        F: FnMut(&str) -> bool,
    {
        let lhs = self.parse_sum(r)?;
        self.skip_ws();
        let op = match self.iter.peek().copied() {
            Some((_, '<')) | Some((_, '>')) | Some((_, '=')) => {
                let mut op = String::new();
                op.push(self.iter.next().unwrap().1);
                if let Some(&(_, c)) = self.iter.peek() {
                    if c == '=' {
                        op.push('=');
                        self.iter.next();
                    }
                }
                op
            }
            _ => return Ok(lhs),
        };
        let rhs = self.parse_sum(r)?;
        let v = match op.as_str() {
            "<" => lhs < rhs,
            "<=" => lhs <= rhs,
            ">" => lhs > rhs,
            ">=" => lhs >= rhs,
            "==" | "=" => (lhs - rhs).abs() < f64::EPSILON,
            _ => unreachable!(),
        };
        Ok(if v { 1.0 } else { 0.0 })
    }

    fn parse_sum<F>(&mut self, r: &mut F) -> Result<f64, BoolExprError>
    where
        F: FnMut(&str) -> bool,
    {
        let mut v = self.parse_or(r)?;
        loop {
            self.skip_ws();
            if self.iter.peek().map(|&(_, c)| c) != Some('+') {
                break;
            }
            self.iter.next();
            let rhs = self.parse_or(r)?;
            v += rhs;
        }
        Ok(v)
    }

    fn parse_or<F>(&mut self, r: &mut F) -> Result<f64, BoolExprError>
    where
        F: FnMut(&str) -> bool,
    {
        let mut v = self.parse_and(r)?;
        loop {
            self.skip_ws();
            if self.iter.peek().map(|&(_, c)| c) != Some('|') {
                break;
            }
            self.iter.next();
            let rhs = self.parse_and(r)?;
            v = if v.abs() > 0.5 || rhs.abs() > 0.5 { 1.0 } else { 0.0 };
        }
        Ok(v)
    }

    fn parse_and<F>(&mut self, r: &mut F) -> Result<f64, BoolExprError>
    where
        F: FnMut(&str) -> bool,
    {
        let mut v = self.parse_not(r)?;
        loop {
            self.skip_ws();
            if self.iter.peek().map(|&(_, c)| c) != Some('&') {
                break;
            }
            self.iter.next();
            let rhs = self.parse_not(r)?;
            v = if v.abs() > 0.5 && rhs.abs() > 0.5 { 1.0 } else { 0.0 };
        }
        Ok(v)
    }

    fn parse_not<F>(&mut self, r: &mut F) -> Result<f64, BoolExprError>
    where
        F: FnMut(&str) -> bool,
    {
        self.skip_ws();
        if self.iter.peek().map(|&(_, c)| c) == Some('!') {
            self.iter.next();
            let v = self.parse_not(r)?;
            return Ok(if v.abs() > 0.5 { 0.0 } else { 1.0 });
        }
        self.parse_atom(r)
    }

    fn parse_atom<F>(&mut self, r: &mut F) -> Result<f64, BoolExprError>
    where
        F: FnMut(&str) -> bool,
    {
        self.skip_ws();
        let (pos, ch) = match self.iter.peek().copied() {
            Some(p) => p,
            None => {
                return Err(BoolExprError::UnexpectedToken {
                    pos: self.src.len(),
                    message: "expected atom, got end of input".into(),
                });
            }
        };
        if ch == '(' {
            self.iter.next();
            let v = self.parse_compare(r)?;
            self.skip_ws();
            match self.iter.next() {
                Some((_, ')')) => Ok(v),
                Some((p, c)) => Err(BoolExprError::UnexpectedToken {
                    pos: p,
                    message: format!("expected ')', got {c:?}"),
                }),
                None => Err(BoolExprError::UnexpectedToken {
                    pos: self.src.len(),
                    message: "expected ')', got end of input".into(),
                }),
            }
        } else if ch == '"' {
            self.iter.next();
            let start = pos + 1;
            let end;
            loop {
                match self.iter.next() {
                    Some((p, '"')) => {
                        end = p;
                        break;
                    }
                    Some(_) => continue,
                    None => {
                        return Err(BoolExprError::UnterminatedQuote { pos });
                    }
                }
            }
            let tok = &self.src[start..end];
            Ok(if r(tok) { 1.0 } else { 0.0 })
        } else if ch.is_ascii_digit() {
            let start = pos;
            let mut end = pos;
            while let Some(&(p, c)) = self.iter.peek() {
                if c.is_ascii_digit() || c == '.' {
                    end = p + c.len_utf8();
                    self.iter.next();
                } else {
                    break;
                }
            }
            self.src[start..end]
                .parse::<f64>()
                .map_err(|_| BoolExprError::UnexpectedToken {
                    pos,
                    message: format!("bad number `{}`", &self.src[start..end]),
                })
        } else if ch.is_alphanumeric() || ch == '_' || ch == '-' {
            let start = pos;
            let mut end = pos;
            while let Some(&(p, c)) = self.iter.peek() {
                if c.is_alphanumeric() || c == '_' || c == '-' {
                    end = p + c.len_utf8();
                    self.iter.next();
                } else {
                    break;
                }
            }
            let tok = &self.src[start..end];
            match tok {
                "TRUE" | "True" | "true" => Ok(1.0),
                "FALSE" | "False" | "false" => Ok(0.0),
                other => Ok(if r(other) { 1.0 } else { 0.0 }),
            }
        } else {
            Err(BoolExprError::UnexpectedToken {
                pos,
                message: format!("unexpected character {ch:?}"),
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn literal_true() {
        assert!(eval("TRUE", |_| false).unwrap());
    }

    #[test]
    fn single_ident_lookup() {
        let set: HashSet<&str> = ["cpd00027"].into_iter().collect();
        assert!(eval("cpd00027", |t| set.contains(t)).unwrap());
        assert!(!eval("cpd99999", |t| set.contains(t)).unwrap());
    }

    #[test]
    fn quoted_pathway() {
        let set: HashSet<&str> = ["N2FIX-PWY"].into_iter().collect();
        assert!(eval(r#""N2FIX-PWY""#, |t| set.contains(t)).unwrap());
        assert!(!eval(r#""MISSING-PWY""#, |t| set.contains(t)).unwrap());
    }

    #[test]
    fn or_and_operators() {
        let set: HashSet<&str> = ["a", "c"].into_iter().collect();
        assert!(eval("a | b", |t| set.contains(t)).unwrap());
        assert!(!eval("b & c", |t| set.contains(t)).unwrap());
        assert!(eval("a & c", |t| set.contains(t)).unwrap());
    }

    #[test]
    fn not_operator() {
        let set: HashSet<&str> = ["a"].into_iter().collect();
        assert!(!eval("!a", |t| set.contains(t)).unwrap());
        assert!(eval("!b", |t| set.contains(t)).unwrap());
        assert!(eval(r#"!("N2FIX-PWY" | "PWY-7576")"#, |t| set.contains(t)).unwrap());
    }

    #[test]
    fn real_gapseq_rule() {
        // Rule for aerobic respiration — O2 gated on any of 4 pathways.
        let rule = r#""PWY-3781" | "PWY1YI0-3" | "PWY-7279" | "PWY1YI0-6""#;
        let set: HashSet<&str> = ["PWY-7279"].into_iter().collect();
        assert!(eval(rule, |t| set.contains(t)).unwrap());
        let nope: HashSet<&str> = HashSet::new();
        assert!(!eval(rule, |t| nope.contains(t)).unwrap());
    }

    #[test]
    fn rxn_mixed_rule() {
        // (rxn05747 | rxn09272) & rxn12494
        let set: HashSet<&str> = ["rxn05747", "rxn12494"].into_iter().collect();
        assert!(eval("(rxn05747 | rxn09272) & rxn12494", |t| set.contains(t)).unwrap());
        let only_one: HashSet<&str> = ["rxn05747"].into_iter().collect();
        assert!(!eval("(rxn05747 | rxn09272) & rxn12494", |t| only_one.contains(t)).unwrap());
    }

    #[test]
    fn counting_rule_with_less_than() {
        // This is the actual rxn02213 rule from the real gapseq table.
        let rule = "rxn02213 + (rxn01740 | rxn01741) + (rxn01739 | rxn38708) + rxn02476 + rxn01255 < 3";
        // 2 TRUE → 2 < 3 → rule fires.
        let two: HashSet<&str> = ["rxn02213", "rxn01740"].into_iter().collect();
        assert!(eval(rule, |t| two.contains(t)).unwrap());
        // 3 TRUE → 3 < 3 → false.
        let three: HashSet<&str> = ["rxn02213", "rxn01740", "rxn02476"].into_iter().collect();
        assert!(!eval(rule, |t| three.contains(t)).unwrap());
        // 4 TRUE → false.
        let four: HashSet<&str> =
            ["rxn02213", "rxn01740", "rxn01739", "rxn02476"].into_iter().collect();
        assert!(!eval(rule, |t| four.contains(t)).unwrap());
    }
}
