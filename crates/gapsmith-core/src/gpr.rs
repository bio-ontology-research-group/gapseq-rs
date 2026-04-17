//! Gene-protein-reaction (GPR) Boolean expressions.
//!
//! Syntax accepted (matching cobrar / COBRApy conventions):
//!
//! - Leaves: bare identifiers `[A-Za-z_][A-Za-z0-9_\.\-]*` or quoted strings.
//! - `and` (case-insensitive) or `&` — logical AND.
//! - `or` (case-insensitive) or `|` — logical OR.
//! - Parentheses for grouping.
//!
//! AND binds tighter than OR (textbook precedence). No NOT operator.
//!
//! A leaf is a [`GeneId`] reference; compound complexes become `And` nodes.
//! Isoforms become `Or` nodes.

use crate::GeneId;
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "lowercase")]
pub enum Gpr {
    Gene { id: GeneId },
    And { operands: Vec<Gpr> },
    Or { operands: Vec<Gpr> },
}

impl Gpr {
    pub fn gene(id: impl Into<GeneId>) -> Self {
        Gpr::Gene { id: id.into() }
    }

    /// Flatten nested same-operator nodes and drop unary/empty collections.
    pub fn normalize(self) -> Self {
        match self {
            g @ Gpr::Gene { .. } => g,
            Gpr::And { operands } => {
                let mut out = Vec::with_capacity(operands.len());
                for op in operands {
                    match op.normalize() {
                        Gpr::And { operands: inner } => out.extend(inner),
                        other => out.push(other),
                    }
                }
                match out.len() {
                    0 => Gpr::And { operands: vec![] },
                    1 => out.into_iter().next().unwrap(),
                    _ => Gpr::And { operands: out },
                }
            }
            Gpr::Or { operands } => {
                let mut out = Vec::with_capacity(operands.len());
                for op in operands {
                    match op.normalize() {
                        Gpr::Or { operands: inner } => out.extend(inner),
                        other => out.push(other),
                    }
                }
                match out.len() {
                    0 => Gpr::Or { operands: vec![] },
                    1 => out.into_iter().next().unwrap(),
                    _ => Gpr::Or { operands: out },
                }
            }
        }
    }

    /// Collect all distinct gene ids referenced anywhere in the tree.
    pub fn collect_genes(&self, out: &mut Vec<GeneId>) {
        match self {
            Gpr::Gene { id } => {
                if !out.contains(id) {
                    out.push(id.clone());
                }
            }
            Gpr::And { operands } | Gpr::Or { operands } => {
                for op in operands {
                    op.collect_genes(out);
                }
            }
        }
    }
}

impl fmt::Display for Gpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Gpr::Gene { id } => f.write_str(id.as_str()),
            Gpr::And { operands } => {
                for (i, op) in operands.iter().enumerate() {
                    if i > 0 {
                        f.write_str(" and ")?;
                    }
                    write_grouped(f, op, matches!(op, Gpr::Or { .. }))?;
                }
                Ok(())
            }
            Gpr::Or { operands } => {
                for (i, op) in operands.iter().enumerate() {
                    if i > 0 {
                        f.write_str(" or ")?;
                    }
                    write_grouped(f, op, false)?;
                }
                Ok(())
            }
        }
    }
}

fn write_grouped(f: &mut fmt::Formatter<'_>, g: &Gpr, paren: bool) -> fmt::Result {
    if paren {
        write!(f, "({})", g)
    } else {
        write!(f, "{}", g)
    }
}

#[derive(Debug, thiserror::Error)]
pub enum GprParseError {
    #[error("unexpected character '{ch}' at position {pos}")]
    UnexpectedChar { ch: char, pos: usize },
    #[error("unexpected end of input at position {pos}")]
    UnexpectedEnd { pos: usize },
    #[error("unclosed parenthesis opened at position {pos}")]
    UnclosedParen { pos: usize },
    #[error("empty expression")]
    Empty,
}

// -- Parser --------------------------------------------------------------
// Precedence: OR < AND. Right-associative AND/OR are fine because of flattening.
// Grammar:
//   expr   := and_expr ('or' and_expr)*
//   and_expr := atom ('and' atom)*
//   atom   := '(' expr ')' | IDENT | QUOTED
//
// Tokens: `(`, `)`, AND, OR, IDENT, QUOTED_IDENT.

impl std::str::FromStr for Gpr {
    type Err = GprParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let tokens = tokenize(s)?;
        if tokens.is_empty() {
            return Err(GprParseError::Empty);
        }
        let mut parser = Parser { tokens: &tokens, pos: 0 };
        let expr = parser.parse_or()?;
        if parser.pos != tokens.len() {
            let (bad_pos, _) = tokens[parser.pos];
            return Err(GprParseError::UnexpectedEnd { pos: bad_pos });
        }
        Ok(expr.normalize())
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum Tok {
    LParen,
    RParen,
    And,
    Or,
    Ident(String),
}

fn tokenize(s: &str) -> Result<Vec<(usize, Tok)>, GprParseError> {
    let mut out = Vec::new();
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i] as char;
        match c {
            ' ' | '\t' | '\n' | '\r' => {
                i += 1;
            }
            '(' => {
                out.push((i, Tok::LParen));
                i += 1;
            }
            ')' => {
                out.push((i, Tok::RParen));
                i += 1;
            }
            '&' => {
                out.push((i, Tok::And));
                i += 1;
            }
            '|' => {
                out.push((i, Tok::Or));
                i += 1;
            }
            '\'' | '"' => {
                let quote = c;
                let start = i + 1;
                i += 1;
                while i < bytes.len() && (bytes[i] as char) != quote {
                    i += 1;
                }
                if i >= bytes.len() {
                    return Err(GprParseError::UnexpectedEnd { pos: start });
                }
                let id: String = s[start..i].to_string();
                out.push((start, Tok::Ident(id)));
                i += 1; // closing quote
            }
            ch if is_ident_start(ch) => {
                let start = i;
                while i < bytes.len() && is_ident_cont(bytes[i] as char) {
                    i += 1;
                }
                let slice = &s[start..i];
                let lower = slice.to_ascii_lowercase();
                let tok = match lower.as_str() {
                    "and" => Tok::And,
                    "or" => Tok::Or,
                    _ => Tok::Ident(slice.to_string()),
                };
                out.push((start, tok));
            }
            other => return Err(GprParseError::UnexpectedChar { ch: other, pos: i }),
        }
    }
    Ok(out)
}

fn is_ident_start(c: char) -> bool {
    c.is_ascii_alphabetic() || c == '_'
}

fn is_ident_cont(c: char) -> bool {
    c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-' | ':' | '+')
}

struct Parser<'a> {
    tokens: &'a [(usize, Tok)],
    pos: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&'a Tok> {
        self.tokens.get(self.pos).map(|(_, t)| t)
    }
    fn peek_pos(&self) -> usize {
        self.tokens
            .get(self.pos)
            .map(|(p, _)| *p)
            .unwrap_or_else(|| self.tokens.last().map(|(p, _)| *p + 1).unwrap_or(0))
    }
    fn bump(&mut self) -> &'a Tok {
        let t = &self.tokens[self.pos].1;
        self.pos += 1;
        t
    }

    fn parse_or(&mut self) -> Result<Gpr, GprParseError> {
        let mut left = self.parse_and()?;
        while matches!(self.peek(), Some(Tok::Or)) {
            self.bump();
            let right = self.parse_and()?;
            left = match left {
                Gpr::Or { mut operands } => {
                    operands.push(right);
                    Gpr::Or { operands }
                }
                other => Gpr::Or { operands: vec![other, right] },
            };
        }
        Ok(left)
    }

    fn parse_and(&mut self) -> Result<Gpr, GprParseError> {
        let mut left = self.parse_atom()?;
        while matches!(self.peek(), Some(Tok::And)) {
            self.bump();
            let right = self.parse_atom()?;
            left = match left {
                Gpr::And { mut operands } => {
                    operands.push(right);
                    Gpr::And { operands }
                }
                other => Gpr::And { operands: vec![other, right] },
            };
        }
        Ok(left)
    }

    fn parse_atom(&mut self) -> Result<Gpr, GprParseError> {
        match self.peek() {
            Some(Tok::LParen) => {
                let open_pos = self.peek_pos();
                self.bump();
                let expr = self.parse_or()?;
                match self.peek() {
                    Some(Tok::RParen) => {
                        self.bump();
                        Ok(expr)
                    }
                    _ => Err(GprParseError::UnclosedParen { pos: open_pos }),
                }
            }
            Some(Tok::Ident(s)) => {
                let id = s.clone();
                self.bump();
                Ok(Gpr::gene(id))
            }
            Some(Tok::And) | Some(Tok::Or) | Some(Tok::RParen) => {
                Err(GprParseError::UnexpectedChar { ch: '?', pos: self.peek_pos() })
            }
            None => Err(GprParseError::UnexpectedEnd { pos: self.peek_pos() }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn single_gene() {
        let g: Gpr = "b0001".parse().unwrap();
        assert_eq!(g, Gpr::gene("b0001"));
    }

    #[test]
    fn precedence_and_binds_tighter() {
        let g: Gpr = "a and b or c and d".parse().unwrap();
        assert_eq!(
            g,
            Gpr::Or {
                operands: vec![
                    Gpr::And { operands: vec![Gpr::gene("a"), Gpr::gene("b")] },
                    Gpr::And { operands: vec![Gpr::gene("c"), Gpr::gene("d")] },
                ]
            }
        );
    }

    #[test]
    fn parens_override_precedence() {
        let g: Gpr = "(a or b) and (c or d)".parse().unwrap();
        assert_eq!(
            g,
            Gpr::And {
                operands: vec![
                    Gpr::Or { operands: vec![Gpr::gene("a"), Gpr::gene("b")] },
                    Gpr::Or { operands: vec![Gpr::gene("c"), Gpr::gene("d")] },
                ]
            }
        );
    }

    #[test]
    fn symbolic_operators() {
        let g: Gpr = "a & (b | c)".parse().unwrap();
        assert_eq!(
            g,
            Gpr::And {
                operands: vec![
                    Gpr::gene("a"),
                    Gpr::Or { operands: vec![Gpr::gene("b"), Gpr::gene("c")] },
                ]
            }
        );
    }

    #[test]
    fn collect_genes_unique() {
        let g: Gpr = "a and (b or a or c)".parse().unwrap();
        let mut v = vec![];
        g.collect_genes(&mut v);
        assert_eq!(
            v.iter().map(|g| g.as_str().to_string()).collect::<Vec<_>>(),
            vec!["a", "b", "c"]
        );
    }

    #[test]
    fn display_roundtrip() {
        let g: Gpr = "a and (b or c)".parse().unwrap();
        let s = g.to_string();
        let g2: Gpr = s.parse().unwrap();
        assert_eq!(g, g2);
    }

    #[test]
    fn empty_is_error() {
        assert!(matches!(Gpr::from_str("   ").unwrap_err(), GprParseError::Empty));
    }

    #[test]
    fn unclosed_paren_is_error() {
        assert!(matches!(
            Gpr::from_str("a and (b or c").unwrap_err(),
            GprParseError::UnclosedParen { .. }
        ));
    }
}
