//! TC (Transporter Classification) helpers.
//!
//! The 5-level TC identifier pattern is `\d\.[A-Z]\.\d+\.\d+\.\d+` —
//! e.g. `1.A.11.4.1`. `analyse_alignments_transport.R:70` extracts
//! the first match from the BLAST `qseqid` column and then indexes
//! the first digit into [`TC_TYPES`] to get a human label.

use regex::Regex;
use std::sync::OnceLock;

pub const TC_TYPES: &[&str] = &[
    "", // placeholder so we can index by TC number (1..=4) directly
    "1.Channels and pores",
    "2.Electrochemical potential-driven transporters",
    "3.Primary active transporters",
    "4.Group translocators",
];

/// Return the `type_index` (1–4) and the canonical type string for a
/// given TC identifier. Types 5/6/8/9 are intentionally filtered out in
/// the R original; we do the same by returning `None`.
pub fn type_of(tc: &str) -> Option<(u8, &'static str)> {
    let first = tc.chars().next()?;
    let n = first.to_digit(10)? as usize;
    if !(1..=4).contains(&n) {
        return None;
    }
    Some((n as u8, TC_TYPES[n]))
}

/// Pull the canonical TC id out of a BLAST `qseqid`. Returns `None` if
/// the query has no matching pattern (e.g. legacy headers).
pub fn extract_tc_id(qseqid: &str) -> Option<&str> {
    static RE: OnceLock<Regex> = OnceLock::new();
    let re = RE.get_or_init(|| {
        Regex::new(r"([1-4]\.[A-Za-z]\.[0-9]+\.[0-9]+\.[0-9]+)")
            .expect("TC regex compiles")
    });
    let m = re.find(qseqid)?;
    Some(m.as_str())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_tc_from_typical_headers() {
        assert_eq!(extract_tc_id("gnl|TC-DB|2.A.1.28.4 metadata"), Some("2.A.1.28.4"));
        assert_eq!(extract_tc_id("random 1.A.11.4.1 proto"), Some("1.A.11.4.1"));
        assert_eq!(extract_tc_id("no tc here"), None);
    }

    #[test]
    fn type_of_returns_name_and_index() {
        assert_eq!(type_of("1.A.11.4.1"), Some((1, "1.Channels and pores")));
        assert_eq!(type_of("4.A.1.1.1"), Some((4, "4.Group translocators")));
        assert!(type_of("5.A.1.1.1").is_none());
        assert!(type_of("bad").is_none());
    }
}
