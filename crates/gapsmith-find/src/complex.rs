//! Port of `src/complex_detection.R`.
//!
//! Input: a vector of FASTA sequence descriptors (the full header text after
//! `>`) taken from one reaction's reference files, plus the reaction id (so
//! per-reaction subunit dictionaries apply).
//!
//! Output: a parallel vector of `Option<String>` giving each sequence's
//! canonical subunit assignment, or `None` when the extractor decided the
//! sequence is unlabeled.
//!
//! Algorithm (faithful to the R original, line numbers from
//! `src/complex_detection.R`):
//!
//! 1. Apply 8 alternating regex patterns (com.pat1..com.pat8) to extract a
//!    raw subunit phrase.
//! 2. Normalize "subunit/chain/polypeptide/component" → "Subunit".
//! 3. Reorder so the label comes after `Subunit` (`alpha Subunit` →
//!    `Subunit alpha`).
//! 4. Subunit-dict translation (`dat/complex_subunit_dict.tsv`): `Subunit
//!    <synonym>` → `Subunit <canonical>`.
//! 5. Numeral mapping: Latin I..XV → 1..15, single letters A..Z → 1..26,
//!    greek alpha..sigma → 1..17, small/medium/large → 1..3.
//! 6. Strip any trailing `[A-z]` on `Subunit N<letter>` (R drops
//!    sub-sub-complexes).
//! 7. Low-count filter: if mean count ≥ 10, drop subunits with count < 5.
//! 8. High-quality selection: if ≥ 66% of hits are numbered subunits,
//!    drop everything else.
//! 9. Global cutoff: if ≤ 20% of inputs produced any hit at all, blank the
//!    whole vector.

use gapsmith_db::ComplexSubunitTable;
use regex::Regex;
use std::collections::HashMap;
use std::sync::OnceLock;

/// Detect one subunit per sequence descriptor for a given reaction. A
/// `None` entry means "no subunit recognized".
pub fn detect_subunits(
    rxn_id: &str,
    descriptors: &[&str],
    dict: &ComplexSubunitTable,
) -> Vec<Option<String>> {
    if descriptors.is_empty() {
        return Vec::new();
    }

    // Step 1: regex extraction.
    let mut hits: Vec<Option<String>> = descriptors
        .iter()
        .map(|s| extract_subunit_phrase(s))
        .collect();

    // Steps 2-3: canonicalization.
    for h in hits.iter_mut().flatten() {
        *h = canonicalize_subunit(h);
    }

    // Step 4: per-reaction dictionary translation.
    if let Some(entries) = dict.for_rxn(rxn_id) {
        let lookup: HashMap<String, String> = entries
            .iter()
            .map(|(synonym, canonical)| {
                (format!("Subunit {synonym}"), format!("Subunit {canonical}"))
            })
            .collect();
        for h in hits.iter_mut().flatten() {
            if let Some(repl) = lookup.get(h.as_str()) {
                *h = repl.clone();
            }
        }
    }

    // Step 5: numeral mapping.
    for h in hits.iter_mut().flatten() {
        *h = apply_numeral_maps(h);
    }

    // Step 6: strip trailing letter on numeric subunits.
    let re_sub_num_letter = sub_num_letter_re();
    for h in hits.iter_mut().flatten() {
        if let Some(caps) = re_sub_num_letter.captures(h) {
            *h = caps.get(1).unwrap().as_str().to_string();
        }
    }

    // Step 7: low-count filter.
    apply_low_count_filter(&mut hits);

    // Step 8: high-quality numbered-subunit selection.
    apply_numbered_quality_filter(&mut hits);

    // Step 9: ≤ 20% coverage ⇒ blank everything.
    let n = hits.len();
    let any = hits.iter().filter(|h| h.is_some()).count();
    if any * 5 <= n {
        for h in &mut hits {
            *h = None;
        }
    }

    hits
}

/// Apply the union of the 8 patterns from `complex_detection.R:11-19` and
/// return the first match.
fn extract_subunit_phrase(descriptor: &str) -> Option<String> {
    // The R original uses case-sensitive matches — subunit words are always
    // lowercase in UniProt descriptors. We replicate that here.
    static COMPILED: OnceLock<Regex> = OnceLock::new();
    let re = COMPILED.get_or_init(|| {
        // synonymes = (subunit|chain|polypeptide|component)
        let syn = r"(?:subunit|chain|polypeptide|component)";
        let numeric = format!(r"{syn} [1-9]+(?:[A-Z])?\b");
        let caps_letters = format!(r"{syn} [A-Z]+\b");
        let pre_caps = format!(r"\b[A-Z]+ {syn}");
        let greek_pre = format!(
            r"(?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) {syn}"
        );
        let greek_post = format!(
            r"{syn} (?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma)"
        );
        let titled = format!(r"{syn} [A-Z][A-Za-z]+\b");
        let size_pre = format!(r"(?:large|medium|small) {syn}");
        let greek_dash = format!(
            r"(?:alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma)-{syn}"
        );
        // Order matters only marginally because `|` picks the leftmost; we
        // keep the same order as the R file.
        let combined = format!(
            r"({numeric}|{caps_letters}|{pre_caps}|{greek_pre}|{greek_post}|{titled}|{size_pre}|{greek_dash})"
        );
        Regex::new(&combined).expect("complex pattern should compile")
    });
    re.find(descriptor).map(|m| m.as_str().to_string())
}

/// Rewrite synonyms (subunit/chain/polypeptide/component) to the canonical
/// "Subunit" and reorder so the tag comes first.
fn canonicalize_subunit(raw: &str) -> String {
    static SYN: OnceLock<Regex> = OnceLock::new();
    let re_syn = SYN.get_or_init(|| {
        Regex::new(r"subunit|chain|polypeptide|component").expect("syn regex")
    });
    let mut out = re_syn.replace_all(raw, "Subunit").to_string();

    // `X Subunit` → `Subunit X` (for ALL-CAPS, greek, small/medium/large).
    static PRE_CAPS: OnceLock<Regex> = OnceLock::new();
    let re_cap = PRE_CAPS
        .get_or_init(|| Regex::new(r"([A-Z]+) Subunit").expect("pre-caps regex"));
    out = re_cap.replace_all(&out, "Subunit $1").to_string();

    static PRE_GREEK: OnceLock<Regex> = OnceLock::new();
    let re_greek = PRE_GREEK.get_or_init(|| {
        Regex::new(
            r"(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) Subunit",
        )
        .expect("greek pre-regex")
    });
    out = re_greek.replace_all(&out, "Subunit $1").to_string();

    static DASH_GREEK: OnceLock<Regex> = OnceLock::new();
    let re_dash = DASH_GREEK.get_or_init(|| {
        Regex::new(
            r"(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma)-Subunit",
        )
        .expect("dash greek regex")
    });
    out = re_dash.replace_all(&out, "Subunit $1").to_string();

    static PRE_SIZE: OnceLock<Regex> = OnceLock::new();
    let re_size =
        PRE_SIZE.get_or_init(|| Regex::new(r"(small|medium|large) Subunit").expect("size regex"));
    out = re_size.replace_all(&out, "Subunit $1").to_string();

    out
}

fn apply_numeral_maps(raw: &str) -> String {
    // Order matters: latin numerals first, then single letters, then greek,
    // then size words (largest-first to avoid substring shadowing — e.g.
    // `XIII` should match before `XII` / `XI`).

    // Latin I..XV → 1..15 (case-insensitive, word-boundary only on the
    // leading side).
    const LATIN: &[(&str, &str)] = &[
        ("XV", "15"),
        ("XIV", "14"),
        ("XIII", "13"),
        ("XII", "12"),
        ("XI", "11"),
        ("X", "10"),
        ("IX", "9"),
        ("VIII", "8"),
        ("VII", "7"),
        ("VI", "6"),
        ("V", "5"),
        ("IV", "4"),
        ("III", "3"),
        ("II", "2"),
        ("I", "1"),
    ];
    static LATIN_RES: OnceLock<Vec<(Regex, String)>> = OnceLock::new();
    let latin_rs = LATIN_RES.get_or_init(|| {
        LATIN
            .iter()
            .map(|(pat, rep)| {
                (Regex::new(&format!(r"(?i)\b{pat}")).expect("latin re"), (*rep).to_string())
            })
            .collect()
    });
    let mut out = raw.to_string();
    for (re, rep) in latin_rs {
        out = re.replace_all(&out, rep).to_string();
    }

    // A..Z → 1..26 (case-insensitive, full word).
    static LETTER_RES: OnceLock<Vec<(Regex, String)>> = OnceLock::new();
    let letter_rs = LETTER_RES.get_or_init(|| {
        (b'A'..=b'Z')
            .enumerate()
            .map(|(i, c)| {
                (
                    Regex::new(&format!(r"(?i)\b{}\b", c as char)).expect("letter re"),
                    (i + 1).to_string(),
                )
            })
            .collect()
    });
    for (re, rep) in letter_rs {
        out = re.replace_all(&out, rep).to_string();
    }

    // Greek alpha..sigma → 1..17 (case-insensitive, full word).
    const GREEK: &[&str] = &[
        "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
        "lambda", "my", "ny", "omikron", "pi", "rho", "sigma",
    ];
    static GREEK_RES: OnceLock<Vec<(Regex, String)>> = OnceLock::new();
    let greek_rs = GREEK_RES.get_or_init(|| {
        GREEK
            .iter()
            .enumerate()
            .map(|(i, g)| {
                (
                    Regex::new(&format!(r"(?i)\b{g}\b")).expect("greek re"),
                    (i + 1).to_string(),
                )
            })
            .collect()
    });
    for (re, rep) in greek_rs {
        out = re.replace_all(&out, rep).to_string();
    }

    // small/medium/large → 1/2/3. Note: R uses case-sensitive substring
    // replacement here (no word boundaries) — we mirror that.
    out = out.replace("small", "1").replace("medium", "2").replace("large", "3");
    out
}

fn sub_num_letter_re() -> &'static Regex {
    static RE: OnceLock<Regex> = OnceLock::new();
    RE.get_or_init(|| Regex::new(r"^(Subunit [0-9]+)[A-Za-z]$").expect("sub-num-letter re"))
}

fn apply_low_count_filter(hits: &mut [Option<String>]) {
    let mut tally: HashMap<&str, usize> = HashMap::new();
    for h in hits.iter().flatten() {
        *tally.entry(h.as_str()).or_insert(0) += 1;
    }
    if tally.is_empty() {
        return;
    }
    let sum: usize = tally.values().sum();
    let mean = sum as f64 / tally.len() as f64;
    if mean < 10.0 {
        return;
    }
    // Snapshot the low-count keys as owned strings to avoid borrow conflict.
    let low: Vec<String> =
        tally.iter().filter(|(_, c)| **c < 5).map(|(k, _)| (*k).to_string()).collect();
    for h in hits.iter_mut() {
        if let Some(inner) = h {
            if low.iter().any(|lo| lo == inner) {
                *h = None;
            }
        }
    }
}

fn apply_numbered_quality_filter(hits: &mut [Option<String>]) {
    static NUMBERED: OnceLock<Regex> = OnceLock::new();
    let re = NUMBERED.get_or_init(|| Regex::new(r"^Subunit [0-9]+").expect("numbered re"));

    let numbered_count = hits
        .iter()
        .flatten()
        .filter(|s| re.is_match(s))
        .count();
    if hits.is_empty() {
        return;
    }
    let cov = numbered_count as f64 / hits.len() as f64;
    if cov < 0.66 {
        return;
    }
    for h in hits.iter_mut() {
        if let Some(inner) = h {
            if !re.is_match(inner) {
                *h = None;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_dict() -> ComplexSubunitTable {
        ComplexSubunitTable::default()
    }

    #[test]
    fn extract_numeric_subunit() {
        let r = extract_subunit_phrase("DNA-directed RNA polymerase subunit 2");
        assert_eq!(r.as_deref(), Some("subunit 2"));
    }

    #[test]
    fn extract_greek_subunit() {
        let r = extract_subunit_phrase("ATP synthase alpha chain");
        assert_eq!(r.as_deref(), Some("alpha chain"));
    }

    #[test]
    fn extract_size_subunit() {
        let r = extract_subunit_phrase("Ribulose bisphosphate carboxylase large chain");
        assert_eq!(r.as_deref(), Some("large chain"));
    }

    #[test]
    fn no_match() {
        assert_eq!(extract_subunit_phrase("Acetaldehyde dehydrogenase"), None);
    }

    #[test]
    fn canonicalize_reorders() {
        assert_eq!(canonicalize_subunit("alpha chain"), "Subunit alpha");
        assert_eq!(canonicalize_subunit("small subunit"), "Subunit small");
        assert_eq!(canonicalize_subunit("ABC component"), "Subunit ABC");
    }

    #[test]
    fn numeral_mapping_greek() {
        assert_eq!(apply_numeral_maps("Subunit alpha"), "Subunit 1");
        assert_eq!(apply_numeral_maps("Subunit sigma"), "Subunit 17");
    }

    #[test]
    fn numeral_mapping_latin() {
        assert_eq!(apply_numeral_maps("Subunit XIII"), "Subunit 13");
        assert_eq!(apply_numeral_maps("Subunit IV"), "Subunit 4");
    }

    #[test]
    fn numeral_mapping_letters() {
        assert_eq!(apply_numeral_maps("Subunit A"), "Subunit 1");
        assert_eq!(apply_numeral_maps("Subunit C"), "Subunit 3");
    }

    #[test]
    fn numeral_mapping_size() {
        assert_eq!(apply_numeral_maps("Subunit small"), "Subunit 1");
        assert_eq!(apply_numeral_maps("Subunit large"), "Subunit 3");
    }

    #[test]
    fn full_pipeline_on_mixed_headers() {
        let dict = empty_dict();
        let descriptors = &[
            "DNA polymerase subunit alpha OS=...",
            "DNA polymerase subunit beta OS=...",
            "DNA polymerase subunit gamma OS=...",
            "Unrelated protein OS=...",
        ];
        let out = detect_subunits("RXN-TEST", descriptors, &dict);
        assert_eq!(out.len(), 4);
        // With only 4 inputs and 3 matched, we're well above the 20%
        // threshold — expect all three labeled, fourth None.
        assert_eq!(out[0].as_deref(), Some("Subunit 1"));
        assert_eq!(out[1].as_deref(), Some("Subunit 2"));
        assert_eq!(out[2].as_deref(), Some("Subunit 3"));
        assert_eq!(out[3], None);
    }

    #[test]
    fn coverage_under_20_percent_blanks_all() {
        let dict = empty_dict();
        // 1 in 10 has subunit info (10%). Should blank.
        let mut headers = vec!["random protein"; 9];
        headers.push("DNA pol subunit alpha");
        let out = detect_subunits("RXN", &headers, &dict);
        assert!(out.iter().all(|h| h.is_none()));
    }

    #[test]
    fn numbered_quality_filter_drops_non_numbered() {
        let dict = empty_dict();
        // 7 numbered subunits (1..7) + 1 greek + 1 random. 7/9 = 77% > 66%
        // ⇒ non-numbered should be filtered out by step 8.
        let hs = (1..=7)
            .map(|i| format!("DNA pol subunit {i}"))
            .collect::<Vec<_>>();
        let mut input: Vec<String> = hs.clone();
        input.push("Ribonuclease subunit alpha".into());
        input.push("random".into());
        let refs: Vec<&str> = input.iter().map(|s| s.as_str()).collect();
        let out = detect_subunits("RXN", &refs, &dict);
        // First 7 should be numbered.
        for (i, h) in out.iter().take(7).enumerate() {
            assert!(h.is_some(), "numbered input {} dropped: {h:?}", i + 1);
            let s = h.as_ref().unwrap();
            assert!(s.starts_with("Subunit "), "expected Subunit prefix: {s}");
        }
        // alpha becomes `Subunit 1` after numeral mapping, so it survives
        // the numbered-quality filter (matches R's behavior exactly).
        assert_eq!(out[7].as_deref(), Some("Subunit 1"));
        // "random" has no subunit token at all.
        assert_eq!(out[8], None);
    }
}
