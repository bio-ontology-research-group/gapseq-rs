//! Build a substrate-filtered `small.fasta` from the TCDB + SEED
//! transporter reference set. Mirrors the shell block at
//! `src/transporter.sh:222-238`:
//!
//! ```text
//! cat $tcdb $otherDB > all.fasta
//! grep ">" all.fasta | sort > fasta_header
//! ...
//! grep -Fiwf SUBkey tcdb_all | cut -f1 > TCkey   # TC ids with a matching substrate
//! grep -Fivf TCkey fasta_header > fasta_header.noTCkey
//! grep -Fivf SUBkey fasta_header.noTCkey > fasta_header.noKey
//! comm -23 fasta_header fasta_header.noKey > fasta_header.small
//! awk 'BEGIN{while((getline<"fasta_header.small")>0)l[$1]=1}/^>/{f=l[$1]}f' all.fasta > small.fasta
//! ```
//!
//! In Rust: walk the union of reference FASTAs, decide for each header
//! whether to keep the sequence, and write the survivors to a single
//! output file.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::data::SeedTransporterRow;

pub struct BuildSmallResult {
    /// Full path of the emitted `small.fasta`.
    pub small_fasta: std::path::PathBuf,
    /// Full header lines (starting with `>`) of every kept sequence.
    pub fasta_header_small: Vec<String>,
    /// Substrate keys with original casing preserved (from `subex.tbl`).
    /// Used as-is when stamping the `sub` column of the output table, so
    /// downstream tooling sees the same casing gapseq does.
    pub sub_key_orig: Vec<String>,
    /// Lower-case view of `sub_key_orig` used for case-insensitive match.
    pub sub_key: HashSet<String>,
    /// `(substrate_lc, exchange_id)` tuples from subex.
    pub subid: Vec<(String, String)>,
}

/// Write `small.fasta` and report the survivor headers.
///
/// `reference_fastas`: concatenation sources — typically
/// `[dat/seq/tcdb.fasta, dat/seq/transporter.fasta]`.
/// `subex_rows`: rows of `dat/subex.tbl`. Only rows with a non-empty
/// `seed` column contribute to the filter (matches the `awk '$4 != ""'`
/// guard in the shell script).
/// `tcdb_all`: TC id → cleaned substrate string (from
/// [`crate::data::load_tcdb_all`]).
/// `only_met`: optional substrate keyword restricting the filter.
pub fn build_small_fasta(
    reference_fastas: &[&Path],
    subex_rows: &[gapsmith_db::SubexRow],
    tcdb_all: &HashMap<String, String>,
    only_met: Option<&str>,
    out_dir: &Path,
) -> std::io::Result<BuildSmallResult> {
    // Step 1: build SUBid + SUBkey from subex.
    // SUBid rows: `(substrate_name_lc, exchange_id)`. A subex row's
    // `name` column can list multiple synonyms separated by `;` — we
    // emit one SUBid row per synonym (`for(i=0;++i<=n;){print a[i]"\t"$4}`).
    let only_met_lc = only_met.map(|s| s.to_ascii_lowercase());
    let mut subid: Vec<(String, String)> = Vec::new();
    let mut sub_key: HashSet<String> = HashSet::new();
    let mut sub_key_orig_set: HashSet<String> = HashSet::new();
    let mut sub_key_orig: Vec<String> = Vec::new();
    for row in subex_rows {
        let seed_id = row.seed.trim();
        if seed_id.is_empty() {
            continue;
        }
        if let Some(ref m) = only_met_lc {
            if !word_contains_ignore_case(&row.name, m) {
                continue;
            }
        }
        for syn in row.name.split(';').map(str::trim).filter(|s| !s.is_empty()) {
            let syn_lc = syn.to_ascii_lowercase();
            subid.push((syn_lc.clone(), seed_id.to_string()));
            sub_key.insert(syn_lc);
            if sub_key_orig_set.insert(syn.to_string()) {
                sub_key_orig.push(syn.to_string());
            }
        }
    }

    // Step 2: build the set of TC ids whose substrate list matches at
    // least one SUBkey (mirrors `grep -Fiwf SUBkey tcdb_all`).
    let tc_keep: HashSet<&str> = tcdb_all
        .iter()
        .filter(|(_, subs)| sub_key_matches(&sub_key, subs))
        .map(|(tc, _)| tc.as_str())
        .collect();

    // Step 3: walk each reference fasta, keep sequences whose header
    // either mentions a TC id in `tc_keep` OR contains a SUBkey word.
    std::fs::create_dir_all(out_dir)?;
    let small_path = out_dir.join("small.fasta");
    let out_file = File::create(&small_path)?;
    let mut w = BufWriter::new(out_file);
    let mut headers_out: Vec<String> = Vec::new();

    let mut current_keep: bool = false;
    for path in reference_fastas {
        let f = File::open(path)?;
        let r = BufReader::new(f);
        for line in r.lines() {
            let line = line?;
            if let Some(rest) = line.strip_prefix('>') {
                let corrected = correct_header_quirks(rest);
                current_keep = header_matches(&corrected, &tc_keep, &sub_key);
                if current_keep {
                    headers_out.push(format!(">{corrected}"));
                    w.write_all(b">")?;
                    w.write_all(corrected.as_bytes())?;
                    w.write_all(b"\n")?;
                }
            } else if current_keep {
                w.write_all(line.as_bytes())?;
                w.write_all(b"\n")?;
            }
        }
    }
    w.flush()?;

    tracing::info!(
        small_fasta = %small_path.display(),
        headers = headers_out.len(),
        sub_keys = sub_key.len(),
        tc_keep = tc_keep.len(),
        "small.fasta built"
    );

    let _ = seed_rows_sentinel();
    Ok(BuildSmallResult {
        small_fasta: small_path,
        fasta_header_small: headers_out,
        sub_key_orig,
        sub_key,
        subid,
    })
}

fn seed_rows_sentinel() -> Option<Vec<SeedTransporterRow>> {
    None
}

/// Apply the two header-text fixups done in `transporter.sh:236-237`
/// (`amino butyrate` → `aminobutyrate`, `GABA butyrate` → `GABA`). We
/// apply them in-memory so downstream regex/word matches see the
/// corrected form.
fn correct_header_quirks(s: &str) -> String {
    s.replace("amino butyrate", "aminobutyrate").replace("GABA butyrate", "GABA")
}

/// `grep -Fiwf SUBkey tcdb_all` — substring search, word-boundary,
/// case-insensitive, against the substrate list column.
fn sub_key_matches(keys: &HashSet<String>, subs: &str) -> bool {
    let subs_lc = subs.to_ascii_lowercase();
    keys.iter().any(|k| word_contains(&subs_lc, k))
}

/// Test whether a header (original casing) should be kept. Mirrors
/// gapseq's (buggy-but-real) shell semantics:
///
/// ```text
/// grep -Fivf TCkey  fasta_header > noTCkey   # substring match, no -w
/// grep -Fivf SUBkey noTCkey      > noKey
/// comm -23 fasta_header noKey    > fasta_header.small
/// ```
///
/// `small` = headers that matched any TCkey OR SUBkey as a *substring*
/// (case-insensitive). Notice: no word boundaries anywhere — so
/// `4.A.2.1.1` accidentally matches a line containing `4.A.2.1.11`.
/// We preserve that behavior so the filter's output matches gapseq's
/// byte-for-byte.
fn header_matches(
    header: &str,
    tc_keep: &HashSet<&str>,
    sub_key: &HashSet<String>,
) -> bool {
    let hdr_lc = header.to_ascii_lowercase();
    // Substring scan of every TCkey entry (already lowercased by caller).
    for tc in tc_keep {
        if hdr_lc.contains(&tc.to_ascii_lowercase()) {
            return true;
        }
    }
    for k in sub_key {
        if hdr_lc.contains(k) {
            return true;
        }
    }
    false
}

/// Case-sensitive word-boundary substring containment. Mirrors
/// `grep -Fw` / `grep -w`: the needle must be bounded on each side by
/// non-word characters (or string ends).
fn word_contains(haystack: &str, needle: &str) -> bool {
    if needle.is_empty() {
        return false;
    }
    let bytes = haystack.as_bytes();
    let nbytes = needle.as_bytes();
    let mut i = 0;
    while i + nbytes.len() <= bytes.len() {
        if &bytes[i..i + nbytes.len()] == nbytes {
            let before_ok = i == 0 || !is_word_byte(bytes[i - 1]);
            let after = i + nbytes.len();
            let after_ok = after == bytes.len() || !is_word_byte(bytes[after]);
            if before_ok && after_ok {
                return true;
            }
        }
        i += 1;
    }
    false
}

fn word_contains_ignore_case(haystack: &str, needle: &str) -> bool {
    word_contains(&haystack.to_ascii_lowercase(), &needle.to_ascii_lowercase())
}

fn is_word_byte(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn word_contains_respects_boundaries() {
        assert!(word_contains("the arginine transporter", "arginine"));
        assert!(!word_contains("linearginine", "arginine"));
        assert!(!word_contains("arginine2", "arginine"));
        assert!(word_contains("arginine", "arginine"));
    }

    #[test]
    fn header_quirks_applied() {
        assert_eq!(correct_header_quirks("amino butyrate xyz"), "aminobutyrate xyz");
        assert_eq!(correct_header_quirks("GABA butyrate"), "GABA");
    }
}
