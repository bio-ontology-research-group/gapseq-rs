//! Reference-FASTA resolver.
//!
//! Given a `(rxn_id, ec, name)` tuple, walk `dat/seq/<taxonomy>/` in the
//! priority order specified by `prepare_batch_alignments.R:150-234, 315-366`:
//!
//! 1. `user/<file>.fasta` — user-supplied overrides win.
//! 2. `rxn/<RXN>.fasta` — direct per-reaction links (MetaCyc style).
//! 3. `rev/<EC>.fasta` — UniProt reviewed clustered at 0.9.
//! 4. `unrev/<EC>.fasta` — UniProt unreviewed clustered at 0.5.
//! 5. `rev/<MD5(name)>.fasta` / `unrev/<MD5(name)>.fasta` — reaction-name
//!    fallback. The MD5 of the bare reaction-name UTF-8 bytes, no trailing
//!    newline, lower-case 32-char hex. Matches `uniprot.sh:155`:
//!    `reaNameHash=$(echo -n "$rea" | md5sum | awk '{print $1}')`.
//!
//! The reaction-name fallback is SKIPPED when the name looks like a
//! synthetic rxn id (`rxn00001`, `Rxn00001`, `RXN-...`) because those
//! are not real enzyme names and gapseq's R code explicitly filters them
//! (`prepare_batch_alignments.R:209`).
//!
//! `SeqfileOptions::seq_src` controls the reviewed-vs-unreviewed policy:
//!
//! - `SeqSrc::Reviewed`    → reviewed only
//! - `SeqSrc::PreferRev`   → reviewed if present, else unreviewed
//! - `SeqSrc::Both`        → reviewed and unreviewed both, concatenated
//! - `SeqSrc::Unreviewed`  → unreviewed only
//!
//! Empty files (0 bytes) are treated as absent — the R pipeline does the
//! same.

use md5;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum SeqSrc {
    Reviewed,
    #[default]
    PreferRev,
    Both,
    Unreviewed,
}

#[derive(Debug, Clone)]
pub struct SeqfileOptions {
    /// Root of the per-taxonomy sequence tree, e.g. `<seq-dir>/Bacteria`.
    /// Used for `rev/`, `unrev/`, `rxn/` resolution.
    pub tax_root: PathBuf,
    /// Separate root for `user/` files — gapseq's R code hard-codes this
    /// to the built-in `<gapseq>/dat/seq/<tax>/` regardless of `-D`
    /// (`prepare_batch_alignments.R:217`). Typically the same as
    /// `tax_root` unless the caller passed an override seqdb.
    pub user_root: PathBuf,
    pub seq_src: SeqSrc,
}

#[derive(Debug, Clone)]
pub struct ResolvedSeq {
    /// Path to the reference FASTA.
    pub path: PathBuf,
    /// Short relative label used as the `file` field in reaction rows
    /// (e.g. `rxn/FOO-RXN.fasta` or `rev/1.1.1.1.fasta`). Matches what
    /// gapseq stamps in the concatenated `query.faa` headers.
    pub label: String,
}

/// Resolve every applicable reference fasta for one reaction. Empty means
/// the reaction has `no_seq_data`.
///
/// Priority order (mirrors `prepare_batch_alignments.R:150-234`):
///
/// 1. `user/<EC>.fasta`, `user/<MD5 name>.fasta`, `user/<rxn>.fasta` —
///    all probed against `opts.user_root` (the built-in gapseq data dir).
/// 2. `rxn/<RXN>.fasta` under `opts.tax_root`.
/// 3. `rev/<EC>.fasta` / `unrev/<EC>.fasta` under `opts.tax_root`.
/// 4. `rev/<MD5 name>.fasta` / `unrev/<MD5 name>.fasta` under `opts.tax_root`.
pub fn resolve_for_reaction(
    opts: &SeqfileOptions,
    rxn_id: &str,
    ec: &str,
    name: &str,
) -> Vec<ResolvedSeq> {
    let mut out: Vec<ResolvedSeq> = Vec::new();

    // 1a. user/<EC>.fasta — one probe per top-level EC (handles `/`-sep).
    for one_ec in ec.split('/').map(str::trim).filter(|s| !s.is_empty()) {
        if let Some(r) = probe(&opts.user_root, "user", &format!("{one_ec}.fasta")) {
            out.push(r);
        }
    }
    // 1b. user/<MD5(name)>.fasta if the name is a real enzyme label.
    let name_trim = name.trim();
    if !name_trim.is_empty() && !looks_like_rxn_id(name_trim) {
        let hash = md5_hex(name_trim);
        if let Some(r) = probe(&opts.user_root, "user", &format!("{hash}.fasta")) {
            out.push(r);
        }
    }
    // 1c. user/<rxn>.fasta.
    if !rxn_id.is_empty() {
        if let Some(r) = probe(&opts.user_root, "user", &format!("{rxn_id}.fasta")) {
            out.push(r);
        }
    }
    if !out.is_empty() {
        return out;
    }

    // 2. Direct MetaCyc-style link.
    if let Some(r) = probe(&opts.tax_root, "rxn", &format!("{rxn_id}.fasta")) {
        return vec![r];
    }

    // 3. EC-based resolution (reviewed / unreviewed per seq_src).
    let ec_trim = ec.trim();
    if !ec_trim.is_empty() {
        let mut ec_paths = Vec::new();
        for one_ec in ec_trim.split('/').map(str::trim).filter(|s| !s.is_empty()) {
            ec_paths.extend(resolve_by_stem(opts, one_ec));
        }
        if !ec_paths.is_empty() {
            return ec_paths;
        }
    }

    // 4. Reaction-name fallback via MD5 hash of the name bytes.
    if !name_trim.is_empty() && !looks_like_rxn_id(name_trim) {
        let hash = md5_hex(name_trim);
        let hits = resolve_by_stem(opts, &hash);
        if !hits.is_empty() {
            return hits;
        }
    }

    Vec::new()
}

/// `true` when a string looks like a synthetic reaction id rather than a
/// human enzyme name. Mirrors the filter in
/// `src/prepare_batch_alignments.R:209`:
/// `^rxn[0-9]+$ | ^Rxn[0-9]+$ | ^RXN`.
pub fn looks_like_rxn_id(name: &str) -> bool {
    if let Some(rest) = name.strip_prefix("rxn") {
        return !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit());
    }
    if let Some(rest) = name.strip_prefix("Rxn") {
        return !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit());
    }
    name.starts_with("RXN")
}

/// MD5 hex digest of the UTF-8 bytes of `s` (no trailing newline).
/// Matches `echo -n "$s" | md5sum`.
pub fn md5_hex(s: &str) -> String {
    use md5::{Digest, Md5};
    let mut h = Md5::new();
    h.update(s.as_bytes());
    let out = h.finalize();
    let mut hex = String::with_capacity(32);
    for b in out.iter() {
        hex.push_str(&format!("{b:02x}"));
    }
    hex
}

fn resolve_by_stem(opts: &SeqfileOptions, ec: &str) -> Vec<ResolvedSeq> {
    let f = format!("{ec}.fasta");
    let rev = probe(&opts.tax_root, "rev", &f);
    let unrev = probe(&opts.tax_root, "unrev", &f);
    match opts.seq_src {
        SeqSrc::Reviewed => rev.into_iter().collect(),
        SeqSrc::PreferRev => rev.or(unrev).into_iter().collect(),
        SeqSrc::Both => vec![rev, unrev].into_iter().flatten().collect(),
        SeqSrc::Unreviewed => unrev.into_iter().collect(),
    }
}

fn probe(tax_root: &Path, dir: &str, filename: &str) -> Option<ResolvedSeq> {
    let p = tax_root.join(dir).join(filename);
    if !p.is_file() {
        return None;
    }
    let meta = std::fs::metadata(&p).ok()?;
    if meta.len() == 0 {
        return None;
    }
    Some(ResolvedSeq { path: p, label: format!("{dir}/{filename}") })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn prefers_rxn_dir() {
        let d = tempfile::tempdir().unwrap();
        let root = d.path();
        for sub in ["rxn", "rev", "unrev", "user"] {
            fs::create_dir_all(root.join(sub)).unwrap();
        }
        fs::write(root.join("rxn/RXN-1.fasta"), ">a\nA\n").unwrap();
        fs::write(root.join("rev/1.1.1.1.fasta"), ">b\nB\n").unwrap();

        let opts = SeqfileOptions {
            tax_root: root.into(),
            user_root: root.into(),
            seq_src: SeqSrc::PreferRev,
        };
        let r = resolve_for_reaction(&opts, "RXN-1", "1.1.1.1", "");
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].label, "rxn/RXN-1.fasta");
    }

    #[test]
    fn falls_back_to_reviewed_ec() {
        let d = tempfile::tempdir().unwrap();
        let root = d.path();
        fs::create_dir_all(root.join("rev")).unwrap();
        fs::write(root.join("rev/1.1.1.1.fasta"), ">b\nB\n").unwrap();
        let opts = SeqfileOptions {
            tax_root: root.into(),
            user_root: root.into(),
            seq_src: SeqSrc::PreferRev,
        };
        let r = resolve_for_reaction(&opts, "RXN-missing", "1.1.1.1", "");
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].label, "rev/1.1.1.1.fasta");
    }

    #[test]
    fn both_mode_returns_rev_and_unrev() {
        let d = tempfile::tempdir().unwrap();
        let root = d.path();
        for sub in ["rev", "unrev"] {
            fs::create_dir_all(root.join(sub)).unwrap();
        }
        fs::write(root.join("rev/1.1.1.1.fasta"), ">r\nR\n").unwrap();
        fs::write(root.join("unrev/1.1.1.1.fasta"), ">u\nU\n").unwrap();
        let opts = SeqfileOptions {
            tax_root: root.into(),
            user_root: root.into(),
            seq_src: SeqSrc::Both,
        };
        let r = resolve_for_reaction(&opts, "RXN-missing", "1.1.1.1", "");
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn user_override_wins() {
        let d = tempfile::tempdir().unwrap();
        let root = d.path();
        for sub in ["rxn", "rev", "user"] {
            fs::create_dir_all(root.join(sub)).unwrap();
        }
        fs::write(root.join("rxn/RXN-1.fasta"), ">a\n").unwrap();
        fs::write(root.join("rev/1.1.1.1.fasta"), ">b\n").unwrap();
        fs::write(root.join("user/RXN-1.fasta"), ">c\nC\n").unwrap();
        let opts = SeqfileOptions {
            tax_root: root.into(),
            user_root: root.into(),
            seq_src: SeqSrc::PreferRev,
        };
        let r = resolve_for_reaction(&opts, "RXN-1", "1.1.1.1", "");
        assert_eq!(r[0].label, "user/RXN-1.fasta");
    }

    #[test]
    fn md5_hex_matches_gnu_md5sum() {
        // Verified with `echo -n "Acetaldehyde dehydrogenase" | md5sum`.
        assert_eq!(md5_hex("Acetaldehyde dehydrogenase"), "1390704749ddc17f2b61599cf204ac4a");
        // Empty input — well-known.
        assert_eq!(md5_hex(""), "d41d8cd98f00b204e9800998ecf8427e");
    }

    #[test]
    fn looks_like_rxn_id_filter() {
        assert!(looks_like_rxn_id("rxn00001"));
        assert!(looks_like_rxn_id("Rxn00001"));
        assert!(looks_like_rxn_id("RXN-8099"));
        assert!(!looks_like_rxn_id("Acetaldehyde dehydrogenase"));
        assert!(!looks_like_rxn_id("ATP synthase"));
        assert!(!looks_like_rxn_id(""));
    }

    #[test]
    fn reaname_md5_fallback() {
        let d = tempfile::tempdir().unwrap();
        let root = d.path();
        fs::create_dir_all(root.join("rev")).unwrap();
        let hash = md5_hex("Acetaldehyde dehydrogenase");
        fs::write(root.join(format!("rev/{hash}.fasta")), ">x\nX\n").unwrap();
        let opts = SeqfileOptions {
            tax_root: root.into(),
            user_root: root.into(),
            seq_src: SeqSrc::PreferRev,
        };
        let r = resolve_for_reaction(&opts, "RXN-missing", "", "Acetaldehyde dehydrogenase");
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].label, format!("rev/{hash}.fasta"));
    }

    #[test]
    fn reaname_fallback_skips_rxn_id_names() {
        let d = tempfile::tempdir().unwrap();
        let root = d.path();
        fs::create_dir_all(root.join("rev")).unwrap();
        let hash = md5_hex("rxn00001");
        fs::write(root.join(format!("rev/{hash}.fasta")), ">x\nX\n").unwrap();
        let opts = SeqfileOptions {
            tax_root: root.into(),
            user_root: root.into(),
            seq_src: SeqSrc::PreferRev,
        };
        let r = resolve_for_reaction(&opts, "", "", "rxn00001");
        // Even though the md5 file exists, we shouldn't pick it — "rxn00001"
        // is a synthetic id.
        assert!(r.is_empty());
    }

    #[test]
    fn empty_file_treated_as_missing() {
        let d = tempfile::tempdir().unwrap();
        let root = d.path();
        fs::create_dir_all(root.join("rxn")).unwrap();
        fs::write(root.join("rxn/RXN-1.fasta"), "").unwrap();
        let opts = SeqfileOptions {
            tax_root: root.into(),
            user_root: root.into(),
            seq_src: SeqSrc::PreferRev,
        };
        let r = resolve_for_reaction(&opts, "RXN-1", "", "");
        assert!(r.is_empty());
    }
}
