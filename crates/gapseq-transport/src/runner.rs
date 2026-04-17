//! End-to-end driver for `gapseq find-transport`.
//!
//! Port of `src/analyse_alignments_transport.R` plus the shell
//! orchestration from `src/transporter.sh`.

use crate::data::{group_seed_by_type_met, load_seed_transporter, load_tcdb_all};
use crate::filter::{build_small_fasta, BuildSmallResult};
use crate::output::TransporterRow;
use crate::tc::{extract_tc_id, type_of};
use gapseq_align::{AlignOpts, Aligner, Hit};
use gapseq_db::SubexRow;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

#[derive(Debug, thiserror::Error)]
pub enum TransportError {
    #[error("alignment error: {0}")]
    Align(#[from] gapseq_align::AlignError),
    #[error("i/o error on `{path}`: {source}")]
    Io { path: PathBuf, #[source] source: std::io::Error },
    #[error("required reference FASTA missing: `{0}`")]
    MissingReference(PathBuf),
    #[error("data error: {0}")]
    Data(String),
}

#[derive(Debug, Clone)]
pub struct TransportOptions<'a> {
    pub bitcutoff: f32,
    pub identcutoff: f32,
    pub coverage_pct: u32,
    /// When true, skip the "alt-transporter" fallback that injects SEED
    /// reactions from another TC type when the exact type has no match.
    pub nouse_alternatives: bool,
    /// Optional single-substrate restriction (matches the `-m` flag).
    pub only_met: Option<&'a str>,
}

impl<'a> Default for TransportOptions<'a> {
    fn default() -> Self {
        Self {
            bitcutoff: 50.0,
            identcutoff: 0.0,
            coverage_pct: 75,
            nouse_alternatives: false,
            only_met: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TransportReport {
    pub rows: Vec<TransporterRow>,
}

/// Re-export for external callers.
pub type Row = TransporterRow;

#[allow(clippy::too_many_arguments)]
pub fn run(
    reference_fastas: &[&Path],
    seed_transporter_path: &Path,
    seed_transporter_custom_path: &Path,
    tcdb_substrates_path: &Path,
    tcdb_custom_path: &Path,
    subex_rows: &[SubexRow],
    metabolite_name_by_id: &HashMap<String, String>,
    genome_fasta: &Path,
    aligner: &dyn Aligner,
    align_opts: &AlignOpts,
    opts: &TransportOptions<'_>,
    workdir: &Path,
) -> Result<TransportReport, TransportError> {
    std::fs::create_dir_all(workdir).map_err(|e| TransportError::Io {
        path: workdir.to_path_buf(),
        source: e,
    })?;

    // Step 1: load the tables we'll need for joining.
    let seed_rows = load_seed_transporter(seed_transporter_path, seed_transporter_custom_path)
        .map_err(|e| TransportError::Io {
            path: seed_transporter_path.to_path_buf(),
            source: e,
        })?;
    let seed_group = group_seed_by_type_met(&seed_rows);
    // metid → all reaction ids across every type (for alt-transporter).
    let mut seed_group_any_type: HashMap<String, Vec<String>> = HashMap::new();
    for ((_t, met), rxns) in &seed_group {
        let entry = seed_group_any_type.entry(met.clone()).or_default();
        for r in rxns {
            if !entry.contains(r) {
                entry.push(r.clone());
            }
        }
    }
    for v in seed_group_any_type.values_mut() {
        v.sort();
    }

    let tcdb_all = load_tcdb_all(tcdb_substrates_path, tcdb_custom_path).map_err(|e| {
        TransportError::Io { path: tcdb_substrates_path.to_path_buf(), source: e }
    })?;

    // Step 2: small.fasta.
    for p in reference_fastas {
        if !p.is_file() {
            return Err(TransportError::MissingReference(p.to_path_buf()));
        }
    }
    let small = build_small_fasta(
        reference_fastas,
        subex_rows,
        &tcdb_all,
        opts.only_met,
        workdir,
    )
    .map_err(|e| TransportError::Io { path: workdir.to_path_buf(), source: e })?;

    if small.fasta_header_small.is_empty() {
        tracing::warn!("small.fasta is empty; no transport alignments will run");
        return Ok(TransportReport { rows: Vec::new() });
    }

    // Step 3: alignment.
    let mut align_opts_use = align_opts.clone();
    align_opts_use.coverage_pct = opts.coverage_pct;
    let hits = aligner.align(&small.small_fasta, genome_fasta, &align_opts_use)?;
    tracing::info!(n = hits.len(), "transport alignment done");

    // Step 4: analyse alignments.
    let analysed = analyse_hits(
        &hits,
        &small,
        &tcdb_all,
        &seed_group,
        &seed_group_any_type,
        metabolite_name_by_id,
        opts,
    );

    Ok(TransportReport { rows: analysed })
}

/// Core analysis routine — port of `analyse_alignments_transport.R`.
#[allow(clippy::too_many_arguments)]
fn analyse_hits(
    hits: &[Hit],
    small: &BuildSmallResult,
    tcdb_all: &HashMap<String, String>,
    seed_by_type_met: &HashMap<(String, String), Vec<String>>,
    seed_by_met: &HashMap<String, Vec<String>>,
    metabolite_name_by_id: &HashMap<String, String>,
    opts: &TransportOptions<'_>,
) -> Vec<TransporterRow> {
    // Identity cutoff is applied first.
    let candidates: Vec<&Hit> = hits.iter().filter(|h| h.pident >= opts.identcutoff).collect();

    // Attach a TC id to every candidate — drop the rest.
    let mut tc_hits: Vec<(&Hit, String)> = Vec::new();
    for h in &candidates {
        if let Some(tc) = extract_tc_id(&h.qseqid) {
            tc_hits.push((*h, tc.to_string()));
        }
    }

    // Pre-build a header-substrate index so findsubs' fallback path
    // doesn't recompile regexes per key. The index maps each sub_key to
    // the set of TC ids whose headers mention it.
    let header_sub_index = build_header_sub_index(
        &small.fasta_header_small,
        &small.sub_key_orig,
    );

    // Per-TC substrate resolution. `findsubs()` cache — Arc avoids
    // cloning the substrate vector on every cache hit.
    let mut subs_cache: HashMap<String, Arc<Vec<String>>> = HashMap::new();
    let resolve_subs =
        |tc: &str, cache: &mut HashMap<String, Arc<Vec<String>>>| -> Arc<Vec<String>> {
            if let Some(v) = cache.get(tc) {
                return Arc::clone(v);
            }
            let v = Arc::new(findsubs(tc, tcdb_all, &header_sub_index));
            cache.insert(tc.to_string(), Arc::clone(&v));
            v
        };

    // Build a lookup substrate-name → exchange id. Same key may map to
    // multiple exids — we split rows across them exactly like R's
    // cartesian merge.
    let mut subid_map: HashMap<String, Vec<String>> = HashMap::new();
    for (name_lc, exid) in &small.subid {
        subid_map.entry(name_lc.clone()).or_default().push(exid.clone());
    }

    #[derive(Debug)]
    struct Expanded<'a> {
        hit: &'a Hit,
        tc: String,
        sub: String,
        exid: String,
    }

    let mut expanded: Vec<Expanded> = Vec::new();
    for (h, tc) in &tc_hits {
        let subs = resolve_subs(tc, &mut subs_cache);
        for sub in subs.iter() {
            let sub_lc = sub.to_ascii_lowercase();
            if let Some(exids) = subid_map.get(&sub_lc) {
                for exid in exids {
                    expanded.push(Expanded {
                        hit: h,
                        tc: tc.clone(),
                        sub: sub.clone(),
                        exid: exid.clone(),
                    });
                }
            }
        }
    }

    // Dedup by (qseqid, stitle, exid) keeping highest bitscore
    // (`analyse_alignments_transport.R:105-107`).
    expanded.sort_by(|a, b| {
        a.hit.qseqid
            .cmp(&b.hit.qseqid)
            .then_with(|| a.hit.stitle.cmp(&b.hit.stitle))
            .then_with(|| a.exid.cmp(&b.exid))
            .then_with(|| {
                b.hit
                    .bitscore
                    .partial_cmp(&a.hit.bitscore)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    // Dedup: consecutive runs of identical (qseqid, stitle, exid) after
    // sort. We can't borrow into `retain` across iterations, so use a
    // simple dedup_by on the already-sorted vector.
    expanded.dedup_by(|a, b| {
        a.hit.qseqid == b.hit.qseqid
            && a.hit.stitle == b.hit.stitle
            && a.exid == b.exid
    });
    let mut rows: Vec<TransporterRow> = Vec::with_capacity(expanded.len());
    let mut mets_with_rxn_highbit: HashSet<String> = HashSet::new();
    let mut mets_with_rxn: HashSet<String> = HashSet::new();

    // First pass: primary transporter assignment.
    for e in &expanded {
        let metid = strip_exchange_affixes(&e.exid);
        let (_, type_str) = match type_of(&e.tc) {
            Some(x) => x,
            None => continue,
        };
        let rxn = seed_by_type_met
            .get(&(type_str.to_string(), metid.clone()))
            .map(|v| v.join(","));
        let comment = if rxn.is_some() { Some("transporter") } else { None };
        if comment.is_some() && e.hit.bitscore >= opts.bitcutoff {
            mets_with_rxn_highbit.insert(metid.clone());
        }
        if comment.is_some() {
            mets_with_rxn.insert(metid.clone());
        }
        let sub_gapseq =
            metabolite_name_by_id.get(&metid).cloned().unwrap_or_default();
        rows.push(TransporterRow {
            id: e.hit.qseqid.clone(),
            tc: e.tc.clone(),
            sub: e.sub.clone(),
            sub_gapseq,
            exid: e.exid.clone(),
            rea: rxn.unwrap_or_default(),
            qseqid: e.hit.qseqid.clone(),
            pident: e.hit.pident,
            evalue: e.hit.evalue,
            bitscore: e.hit.bitscore,
            qcov: e.hit.qcov,
            stitle: e.hit.stitle.clone(),
            sstart: e.hit.sstart,
            send: e.hit.send,
            comment: comment.map(|s| s.to_string()),
            metid,
        });
    }

    // Second pass: alt-transporter fallback.
    // Any metid without a SEED rxn AND not already accounted for in the
    // high-bit set should inherit reactions from seed_by_met (union
    // across TC types).
    if !opts.nouse_alternatives {
        for r in rows.iter_mut() {
            if r.rea.is_empty() && !mets_with_rxn.contains(&r.metid) {
                if let Some(alts) = seed_by_met.get(&r.metid) {
                    r.rea = alts.join(",");
                    r.comment = Some("alt-transporter".into());
                }
            }
        }
    }

    // Ordering: comment (transporter < alt-transporter), tc, -bitscore
    // (`analyse_alignments_transport.R:130`).
    rows.sort_by(|a, b| {
        comment_rank(&a.comment)
            .cmp(&comment_rank(&b.comment))
            .then_with(|| a.tc.cmp(&b.tc))
            .then_with(|| {
                b.bitscore.partial_cmp(&a.bitscore).unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    rows
}

fn comment_rank(c: &Option<String>) -> u8 {
    match c.as_deref() {
        Some("transporter") => 0,
        Some("alt-transporter") => 1,
        _ => 2,
    }
}

/// Strip `EX_` prefix and `_eN` suffix from an exchange-reaction id to
/// get the bare metabolite id. Replaces the previous `Regex::new(r"^EX_|_e.$")`.
fn strip_exchange_affixes(exid: &str) -> String {
    let s = exid.strip_prefix("EX_").unwrap_or(exid);
    if s.len() >= 3 {
        let tail = &s[s.len() - 3..];
        if tail.starts_with("_e") && tail.as_bytes()[2].is_ascii_digit() {
            return s[..s.len() - 3].to_string();
        }
    }
    s.to_string()
}

/// Pre-computed index mapping each substrate key (original case) to the
/// set of TC ids whose FASTA headers mention it. Built once, queried
/// from every `findsubs()` fallback. Replaces the per-TC-id regex
/// compilation loop.
struct HeaderSubIndex {
    /// sub_key_orig → set of TC ids whose headers contain this substrate.
    by_key: HashMap<String, HashSet<String>>,
}

fn build_header_sub_index(
    fasta_headers: &[String],
    sub_key_orig: &[String],
) -> HeaderSubIndex {
    let headers_lc: Vec<String> = fasta_headers
        .iter()
        .map(|h| h.to_ascii_lowercase())
        .collect();

    let mut by_key: HashMap<String, HashSet<String>> = HashMap::new();
    for k in sub_key_orig {
        let k_lc = k.to_ascii_lowercase();
        let mut tcs: HashSet<String> = HashSet::new();
        for (header_lc, header_orig) in headers_lc.iter().zip(fasta_headers.iter()) {
            if !header_lc.contains(&k_lc) {
                continue;
            }
            if let Some(tc) = extract_tc_id(header_orig) {
                tcs.insert(tc.to_string());
            }
        }
        if !tcs.is_empty() {
            by_key.insert(k.clone(), tcs);
        }
    }
    HeaderSubIndex { by_key }
}

/// Port of `findsubs(tc)` — resolve substrate names either via tcdb_all
/// or via the pre-built header index. Returns the `sub_key_orig` entries
/// (original casing) that match.
fn findsubs(
    tc: &str,
    tcdb_all: &HashMap<String, String>,
    header_index: &HeaderSubIndex,
) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    let mut out_seen: HashSet<String> = HashSet::new();

    // (a) Via tcdb_all — O(K + S) with HashSet.
    if let Some(subs) = tcdb_all.get(tc) {
        let candidates_lc: HashSet<String> = subs
            .split(';')
            .map(|s| s.trim().to_ascii_lowercase())
            .filter(|s| !s.is_empty())
            .collect();
        for (k, tcs) in &header_index.by_key {
            let _ = tcs; // used only in path (b)
            let k_lc = k.to_ascii_lowercase();
            if candidates_lc.contains(&k_lc) && out_seen.insert(k.clone()) {
                out.push(k.clone());
            }
        }
        if !out.is_empty() {
            return out;
        }
    }

    // (b) Fallback: check if this TC id appears in any header that
    // mentions a substrate key. O(K) lookups against the pre-built index.
    for (k, tcs) in &header_index.by_key {
        if tcs.contains(tc) && out_seen.insert(k.clone()) {
            out.push(k.clone());
        }
    }
    out
}
