//! End-to-end driver for `gapsmith find`.
//!
//! Orchestrates pathway expansion, seqfile resolution, one concatenated
//! alignment, per-reaction classification, and pathway completeness
//! scoring. Returns a [`FindReport`] that downstream code (CLI + output
//! writers) consumes.
//!
//! Complex / subunit detection lives in [`crate::complex`]; this runner
//! invokes it per-reaction to populate `is_complex` / `complex_status`.

use crate::classify::ClassifyOptions;
use crate::complex::detect_subunits;
use crate::dbhit::DbhitIndex;
use crate::pathways::{self, ExpandedReaction, MatchMode, PathwaySelectOptions};
use crate::seqfile::{self, SeqfileOptions};
use crate::types::{HitStatus, PathwayResult, PwyStatus, ReactionHit};
use gapsmith_align::{AlignOpts, Aligner, Hit};
use gapsmith_db::{ComplexSubunitTable, ExceptionRow, PathwayTable};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum FindError {
    #[error("alignment error: {0}")]
    Align(#[from] gapsmith_align::AlignError),
    #[error("i/o error on `{path}`: {source}")]
    Io { path: std::path::PathBuf, #[source] source: std::io::Error },
    #[error("no reference FASTAs resolved — check --seq-dir and --taxonomy")]
    NoReferences,
}

#[derive(Debug, Clone)]
pub struct FindOptions<'a> {
    pub keyword: &'a str,
    pub match_mode: MatchMode,
    pub exclude_superpathways: bool,
    pub only_active: bool,
    pub bitcutoff: f32,
    pub identcutoff: f32,
    pub ident_exception: f32,
    /// Coverage cutoff passed through to the aligner. Not re-checked in
    /// classification — we trust the aligner to filter.
    pub coverage_pct: u32,
    /// Fraction-of-reactions threshold above which vague (no_seq_data)
    /// reactions dilute completeness. Matches the `vagueCutoff` arg in
    /// analyse_alignments.R.
    pub vague_cutoff: f32,
    /// Main completeness threshold, 0–1 (e.g. 0.80).
    pub completeness_hint_off: f32,
    /// Relaxed threshold applied when all key reactions are present, 0–1.
    pub completeness_hint_on: f32,
    pub strict_candidates: bool,
    /// Fraction of known subunits needed for a complex to count as
    /// present (0–1). Default 0.5.
    pub subunit_cutoff: f32,
    /// Tax IDs kept when filtering by `taxrange`. Empty → no filter.
    pub valid_tax_ids: &'a [String],
}

impl<'a> Default for FindOptions<'a> {
    fn default() -> Self {
        Self {
            keyword: "all",
            match_mode: MatchMode::Regex,
            exclude_superpathways: true,
            only_active: true,
            bitcutoff: 200.0,
            identcutoff: 0.0,
            ident_exception: 70.0,
            coverage_pct: 75,
            vague_cutoff: 0.3,
            completeness_hint_off: 0.80,
            completeness_hint_on: 0.66,
            strict_candidates: false,
            subunit_cutoff: 0.5,
            valid_tax_ids: &[],
        }
    }
}

#[derive(Debug, Clone)]
pub struct FindReport {
    pub reactions: Vec<ReactionHit>,
    pub pathways: Vec<PathwayResult>,
}

/// Run the full find pipeline.
#[allow(clippy::too_many_arguments)]
pub fn run_find(
    table: &PathwayTable,
    exceptions: &[ExceptionRow],
    subunit_dict: &ComplexSubunitTable,
    dbhit_index: &DbhitIndex,
    seq_opts: &SeqfileOptions,
    genome_fasta: &Path,
    aligner: &dyn Aligner,
    align_opts: &AlignOpts,
    opts: &FindOptions<'_>,
    workdir: &Path,
) -> Result<FindReport, FindError> {
    std::fs::create_dir_all(workdir).map_err(|e| FindError::Io {
        path: workdir.to_path_buf(),
        source: e,
    })?;

    // 1. Pick pathways and expand to reactions.
    let expanded = pathways::select(
        table,
        &PathwaySelectOptions {
            keyword: opts.keyword,
            mode: opts.match_mode,
            exclude_superpathways: opts.exclude_superpathways,
            only_active: opts.only_active,
            valid_tax_ids: opts.valid_tax_ids,
        },
    );
    tracing::info!(reactions = expanded.len(), "pathways expanded to reactions");
    if expanded.is_empty() {
        return Ok(FindReport { reactions: Vec::new(), pathways: Vec::new() });
    }

    // 2. For each unique reaction, resolve reference FASTAs and scan each
    //    one for subunit assignments.
    //
    // Multiple reactions frequently share the same reference FASTA (two
    // reactions with the same EC → same `rev/<EC>.fasta`). Cache the
    // parsed header list keyed by absolute path so each file is read
    // and parsed exactly once per `find` invocation. Entries are wrapped
    // in `Arc` so the cache holds ownership while callers borrow the
    // parsed headers freely.
    let unique_reactions = dedup_reactions(&expanded);
    let mut seq_by_key: HashMap<ReactionKey, Vec<seqfile::ResolvedSeq>> = HashMap::new();
    let mut subunit_by_qseqid: HashMap<(String, String), Option<String>> = HashMap::new();
    let mut subunit_count_by_key: HashMap<ReactionKey, (bool, u32, String)> = HashMap::new();
    let mut fasta_header_cache: HashMap<std::path::PathBuf, std::sync::Arc<Vec<(String, String)>>> =
        HashMap::new();
    for k in &unique_reactions {
        let resolved = seqfile::resolve_for_reaction(seq_opts, &k.rxn, &k.ec, &k.name);
        // Scan every resolved fasta's headers for subunit clues.
        let mut desc_owned: Vec<(String, String, String)> = Vec::new();
        for r in &resolved {
            let entries = fasta_header_cache
                .entry(r.path.clone())
                .or_insert_with(|| {
                    std::sync::Arc::new(read_fasta_headers(&r.path).unwrap_or_default())
                })
                .clone();
            for (acc, full) in entries.iter() {
                desc_owned.push((r.label.clone(), acc.clone(), full.clone()));
            }
        }
        // Run detect_subunits across the combined descriptor list (gapseq's
        // R code operates on the unioned set per reaction).
        let refs: Vec<&str> =
            desc_owned.iter().map(|(_, _, full)| full.as_str()).collect();
        let assigns = detect_subunits(&k.rxn, &refs, subunit_dict);
        let mut distinct: HashSet<String> = HashSet::new();
        for ((label, acc, _), su) in desc_owned.iter().zip(assigns.iter()) {
            subunit_by_qseqid.insert((label.clone(), acc.clone()), su.clone());
            if let Some(name) = su {
                if name != "Subunit undefined" {
                    distinct.insert(name.clone());
                }
            }
        }
        let subunit_count = distinct.len() as u32;
        let is_complex = subunit_count >= 2;
        let mut list: Vec<String> = distinct.into_iter().collect();
        list.sort();
        subunit_count_by_key
            .insert(k.clone(), (is_complex, subunit_count, list.join(",")));
        seq_by_key.insert(k.clone(), resolved);
    }

    // 3. Concatenate all reference fastas into one query.faa.
    let query_path = workdir.join("query.faa");
    let n_concat = concat_refs(&query_path, &seq_by_key)?;
    tracing::info!(seqfiles = n_concat, "concatenated reference fastas");

    // 4. If nothing to align, emit empty results.
    let hits_by_file: HashMap<String, Vec<Hit>> = if n_concat > 0 {
        let raw_hits = aligner.align(&query_path, genome_fasta, align_opts)?;
        tracing::info!(hits = raw_hits.len(), "alignment complete");
        group_by_file(raw_hits)
    } else {
        HashMap::new()
    };

    // 5. Build per-reaction classified results.
    let exception_set: HashSet<String> = exceptions.iter().map(|e| e.id.clone()).collect();
    let classify_opts = ClassifyOptions {
        bitcutoff: opts.bitcutoff,
        identcutoff: opts.identcutoff,
        ident_exception: opts.ident_exception,
        exception_ecs: &exception_set,
    };

    let mut reactions: Vec<ReactionHit> = Vec::with_capacity(expanded.len());
    for e in &expanded {
        // R emits one row per unique (rxn, name, ec, stitle, complex). We
        // reproduce that here: if the reaction has ≥ 1 hit, emit one row
        // per hit; otherwise emit one placeholder row.
        let rows = build_reaction_rows(
            e,
            &seq_by_key,
            &hits_by_file,
            &subunit_by_qseqid,
            &subunit_count_by_key,
            dbhit_index,
            &classify_opts,
            opts.subunit_cutoff,
        );
        reactions.extend(rows);
    }

    // 5b. Sort reactions by (pathway, rxn, complex, -bitscore). Matches
    //     R's `rxndt <- rxndt[order(pathway, rxn, complex, -bitscore)]`.
    reactions.sort_by(|a, b| {
        a.pathway
            .cmp(&b.pathway)
            .then_with(|| a.rxn.cmp(&b.rxn))
            .then_with(|| {
                a.complex
                    .as_deref()
                    .unwrap_or("")
                    .cmp(b.complex.as_deref().unwrap_or(""))
            })
            .then_with(|| {
                b.bitscore
                    .partial_cmp(&a.bitscore)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    // 6. Pathway scoring.
    let pathway_results = score_pathways(&reactions, opts);

    // 7. Stamp each reaction row with its pathway status and `has_dbhit`.
    let pwy_status_map: HashMap<&str, Option<PwyStatus>> =
        pathway_results.iter().map(|p| (p.id.as_str(), p.status)).collect();
    for r in &mut reactions {
        r.pathway_status = pwy_status_map.get(r.pathway.as_str()).copied().unwrap_or(None);
        r.has_dbhit = !r.dbhit.is_empty();
    }

    Ok(FindReport { reactions, pathways: pathway_results })
}

// -- Internal plumbing --

#[derive(Clone, Hash, PartialEq, Eq)]
struct ReactionKey {
    rxn: String,
    name: String,
    ec: String,
}

fn dedup_reactions(expanded: &[ExpandedReaction]) -> Vec<ReactionKey> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for e in expanded {
        let k = ReactionKey { rxn: e.rxn.clone(), name: e.name.clone(), ec: e.ec.clone() };
        if seen.insert(k.clone()) {
            out.push(k);
        }
    }
    out
}

/// Concatenate every resolved reference FASTA into one big `query.faa`,
/// prepending `<label>|` to each sequence header so the downstream hit
/// parser can trace a hit back to its originating reference file (matching
/// `prepare_batch_alignments.R:374-407`).
fn concat_refs(
    out: &Path,
    seq_by_key: &HashMap<ReactionKey, Vec<seqfile::ResolvedSeq>>,
) -> Result<usize, FindError> {
    let f = File::create(out).map_err(|e| FindError::Io {
        path: out.to_path_buf(),
        source: e,
    })?;
    let mut w = BufWriter::new(f);
    let mut emitted: HashSet<&str> = HashSet::new();
    let mut n = 0usize;
    for list in seq_by_key.values() {
        for r in list {
            if !emitted.insert(r.label.as_str()) {
                continue;
            }
            let mut buf = Vec::new();
            File::open(&r.path)
                .and_then(|mut f| f.read_to_end(&mut buf))
                .map_err(|e| FindError::Io { path: r.path.clone(), source: e })?;
            // rewrite headers `>...` to `>LABEL|...`
            for line in buf.split(|b| *b == b'\n') {
                if line.is_empty() {
                    continue;
                }
                if line.first() == Some(&b'>') {
                    write!(w, ">{}|", r.label).map_err(|e| FindError::Io {
                        path: out.to_path_buf(),
                        source: e,
                    })?;
                    w.write_all(&line[1..]).map_err(|e| FindError::Io {
                        path: out.to_path_buf(),
                        source: e,
                    })?;
                    writeln!(w).map_err(|e| FindError::Io {
                        path: out.to_path_buf(),
                        source: e,
                    })?;
                } else {
                    w.write_all(line).map_err(|e| FindError::Io {
                        path: out.to_path_buf(),
                        source: e,
                    })?;
                    writeln!(w).map_err(|e| FindError::Io {
                        path: out.to_path_buf(),
                        source: e,
                    })?;
                }
            }
            n += 1;
        }
    }
    w.flush().map_err(|e| FindError::Io {
        path: out.to_path_buf(),
        source: e,
    })?;
    Ok(n)
}

fn group_by_file(hits: Vec<Hit>) -> HashMap<String, Vec<Hit>> {
    let mut map: HashMap<String, Vec<Hit>> = HashMap::new();
    for h in hits {
        // The qseqid is the full FASTA header first-token, i.e.
        // `<label>|<orig_acc>` (BLAST/diamond preserve first-token only).
        let (file, orig) = match h.qseqid.split_once('|') {
            Some((f, rest)) => (f.to_string(), rest.to_string()),
            None => (String::new(), h.qseqid.clone()),
        };
        let mut hh = h;
        hh.qseqid = orig;
        map.entry(file).or_default().push(hh);
    }
    map
}

#[allow(clippy::too_many_arguments)]
fn build_reaction_rows(
    e: &ExpandedReaction,
    seq_by_key: &HashMap<ReactionKey, Vec<seqfile::ResolvedSeq>>,
    hits_by_file: &HashMap<String, Vec<Hit>>,
    subunit_by_qseqid: &HashMap<(String, String), Option<String>>,
    subunit_count_by_key: &HashMap<ReactionKey, (bool, u32, String)>,
    dbhit_index: &DbhitIndex,
    classify_opts: &ClassifyOptions<'_>,
    subunit_cutoff: f32,
) -> Vec<ReactionHit> {
    let key = ReactionKey { rxn: e.rxn.clone(), name: e.name.clone(), ec: e.ec.clone() };
    let resolved = seq_by_key.get(&key);
    let files: Vec<String> = resolved
        .map(|v| v.iter().map(|r| r.label.clone()).collect())
        .unwrap_or_default();
    let has_seq_data = !files.is_empty();

    // Collect per-file hit list (preserves file attribution).
    let mut all_hits: Vec<(String, Hit)> = Vec::new();
    for f in &files {
        if let Some(v) = hits_by_file.get(f) {
            for h in v {
                all_hits.push((f.clone(), h.clone()));
            }
        }
    }

    // Pre-compute dbhit once per reaction.
    let dbhit = dbhit_index.lookup(&e.rxn, &e.name, &e.ec);

    // Subunit metadata (per-reaction, applied to every emitted row).
    // `complex_scanned` is true iff we had at least one reference FASTA
    // header to analyse — matches gapseq's convention of printing `NA`
    // only when complex_detection actually ran.
    let (is_complex, subunit_count, subunits_str) = subunit_count_by_key
        .get(&key)
        .cloned()
        .unwrap_or((false, 0, String::new()));
    let complex_scanned = has_seq_data;

    // Per-reaction exception flag — computed once so every emitted row
    // carries the same value.
    let exception = crate::classify::ec_is_exception(&e.ec, classify_opts.exception_ecs);

    // Complex-status computation across the full hit list.
    let (complex_subunits_found, complex_undef_found, complex_status) = if is_complex {
        let mut found: HashSet<String> = HashSet::new();
        let mut undef = false;
        for (f, h) in &all_hits {
            if h.bitscore < classify_opts.bitcutoff || h.pident < classify_opts.identcutoff {
                continue;
            }
            if exception && h.pident < classify_opts.ident_exception {
                continue;
            }
            if let Some(Some(sub)) = subunit_by_qseqid.get(&(f.clone(), h.qseqid.clone())) {
                if sub == "Subunit undefined" {
                    undef = true;
                } else {
                    found.insert(sub.clone());
                }
            }
        }
        let sf = found.len() as u32;
        let ratio = if subunit_count == 0 { 0.0 } else { sf as f32 / subunit_count as f32 };
        let st = if ratio > subunit_cutoff || (ratio == subunit_cutoff && undef) {
            Some(1u8)
        } else {
            None
        };
        (Some(sf), Some(undef), st)
    } else {
        (None, None, None)
    };

    // Emit one row per hit. When there are no hits, emit a single
    // placeholder row with status ∈ {NoBlast, NoSeqData, Spontaneous}.
    if all_hits.is_empty() {
        let status = if !has_seq_data && e.spont {
            HitStatus::Spontaneous
        } else if has_seq_data {
            HitStatus::NoBlast
        } else {
            HitStatus::NoSeqData
        };
        let label_for_meta = files.first().cloned().unwrap_or_default();
        let (src, reftype) = derive_src_and_type(&label_for_meta);
        return vec![ReactionHit {
            pathway: e.pathway.clone(),
            pathway_status: None,
            rxn: e.rxn.clone(),
            name: e.name.clone(),
            ec: e.ec.clone(),
            keyrea: e.keyrea,
            spont: e.spont,
            is_complex,
            subunit_count,
            // gapseq leaves the column blank when complex detection didn't
            // run (no reference fastas). The output layer turns a blank
            // with complex_scanned=true into `NA`; keep blank otherwise.
            subunits: if complex_scanned { subunits_str.clone() } else { String::new() },
            complex: None,
            subunits_found: if complex_scanned { complex_subunits_found.or(Some(0)) } else { None },
            subunit_undefined_found: complex_undef_found,
            complex_status,
            file: files.first().cloned(),
            dbhit,
            has_dbhit: false,
            src,
            reftype,
            qseqid: None,
            pident: None,
            evalue: None,
            bitscore: None,
            qcov: None,
            stitle: None,
            sstart: None,
            send: None,
            exception,
            status,
        }];
    }

    // Per-hit rows, sorted by (-bitscore, stitle) to mirror R's
    // `order(rxn, name, ec, stitle, complex, -bitscore)` then unique().
    let mut rows: Vec<ReactionHit> = Vec::with_capacity(all_hits.len());
    for (f, h) in &all_hits {
        let pass_main = h.bitscore >= classify_opts.bitcutoff && h.pident >= classify_opts.identcutoff;
        let pass_exc = if exception { h.pident >= classify_opts.ident_exception } else { true };
        let status = if pass_main && pass_exc {
            HitStatus::GoodBlast
        } else {
            HitStatus::BadBlast
        };
        let complex = subunit_by_qseqid.get(&(f.clone(), h.qseqid.clone())).cloned().flatten();
        let (src, reftype) = derive_src_and_type(f);
        rows.push(ReactionHit {
            pathway: e.pathway.clone(),
            pathway_status: None,
            rxn: e.rxn.clone(),
            name: e.name.clone(),
            ec: e.ec.clone(),
            keyrea: e.keyrea,
            spont: e.spont,
            is_complex,
            subunit_count,
            subunits: subunits_str.clone(),
            complex,
            subunits_found: complex_subunits_found,
            subunit_undefined_found: complex_undef_found,
            complex_status,
            file: Some(f.clone()),
            dbhit: dbhit.clone(),
            has_dbhit: false,
            src,
            reftype,
            qseqid: Some(h.qseqid.clone()),
            pident: Some(h.pident),
            evalue: Some(h.evalue),
            bitscore: Some(h.bitscore),
            qcov: Some(h.qcov),
            stitle: Some(h.stitle.clone()),
            sstart: Some(h.sstart),
            send: Some(h.send),
            exception,
            status,
        });
    }

    // R dedups by (rxn, name, ec, stitle, complex) keeping the row with
    // the highest bitscore (because the preceding sort puts -bitscore last).
    rows.sort_by(|a, b| {
        let key_a = (a.stitle.as_deref().unwrap_or(""), a.complex.as_deref().unwrap_or(""));
        let key_b = (b.stitle.as_deref().unwrap_or(""), b.complex.as_deref().unwrap_or(""));
        key_a
            .cmp(&key_b)
            .then_with(|| b.bitscore.partial_cmp(&a.bitscore).unwrap_or(std::cmp::Ordering::Equal))
    });
    let mut seen: HashSet<(String, String)> = HashSet::new();
    rows.retain(|r| {
        seen.insert((
            r.stitle.clone().unwrap_or_default(),
            r.complex.clone().unwrap_or_default(),
        ))
    });
    rows
}

/// Map a file label like `rev/1.1.1.1.fasta` into (`src`, `type`).
/// - `rxn/*.fasta` → src=rxn, type=metacyc
/// - `rev/*.fasta` / `unrev/*.fasta` with `<EC>.fasta` → type=EC
/// - `rev/*.fasta` / `unrev/*.fasta` with MD5 stem → type=reaName
/// - `user/*.fasta` → src=user, type preserved from filename if distinguishable
fn derive_src_and_type(label: &str) -> (String, String) {
    if label.is_empty() {
        return (String::new(), String::new());
    }
    let (src, rest) = match label.split_once('/') {
        Some(p) => p,
        None => return (String::new(), String::new()),
    };
    let stem = rest.trim_end_matches(".fasta");
    let rtype = match src {
        "rxn" => "metacyc",
        "user" => "user",
        _ => {
            if is_ec_stem(stem) {
                "EC"
            } else {
                "reaName"
            }
        }
    };
    (src.to_string(), rtype.to_string())
}

/// True if `stem` looks like an EC number: `N.N.N.N` with digits and
/// optional `n` in each slot.
fn is_ec_stem(stem: &str) -> bool {
    let parts: Vec<&str> = stem.split('.').collect();
    if parts.len() != 4 {
        return false;
    }
    parts.iter().all(|p| !p.is_empty() && p.chars().all(|c| c.is_ascii_digit() || c == 'n'))
}

/// Scan a FASTA file and return `(accession, full_descriptor)` per header.
/// `accession` is the first whitespace-delimited token of the header;
/// `full_descriptor` is the entire header line after `>`.
fn read_fasta_headers(path: &Path) -> Result<Vec<(String, String)>, std::io::Error> {
    let f = File::open(path)?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    for line in r.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            let acc = rest
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            out.push((acc, rest.to_string()));
        }
    }
    Ok(out)
}

fn score_pathways(reactions: &[ReactionHit], opts: &FindOptions<'_>) -> Vec<PathwayResult> {
    // One entry per pathway; reactions are already one-per-(pathway,rxn).
    let mut by_pwy: BTreeMap<String, (String, Vec<&ReactionHit>)> = BTreeMap::new();
    for r in reactions {
        let entry = by_pwy.entry(r.pathway.clone()).or_insert_with(|| (String::new(), Vec::new()));
        entry.1.push(r);
    }

    // Reactions were emitted with pathway_name on the ExpandedReaction,
    // but not stored on ReactionHit. We infer from the first pass by
    // keeping a separate lookup elsewhere; here we leave blank and let
    // the caller supply it. Simplification: we stash nothing and the
    // caller downstream can join with the pathway table by id.
    let _ = &mut by_pwy; // placeholder for name injection below

    let mut out = Vec::new();
    for (pwy_id, (_name, rxns)) in &by_pwy {
        // Deduplicate by reaction id within a pathway (same rxn may
        // appear twice if reaEc has duplicates). Sort alphabetically by
        // rxn id to mirror `analyse_alignments.R:143` which runs
        // `order(pathway, rxn, complex, -bitscore)` before aggregation.
        let mut seen = HashSet::new();
        let mut rxns: Vec<&ReactionHit> = rxns
            .iter()
            .copied()
            .filter(|r| seen.insert(r.rxn.clone()))
            .collect();
        rxns.sort_by(|a, b| a.rxn.cmp(&b.rxn));

        let n_reaction = rxns.len() as u32;
        let n_spont = rxns.iter().filter(|r| r.spont).count() as u32;
        let n_vague = rxns.iter().filter(|r| r.status == HitStatus::NoSeqData).count() as u32;
        let n_key = rxns
            .iter()
            .filter(|r| r.keyrea && r.status != HitStatus::NoSeqData)
            .count() as u32;
        let reaction_found = |r: &&ReactionHit| {
            (!r.is_complex && r.status == HitStatus::GoodBlast)
                || (r.is_complex && r.complex_status.is_some())
        };
        let n_found = rxns.iter().filter(|r| reaction_found(r)).count() as u32;
        let n_key_found = rxns
            .iter()
            .filter(|r| reaction_found(r) && r.keyrea)
            .count() as u32;

        let found_ids: Vec<String> = rxns
            .iter()
            .filter(|r| reaction_found(r))
            .map(|r| r.rxn.clone())
            .collect();
        let spont_ids: Vec<String> =
            rxns.iter().filter(|r| r.spont).map(|r| r.rxn.clone()).collect();
        let key_ids: Vec<String> =
            rxns.iter().filter(|r| r.keyrea).map(|r| r.rxn.clone()).collect();

        // Match R's analyse_alignments.R:164-165 exactly — R uses the ratio
        // `NrVague/(NrReaction-NrSpontaneous)` with 0/0 guarded by data.table
        // (coerces to NaN which is then neither < nor >=). We handle that
        // explicitly: if the denominator is 0, force completeness to 0.
        let denom_no_spont = n_reaction.saturating_sub(n_spont);
        let hint_off = opts.completeness_hint_off as f64;
        let hint_on = opts.completeness_hint_on as f64;
        let completeness: f64 = if denom_no_spont == 0 {
            0.0
        } else {
            let vague_frac = n_vague as f64 / denom_no_spont as f64;
            if vague_frac < opts.vague_cutoff as f64 {
                let denom = n_reaction.saturating_sub(n_vague).saturating_sub(n_spont);
                if denom == 0 {
                    0.0
                } else {
                    n_found as f64 / denom as f64 * 100.0
                }
            } else {
                n_found as f64 / denom_no_spont as f64 * 100.0
            }
        };

        let all_key_found = n_key == n_key_found;
        let mut prediction = if !opts.strict_candidates {
            completeness >= hint_off * 100.0 && all_key_found
        } else {
            completeness >= hint_off * 100.0
        };
        if !opts.strict_candidates
            && n_key > 0
            && all_key_found
            && completeness >= hint_on * 100.0
        {
            prediction = true;
        }
        if n_reaction == n_vague + n_spont {
            prediction = false;
        }

        let status = if (completeness - 100.0).abs() < 1e-9 {
            Some(PwyStatus::Full)
        } else if prediction && all_key_found {
            Some(PwyStatus::Threshold)
        } else if prediction && completeness < hint_off * 100.0 {
            Some(PwyStatus::Keyenzyme)
        } else {
            None
        };

        out.push(PathwayResult {
            id: pwy_id.clone(),
            name: String::new(), // stamped below
            prediction,
            completeness,
            status,
            n_reaction,
            n_spontaneous: n_spont,
            n_vague,
            n_key_reaction: n_key,
            n_reaction_found: n_found,
            n_key_reaction_found: n_key_found,
            reactions_found: found_ids,
            spontaneous_reactions: spont_ids,
            key_reactions: key_ids,
        });
    }
    out
}
