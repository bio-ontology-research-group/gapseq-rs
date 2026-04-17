//! Candidate reaction selection — simplified port of
//! `src/prepare_candidate_reaction_tables.R`.
//!
//! The R original runs EC/TC conflict-overlap math and complex-subunit
//! deduplication; this port focuses on the core flow needed to select
//! "draft" reactions:
//!
//! 1. For each reaction row, explode `dbhit` (space-separated SEED
//!    reaction ids) into individual candidates.
//! 2. For each transport row, explode `rea` (comma-separated) similarly.
//! 3. Attach per-candidate `bitscore`, `pathway_status`, `status`,
//!    `is_complex`, `complex_status`.
//! 4. Compute an "effective" bitscore (max across all rows contributing
//!    to the same SEED id).
//! 5. Report a `CandidateTable` sorted by best-bitscore.
//!
//! Known gaps vs R:
//!
//! - `resolve_common_EC_conflicts` / `resolve_common_TC_conflicts`
//!   (IRanges overlap filter) is deferred.
//! - Topology-based pathway-status promotion (setting bitscore to
//!   `high_evi_rxn_bs` when the status is full/threshold/keyenzyme) is
//!   preserved.

use crate::reactions_tbl::{ReactionRow, TransporterRow};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct CandidateReaction {
    pub seed: String,
    pub best_bitscore: Option<f32>,
    pub best_stitle: String,
    pub best_qseqid: String,
    pub best_pident: Option<f32>,
    pub best_ec: String,
    pub best_tc: String,
    pub best_complex: String,
    pub best_pathway: String,
    pub pathway_status: String,
    pub status: String,
    pub is_complex: bool,
    pub complex_status: Option<u8>,
    pub exception: bool,
    pub weight: f32,
    pub from_transport: bool,
}

#[derive(Debug, Clone, Default)]
pub struct CandidateTable {
    pub rows: Vec<CandidateReaction>,
}

pub struct CandidateOptions {
    /// Reactions with a hit above this bitscore go to the draft model as
    /// "core" reactions.
    pub high_evi_rxn_bs: f32,
    /// Reactions below this bitscore are treated as "no hit".
    pub min_bs_for_core: f32,
}

impl Default for CandidateOptions {
    fn default() -> Self {
        Self { high_evi_rxn_bs: 200.0, min_bs_for_core: 50.0 }
    }
}

/// Assign a priority to a `pathway.status` string so we can always prefer
/// the "more supported" one across rows. Any of `full`, `threshold`,
/// `keyenzyme` trumps empty / other values; within the top tier,
/// `keyenzyme > full > threshold` matches gapseq's own scoring.
fn status_rank(s: &str) -> u8 {
    match s {
        "keyenzyme" => 3,
        "full" => 2,
        "threshold" => 1,
        _ => 0,
    }
}

fn promote_status(incoming: &str, current: &str) -> bool {
    status_rank(incoming) > status_rank(current)
}

pub fn build_candidates(
    reactions: &[ReactionRow],
    transporters: &[TransporterRow],
    opts: &CandidateOptions,
) -> CandidateTable {
    // Accumulator keyed by SEED id; track the best row so far.
    let mut acc: HashMap<String, CandidateReaction> = HashMap::new();

    let promote = |bs: Option<f32>, pwy_status: &str| -> Option<f32> {
        if matches!(pwy_status, "full" | "threshold" | "keyenzyme") {
            match bs {
                Some(v) => Some(v.max(opts.high_evi_rxn_bs)),
                None => Some(opts.high_evi_rxn_bs),
            }
        } else {
            bs
        }
    };

    // -- Reaction rows --
    for r in reactions {
        for seed in r.dbhit.split_whitespace().filter(|s| !s.is_empty()) {
            let bs = promote(r.bitscore, &r.pathway_status);
            let entry = acc.entry(seed.to_string()).or_insert_with(|| CandidateReaction {
                seed: seed.to_string(),
                best_bitscore: None,
                best_stitle: String::new(),
                best_qseqid: String::new(),
                best_pident: None,
                best_ec: String::new(),
                best_tc: String::new(),
                best_complex: String::new(),
                best_pathway: String::new(),
                pathway_status: String::new(),
                status: String::new(),
                is_complex: false,
                complex_status: None,
                exception: false,
                weight: 0.0,
                from_transport: false,
            });
            let beats = match (bs, entry.best_bitscore) {
                (Some(a), Some(b)) => a > b,
                (Some(_), None) => true,
                _ => false,
            };
            if beats || entry.best_bitscore.is_none() {
                entry.best_bitscore = bs;
                entry.best_stitle = r.stitle.clone();
                entry.best_qseqid = r.qseqid.clone();
                entry.best_pident = r.pident;
                entry.best_ec = r.ec.clone();
                entry.best_complex = r.complex.clone();
                entry.best_pathway = r.pathway.clone();
                entry.status = r.status.clone();
                entry.exception = r.exception;
            }
            // Best-pathway-status across all rows for this seed. A reaction
            // is "supported" if ANY of its rows has pathway.status in
            // full/threshold/keyenzyme. gapseq's R code takes an OR across
            // rows rather than a strict "best-row" copy.
            if promote_status(&r.pathway_status, &entry.pathway_status) {
                entry.pathway_status = r.pathway_status.clone();
            }
            // is_complex is TRUE if any row marks this seed as complex.
            if r.is_complex {
                entry.is_complex = true;
            }
            // complex_status is the max observed across rows (so a
            // partial-subunit row "1" wins over NA).
            if let Some(cs) = r.complex_status {
                entry.complex_status = match entry.complex_status {
                    None => Some(cs),
                    Some(old) => Some(old.max(cs)),
                };
            }
        }
    }

    // -- Transporter rows --
    for t in transporters {
        for seed in t.rea.split(',').filter(|s| !s.is_empty()) {
            let status = if t.bitscore.unwrap_or(0.0) >= opts.high_evi_rxn_bs {
                "good_blast"
            } else {
                "bad_blast"
            };
            let entry = acc.entry(seed.to_string()).or_insert_with(|| CandidateReaction {
                seed: seed.to_string(),
                best_bitscore: None,
                best_stitle: String::new(),
                best_qseqid: String::new(),
                best_pident: None,
                best_ec: String::new(),
                best_tc: String::new(),
                best_complex: String::new(),
                best_pathway: String::new(),
                pathway_status: String::new(),
                status: String::new(),
                is_complex: false,
                complex_status: None,
                exception: false,
                weight: 0.0,
                from_transport: true,
            });
            let beats = match (t.bitscore, entry.best_bitscore) {
                (Some(a), Some(b)) => a > b,
                (Some(_), None) => true,
                _ => false,
            };
            if beats || entry.best_bitscore.is_none() {
                entry.best_bitscore = t.bitscore;
                entry.best_stitle = t.stitle.clone();
                entry.best_qseqid = t.qseqid.clone();
                entry.best_pident = t.pident;
                entry.best_tc = t.tc.clone();
                entry.status = status.into();
                entry.from_transport = true;
            }
        }
    }

    // Compute weights (port of R's `weight := (bs - high_evi) * ((0.005 -
    // dummy) / (high_evi - min_bs_for_core)) + 0.005`).
    let dummy_weight = 100.0_f32;
    let slope = (0.005 - dummy_weight) / (opts.high_evi_rxn_bs - opts.min_bs_for_core);
    for cand in acc.values_mut() {
        let mut bs = cand.best_bitscore.unwrap_or(0.0);
        if bs > opts.high_evi_rxn_bs {
            bs = opts.high_evi_rxn_bs;
        }
        let mut w = (bs - opts.high_evi_rxn_bs) * slope + 0.005;
        if w < 0.005 {
            w = 0.005;
        }
        if w > dummy_weight {
            w = dummy_weight;
        }
        cand.weight = w;
    }

    let mut rows: Vec<CandidateReaction> = acc.into_values().collect();
    rows.sort_by(|a, b| {
        b.best_bitscore
            .partial_cmp(&a.best_bitscore)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.seed.cmp(&b.seed))
    });
    CandidateTable { rows }
}
