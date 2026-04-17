//! Candidate-reaction pool — the "full model" from which the gap-filler
//! pulls reactions. Port of `gapfill4.R:12-56` and gapseq's
//! `construct_full_model.R`.
//!
//! The pool is the union of the draft model's reactions plus every
//! gapseq-approved SEED reaction not already present, deduped by
//! stoichiometric hash. Each candidate has a pFBA weight derived from the
//! bitscore of its best homology hit (if any) — see [`rxn_weight`].

use crate::medium::{apply_medium, MediumEntry};
use gapsmith_core::{Model, RxnId};
use gapsmith_db::SeedRxnRow;
use gapsmith_draft::rxn_stoich_hash;
use std::collections::{HashMap, HashSet};

/// Bitscore threshold above which a reaction counts as "core" for
/// gap-filling. Matches gapseq `-b`, default `50`.
pub const DEFAULT_BCORE: f64 = 50.0;
/// Bitscore at which weight saturates at `0.005`. Matches gapseq's
/// `high.evi.rxn.BS`, default `200`.
pub const DEFAULT_HIGH_EVI: f64 = 200.0;
/// Weight assigned to reactions with no homology hit. Matches gapseq's
/// `dummy.weight`, default `100`.
pub const DEFAULT_DUMMY_WEIGHT: f64 = 100.0;

/// Per-reaction hit record produced by `find` / `find-transport`. Parsed
/// from a Reactions.tbl / Transporter.tbl's `(seed, bitscore)` pair.
#[derive(Debug, Clone)]
pub struct HitRecord {
    pub seed: String,
    pub bitscore: f64,
}

/// Lookup table mapping SEED rxn id → best bitscore observed.
#[derive(Debug, Clone, Default)]
pub struct RxnWeights {
    pub by_seed: HashMap<String, f64>,
    /// Threshold below which a hit is ignored as noise (the `bcore` flag).
    pub bcore: f64,
    /// Upper calibration anchor (the `high.evi.rxn.BS` flag).
    pub high_evi: f64,
    /// Weight when no hit is recorded.
    pub dummy_weight: f64,
}

impl RxnWeights {
    pub fn new() -> Self {
        Self {
            by_seed: HashMap::new(),
            bcore: DEFAULT_BCORE,
            high_evi: DEFAULT_HIGH_EVI,
            dummy_weight: DEFAULT_DUMMY_WEIGHT,
        }
    }

    /// Insert `seed → bitscore`, keeping the maximum on collisions.
    pub fn update(&mut self, seed: &str, bitscore: f64) {
        let cur = self.by_seed.entry(seed.to_string()).or_insert(f64::MIN);
        if bitscore > *cur {
            *cur = bitscore;
        }
    }

    pub fn bitscore(&self, seed: &str) -> Option<f64> {
        self.by_seed.get(seed).copied()
    }

    /// Port of `prepare_candidate_reaction_tables.R:223-228`.
    pub fn weight(&self, seed: &str) -> f64 {
        let bs = match self.by_seed.get(seed) {
            Some(&b) if b >= self.bcore => b,
            _ => return self.dummy_weight,
        };
        rxn_weight(bs, self.high_evi, self.bcore, self.dummy_weight)
    }

    /// Is this SEED reaction "core" for gap-fill purposes? Core == bitscore
    /// at least `bcore`.
    pub fn is_core(&self, seed: &str) -> bool {
        self.by_seed.get(seed).is_some_and(|&b| b >= self.bcore)
    }
}

/// Slope-scaled weight: bitscore in `[bcore, high_evi]` maps linearly to
/// `[dummy_weight, 0.005]`, clamped outside. Higher bitscore → smaller
/// weight → pFBA prefers this reaction.
pub fn rxn_weight(bitscore: f64, high_evi: f64, bcore: f64, dummy_weight: f64) -> f64 {
    let bs_clamped = bitscore.min(high_evi);
    let slope = (0.005 - dummy_weight) / (high_evi - bcore);
    let w = (bs_clamped - high_evi) * slope + 0.005;
    w.clamp(0.005, dummy_weight)
}

/// Bulk-load a Reactions.tbl / Transporter.tbl-style TSV and record the
/// best bitscore per SEED id.
///
/// We take the `dbhit` column as a space-separated list of SEED ids and
/// attribute the row's bitscore to each. This matches how `gapfill4.R`
/// intersects its `rxn.weights` lookup.
pub fn read_weights_from_reactions_tbl(
    path: &std::path::Path,
    bcore: f64,
    high_evi: f64,
    dummy_weight: f64,
) -> std::io::Result<RxnWeights> {
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(path)?;
    let r = BufReader::new(f);
    let mut lines = Vec::<String>::new();
    for line in r.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        lines.push(line);
    }
    let mut w = RxnWeights { bcore, high_evi, dummy_weight, ..Default::default() };
    if lines.is_empty() {
        return Ok(w);
    }
    let header: Vec<&str> = lines[0].split('\t').collect();
    let idx = |name: &str| header.iter().position(|c| *c == name);
    let i_dbhit = idx("dbhit");
    let i_rxn = idx("rxn");
    let i_bitscore = idx("bitscore");
    let (Some(i_dbhit), Some(i_bitscore)) = (i_dbhit, i_bitscore) else {
        return Ok(w);
    };
    for row in &lines[1..] {
        let cols: Vec<&str> = row.split('\t').collect();
        if cols.len() <= i_dbhit.max(i_bitscore) {
            continue;
        }
        let bs: f64 = match cols[i_bitscore].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        for seed in cols[i_dbhit].split_whitespace().filter(|s| !s.is_empty()) {
            w.update(seed, bs);
        }
        // If there's no dbhit but the `rxn` column is itself a SEED id
        // (transporter tables sometimes use it that way), record under
        // that key.
        if cols[i_dbhit].trim().is_empty() {
            if let Some(i_rxn) = i_rxn {
                if let Some(rxn) = cols.get(i_rxn) {
                    if rxn.starts_with("rxn") {
                        w.update(rxn, bs);
                    }
                }
            }
        }
    }
    Ok(w)
}

/// Build the "full model" = draft + every approved SEED candidate not
/// already present. Candidates are deduped against each other and against
/// the draft by stoichiometric hash.
///
/// Returns the extended model alongside the set of seed ids that ended up
/// in the pool (in insertion order) — the gapfiller uses this to decide
/// which columns are candidates vs pre-existing.
pub fn build_full_model(
    draft: &Model,
    seed_rxns: &[SeedRxnRow],
    weights: &RxnWeights,
) -> Result<(Model, Vec<String>), String> {
    let mut full = draft.clone();

    // Collect seed ids already in the draft; their stoich hashes too.
    let present_seeds: HashSet<String> = full
        .rxns
        .iter()
        .map(|r| strip_compartment(r.id.as_str()).to_string())
        .collect();
    let mut taken_hashes: HashSet<String> = HashSet::new();
    for r in seed_rxns
        .iter()
        .filter(|r| present_seeds.contains(r.id.as_str()))
    {
        // Prefer the hash cached at load time; fall back to recompute
        // for rows constructed outside of `load_seed_reactions` (tests).
        let h = match r.stoich_hash.as_deref() {
            Some(h) => h.to_string(),
            None => rxn_stoich_hash(&r.stoichiometry, r.reversibility.as_str())
                .map_err(|e| e.to_string())?,
        };
        taken_hashes.insert(h);
    }

    // Filter the approved SEED rows to ones not in the draft, dedup by
    // hash with core preference (core reactions win ties).
    let mut candidates: Vec<&SeedRxnRow> = seed_rxns
        .iter()
        .filter(|r| r.gapseq_status.is_usable())
        .filter(|r| !present_seeds.contains(r.id.as_str()))
        .collect();

    // Sort: core first (higher priority on hash ties), then by id for
    // stability.
    candidates.sort_by(|a, b| {
        let ac = weights.is_core(a.id.as_str());
        let bc = weights.is_core(b.id.as_str());
        bc.cmp(&ac).then_with(|| a.id.as_str().cmp(b.id.as_str()))
    });

    let mut added_seeds = Vec::<String>::new();
    for row in candidates {
        let h = match row.stoich_hash.as_deref() {
            Some(h) => h.to_string(),
            None => match rxn_stoich_hash(&row.stoichiometry, row.reversibility.as_str()) {
                Ok(h) => h,
                Err(_) => continue,
            },
        };
        if !taken_hashes.insert(h) {
            continue;
        }
        gapsmith_draft::builder::add_seed_reaction(&mut full, row, None);
        added_seeds.push(row.id.as_str().to_string());
    }
    gapsmith_draft::builder::rebuild_s_matrix(&mut full);

    Ok((full, added_seeds))
}

/// Strip the `_c0` / `_e0` / `_p0` compartment suffix from a reaction id.
/// Matches `gsub("_.*","", react_id)` used throughout the R code.
pub fn strip_compartment(rxn_id: &str) -> &str {
    match rxn_id.rfind('_') {
        Some(i) if rxn_id[i + 1..].starts_with(['c', 'e', 'p']) => &rxn_id[..i],
        _ => rxn_id,
    }
}

/// Synchronize the draft's exchange bounds onto the full model. The
/// gapfill4 R code calls `sync_full_mod` (lines 305-344) which does two
/// things we handle here: (1) every reaction unique to the draft is
/// already in `full` because we started by cloning the draft; and (2)
/// exchange bounds from the draft — which the caller may have already
/// mutated via `apply_medium` — need to be mirrored onto `full`.
///
/// Because `full` is built by cloning the draft and only APPENDING
/// candidate reactions, point (1) is free and point (2) just means
/// `apply_medium` on `full` after construction sets the same constraints.
pub fn apply_medium_to_full(full: &mut Model, medium: &[MediumEntry]) {
    apply_medium(full, medium, 1.0, 1000.0);
}

/// Build a parallel `Vec<f64>` of pFBA weights, one per reaction in
/// `full.rxns`. Matches the cost-coefficient table assembled by
/// `gapfill4.R:82-86`:
///
/// - Non-candidate and exchange / demand / biomass reactions get `pres_w`
///   (default `1e-5`) — effectively free to use.
/// - Reactions also present in the draft get `pres_w` as well.
/// - Candidate reactions get the slope-scaled weight from [`RxnWeights`].
///
/// `draft_rxn_ids` is the set of reaction ids already in the draft *before*
/// the pool expansion; the caller typically passes
/// `draft.rxns.iter().map(|r| r.id.to_string()).collect()`.
pub fn pfba_weights(
    full: &Model,
    draft_rxn_ids: &HashSet<RxnId>,
    weights: &RxnWeights,
    pres_w: f64,
) -> Vec<f64> {
    full.rxns
        .iter()
        .map(|r| {
            let id_str = r.id.as_str();
            if draft_rxn_ids.contains(&r.id)
                || r.is_exchange
                || id_str.starts_with("EX_")
                || id_str.starts_with("DM_")
                || r.is_biomass
                || id_str.starts_with("bio")
            {
                return pres_w;
            }
            let seed = strip_compartment(id_str);
            weights.weight(seed)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weight_slope_matches_r() {
        // With defaults bcore=50, high_evi=200, dummy=100:
        assert!((rxn_weight(200.0, 200.0, 50.0, 100.0) - 0.005).abs() < 1e-9);
        assert!((rxn_weight(50.0, 200.0, 50.0, 100.0) - 100.0).abs() < 1e-3);
        // Below bcore → clamp to dummy.
        assert_eq!(rxn_weight(0.0, 200.0, 50.0, 100.0), 100.0);
    }

    #[test]
    fn strip_compartment_cases() {
        assert_eq!(strip_compartment("rxn00001_c0"), "rxn00001");
        assert_eq!(strip_compartment("EX_cpd00027_e0"), "EX_cpd00027");
        assert_eq!(strip_compartment("bio1"), "bio1");
    }
}
