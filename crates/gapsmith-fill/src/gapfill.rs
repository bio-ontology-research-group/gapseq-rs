//! Single-iteration gap-fill. Port of `src/gapfill4.R:1–303`.
//!
//! Workflow, given a draft model + a pre-assembled full model (draft
//! + every candidate SEED reaction):
//!
//! 1. Set the biomass-reaction lower bound of `full` to `min_growth` so
//!    pFBA is forced to find a feasible pathway.
//! 2. Run [`crate::pfba_heuristic`] with per-reaction weights.
//! 3. Extract candidate reactions (i.e. reactions in `full` but not in
//!    the draft) whose flux is non-zero at the optimum.
//! 4. Add them to the draft; rebuild the S matrix.
//! 5. KO loop — port of `gapfill4.R:247–280`. For each non-core added
//!    reaction, zero its bounds, solve FBA, and drop the reaction if
//!    growth remains ≥ `min_growth` (i.e. reaction was non-essential).

use crate::error::FillError;
use crate::fba::{fba, FbaOptions};
use crate::pfba::{pfba_heuristic, PfbaHeuristicOptions};
use crate::pool::{pfba_weights, strip_compartment, RxnWeights};
use crate::SolveStatus;
use gapsmith_core::{Model, RxnId};
use gapsmith_db::SeedRxnRow;
use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct GapfillOptions {
    pub min_growth: f64,
    /// Weight for reactions already present in the draft + exchanges /
    /// biomass. Matches gapseq's `presRxnWeight` = `1e-5`.
    pub pres_weight: f64,
    /// gs.origin marker stamped on added reactions (port of the
    /// `gs.origin` R parameter).
    pub gs_origin: Option<i8>,
    pub pfba_heuristic: PfbaHeuristicOptions,
}

impl GapfillOptions {
    pub fn new(min_growth: f64, n_rxns_full: usize) -> Self {
        // Weights filled in by the caller after knowing the full model.
        Self {
            min_growth,
            pres_weight: 1e-5,
            gs_origin: Some(1),
            pfba_heuristic: PfbaHeuristicOptions::new(vec![1.0; n_rxns_full], min_growth),
        }
    }
}

#[derive(Debug, Clone)]
pub struct GapfillReport {
    pub model: Model,
    pub rxns_added: Vec<String>,
    pub rxns_added_core: Vec<String>,
    pub growth_rate: f64,
    pub status: SolveStatus,
}

/// Run a single gap-fill iteration.
///
/// `draft` is the input model (will be modified in the returned report).
/// `full` is the pre-expanded candidate-pool model from
/// [`crate::build_full_model`]. Both models must carry consistent biomass
/// and exchange reactions — apply the same medium to both before calling.
/// `weights` maps SEED id → best bitscore for core/weight classification.
/// `seed_rxns` is the full SEED DB used to look up stoichiometry for
/// candidates being added to the draft.
pub fn gapfill4(
    draft: &Model,
    full: &Model,
    weights: &RxnWeights,
    seed_rxns: &[SeedRxnRow],
    opts: &GapfillOptions,
) -> Result<GapfillReport, FillError> {
    // Build the pFBA weight vector aligned to `full.rxns`.
    let draft_ids: HashSet<RxnId> = draft.rxns.iter().map(|r| r.id.clone()).collect();
    let mut weight_vec = pfba_weights(full, &draft_ids, weights, opts.pres_weight);

    // Enforce min growth by pinning the biomass lb on `full`.
    let mut full_forced = full.clone();
    for r in &mut full_forced.rxns {
        if r.obj_coef.abs() >= f64::EPSILON {
            r.lb = opts.min_growth.max(r.lb);
        }
    }

    // Compute pFBA.
    let mut heur = opts.pfba_heuristic.clone();
    heur.weights = std::mem::take(&mut weight_vec);
    heur.min_growth = opts.min_growth;
    let sol = pfba_heuristic(&full_forced, &heur)?;
    if !matches!(sol.status, SolveStatus::Optimal) {
        return Ok(GapfillReport {
            model: draft.clone(),
            rxns_added: Vec::new(),
            rxns_added_core: Vec::new(),
            growth_rate: 0.0,
            status: sol.status,
        });
    }

    // Identify utilized candidate reactions: rxns in `full` but not in the
    // draft, with non-zero flux.
    let tol = 1e-6;
    let utilized: Vec<(String, &SeedRxnRow)> = full
        .rxns
        .iter()
        .zip(sol.fluxes.iter())
        .filter(|(r, _)| !draft_ids.contains(&r.id))
        .filter(|(_, v)| v.abs() > tol)
        .filter_map(|(r, _)| {
            let seed = strip_compartment(r.id.as_str()).to_string();
            seed_rxns
                .iter()
                .find(|sr| sr.id.as_str() == seed)
                .map(|sr| (seed, sr))
        })
        .collect();

    if utilized.is_empty() {
        // pFBA feasible but no candidate reactions used — draft already
        // grows on this medium. Return it unchanged.
        return Ok(GapfillReport {
            model: draft.clone(),
            rxns_added: Vec::new(),
            rxns_added_core: Vec::new(),
            growth_rate: sol.growth,
            status: SolveStatus::Optimal,
        });
    }

    // Add utilized candidates to a clone of the draft.
    let mut filled = draft.clone();
    let mut added_ids = Vec::<String>::new();
    let mut added_core = Vec::<String>::new();
    for (seed, row) in &utilized {
        gapsmith_draft::builder::add_seed_reaction(&mut filled, row, opts.gs_origin);
        added_ids.push(format!("{seed}_c0"));
        if weights.is_core(seed) {
            added_core.push(seed.clone());
        }
    }
    gapsmith_draft::builder::rebuild_s_matrix(&mut filled);

    // Verify the filled draft meets min_growth via plain FBA.
    let fba_sol = fba(
        &filled,
        &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
    )?;
    if !matches!(fba_sol.status, SolveStatus::Optimal) || fba_sol.objective < opts.min_growth {
        tracing::warn!(
            growth = fba_sol.objective,
            min_growth = opts.min_growth,
            "post-fill FBA failed to hit growth floor — returning unfilled draft"
        );
        return Ok(GapfillReport {
            model: draft.clone(),
            rxns_added: Vec::new(),
            rxns_added_core: Vec::new(),
            growth_rate: 0.0,
            status: fba_sol.status,
        });
    }

    // KO loop: for each non-core added reaction, zero its bounds and
    // check whether growth remains ≥ min_growth. Drop if non-essential.
    //
    // R code orders `ko.dt[order(core, -dummy.weight)]` — non-core first,
    // highest-weight first (most "expensive" to keep).
    let mut ko_order: Vec<usize> = (0..added_ids.len())
        .filter(|i| {
            let seed = strip_compartment(&added_ids[*i]);
            !weights.is_core(seed)
        })
        .collect();
    ko_order.sort_by(|&a, &b| {
        let wa = weights.weight(strip_compartment(&added_ids[a]));
        let wb = weights.weight(strip_compartment(&added_ids[b]));
        wb.partial_cmp(&wa).unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut removed: HashSet<String> = HashSet::new();
    for idx in ko_order {
        let rxn_id = added_ids[idx].clone();
        if removed.contains(&rxn_id) {
            continue;
        }
        // Find the index in `filled.rxns`.
        let pos = match filled.rxns.iter().position(|r| r.id.as_str() == rxn_id) {
            Some(p) => p,
            None => continue,
        };
        let bu_lb = filled.rxns[pos].lb;
        let bu_ub = filled.rxns[pos].ub;
        filled.rxns[pos].lb = 0.0;
        filled.rxns[pos].ub = 0.0;
        let probe = fba(
            &filled,
            &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
        )?;
        let still_grows = matches!(probe.status, SolveStatus::Optimal)
            && probe.objective >= opts.min_growth;
        if still_grows {
            removed.insert(rxn_id);
            // Keep the zeroed bounds — reaction is effectively gone. We'll
            // delete it from the model after the loop.
        } else {
            filled.rxns[pos].lb = bu_lb;
            filled.rxns[pos].ub = bu_ub;
        }
    }

    if !removed.is_empty() {
        drop_reactions(&mut filled, &removed);
    }

    let kept_ids: Vec<String> = added_ids
        .into_iter()
        .filter(|id| !removed.contains(id))
        .collect();

    let final_fba = fba(
        &filled,
        &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
    )?;

    Ok(GapfillReport {
        model: filled,
        rxns_added: kept_ids,
        rxns_added_core: added_core,
        growth_rate: final_fba.objective,
        status: final_fba.status,
    })
}

/// Delete reactions by id from a model, preserving the S matrix.
pub fn drop_reactions(model: &mut Model, rxn_ids: &HashSet<String>) {
    if rxn_ids.is_empty() {
        return;
    }
    let keep_mask: Vec<bool> = model
        .rxns
        .iter()
        .map(|r| !rxn_ids.contains(r.id.as_str()))
        .collect();

    // Build new S from triplets, remapping column indices.
    let n_mets = model.mets.len();
    let n_rxns_old = model.rxns.len();
    let mut new_col_ix: Vec<Option<usize>> = vec![None; n_rxns_old];
    let mut next = 0usize;
    for (i, keep) in keep_mask.iter().enumerate() {
        if *keep {
            new_col_ix[i] = Some(next);
            next += 1;
        }
    }
    let n_rxns_new = next;

    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(model.s.nnz());
    for c in 0..n_rxns_old {
        if let Some(new_c) = new_col_ix[c] {
            for (row, v) in model.s.column(c) {
                triplets.push((row, new_c, v));
            }
        }
    }

    // Filter reactions.
    let mut new_rxns = Vec::with_capacity(n_rxns_new);
    for (i, r) in model.rxns.iter().enumerate() {
        if keep_mask[i] {
            new_rxns.push(r.clone());
        }
    }
    model.rxns = new_rxns;
    model.s = gapsmith_core::StoichMatrix::from_triplets(n_mets, n_rxns_new, triplets);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drop_reactions_preserves_nnz() {
        use gapsmith_core::{CompartmentId, Metabolite, Reaction, StoichMatrix};
        let mut m = Model::new("x");
        m.mets.push(Metabolite::new("A", "A", CompartmentId::CYTOSOL));
        m.mets.push(Metabolite::new("B", "B", CompartmentId::CYTOSOL));
        m.rxns.push(Reaction::new("r1", "", 0.0, 1.0));
        m.rxns.push(Reaction::new("r2", "", 0.0, 1.0));
        m.rxns.push(Reaction::new("r3", "", 0.0, 1.0));
        m.s = StoichMatrix::from_triplets(
            2,
            3,
            vec![(0, 0, -1.0), (1, 1, 1.0), (0, 2, -1.0), (1, 2, 1.0)],
        );
        drop_reactions(&mut m, &HashSet::from(["r2".to_string()]));
        assert_eq!(m.rxn_count(), 2);
        assert_eq!(m.rxns[0].id.as_str(), "r1");
        assert_eq!(m.rxns[1].id.as_str(), "r3");
        // nnz was 4, removed 1 (from r2). 3 remaining.
        assert_eq!(m.s.nnz(), 3);
    }
}
