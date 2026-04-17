//! Futile-cycle detection for the candidate pool.
//!
//! A reaction participates in a thermodynamically infeasible cycle when it
//! can carry flux at capacity even with every exchange / demand reaction
//! closed. Such reactions, if added during gap-filling, would inflate the
//! flux solution without any real biochemistry behind them.
//!
//! Port of the recent upstream commit `cccbb6f0`: for each candidate
//! reaction, solve `max ±v[r]` on the full model with all EX/DM zeroed.
//! If the optimum saturates at `0.99 · max_flux`, drop the reaction.
//!
//! Parallelism: rayon `par_iter` over the candidate list. Each worker
//! builds its own tiny LP (one optimisation), so there's no shared state.

use crate::error::FillError;
use crate::fba::{fba, FbaOptions};
use crate::SolveStatus;
use gapsmith_core::Model;
use rayon::prelude::*;
use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct FutileOptions {
    /// Upper bound on net flux used by the LP. Matches `max_flux` = 1000
    /// in gapseq.
    pub max_flux: f64,
    /// Flux-saturation threshold (fraction of `max_flux`). Default `0.99`.
    pub saturation: f64,
}

impl Default for FutileOptions {
    fn default() -> Self {
        Self { max_flux: 1000.0, saturation: 0.99 }
    }
}

/// Detect futile-cycle reactions in `candidate_rxn_ids` on `full`.
///
/// Returns the subset of ids that saturate `sat · max_flux` in either the
/// forward or backward direction with every boundary reaction closed.
pub fn detect_futile_cycles(
    full: &Model,
    candidate_rxn_ids: &[String],
    opts: &FutileOptions,
) -> Result<HashSet<String>, FillError> {
    // Build a closed-boundary variant once, cloned by workers.
    let closed = close_boundaries(full);
    let threshold = opts.max_flux * opts.saturation;

    let bad: Vec<Option<String>> = candidate_rxn_ids
        .par_iter()
        .map(|rxn_id| {
            let idx = match closed.rxns.iter().position(|r| r.id.as_str() == rxn_id) {
                Some(i) => i,
                None => return None,
            };
            // Maximise v[idx] (net flux).
            let mut obj = vec![0.0f64; closed.rxn_count()];
            obj[idx] = 1.0;
            let fwd = fba(
                &closed,
                &FbaOptions { objective: Some(obj.clone()), maximise: true, max_flux: opts.max_flux },
            );
            let mut caught = matches!(&fwd, Ok(s) if matches!(s.status, SolveStatus::Optimal) && s.objective >= threshold);
            if !caught {
                // Try the reverse direction: maximise -v[idx].
                let mut obj_neg = obj.clone();
                obj_neg[idx] = -1.0;
                let bwd = fba(
                    &closed,
                    &FbaOptions { objective: Some(obj_neg), maximise: true, max_flux: opts.max_flux },
                );
                caught = matches!(&bwd, Ok(s) if matches!(s.status, SolveStatus::Optimal) && s.objective >= threshold);
            }
            caught.then(|| rxn_id.clone())
        })
        .collect();

    Ok(bad.into_iter().flatten().collect())
}

/// Clone `full` with every `is_exchange` / `EX_*` / `DM_*` reaction's bounds
/// pinned to zero. The remaining network must then carry flux only through
/// genuine metabolic recycling — any flux at capacity is a cycle.
fn close_boundaries(full: &Model) -> Model {
    let mut m = full.clone();
    for r in &mut m.rxns {
        let id = r.id.as_str();
        if r.is_exchange
            || r.is_biomass
            || id.starts_with("EX_")
            || id.starts_with("DM_")
            || id == "bio1"
            || id.starts_with("bio")
        {
            r.lb = 0.0;
            r.ub = 0.0;
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    use gapsmith_core::{CompartmentId, Metabolite, Reaction, StoichMatrix};

    #[test]
    fn detects_a_cycle() {
        // Two mets A, B and two reactions that together cycle: R1: A -> B,
        // R2: B -> A. Both reversible, no exchanges. This is a futile
        // cycle — running both at flux 1000 violates thermodynamics but
        // LP-feasible.
        let mut m = Model::new("cycle");
        m.mets.push(Metabolite::new("A", "A", CompartmentId::CYTOSOL));
        m.mets.push(Metabolite::new("B", "B", CompartmentId::CYTOSOL));
        m.rxns.push(Reaction::new("R1", "A<->B", -1000.0, 1000.0));
        m.rxns.push(Reaction::new("R2", "B<->A", -1000.0, 1000.0));
        m.s = StoichMatrix::from_triplets(
            2,
            2,
            vec![(0, 0, -1.0), (1, 0, 1.0), (1, 1, -1.0), (0, 1, 1.0)],
        );
        let bad = detect_futile_cycles(
            &m,
            &vec!["R1".into(), "R2".into()],
            &FutileOptions::default(),
        )
        .unwrap();
        assert_eq!(bad.len(), 2);
    }

    #[test]
    fn does_not_flag_non_cycled() {
        // A -> B with only the forward direction and an exchange. With
        // the exchange closed, no flux is possible — not a cycle.
        let mut m = Model::new("linear");
        m.mets.push(Metabolite::new("A", "A", CompartmentId::CYTOSOL));
        m.mets.push(Metabolite::new("B", "B", CompartmentId::CYTOSOL));
        let mut ex = Reaction::new("EX_A", "", -1000.0, 1000.0);
        ex.is_exchange = true;
        m.rxns.push(ex);
        m.rxns.push(Reaction::new("R1", "A->B", 0.0, 1000.0));
        m.s = StoichMatrix::from_triplets(
            2,
            2,
            vec![(0, 0, -1.0), (0, 1, -1.0), (1, 1, 1.0)],
        );
        let bad = detect_futile_cycles(
            &m,
            &vec!["R1".into()],
            &FutileOptions::default(),
        )
        .unwrap();
        assert!(bad.is_empty(), "got unexpected cycle: {bad:?}");
    }
}
