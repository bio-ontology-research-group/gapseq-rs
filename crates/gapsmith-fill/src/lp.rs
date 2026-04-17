//! Split-flux LP encoding used by FBA, pFBA, and the gap-filler.
//!
//! # Why split-flux?
//!
//! Cobrar's `pfbaHeuristic` minimises `Σ w_r · |v_r|`. Because `|v_r|` is
//! non-linear, we split each reaction into a forward and backward part:
//!
//! ```text
//!   v_r   = vp_r − vn_r           (net flux)
//!   |v_r| = vp_r + vn_r            (when at most one is non-zero at optimum)
//! ```
//!
//! Both `vp_r, vn_r ≥ 0`. This yields a pure LP (no integer / quadratic
//! terms) that HiGHS / CBC handle efficiently.
//!
//! # Bound translation
//!
//! Reaction bounds `[lb, ub]` on the net flux become split-variable bounds:
//!
//! | lb, ub                | vp upper bound   | vn upper bound   |
//! |-----------------------|------------------|------------------|
//! | `lb ≥ 0, ub ≥ 0`      | `ub`             | `0`              |
//! | `lb ≤ 0, ub ≥ 0`      | `ub`             | `−lb`            |
//! | `lb ≤ 0, ub ≤ 0`      | `0`              | `−lb`            |
//!
//! For `lb > 0` we additionally require `vp_r ≥ lb`; for `ub < 0` we require
//! `vn_r ≥ −ub`. In both corner cases the direction is forced by the
//! corresponding opposite variable being fixed to zero.

use gapsmith_core::Model;

/// A split-flux LP skeleton carrying the per-reaction variable bounds.
///
/// The actual variables / constraints are created lazily inside each solver
/// entry point; this struct just holds the bound translation so the FBA and
/// pFBA builders can share the logic.
#[derive(Debug, Clone)]
pub struct SplitFluxLp<'m> {
    pub model: &'m Model,
    /// Upper bound on `vp_r` for each reaction, in column order.
    pub vp_ub: Vec<f64>,
    /// Lower bound on `vp_r` (non-zero only when `lb > 0`).
    pub vp_lb: Vec<f64>,
    /// Upper bound on `vn_r`.
    pub vn_ub: Vec<f64>,
    /// Lower bound on `vn_r` (non-zero only when `ub < 0`).
    pub vn_lb: Vec<f64>,
}

impl<'m> SplitFluxLp<'m> {
    pub fn from_model(model: &'m Model) -> Self {
        let n = model.rxn_count();
        let mut vp_ub = Vec::with_capacity(n);
        let mut vp_lb = Vec::with_capacity(n);
        let mut vn_ub = Vec::with_capacity(n);
        let mut vn_lb = Vec::with_capacity(n);

        for r in &model.rxns {
            let lb = r.lb;
            let ub = r.ub;

            // vp covers the positive half; vn covers the negative half.
            vp_ub.push(ub.max(0.0));
            vp_lb.push(lb.max(0.0));
            vn_ub.push((-lb).max(0.0));
            vn_lb.push((-ub).max(0.0));
        }

        Self { model, vp_ub, vp_lb, vn_ub, vn_lb }
    }

    /// Reassemble net fluxes `v_r = vp_r − vn_r` from separate variable value
    /// vectors.
    pub fn net_flux(&self, vp_vals: &[f64], vn_vals: &[f64]) -> Vec<f64> {
        assert_eq!(vp_vals.len(), vn_vals.len());
        vp_vals.iter().zip(vn_vals).map(|(p, n)| p - n).collect()
    }
}

#[cfg(test)]
pub(crate) fn toy_model() -> Model {
    use gapsmith_core::{CompartmentId, Metabolite, Reaction, StoichMatrix};
    let mut m = Model::new("toy");
    m.mets.push(Metabolite::new("cpdA", "A", CompartmentId::CYTOSOL));
    m.mets.push(Metabolite::new("cpdB", "B", CompartmentId::CYTOSOL));
    m.mets.push(Metabolite::new("cpdC", "C", CompartmentId::CYTOSOL));

    let mut ex_a = Reaction::new("EX_A", "EX A", -10.0, 1000.0);
    ex_a.is_exchange = true;
    m.rxns.push(ex_a);
    m.rxns.push(Reaction::new("R_AB", "A -> B", 0.0, 1000.0));
    let mut r_bc = Reaction::new("R_BC", "B -> C", 0.0, 1000.0);
    r_bc.obj_coef = 1.0;
    m.rxns.push(r_bc);
    let mut ex_c = Reaction::new("EX_C", "EX C", -1000.0, 1000.0);
    ex_c.is_exchange = true;
    m.rxns.push(ex_c);

    m.s = StoichMatrix::from_triplets(
        3,
        4,
        vec![
            (0, 0, -1.0),
            (0, 1, -1.0),
            (1, 1, 1.0),
            (1, 2, -1.0),
            (2, 2, 1.0),
            (2, 3, -1.0),
        ],
    );
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bound_translation_forward() {
        let m = toy_model();
        let lp = SplitFluxLp::from_model(&m);
        // EX_A has lb=-10, ub=1000
        assert!((lp.vp_ub[0] - 1000.0).abs() < 1e-9);
        assert!((lp.vn_ub[0] - 10.0).abs() < 1e-9);
        // R_AB has lb=0, ub=1000
        assert!((lp.vp_ub[1] - 1000.0).abs() < 1e-9);
        assert!(lp.vn_ub[1].abs() < 1e-9);
    }
}
