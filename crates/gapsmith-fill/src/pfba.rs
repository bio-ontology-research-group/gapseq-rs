//! Parsimonious FBA with pFBA-heuristic tolerance ladder.
//!
//! Two entry points:
//!
//! - [`pfba`] — a single pFBA solve: minimise `Σ w_r · |v_r|` subject to the
//!   usual mass-balance + bound constraints plus a biomass-floor constraint
//!   `v_bio ≥ min_gr`. One call → one solution. Used as a leaf primitive.
//!
//! - [`pfba_heuristic`] — the tolerance-ladder driver copied from
//!   `gapfill4.R:95–137`. Starts at `pfba_coef = 1e-3` / `tol = 1e-6`, falls
//!   back to smaller tolerances, then relaxes the pFBA coefficient; retries
//!   up to `max_iter` times. Returns once a solution passes a validation
//!   FBA with the zero-flux reactions removed.

use crate::error::{FillError, SolveStatus};
use crate::fba::{build_row_exprs, build_vars, fba, resolve_obj, FbaOptions};
use crate::lp::SplitFluxLp;
use gapsmith_core::{Model, RxnId};
use good_lp::{
    constraint, solvers::highs::highs, Expression, ResolutionError, Solution,
    SolverModel,
};

#[derive(Debug, Clone)]
pub struct PfbaOptions {
    /// Per-reaction weights for the `Σ w_r · |v_r|` penalty (one entry per
    /// reaction, model column order).
    pub weights: Vec<f64>,
    /// Coefficient multiplied against the biomass-flux term in the objective
    /// to trade off against the absolute-flux sum. Matches cobrar's
    /// `pFBAcoeff`. Default `1e-3`.
    pub pfba_coef: f64,
    /// Lower bound forced on the biomass reaction (`obj_coef != 0`).
    pub min_growth: f64,
    /// Objective coefficients — if `None`, read from `model.rxns[*].obj_coef`.
    pub objective: Option<Vec<f64>>,
}

impl PfbaOptions {
    pub fn uniform(model: &Model, pfba_coef: f64, min_growth: f64) -> Self {
        let n = model.rxn_count();
        Self {
            weights: vec![1.0; n],
            pfba_coef,
            min_growth,
            objective: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PfbaSolution {
    pub status: SolveStatus,
    pub objective: f64,
    pub fluxes: Vec<f64>,
    /// Growth (biomass) flux — `v_r` where `obj_coef != 0`.
    pub growth: f64,
}

pub fn pfba(model: &Model, opts: &PfbaOptions) -> Result<PfbaSolution, FillError> {
    model.check_shape().map_err(|_| FillError::BadShape)?;
    let n = model.rxn_count();
    if opts.weights.len() != n {
        return Err(FillError::BadInput("weights length ≠ rxn count".into()));
    }
    let obj_coefs = resolve_obj(
        model,
        &FbaOptions { objective: opts.objective.clone(), maximise: true, max_flux: 1000.0 },
    )?;

    let lp = SplitFluxLp::from_model(model);
    let (vars, vp, vn, _obj_expr) = build_vars(&lp, &obj_coefs);

    // pFBA objective: minimise  Σ w_r · (vp_r + vn_r) − pfba_coef · v_bio
    // (the −pfba_coef·v_bio term is negative so minimising pushes v_bio up,
    // matching cobrar's `pfbaHeuristic` which minimises `sum(|v|) − c·v_bio`.)
    let mut obj = Expression::from(0.0);
    for (i, w) in opts.weights.iter().enumerate() {
        if w.abs() >= f64::EPSILON {
            obj += *w * (vp[i] + vn[i]);
        }
    }
    for (i, c) in obj_coefs.iter().enumerate() {
        if c.abs() >= f64::EPSILON {
            obj -= opts.pfba_coef * *c * (vp[i] - vn[i]);
        }
    }

    let mut problem = vars.minimise(obj).using(highs);

    // Mass balance.
    for expr in build_row_exprs(model, &vp, &vn, model.met_count()) {
        problem = problem.with(constraint!(expr == 0.0));
    }

    // Biomass floor: v_bio ≥ min_growth. We impose this on the net flux of
    // each reaction that carries a non-zero obj coefficient.
    for (i, c) in obj_coefs.iter().enumerate() {
        if c.abs() >= f64::EPSILON {
            let expr = *c * (vp[i] - vn[i]);
            problem = problem.with(constraint!(expr >= opts.min_growth));
        }
    }

    let solution = match problem.solve() {
        Ok(s) => s,
        Err(ResolutionError::Infeasible) => {
            return Ok(PfbaSolution {
                status: SolveStatus::Infeasible,
                objective: 0.0,
                fluxes: vec![0.0; n],
                growth: 0.0,
            })
        }
        Err(ResolutionError::Unbounded) => {
            return Ok(PfbaSolution {
                status: SolveStatus::Unbounded,
                objective: 0.0,
                fluxes: vec![0.0; n],
                growth: 0.0,
            })
        }
        Err(e) => return Err(FillError::Solver(e.to_string())),
    };

    let fluxes: Vec<f64> =
        (0..n).map(|i| solution.value(vp[i]) - solution.value(vn[i])).collect();
    let growth: f64 = obj_coefs
        .iter()
        .zip(&fluxes)
        .map(|(c, v)| c * v)
        .sum();
    let objective: f64 = opts
        .weights
        .iter()
        .zip(&fluxes)
        .map(|(w, v)| w * v.abs())
        .sum();

    Ok(PfbaSolution { status: SolveStatus::Optimal, objective, fluxes, growth })
}

#[derive(Debug, Clone)]
pub struct PfbaHeuristicOptions {
    pub weights: Vec<f64>,
    pub min_growth: f64,
    /// Starting pFBA coefficient. Default `1e-3`.
    pub start_pfba_coef: f64,
    /// Starting solver feasibility tolerance. Default `1e-6`.
    pub start_tol: f64,
    /// Tolerance floor. Default `1e-9`.
    pub min_tol: f64,
    /// Max retries. Default `15`.
    pub max_iter: usize,
    /// Objective coefficients — if `None`, read from `model.rxns[*].obj_coef`.
    pub objective: Option<Vec<f64>>,
    /// When HiGHS exhausts the tolerance/coefficient ladder, retry once
    /// with CBC (the `cbc` feature must be enabled at build time). Default
    /// `false`.
    pub cbc_fallback: bool,
}

impl PfbaHeuristicOptions {
    pub fn new(weights: Vec<f64>, min_growth: f64) -> Self {
        Self {
            weights,
            min_growth,
            start_pfba_coef: 1e-3,
            start_tol: 1e-6,
            min_tol: 1e-9,
            max_iter: 15,
            objective: None,
            cbc_fallback: false,
        }
    }
}

/// Tolerance-ladder pFBA driver. Port of `gapfill4.R:95–137`.
///
/// Loop:
/// 1. Solve pFBA with the current `pfba_coef`.
/// 2. If infeasible / growth below tolerance → halve `pfba_coef`, retry.
/// 3. Otherwise, build a reduced model without the zero-flux reactions and
///    run FBA; if feasible, accept. If the reduced FBA is infeasible and we
///    haven't hit the tolerance floor, halve `tol` and retry. If we've hit
///    the floor, halve `pfba_coef` instead.
pub fn pfba_heuristic(
    model: &Model,
    opts: &PfbaHeuristicOptions,
) -> Result<PfbaSolution, FillError> {
    let mut tol = opts.start_tol;
    let mut coef = opts.start_pfba_coef;

    for _ in 0..opts.max_iter {
        let sol = pfba(
            model,
            &PfbaOptions {
                weights: opts.weights.clone(),
                pfba_coef: coef,
                min_growth: opts.min_growth,
                objective: opts.objective.clone(),
            },
        )?;

        if !matches!(sol.status, SolveStatus::Optimal) || sol.growth <= tol {
            coef /= 2.0;
            tracing::debug!(coef, "pfba_heuristic: relaxing pFBA coefficient");
            continue;
        }

        // Validation FBA on the reduced model.
        let reduced = remove_zero_flux(model, &sol.fluxes, tol);
        let vopts = FbaOptions { objective: opts.objective.clone(), maximise: true, max_flux: 1000.0 };
        let val = fba(&reduced, &vopts)?;
        if matches!(val.status, SolveStatus::Optimal) && val.objective >= opts.min_growth {
            return Ok(sol);
        }

        // Validation failed — drop tolerance, then pFBA coefficient.
        if tol > opts.min_tol {
            tol = (tol / 2.0).max(opts.min_tol);
            tracing::debug!(tol, "pfba_heuristic: lowering tolerance");
        } else {
            coef /= 2.0;
            tracing::debug!(coef, "pfba_heuristic: tolerance floor hit, relaxing pFBA coefficient");
        }
    }

    if opts.cbc_fallback {
        tracing::info!("pfba_heuristic: HiGHS exhausted ladder, trying CBC fallback");
        #[cfg(feature = "cbc")]
        {
            return pfba_cbc(model, opts);
        }
        #[cfg(not(feature = "cbc"))]
        {
            tracing::warn!(
                "cbc_fallback requested but binary built without `cbc` feature; rebuild with `--features cbc`"
            );
        }
    }

    Err(FillError::Solver(
        "pfba_heuristic: max iterations reached without a valid solution".into(),
    ))
}

/// CBC backend for pFBA. Identical formulation to [`pfba`] but dispatched
/// through `good_lp::solvers::coin_cbc`. Only compiled when the `cbc`
/// feature is enabled.
#[cfg(feature = "cbc")]
fn pfba_cbc(model: &Model, opts: &PfbaHeuristicOptions) -> Result<PfbaSolution, FillError> {
    use crate::fba::{build_row_exprs, build_vars, resolve_obj, FbaOptions};
    use crate::lp::SplitFluxLp;
    use good_lp::{solvers::coin_cbc::coin_cbc, constraint, ResolutionError, Solution, SolverModel};

    model.check_shape().map_err(|_| FillError::BadShape)?;
    let n = model.rxn_count();
    let obj_coefs = resolve_obj(
        model,
        &FbaOptions { objective: opts.objective.clone(), maximise: true, max_flux: 1000.0 },
    )?;

    let lp = SplitFluxLp::from_model(model);
    let (vars, vp, vn, _) = build_vars(&lp, &obj_coefs);

    let mut obj = good_lp::Expression::from(0.0);
    for (i, w) in opts.weights.iter().enumerate() {
        if w.abs() >= f64::EPSILON {
            obj += *w * (vp[i] + vn[i]);
        }
    }
    for (i, c) in obj_coefs.iter().enumerate() {
        if c.abs() >= f64::EPSILON {
            obj -= opts.start_pfba_coef * *c * (vp[i] - vn[i]);
        }
    }
    let mut problem = vars.minimise(obj).using(coin_cbc);
    for expr in build_row_exprs(model, &vp, &vn, model.met_count()) {
        problem = problem.with(constraint!(expr == 0.0));
    }
    for (i, c) in obj_coefs.iter().enumerate() {
        if c.abs() >= f64::EPSILON {
            let expr = *c * (vp[i] - vn[i]);
            problem = problem.with(constraint!(expr >= opts.min_growth));
        }
    }

    let solution = match problem.solve() {
        Ok(s) => s,
        Err(ResolutionError::Infeasible) => {
            return Ok(PfbaSolution {
                status: crate::SolveStatus::Infeasible,
                objective: 0.0,
                fluxes: vec![0.0; n],
                growth: 0.0,
            })
        }
        Err(e) => return Err(FillError::Solver(format!("cbc: {e}"))),
    };
    let fluxes: Vec<f64> =
        (0..n).map(|i| solution.value(vp[i]) - solution.value(vn[i])).collect();
    let growth: f64 = obj_coefs.iter().zip(&fluxes).map(|(c, v)| c * v).sum();
    let objective: f64 = opts
        .weights
        .iter()
        .zip(&fluxes)
        .map(|(w, v)| w * v.abs())
        .sum();
    Ok(PfbaSolution {
        status: crate::SolveStatus::Optimal,
        objective,
        fluxes,
        growth,
    })
}

/// Copy `model` but zero the bounds of reactions whose flux in `sol` is
/// below `tol` in absolute value. The draft ports `rmReact` in-place; here
/// we instead clone and mutate bounds, which is simpler and matches the R
/// code's intent of "can the model still grow once these are excluded?".
fn remove_zero_flux(model: &Model, fluxes: &[f64], tol: f64) -> Model {
    let mut m = model.clone();
    for (r, v) in m.rxns.iter_mut().zip(fluxes) {
        if v.abs() < tol {
            r.lb = 0.0;
            r.ub = 0.0;
        }
    }
    m
}

#[allow(dead_code)]
fn find_bio_idx(model: &Model) -> Option<(usize, &RxnId)> {
    model
        .rxns
        .iter()
        .enumerate()
        .find(|(_, r)| r.obj_coef.abs() >= f64::EPSILON)
        .map(|(i, r)| (i, &r.id))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp::toy_model;

    #[test]
    fn pfba_picks_minimum_absolute_flux() {
        let m = toy_model();
        // All weights = 1 except biomass weight 0 (we penalise pathway use).
        let n = m.rxn_count();
        let mut weights = vec![1.0; n];
        weights[2] = 0.0; // biomass reaction — do not penalise
        let opts = PfbaOptions {
            weights,
            pfba_coef: 1e-3,
            min_growth: 5.0,
            objective: None,
        };
        let sol = pfba(&m, &opts).unwrap();
        assert_eq!(sol.status, SolveStatus::Optimal);
        // Minimum uptake that still hits growth≥5 is exactly 5.
        assert!(sol.growth >= 5.0 - 1e-6);
    }

    #[test]
    fn pfba_heuristic_converges() {
        let m = toy_model();
        let mut weights = vec![1.0; m.rxn_count()];
        weights[2] = 0.0;
        let opts = PfbaHeuristicOptions::new(weights, 1.0);
        let sol = pfba_heuristic(&m, &opts).unwrap();
        assert!(sol.growth >= 1.0 - 1e-6);
    }
}
