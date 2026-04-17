//! Flux-balance analysis — maximise `Σ c_r · v_r` s.t. steady state + bounds.
//!
//! Port of cobrar's `fba()` as used across gapseq's fill / doall pipeline.
//! We translate each reaction into a (vp, vn) pair so the split-flux LP from
//! [`crate::lp`] can later drive pFBA under the same encoding.

use crate::error::{FillError, SolveStatus};
use crate::lp::SplitFluxLp;
use gapsmith_core::Model;
use good_lp::{
    constraint, solvers::highs::highs, variable, Expression, ProblemVariables,
    ResolutionError, Solution, SolverModel, Variable,
};

#[derive(Debug, Clone)]
pub struct FbaOptions {
    /// Explicit objective coefficients. If `None`, the `obj_coef` field on
    /// each reaction is used.
    pub objective: Option<Vec<f64>>,
    /// Direction: `true` → maximise (default), `false` → minimise.
    pub maximise: bool,
    /// Bound applied when a reaction declares an unbounded (0 / 1000)
    /// default; anything ≥ this is treated as +∞ internally. Matches
    /// cobrar's `COBRAR_MAX_FLUX`.
    pub max_flux: f64,
}

impl Default for FbaOptions {
    fn default() -> Self {
        Self { objective: None, maximise: true, max_flux: 1000.0 }
    }
}

#[derive(Debug, Clone)]
pub struct FbaSolution {
    pub status: SolveStatus,
    pub objective: f64,
    pub fluxes: Vec<f64>,
}

pub fn fba(model: &Model, opts: &FbaOptions) -> Result<FbaSolution, FillError> {
    model.check_shape().map_err(|_| FillError::BadShape)?;
    let n = model.rxn_count();
    let m = model.met_count();

    let coefs = resolve_obj(model, opts)?;
    let lp = SplitFluxLp::from_model(model);

    let (vars, vp, vn, obj_expr) = build_vars(&lp, &coefs);

    let mut problem = if opts.maximise {
        vars.maximise(obj_expr).using(highs)
    } else {
        vars.minimise(obj_expr).using(highs)
    };

    // Mass balance: Σ_r S[i,r] · (vp[r] - vn[r]) = 0 for each met i.
    // We accumulate per-row expressions by walking the CSC matrix once.
    let row_exprs = build_row_exprs(model, &vp, &vn, m);
    for expr in row_exprs {
        problem = problem.with(constraint!(expr == 0.0));
    }

    let solution = match problem.solve() {
        Ok(s) => s,
        Err(ResolutionError::Infeasible) => {
            return Ok(FbaSolution {
                status: SolveStatus::Infeasible,
                objective: 0.0,
                fluxes: vec![0.0; n],
            })
        }
        Err(ResolutionError::Unbounded) => {
            return Ok(FbaSolution {
                status: SolveStatus::Unbounded,
                objective: 0.0,
                fluxes: vec![0.0; n],
            })
        }
        Err(e) => return Err(FillError::Solver(e.to_string())),
    };

    let fluxes: Vec<f64> =
        (0..n).map(|i| solution.value(vp[i]) - solution.value(vn[i])).collect();
    let objective: f64 = coefs.iter().zip(&fluxes).map(|(c, v)| c * v).sum();

    Ok(FbaSolution { status: SolveStatus::Optimal, objective, fluxes })
}

pub(crate) fn resolve_obj(model: &Model, opts: &FbaOptions) -> Result<Vec<f64>, FillError> {
    let n = model.rxn_count();
    let coefs = match &opts.objective {
        Some(v) if v.len() == n => v.clone(),
        Some(_) => {
            return Err(FillError::BadInput(
                "objective length must equal reaction count".into(),
            ));
        }
        None => model.rxns.iter().map(|r| r.obj_coef).collect(),
    };
    if coefs.iter().all(|c| c.abs() < f64::EPSILON) {
        return Err(FillError::NoObjective);
    }
    Ok(coefs)
}

pub(crate) fn build_vars(
    lp: &SplitFluxLp,
    coefs: &[f64],
) -> (ProblemVariables, Vec<Variable>, Vec<Variable>, Expression) {
    let n = lp.model.rxn_count();
    let mut vars = ProblemVariables::new();

    let vp: Vec<Variable> = (0..n)
        .map(|i| {
            let mut v = variable().min(lp.vp_lb[i]).max(lp.vp_ub[i]);
            if !lp.model.rxns[i].id.as_str().is_empty() {
                v = v.name(format!("vp_{}", lp.model.rxns[i].id.as_str()));
            }
            vars.add(v)
        })
        .collect();

    let vn: Vec<Variable> = (0..n)
        .map(|i| {
            let mut v = variable().min(lp.vn_lb[i]).max(lp.vn_ub[i]);
            if !lp.model.rxns[i].id.as_str().is_empty() {
                v = v.name(format!("vn_{}", lp.model.rxns[i].id.as_str()));
            }
            vars.add(v)
        })
        .collect();

    let obj_expr: Expression = coefs
        .iter()
        .enumerate()
        .filter(|(_, c)| c.abs() >= f64::EPSILON)
        .map(|(i, c)| *c * (vp[i] - vn[i]))
        .sum();

    (vars, vp, vn, obj_expr)
}

/// Build one [`Expression`] per metabolite row — the mass-balance LHS.
///
/// We walk the CSC matrix column-by-column (the cheap direction) and push
/// each non-zero `coef · (vp[c] - vn[c])` onto row `r`'s expression. Result
/// is `O(nnz)` total vs `O(m · n)` from the naive per-row scan.
pub(crate) fn build_row_exprs(
    model: &Model,
    vp: &[Variable],
    vn: &[Variable],
    m: usize,
) -> Vec<Expression> {
    let s = model.s.inner();
    let n = model.rxn_count();
    let mut out: Vec<Expression> = (0..m).map(|_| Expression::from(0.0)).collect();
    if s.is_csc() {
        for c in 0..n {
            if let Some(col) = s.outer_view(c) {
                for (r, &coef) in col.iter() {
                    out[r] = out[r].clone() + coef * (vp[c] - vn[c]);
                }
            }
        }
    } else {
        // Fallback — triplet iteration for unexpected CSR layouts.
        for (val, (r, c)) in s.iter() {
            out[r] = out[r].clone() + *val * (vp[c] - vn[c]);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chain_reaches_uptake_bound() {
        let m = crate::lp::toy_model();
        let sol = fba(&m, &FbaOptions::default()).unwrap();
        assert_eq!(sol.status, SolveStatus::Optimal);
        // With uptake bounded at 10 on EX_A, steady-state says flux through
        // R_BC (biomass) is also 10.
        assert!((sol.objective - 10.0).abs() < 1e-6, "obj = {}", sol.objective);
        assert!((sol.fluxes[2] - 10.0).abs() < 1e-6);
    }

    #[test]
    fn infeasible_when_uptake_blocked() {
        let mut m = crate::lp::toy_model();
        // EX_A closed → no flux can enter.
        m.rxns[0].lb = 0.0;
        m.rxns[0].ub = 0.0;
        let sol = fba(&m, &FbaOptions::default()).unwrap();
        // Still "optimal" but objective is zero.
        assert_eq!(sol.status, SolveStatus::Optimal);
        assert!(sol.objective.abs() < 1e-6);
    }

    #[test]
    fn no_objective_errors() {
        let mut m = crate::lp::toy_model();
        for r in &mut m.rxns {
            r.obj_coef = 0.0;
        }
        let err = fba(&m, &FbaOptions::default()).unwrap_err();
        matches!(err, FillError::NoObjective);
    }
}
