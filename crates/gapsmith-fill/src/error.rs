//! Error and status types for the LP layer.

/// Solver status — collapsed from `good_lp::ResolutionError` so callers can
/// branch without depending on the solver crate.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SolveStatus {
    /// Solver found an optimal primal solution.
    Optimal,
    /// Problem has no feasible solution.
    Infeasible,
    /// Problem is unbounded (the objective can be driven to ±∞).
    Unbounded,
    /// Solver returned a non-optimal status (time-out, numeric trouble, etc.).
    Other,
}

#[derive(Debug, thiserror::Error)]
pub enum FillError {
    #[error("model has no reaction with obj_coef != 0")]
    NoObjective,
    #[error("model stoichiometric matrix shape does not match mets × rxns")]
    BadShape,
    #[error("solver error: {0}")]
    Solver(String),
    #[error("bad input: {0}")]
    BadInput(String),
}
