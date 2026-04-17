//! Sparse stoichiometric matrix.
//!
//! Orientation is CSC (compressed sparse column): one column per reaction,
//! one row per metabolite. FBA-style constraint assembly iterates each
//! reaction column once per solve, so column-major layout keeps cache behavior
//! favorable.
//!
//! For construction we accept `(row, col, coef)` triplets. Internally we build
//! a [`sprs::TriMat`], then freeze to [`sprs::CsMat`] in CSC orientation.

use serde::{Deserialize, Serialize};
use sprs::{CsMat, TriMat};

/// Newtype around [`sprs::CsMat`] that serializes transparently.
///
/// Keeping this as a newtype (rather than a bare `CsMat`) lets us evolve the
/// representation later — e.g. switching to a hand-rolled CSR or adding
/// per-column cached metadata — without touching call sites.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct StoichMatrix(CsMat<f64>);

impl Default for StoichMatrix {
    fn default() -> Self {
        Self::zeros(0, 0)
    }
}

impl StoichMatrix {
    /// Build from `(row, col, value)` triplets. Shape is `(n_rows, n_cols)`.
    ///
    /// Duplicate `(row, col)` entries are summed (matches `sprs` semantics).
    pub fn from_triplets(
        n_rows: usize,
        n_cols: usize,
        triplets: impl IntoIterator<Item = (usize, usize, f64)>,
    ) -> Self {
        let mut tri: TriMat<f64> = TriMat::new((n_rows, n_cols));
        for (r, c, v) in triplets {
            tri.add_triplet(r, c, v);
        }
        Self(tri.to_csc())
    }

    /// Empty matrix of the given shape (CSC orientation).
    pub fn zeros(n_rows: usize, n_cols: usize) -> Self {
        let tri: TriMat<f64> = TriMat::new((n_rows, n_cols));
        Self(tri.to_csc())
    }

    pub fn rows(&self) -> usize {
        self.0.rows()
    }
    pub fn cols(&self) -> usize {
        self.0.cols()
    }
    pub fn nnz(&self) -> usize {
        self.0.nnz()
    }

    pub fn inner(&self) -> &CsMat<f64> {
        &self.0
    }
    pub fn into_inner(self) -> CsMat<f64> {
        self.0
    }

    /// Iterate non-zero entries of a single reaction column as `(met_row, coef)`.
    ///
    /// Assumes CSC orientation (the invariant enforced by our constructors).
    /// If the underlying matrix happens to be CSR, we fall back to a linear
    /// scan of all triplets — correct, just slower.
    pub fn column(&self, col: usize) -> Vec<(usize, f64)> {
        let mut out = Vec::new();
        if self.0.is_csc() {
            if let Some(view) = self.0.outer_view(col) {
                for (row, &val) in view.iter() {
                    out.push((row, val));
                }
            }
        } else {
            for (val, (row, c)) in self.0.iter() {
                if c == col {
                    out.push((row, *val));
                }
            }
        }
        out
    }
}

impl From<CsMat<f64>> for StoichMatrix {
    fn from(m: CsMat<f64>) -> Self {
        Self(m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_and_query() {
        let s = StoichMatrix::from_triplets(3, 2, vec![(0, 0, -1.0), (1, 0, 1.0), (2, 1, 2.0)]);
        assert_eq!(s.rows(), 3);
        assert_eq!(s.cols(), 2);
        assert_eq!(s.nnz(), 3);

        let col0 = s.column(0);
        assert_eq!(col0.len(), 2);
        assert!(col0.contains(&(0, -1.0)));
        assert!(col0.contains(&(1, 1.0)));
    }

    #[test]
    fn serde_json_roundtrip() {
        let s = StoichMatrix::from_triplets(2, 2, vec![(0, 0, -1.0), (1, 1, 1.0)]);
        let j = serde_json::to_string(&s).unwrap();
        let back: StoichMatrix = serde_json::from_str(&j).unwrap();
        assert_eq!(back.rows(), 2);
        assert_eq!(back.cols(), 2);
        assert_eq!(back.nnz(), 2);
    }
}
