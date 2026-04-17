//! Pre-computed alignment-TSV aligner.
//!
//! The user provides a TSV produced by an earlier diamond / blastp / mmseqs2
//! run (or a concatenation thereof — e.g. from the M4.5 batch-cluster mode).
//! The `align` method ignores the `query_fasta` and `target_fasta`
//! arguments and simply returns every parsed hit. Downstream layers
//! (find/transport) are responsible for grouping hits by `qseqid` and
//! applying bitscore / identity / coverage cutoffs.
//!
//! Motivation: batch annotation of many genomes can amortize a single
//! alignment run over a clustered reference set instead of re-running
//! blast/diamond per genome — exactly the use case gapsmith is being
//! built for.

use crate::error::{io_err, AlignError};
use crate::hit::Hit;
use crate::tsv::parse_tsv;
use crate::{AlignOpts, Aligner};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

pub struct PrecomputedTsvAligner {
    path: PathBuf,
    coverage_is_fraction: bool,
    cache: OnceLock<Vec<Hit>>,
}

impl PrecomputedTsvAligner {
    /// TSV where coverage is reported 0–100 (blast / diamond convention).
    pub fn new_percentage(path: impl Into<PathBuf>) -> Self {
        Self { path: path.into(), coverage_is_fraction: false, cache: OnceLock::new() }
    }

    /// TSV where coverage is a 0–1 fraction (mmseqs2 native output).
    pub fn new_fraction(path: impl Into<PathBuf>) -> Self {
        Self { path: path.into(), coverage_is_fraction: true, cache: OnceLock::new() }
    }

    fn load(&self) -> Result<&Vec<Hit>, AlignError> {
        if let Some(hits) = self.cache.get() {
            return Ok(hits);
        }
        let f = File::open(&self.path).map_err(|e| io_err(&self.path, e))?;
        let rdr = BufReader::new(f);
        let hits = parse_tsv(rdr, self.coverage_is_fraction)?;
        let _ = self.cache.set(hits);
        Ok(self.cache.get().unwrap())
    }
}

impl Aligner for PrecomputedTsvAligner {
    fn name(&self) -> &'static str {
        "precomputed"
    }

    fn align(
        &self,
        _query_fasta: &Path,
        _target_fasta: &Path,
        _opts: &AlignOpts,
    ) -> Result<Vec<Hit>, AlignError> {
        let hits = self.load()?;
        Ok(hits.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn loads_and_caches() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("a.tsv");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "q1\t95.5\t1e-100\t312.5\t90\tt1\t1\t200").unwrap();
        writeln!(f, "q2\t80.0\t1e-20\t120.0\t70\tt2 descr\t10\t150").unwrap();
        drop(f);
        let a = PrecomputedTsvAligner::new_percentage(&p);
        let h = a.align(Path::new("unused"), Path::new("unused"), &AlignOpts::default()).unwrap();
        assert_eq!(h.len(), 2);
        assert_eq!(h[0].qseqid, "q1");
        assert_eq!(h[1].stitle, "t2 descr");
        // Second call serves from cache — we can't observe timing but the
        // result stays identical.
        let h2 = a.align(Path::new("x"), Path::new("y"), &AlignOpts::default()).unwrap();
        assert_eq!(h, h2);
    }
}
