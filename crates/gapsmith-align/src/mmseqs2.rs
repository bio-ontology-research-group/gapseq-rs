//! MMseqs2 aligner (protein-protein).
//!
//! Mirrors gapseq's explicit 4-command pipeline from
//! `src/gapseq_find.sh:603–615`:
//!
//! ```text
//! mmseqs createdb TARGET targetDB
//! mmseqs createdb QUERY queryDB
//! mmseqs search queryDB targetDB resultDB tmp --threads N -c 0.C
//! mmseqs convertalis queryDB targetDB resultDB out.tsv
//!     --format-output "qheader,pident,evalue,bits,qcov,theader,tstart,tend"
//! ```
//!
//! Rationale: `mmseqs easy-search` would collapse the four steps into one,
//! but it defaults to `align --alignment-mode 3` (full-alignment identity)
//! whereas `search` reports the k-mer prefilter identity. Identity cutoffs
//! downstream (`analyse_alignments.R`) were calibrated against the latter —
//! so we match that semantics exactly. Native coverage is a 0–1 fraction;
//! we rescale to 0–100 via [`crate::tsv::parse_tsv`].

use crate::blast::which;
use crate::error::{io_err, AlignError};
use crate::hit::Hit;
use crate::tsv::parse_tsv;
use crate::{AlignOpts, Aligner};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::process::Command;

/// Columns in the mmseqs2 native format. Note `bits` (not `bitscore`),
/// `qcov` (fraction), `theader` instead of `stitle`. We use `qheader`
/// (rather than `query`) so the emitted first column matches what diamond
/// and blastp produce from the same FASTA — both of those record the
/// entire first whitespace-delimited token of the header. Post-processing
/// in [`normalize_qheader`] replicates gapseq's sed fixup from
/// `src/gapseq_find.sh`.
const FORMAT: &str = "qheader,pident,evalue,bits,qcov,theader,tstart,tend";

pub struct Mmseqs2Aligner {
    /// Sensitivity (`-s`). Default 5.7 — matches mmseqs's default.
    pub sensitivity: f32,
}

impl Mmseqs2Aligner {
    pub fn new() -> Self {
        Self { sensitivity: 5.7 }
    }
}
impl Default for Mmseqs2Aligner {
    fn default() -> Self {
        Self::new()
    }
}

fn require(tool: &'static str) -> Result<(), AlignError> {
    if which(tool).is_some() {
        Ok(())
    } else {
        Err(AlignError::ToolMissing { tool })
    }
}

impl Aligner for Mmseqs2Aligner {
    fn name(&self) -> &'static str {
        "mmseqs2"
    }

    fn align(
        &self,
        query_fasta: &Path,
        target_fasta: &Path,
        opts: &AlignOpts,
    ) -> Result<Vec<Hit>, AlignError> {
        require("mmseqs")?;
        let tmp = tempfile::tempdir().map_err(|e| io_err(Path::new("/tmp"), e))?;
        let target_db = tmp.path().join("targetDB");
        let query_db = tmp.path().join("queryDB");
        let result_db = tmp.path().join("resultDB");
        let out_tsv = tmp.path().join("aln.tsv");
        let scratch = tmp.path().join("scratch");
        std::fs::create_dir_all(&scratch).map_err(|e| io_err(&scratch, e))?;

        let coverage = opts.coverage_pct as f64 / 100.0;
        let verbosity = if opts.quiet { "1" } else { "3" };

        // 1) createdb for target
        run_mmseqs(&[
            "createdb",
            target_fasta.to_str().unwrap_or(""),
            target_db.to_str().unwrap_or(""),
            "-v",
            verbosity,
        ])?;

        // 2) createdb for query
        run_mmseqs(&[
            "createdb",
            query_fasta.to_str().unwrap_or(""),
            query_db.to_str().unwrap_or(""),
            "-v",
            verbosity,
        ])?;

        // 3) search
        let mut search = vec![
            "search".to_string(),
            query_db.display().to_string(),
            target_db.display().to_string(),
            result_db.display().to_string(),
            scratch.display().to_string(),
            "--threads".to_string(),
            opts.threads.to_string(),
            "-c".to_string(),
            format!("{coverage:.2}"),
            "-s".to_string(),
            self.sensitivity.to_string(),
            "-v".to_string(),
            verbosity.to_string(),
        ];
        if let Some(e) = opts.evalue {
            search.push("-e".into());
            search.push(e.to_string());
        }
        for a in &opts.extra_args {
            search.push(a.clone());
        }
        run_mmseqs(&search.iter().map(|s| s.as_str()).collect::<Vec<_>>())?;

        // 4) convertalis → TSV
        run_mmseqs(&[
            "convertalis",
            query_db.to_str().unwrap_or(""),
            target_db.to_str().unwrap_or(""),
            result_db.to_str().unwrap_or(""),
            out_tsv.to_str().unwrap_or(""),
            "--format-mode",
            "0",
            "--format-output",
            FORMAT,
            "-v",
            verbosity,
        ])?;

        let f = File::open(&out_tsv).map_err(|e| io_err(&out_tsv, e))?;
        let mut hits = parse_tsv(BufReader::new(f), true)?;
        for h in &mut hits {
            h.qseqid = normalize_qheader(&h.qseqid);
        }
        Ok(hits)
    }
}

fn run_mmseqs(args: &[&str]) -> Result<(), AlignError> {
    let mut cmd = Command::new("mmseqs");
    for a in args {
        cmd.arg(a);
    }
    tracing::debug!(?cmd, "running mmseqs");
    let out = cmd
        .output()
        .map_err(|e| io_err(Path::new("mmseqs"), e))?;
    if !out.status.success() {
        return Err(AlignError::ToolFailed {
            tool: "mmseqs2",
            status: out.status,
            stderr: String::from_utf8_lossy(&out.stderr).to_string(),
        });
    }
    Ok(())
}

/// Reduce a mmseqs2 `qheader` field to its first whitespace-delimited token —
/// matching diamond and blastp's `qseqid` output. Mirrors the sed fixup in
/// `src/gapseq_find.sh` after a `convertalis` call.
fn normalize_qheader(h: &str) -> String {
    match h.split_once(|c: char| c.is_whitespace()) {
        Some((first, _)) => first.to_string(),
        None => h.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strip_qheader_first_token() {
        assert_eq!(
            normalize_qheader("sp|P77580|ACDH_ECOLI Acetaldehyde dehydrogenase"),
            "sp|P77580|ACDH_ECOLI"
        );
        assert_eq!(normalize_qheader("noSpaces"), "noSpaces");
        assert_eq!(normalize_qheader("tab\tdelimited"), "tab");
    }
}
