//! Diamond aligner (protein-protein).
//!
//! Command template mirrors `src/gapseq_find.sh:592–600`:
//!
//! ```text
//! diamond makedb --in TARGET -d orgdb
//! diamond blastp -d orgdb.dmnd -q QUERY --threads N --out alignments.tsv
//!                --outfmt 6 qseqid pident evalue bitscore qcovhsp stitle sstart send
//!                --query-cover 75
//! ```
//!
//! By default we pass `--more-sensitive`, matching gapseq's default `aliArgs`
//! for diamond.

use crate::blast::which;
use crate::error::{io_err, AlignError};
use crate::hit::Hit;
use crate::tsv::parse_tsv;
use crate::{AlignOpts, Aligner};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::process::Command;

/// Columns emitted by our diamond invocation. Note `qcovhsp` (per-HSP) rather
/// than `qcovs` — diamond doesn't support `qcovs`. Both are 0–100.
const COLUMNS: &[&str] = &[
    "qseqid", "pident", "evalue", "bitscore", "qcovhsp", "stitle", "sstart", "send",
];

pub struct DiamondAligner {
    /// Pass `--more-sensitive` (default: true). Matches gapseq's aliArgs
    /// default when the user doesn't override `-R`.
    pub more_sensitive: bool,
}

impl DiamondAligner {
    pub fn new() -> Self {
        Self { more_sensitive: true }
    }
}
impl Default for DiamondAligner {
    fn default() -> Self {
        Self::new()
    }
}

impl Aligner for DiamondAligner {
    fn name(&self) -> &'static str {
        "diamond"
    }

    fn align(
        &self,
        query_fasta: &Path,
        target_fasta: &Path,
        opts: &AlignOpts,
    ) -> Result<Vec<Hit>, AlignError> {
        require("diamond")?;
        let tmp = tempfile::tempdir().map_err(|e| io_err(Path::new("/tmp"), e))?;
        let db_prefix = tmp.path().join("orgdb");
        let out_tsv = tmp.path().join("aln.tsv");

        // makedb
        let out = Command::new("diamond")
            .arg("makedb")
            .arg("--in")
            .arg(target_fasta)
            .arg("-d")
            .arg(&db_prefix)
            .arg("--quiet")
            .output()
            .map_err(|e| io_err(Path::new("diamond"), e))?;
        if !out.status.success() {
            return Err(AlignError::ToolFailed {
                tool: "diamond",
                status: out.status,
                stderr: String::from_utf8_lossy(&out.stderr).to_string(),
            });
        }

        // blastp
        let dmnd_db = {
            let mut p = db_prefix.clone();
            p.set_extension("dmnd");
            p
        };
        let mut cmd = Command::new("diamond");
        cmd.arg("blastp")
            .arg("-d")
            .arg(&dmnd_db)
            .arg("-q")
            .arg(query_fasta)
            .arg("--threads")
            .arg(opts.threads.to_string())
            .arg("--out")
            .arg(&out_tsv)
            .arg("--outfmt")
            .arg("6");
        for c in COLUMNS {
            cmd.arg(c);
        }
        cmd.arg("--query-cover").arg(opts.coverage_pct.to_string());
        if self.more_sensitive {
            cmd.arg("--more-sensitive");
        }
        if let Some(e) = opts.evalue {
            cmd.arg("--evalue").arg(e.to_string());
        }
        if opts.quiet {
            cmd.arg("--quiet");
        }
        for a in &opts.extra_args {
            cmd.arg(a);
        }

        tracing::debug!(?cmd, "running diamond");
        let out = cmd
            .output()
            .map_err(|e| io_err(Path::new("diamond"), e))?;
        if !out.status.success() {
            return Err(AlignError::ToolFailed {
                tool: "diamond",
                status: out.status,
                stderr: String::from_utf8_lossy(&out.stderr).to_string(),
            });
        }
        if !opts.quiet && !out.stderr.is_empty() {
            eprintln!("{}", String::from_utf8_lossy(&out.stderr));
        }

        let f = File::open(&out_tsv).map_err(|e| io_err(&out_tsv, e))?;
        parse_tsv(BufReader::new(f), false)
    }
}

fn require(tool: &'static str) -> Result<(), AlignError> {
    if which(tool).is_some() {
        Ok(())
    } else {
        Err(AlignError::ToolMissing { tool })
    }
}
