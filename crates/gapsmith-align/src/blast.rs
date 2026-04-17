//! BLAST+ aligners: `blastp` (protein-protein) and `tblastn` (protein query vs.
//! nucleotide database, translated in all 6 frames).
//!
//! Shell-out target: `makeblastdb` + `blastp` / `tblastn`. Command template
//! mirrors `src/gapseq_find.sh:585–589`:
//!
//! ```text
//! makeblastdb -in TARGET -dbtype prot -out orgdb
//! blastp -db orgdb -query QUERY -qcov_hsp_perc 75 -num_threads N
//!        -outfmt "6 qseqid pident evalue bitscore qcovs stitle sstart send"
//! ```
//!
//! `tblastn` differs only in `-dbtype nucl` and the search executable.

use crate::error::{io_err, AlignError};
use crate::hit::Hit;
use crate::tsv::parse_tsv;
use crate::{AlignOpts, Aligner};
use std::io::Cursor;
use std::path::Path;
use std::process::Command;

const COLUMNS: &str = "qseqid pident evalue bitscore qcovs stitle sstart send";

pub struct BlastpAligner;
pub struct TblastnAligner;

impl BlastpAligner {
    pub fn new() -> Self {
        Self
    }
}
impl Default for BlastpAligner {
    fn default() -> Self {
        Self::new()
    }
}
impl TblastnAligner {
    pub fn new() -> Self {
        Self
    }
}
impl Default for TblastnAligner {
    fn default() -> Self {
        Self::new()
    }
}

fn run_blast(
    search_tool: &'static str,
    dbtype: &'static str,
    query_fasta: &Path,
    target_fasta: &Path,
    opts: &AlignOpts,
) -> Result<Vec<Hit>, AlignError> {
    require_tool("makeblastdb")?;
    require_tool(search_tool)?;

    let tmp = tempfile::tempdir().map_err(|e| io_err(Path::new("/tmp"), e))?;
    let db_prefix = tmp.path().join("orgdb");

    // makeblastdb -in TARGET -dbtype {prot,nucl} -out db_prefix
    let out = Command::new("makeblastdb")
        .arg("-in")
        .arg(target_fasta)
        .arg("-dbtype")
        .arg(dbtype)
        .arg("-out")
        .arg(&db_prefix)
        .output()
        .map_err(|e| io_err(Path::new("makeblastdb"), e))?;
    if !out.status.success() {
        return Err(AlignError::ToolFailed {
            tool: "makeblastdb",
            status: out.status,
            stderr: String::from_utf8_lossy(&out.stderr).to_string(),
        });
    }

    // blastp / tblastn
    let mut cmd = Command::new(search_tool);
    cmd.arg("-db")
        .arg(&db_prefix)
        .arg("-query")
        .arg(query_fasta)
        .arg("-qcov_hsp_perc")
        .arg(opts.coverage_pct.to_string())
        .arg("-num_threads")
        .arg(opts.threads.to_string())
        .arg("-outfmt")
        .arg(format!("6 {COLUMNS}"));
    if let Some(e) = opts.evalue {
        cmd.arg("-evalue").arg(e.to_string());
    }
    for a in &opts.extra_args {
        cmd.arg(a);
    }

    tracing::debug!(?cmd, "running blast");
    let out = cmd
        .output()
        .map_err(|e| io_err(Path::new(search_tool), e))?;
    if !out.status.success() {
        return Err(AlignError::ToolFailed {
            tool: search_tool,
            status: out.status,
            stderr: String::from_utf8_lossy(&out.stderr).to_string(),
        });
    }
    if !opts.quiet && !out.stderr.is_empty() {
        eprintln!("{}", String::from_utf8_lossy(&out.stderr));
    }
    parse_tsv(Cursor::new(out.stdout), false)
}

impl Aligner for BlastpAligner {
    fn name(&self) -> &'static str {
        "blastp"
    }
    fn align(
        &self,
        query_fasta: &Path,
        target_fasta: &Path,
        opts: &AlignOpts,
    ) -> Result<Vec<Hit>, AlignError> {
        run_blast("blastp", "prot", query_fasta, target_fasta, opts)
    }
}

impl Aligner for TblastnAligner {
    fn name(&self) -> &'static str {
        "tblastn"
    }
    fn align(
        &self,
        query_fasta: &Path,
        target_fasta: &Path,
        opts: &AlignOpts,
    ) -> Result<Vec<Hit>, AlignError> {
        run_blast("tblastn", "nucl", query_fasta, target_fasta, opts)
    }
}

fn require_tool(tool: &'static str) -> Result<(), AlignError> {
    if which(tool).is_some() {
        Ok(())
    } else {
        Err(AlignError::ToolMissing { tool })
    }
}

pub(crate) fn which(name: &str) -> Option<std::path::PathBuf> {
    let path = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path) {
        let candidate = dir.join(name);
        if is_executable(&candidate) {
            return Some(candidate);
        }
    }
    None
}

#[cfg(unix)]
pub(crate) fn is_executable(p: &std::path::Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    p.metadata()
        .map(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
        .unwrap_or(false)
}
#[cfg(not(unix))]
pub(crate) fn is_executable(p: &std::path::Path) -> bool {
    p.is_file()
}
