//! Sequence alignment abstraction for gapsmith.
//!
//! Every aligner (blast, diamond, mmseqs2, precomputed TSV) exposes a common
//! [`Aligner`] trait that takes a query FASTA and a target FASTA and returns
//! a vector of [`Hit`]. Internally the shell-out implementations manage
//! their own temp work-directories so callers just see FASTA-in, hits-out.
//!
//! # Example
//!
//! ```ignore
//! use gapsmith_align::{AlignOpts, Aligner, DiamondAligner};
//! use std::path::Path;
//!
//! let aligner = DiamondAligner;
//! let hits = aligner.run(
//!     Path::new("query.faa"),
//!     &[Path::new("reference.faa")],
//!     &AlignOpts::default(),
//! ).unwrap();
//! for h in hits.iter().take(5) {
//!     println!("{}\t{}\t{}", h.qseqid, h.pident, h.bitscore);
//! }
//! ```
//!
//! # Backend selection
//!
//! - [`BlastpAligner`] — protein-vs-protein; always available if NCBI BLAST+
//!   is on `PATH`. Slow on large genomes but the gapseq reference.
//! - [`TblastnAligner`] — protein query vs nucleotide subject (rare; used
//!   for nucleotide-based reference FASTAs).
//! - [`DiamondAligner`] — 5-20× faster than BLASTp on large proteomes;
//!   comparable sensitivity at `--more-sensitive` (which we default on).
//! - [`Mmseqs2Aligner`] — fast k-mer-based alternative; we replicate
//!   gapseq's 4-command pipeline (createdb → search → convertalis) rather
//!   than `easy-search`, because the latter reports full-alignment
//!   identities instead of the k-mer prefilter identities gapseq
//!   calibrates against.
//! - [`PrecomputedTsvAligner`] — skips the aligner entirely; reads a TSV
//!   the caller produced with their own tool. Used by `gapsmith`'s
//!   `--aligner precomputed` mode and by [`BatchClusterAligner`].
//! - [`BatchClusterAligner`] — new in gapsmith. mmseqs2-clusters N
//!   genomes, runs one alignment against the reference, then expands the
//!   cluster membership to per-genome TSVs. Amortises aligner cost over
//!   many genomes.
//!
//! Columns always emitted by our wrappers (matching gapseq's convention):
//!
//! | column  | meaning                                             |
//! |---------|-----------------------------------------------------|
//! | qseqid  | query identifier (full FASTA header, up to a space) |
//! | pident  | percent identity (0–100)                            |
//! | evalue  | BLAST-style e-value                                  |
//! | bitscore| bit score                                           |
//! | qcov    | query coverage (0–100)                              |
//! | stitle  | subject title (may contain spaces)                  |
//! | sstart  | subject start                                       |
//! | send    | subject end                                         |
//!
//! This keeps parity with `src/gapseq_find.sh` lines 249–255.

pub mod batch;
pub mod blast;
pub mod diamond;
pub mod error;
pub mod hit;
pub mod mmseqs2;
pub mod precomputed;
pub mod tsv;

pub use batch::{BatchClusterAligner, ClusterResult, GenomeHitSet, GenomeInput};

// Re-exported for external parity tests that need to reuse the TSV parser.
// Not part of the stable public API; prefer the per-backend aligners.

pub use blast::{BlastpAligner, TblastnAligner};
pub use diamond::DiamondAligner;
pub use error::AlignError;
pub use hit::Hit;
pub use mmseqs2::Mmseqs2Aligner;
pub use precomputed::PrecomputedTsvAligner;

use std::path::Path;

/// Options tuning an alignment run. Sensible gapseq defaults: coverage 75%,
/// use all detected cores, no extra user args.
#[derive(Debug, Clone)]
pub struct AlignOpts {
    /// Worker threads to pass to the underlying binary.
    pub threads: usize,
    /// Minimum query coverage, 0–100 (not enforced by PrecomputedTsvAligner).
    pub coverage_pct: u32,
    /// Optional e-value cutoff. None = leave to the tool's default.
    pub evalue: Option<f64>,
    /// Additional free-form CLI flags appended verbatim. Unsafe if fed
    /// untrusted input — matches gapseq's `-R` flag in `src/gapseq_find.sh`.
    pub extra_args: Vec<String>,
    /// Suppress the tool's stderr when true. Useful in batch pipelines.
    pub quiet: bool,
}

impl Default for AlignOpts {
    fn default() -> Self {
        Self {
            threads: std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            coverage_pct: 75,
            evalue: None,
            extra_args: Vec::new(),
            quiet: false,
        }
    }
}

/// Common trait implemented by every aligner backend.
pub trait Aligner: Send + Sync {
    /// Short label printed in logs (`"blast"`, `"diamond"`, `"mmseqs2"`,
    /// `"tblastn"`, `"precomputed"`).
    fn name(&self) -> &'static str;

    /// Align `query_fasta` against `target_fasta`. The `target_fasta` may be
    /// a protein or nucleotide file depending on the backend — check the
    /// impl's docs.
    fn align(
        &self,
        query_fasta: &Path,
        target_fasta: &Path,
        opts: &AlignOpts,
    ) -> Result<Vec<Hit>, AlignError>;
}
