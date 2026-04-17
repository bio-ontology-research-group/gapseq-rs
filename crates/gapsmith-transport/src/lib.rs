//! Transporter detection.
//!
//! Port of `src/transporter.sh` + `src/analyse_alignments_transport.R`.
//! The pipeline:
//!
//! 1. Union the TCDB + SEED transporter reference FASTAs.
//! 2. Filter by substrate keywords (from `dat/subex.tbl`) to produce a
//!    small query FASTA — only TC entries whose substrates match any
//!    interesting metabolite.
//! 3. Align that query against the genome (via [`gapsmith_align::Aligner`]).
//! 4. Parse the TC id out of every hit's `qseqid`, resolve its substrate
//!    (either via `tcdb_substrates`/`tcdb_custom`, or by scanning the
//!    FASTA header text as fallback), then join with the SEED transporter
//!    DB (`seed_transporter.tbl` + `seed_transporter_custom.tbl`) to
//!    attach a reaction id.
//! 5. Emit a gapseq-compatible `*-Transporter.tbl`.
//!
//! `gapsmith-transport` exposes both a high-level [`run`] entry point and
//! the individual building blocks (substrate filter, TC parse, etc.) so
//! downstream tooling and tests can compose them.

pub mod data;
pub mod filter;
pub mod output;
pub mod runner;
pub mod tc;

pub use data::{load_seed_transporter, load_tcdb_all, SeedTransporterRow};
pub use filter::{build_small_fasta, BuildSmallResult};
pub use output::{write_transporter_tbl, TransporterRow};
pub use runner::{run, TransportError, TransportOptions, TransportReport};
pub use tc::{extract_tc_id, TC_TYPES};
