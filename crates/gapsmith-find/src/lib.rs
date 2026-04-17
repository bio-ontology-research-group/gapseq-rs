//! Pathway / reaction finder.
//!
//! Mirrors the R reference implementation spread across
//! `src/prepare_batch_alignments.R`, `src/analyse_alignments.R`, and
//! `src/gapseq_find.sh`. The pipeline:
//!
//! 1. [`pathways::select`] picks a subset of [`gapsmith_db::PathwayRow`] to
//!    evaluate, based on a user keyword or pattern.
//! 2. [`seqfile::resolve_for_reaction`] walks the `dat/seq/` tree to find
//!    the reference FASTA(s) for each reaction.
//! 3. [`runner::run_find`] builds one concatenated `query.faa`, runs the
//!    alignment (via [`gapsmith_align::Aligner`]) against the input genome,
//!    classifies every hit, and aggregates per-pathway completeness.
//! 4. [`output`] emits `*-Reactions.tbl` and `*-Pathways.tbl` in gapseq's
//!    column order.
//!
//! Complex / subunit detection (`src/complex_detection.R`) lives in
//! [`complex`]; it has point-by-point R-parity on 9 handcrafted cases
//! (greek / latin numerals, size-dict synonyms, coverage edges). See
//! `crates/gapsmith-find/tests/complex_parity.rs`.
//!
//! # End-to-end parity
//!
//! [`runner::run_find`] produces byte-identical `*-Pathways.tbl` output
//! against real gapseq on two test cases (`-p PWY-6587` and `-p amino` on
//! `toy/ecore.faa`). See `crates/gapsmith-find/tests/pipeline_parity.rs`.

pub mod classify;
pub mod complex;
pub mod dbhit;
pub mod output;
pub mod pathways;
pub mod runner;
pub mod seqfile;
pub mod taxonomy;
pub mod types;

pub use classify::{classify_hits, ClassifyOptions};
pub use output::{write_pathways_tbl, write_reactions_tbl};
pub use pathways::{select, ExpandedReaction, PathwaySelectOptions};
pub use runner::{run_find, FindError, FindOptions, FindReport};
pub use seqfile::{resolve_for_reaction, ResolvedSeq, SeqfileOptions};
pub use types::{HitStatus, PathwayResult, PwyStatus, ReactionHit};
