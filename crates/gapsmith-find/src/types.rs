//! Result types produced by the find pipeline.

use serde::{Deserialize, Serialize};

/// One row of the Reactions.tbl output, carrying both the reaction metadata
/// (pathway, name, EC, keyrea, spont) and the best blast hit for the
/// reaction (if any). Mirrors `src/analyse_alignments.R`'s `rxndt`.
///
/// Column order matches the golden output at
/// `toy/ecoli-all-Reactions.tbl` so downstream tools (R's
/// `generate_GSdraft.R`, `predict_medium.R`, etc.) parse it unchanged.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionHit {
    pub pathway: String,
    pub pathway_status: Option<PwyStatus>,
    pub rxn: String,
    pub name: String,
    pub ec: String,
    pub keyrea: bool,
    pub spont: bool,

    // Complex-related fields.
    pub is_complex: bool,
    pub subunit_count: u32,
    pub subunits: String,
    pub complex: Option<String>,
    pub subunits_found: Option<u32>,
    pub subunit_undefined_found: Option<bool>,
    pub complex_status: Option<u8>,

    /// The originating reference fasta path (e.g. `rev/1.1.1.1.fasta`) or
    /// `None` if no fasta was found.
    pub file: Option<String>,
    /// Space-joined, sorted, deduplicated SEED reaction IDs associated with
    /// this reaction via EC / MetaCyc-id / enzyme-name lookups. Matches
    /// gapseq's `dbhit` column exactly (see `src/getDBhit.R`). Empty when
    /// no SEED reaction matches.
    pub dbhit: String,
    /// True when `dbhit` is non-empty. Not written to output; handy for
    /// downstream consumers.
    pub has_dbhit: bool,

    /// Source subdirectory (`rxn`, `rev`, `unrev`, `user`) derived from
    /// `file`. Empty when no fasta was found.
    pub src: String,
    /// Reference-sequence type (`EC`, `metacyc`, `reaName`). Empty when no
    /// fasta was found.
    pub reftype: String,

    pub qseqid: Option<String>,
    pub pident: Option<f32>,
    pub evalue: Option<f64>,
    pub bitscore: Option<f32>,
    pub qcov: Option<f32>,
    pub stitle: Option<String>,
    pub sstart: Option<i32>,
    pub send: Option<i32>,
    pub exception: bool,
    pub status: HitStatus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HitStatus {
    /// Hit passes bitscore + identity cutoffs (and exception-table cutoff
    /// if applicable). Reaction is considered present.
    GoodBlast,
    /// Hit present but fails one of the cutoffs.
    BadBlast,
    /// Reference FASTA exists for this reaction but no hit was produced.
    NoBlast,
    /// No reference FASTA could be found for this reaction.
    NoSeqData,
    /// Reaction is marked as spontaneous; no hit and no sequence data.
    Spontaneous,
}

impl HitStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            HitStatus::GoodBlast => "good_blast",
            HitStatus::BadBlast => "bad_blast",
            HitStatus::NoBlast => "no_blast",
            HitStatus::NoSeqData => "no_seq_data",
            HitStatus::Spontaneous => "spontaneous",
        }
    }
}

/// High-level presence tag attached to every predicted-present pathway.
/// See `src/analyse_alignments.R:180-189`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum PwyStatus {
    /// Completeness exactly 100%.
    Full,
    /// Predicted present via the lower (`completenessCutoff`) threshold,
    /// with all key reactions present.
    Threshold,
    /// Predicted present via key-enzyme hint below the main threshold.
    Keyenzyme,
}

impl PwyStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            PwyStatus::Full => "full",
            PwyStatus::Threshold => "threshold",
            PwyStatus::Keyenzyme => "keyenzyme",
        }
    }
}

/// One row of the Pathways.tbl output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathwayResult {
    pub id: String,
    pub name: String,
    pub prediction: bool,
    pub completeness: f64,
    pub status: Option<PwyStatus>,
    pub n_reaction: u32,
    pub n_spontaneous: u32,
    pub n_vague: u32,
    pub n_key_reaction: u32,
    pub n_reaction_found: u32,
    pub n_key_reaction_found: u32,
    pub reactions_found: Vec<String>,
    pub spontaneous_reactions: Vec<String>,
    pub key_reactions: Vec<String>,
}
