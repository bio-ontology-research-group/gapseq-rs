//! Rule-based medium inference — port of `src/predict_medium.R`.
//!
//! Evaluates the boolean rules in `dat/medium_prediction_rules.tsv`
//! against a draft model and a Pathways.tbl prediction set to produce a
//! medium CSV suitable for passing to `gapsmith fill`.
//!
//! # Submodules
//!
//! - [`boolexpr`] — tiny recursive-descent parser / evaluator for the
//!   rule language (`|`, `&`, `!`, `<`, `>`, `<=`, `>=`, `==`, `+`).
//!   Returns a `bool`; any non-zero numeric result counts as true.
//! - [`rules`] — TSV loader for the rules table. Tolerant of Latin-1
//!   comment columns (the real gapseq file contains non-ASCII glyphs
//!   like β-D-xylose).
//! - [`predict`] — the [`predict_medium`] entry point that ties the
//!   pieces together, dedups by cpd, balances protons, applies manual
//!   overrides.
//!
//! # Example
//!
//! ```no_run
//! use gapsmith_medium::{load_rules, predict_medium, parse_manual_flux};
//! use gapsmith_core::Model;
//! use std::collections::HashSet;
//!
//! let rules = load_rules(std::path::Path::new("dat/medium_prediction_rules.tsv")).unwrap();
//! let predicted_pathways: HashSet<String> = /* from Pathways.tbl */
//! #    HashSet::new();
//! let model: Model = /* from draft .gmod.cbor */
//! #    Model::new("m");
//! let seed_cpds: Vec<gapsmith_db::SeedCpdRow> = /* from seed_metabolites_edited.tsv */
//! #    Vec::new();
//! let manual = parse_manual_flux("cpd00007:0").unwrap();
//! let medium = predict_medium(&model, &predicted_pathways, &rules, &manual, &seed_cpds).unwrap();
//! ```

pub mod boolexpr;
pub mod predict;
pub mod rules;

pub use boolexpr::{eval, BoolExprError};
pub use predict::{
    parse_manual_flux, predict_medium, MediumLine, MediumPredictError, PredictedMedium,
};
pub use rules::{load_rules, MediumRule, RulesError};
