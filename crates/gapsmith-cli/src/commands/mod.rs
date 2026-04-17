//! Subcommand implementations.
//!
//! Each submodule owns the clap `Args` struct plus its `run` function. The
//! binary's `main` only dispatches — the heavy lifting lives here.

pub mod adapt;
pub mod align;
pub mod batch_align;
pub mod community;
pub mod convert;
pub mod db;
pub mod doall;
pub mod doall_batch;
pub mod draft;
pub mod example_model;
pub mod export_sbml;
pub mod fba;
pub mod fill;
pub mod find;
pub mod find_transport;
pub mod medium;
pub mod pan;
pub mod test;
pub mod update_data;
pub mod update_sequences;
