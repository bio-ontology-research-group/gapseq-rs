//! `gapsmith export-sbml` — read a CBOR/JSON model and emit SBML L3V1 + FBC2.
//!
//! Kept as its own subcommand (rather than folded into `convert`) because
//! SBML is write-only here — gapsmith loads its models from CBOR/JSON, not
//! SBML.

use clap::{Parser, ValueEnum};
use gapsmith_io::{read_model_cbor, read_model_json, ModelFormat};
use gapsmith_sbml::{write_sbml, ObjectiveSense, WriteOptions};
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    /// Input model (`.gmod.cbor`, `.cbor`, or `.json`).
    pub input: PathBuf,

    /// Output SBML path (`.xml`).
    pub output: PathBuf,

    /// Objective identifier to emit. Default: `obj`.
    #[arg(long, default_value = "obj")]
    pub objective_id: String,

    /// Objective sense.
    #[arg(long, value_enum, default_value_t = SenseArg::Maximize)]
    pub objective_sense: SenseArg,

    /// Write a compact (non-indented) document.
    #[arg(long)]
    pub compact: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum SenseArg {
    Maximize,
    Minimize,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let model = match ModelFormat::from_path(&args.input) {
        ModelFormat::Cbor => read_model_cbor(&args.input)?,
        ModelFormat::Json => read_model_json(&args.input)?,
    };
    let opts = WriteOptions {
        pretty: !args.compact,
        objective_id: args.objective_id,
        objective_sense: match args.objective_sense {
            SenseArg::Maximize => ObjectiveSense::Maximize,
            SenseArg::Minimize => ObjectiveSense::Minimize,
        },
    };
    write_sbml(&model, &args.output, &opts)?;
    eprintln!(
        "wrote `{}`  ({} mets, {} rxns)",
        args.output.display(),
        model.met_count(),
        model.rxn_count()
    );
    Ok(())
}
