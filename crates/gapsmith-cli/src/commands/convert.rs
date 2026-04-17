//! `gapsmith convert` — round-trip a model between CBOR and JSON.
//!
//! Used for inspection and format migration. Format is inferred from file
//! extension (`.json` → JSON, anything else → CBOR); override with
//! `--to cbor|json`.

use clap::Parser;
use gapsmith_io::{
    read_model_cbor, read_model_json, write_model_cbor, write_model_json, ModelFormat,
};
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    /// Input model file (`.gmod.cbor`, `.cbor`, or `.json`).
    pub input: PathBuf,

    /// Output path. Extension determines format unless `--to` is set.
    pub output: PathBuf,

    /// Force output format.
    #[arg(long, value_enum)]
    pub to: Option<FormatArg>,

    /// Pretty-print JSON output.
    #[arg(long)]
    pub pretty: bool,
}

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum FormatArg {
    Cbor,
    Json,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let in_fmt = ModelFormat::from_path(&args.input);
    let out_fmt = match args.to {
        Some(FormatArg::Cbor) => ModelFormat::Cbor,
        Some(FormatArg::Json) => ModelFormat::Json,
        None => ModelFormat::from_path(&args.output),
    };

    tracing::info!(
        input = %args.input.display(),
        output = %args.output.display(),
        in_fmt = ?in_fmt,
        out_fmt = ?out_fmt,
        "converting model"
    );

    let model = match in_fmt {
        ModelFormat::Cbor => read_model_cbor(&args.input)?,
        ModelFormat::Json => read_model_json(&args.input)?,
    };

    model.check_shape()?;

    match out_fmt {
        ModelFormat::Cbor => write_model_cbor(&model, &args.output)?,
        ModelFormat::Json => write_model_json(&model, &args.output, args.pretty)?,
    }

    eprintln!(
        "wrote `{}`  ({} mets, {} rxns, {} nnz)",
        args.output.display(),
        model.met_count(),
        model.rxn_count(),
        model.s.nnz()
    );
    Ok(())
}
