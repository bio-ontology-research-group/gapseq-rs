//! `gapsmith` command-line entry point.
//!
//! Mirrors the bash-dispatched `gapseq` subcommand structure (upstream). Every subcommand
//! implementation lives in [`commands`] so the binary stays a thin dispatch
//! layer; all logic is library code that can be tested without spawning a
//! process.
//!
//! # Subcommands
//!
//! | Command | Module | Purpose |
//! |---|---|---|
//! | `doall` | [`commands::doall`] | Full pipeline: find → find-transport → draft → medium → fill |
//! | `find` | [`commands::find`] | Pathway / reaction detection |
//! | `find-transport` | [`commands::find_transport`] | Transporter detection |
//! | `draft` | [`commands::draft`] | Build a draft metabolic model |
//! | `medium` | [`commands::medium`] | Rule-based medium inference |
//! | `fill` | [`commands::fill`] | Iterative gap-filling (pFBA + KO) |
//! | `adapt` | [`commands::adapt`] | Edit reactions / force growth |
//! | `pan` | [`commands::pan`] | Pan-draft union from N drafts |
//! | `update-sequences` | [`commands::update_sequences`] | Zenodo seqdb sync |
//! | `fba` | [`commands::fba`] | FBA / pFBA on an existing model |
//! | `align` | [`commands::align`] | Standalone aligner debug wrapper |
//! | `batch-align` | [`commands::batch_align`] | Cluster N genomes + single alignment |
//! | `convert` | [`commands::convert`] | CBOR ↔ JSON round-trip |
//! | `export-sbml` | [`commands::export_sbml`] | Write SBML L3V1+FBC2+groups |
//! | `db inspect` | [`commands::db`] | Reference-data smoke test |
//! | `example-model` | [`commands::example_model`] | Toy model fixture |
//! | `test` | [`commands::test`] | Report resolved paths + tool availability |
//!
//! Each subcommand's module-level docstring contains a runnable example.

mod commands;

use clap::{Parser, Subcommand};
use std::path::PathBuf;

/// Top-level CLI definition.
///
/// New subcommands arrive per milestone — see `gapsmith/README.md`. Anything
/// not yet wired up returns a "not implemented yet" error instead of silently
/// doing nothing.
#[derive(Debug, Parser)]
#[command(name = "gapsmith", version, about, long_about = None)]
pub struct Cli {
    /// Override the reference-data directory (default: auto-detect).
    #[arg(long, global = true)]
    pub data_dir: Option<PathBuf>,

    /// Override the sequence database directory (default: `<data-dir>/seq`).
    #[arg(long, global = true)]
    pub seq_dir: Option<PathBuf>,

    /// Number of worker threads (default: all cores).
    #[arg(long, short = 'K', global = true)]
    pub threads: Option<usize>,

    /// Increase log verbosity (-v info, -vv debug, -vvv trace).
    #[arg(long, short = 'v', global = true, action = clap::ArgAction::Count)]
    pub verbose: u8,

    #[command(subcommand)]
    pub cmd: Cmd,
}

#[derive(Debug, Subcommand)]
pub enum Cmd {
    /// Convert a model file between CBOR and JSON (round-trip utility).
    Convert(commands::convert::Args),

    /// Print resolved paths and tool versions; used by `gapsmith test`.
    Test(commands::test::Args),

    /// Emit a tiny hand-built model (for smoke tests and downstream tooling).
    #[command(name = "example-model")]
    ExampleModel(commands::example_model::Args),

    /// Introspect the reference-data tables under `--data-dir`.
    Db(commands::db::Args),

    /// Write a model as SBML Level 3 V1 + FBC2 + groups.
    #[command(name = "export-sbml")]
    ExportSbml(commands::export_sbml::Args),

    /// Run a sequence aligner (blast/diamond/mmseqs2) or parse a precomputed TSV.
    Align(commands::align::Args),

    // The following subcommands are intentionally scaffolded but not yet
    // implemented. They return `NotImplemented` so users see a clear signal
    // rather than misleading success. Implemented per milestone M2..M10.
    /// Pathway & reaction detection for a genome.
    Find(commands::find::Args),
    /// Transporter detection for a genome.
    #[command(name = "find-transport")]
    FindTransport(commands::find_transport::Args),
    /// Build a draft metabolic model from `find` + `find-transport` output.
    Draft(commands::draft::Args),
    /// Solve FBA / pFBA on a model (useful for validating a draft pre-fill).
    Fba(commands::fba::Args),
    /// Iterative gap-filling via pFBA + KO essentiality.
    Fill(commands::fill::Args),
    /// Rule-based medium inference from a draft + Pathways.tbl.
    Medium(commands::medium::Args),
    /// Chain find → find-transport → draft → medium → fill end-to-end.
    Doall(commands::doall::Args),
    /// Add / remove reactions or force growth on an existing model.
    Adapt(commands::adapt::Args),
    /// Build a pan-draft model from N drafts.
    Pan(commands::pan::Args),
    /// Sync the reference sequence database from Zenodo.
    #[command(name = "update-sequences")]
    UpdateSequences(commands::update_sequences::Args),
    /// Cluster N genomes, align once, expand per-genome TSVs.
    #[command(name = "batch-align")]
    BatchAlign(commands::batch_align::Args),
}

fn init_logging(verbose: u8) {
    let level = match verbose {
        0 => "gapsmith=warn",
        1 => "gapsmith=info",
        2 => "gapsmith=debug",
        _ => "gapsmith=trace",
    };
    let filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new(level));
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .with_writer(std::io::stderr)
        .init();
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose);

    match cli.cmd {
        Cmd::Convert(args) => commands::convert::run(args),
        Cmd::Test(args) => commands::test::run(args, cli.data_dir.as_deref()),
        Cmd::ExampleModel(args) => commands::example_model::run(args),
        Cmd::Db(args) => commands::db::run(args, cli.data_dir.as_deref()),
        Cmd::ExportSbml(args) => commands::export_sbml::run(args),
        Cmd::Align(args) => commands::align::run(args),
        Cmd::BatchAlign(args) => commands::batch_align::run(args),
        Cmd::Find(args) => commands::find::run(
            args,
            cli.data_dir.as_deref(),
            cli.seq_dir.as_deref(),
        ),
        Cmd::FindTransport(args) => commands::find_transport::run_cli(
            args,
            cli.data_dir.as_deref(),
            cli.seq_dir.as_deref(),
        ),
        Cmd::Draft(args) => commands::draft::run_cli(args, cli.data_dir.as_deref()),
        Cmd::Fba(args) => commands::fba::run(args),
        Cmd::Fill(args) => commands::fill::run_cli(args, cli.data_dir.as_deref()),
        Cmd::Medium(args) => commands::medium::run_cli(args, cli.data_dir.as_deref()),
        Cmd::Doall(args) => commands::doall::run_cli(
            args,
            cli.data_dir.as_deref(),
            cli.seq_dir.as_deref(),
        ),
        Cmd::Adapt(args) => commands::adapt::run_cli(args, cli.data_dir.as_deref()),
        Cmd::Pan(args) => commands::pan::run_cli(args),
        Cmd::UpdateSequences(args) => commands::update_sequences::run_cli(
            args,
            cli.data_dir.as_deref(),
            cli.seq_dir.as_deref(),
        ),
        other => {
            anyhow::bail!(
                "subcommand `{}` is scaffolded but not yet implemented in this milestone",
                subcommand_label(&other)
            );
        }
    }
}

fn subcommand_label(c: &Cmd) -> &'static str {
    match c {
        Cmd::Convert(_) => "convert",
        Cmd::Test(_) => "test",
        Cmd::ExampleModel(_) => "example-model",
        Cmd::Db(_) => "db",
        Cmd::ExportSbml(_) => "export-sbml",
        Cmd::Align(_) => "align",
        Cmd::BatchAlign(_) => "batch-align",
        Cmd::Find(_) => "find",
        Cmd::FindTransport(_) => "find-transport",
        Cmd::Draft(_) => "draft",
        Cmd::Fba(_) => "fba",
        Cmd::Fill(_) => "fill",
        Cmd::Medium(_) => "medium",
        Cmd::Doall(_) => "doall",
        Cmd::Adapt(_) => "adapt",
        Cmd::Pan(_) => "pan",
        Cmd::UpdateSequences(_) => "update-sequences",
    }
}
