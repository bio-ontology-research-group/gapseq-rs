//! `gapsmith draft` — build a draft metabolic model from find /
//! find-transport output.
//!
//! # Example
//!
//! ```bash
//! gapsmith draft \
//!     -r ecoli-all-Reactions.tbl \
//!     -t ecoli-Transporter.tbl \
//!     -b neg \
//!     -o drafts/
//! # → drafts/ecoli-draft.gmod.cbor   native CBOR model
//! #   drafts/ecoli-draft.xml         SBML L3V1 + FBC2 + groups
//! ```
//!
//! Biomass flag `-b auto|pos|neg|archaea|<path/to/custom.json>`:
//! `auto` reads the `gram=` line from the Reactions.tbl header.
//! The SBML output validates with 0 libSBML errors / warnings and
//! loads cleanly in COBRApy.

use clap::Parser;
use gapsmith_draft::{run, DraftOptions};
use gapsmith_io::{resolve_data_dir, write_model_cbor};
use gapsmith_sbml::{write_sbml, WriteOptions};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
pub struct Args {
    /// `*-Reactions.tbl` from `gapsmith find`.
    #[arg(long, short = 'r')]
    pub reactions: PathBuf,

    /// `*-Transporter.tbl` from `gapsmith find-transport`.
    #[arg(long, short = 't')]
    pub transporter: PathBuf,

    /// Biomass template: `auto` (default), `pos`, `neg`, `archaea`,
    /// or a path to a custom JSON.
    #[arg(long, short = 'b', default_value = "auto")]
    pub biomass: String,

    /// Model name / id. Default: the basename of `--reactions` with the
    /// trailing `-<suffix>-Reactions.tbl` stripped.
    #[arg(long, short = 'n')]
    pub name: Option<String>,

    #[arg(long, short = 'u', default_value_t = 200.0)]
    pub high_evi_rxn_bs: f32,

    #[arg(long, short = 'l', default_value_t = 50.0)]
    pub min_bs_for_core: f32,

    #[arg(long, short = 'o', default_value = ".")]
    pub out_dir: PathBuf,

    /// Skip SBML output.
    #[arg(long = "no-sbml")]
    pub no_sbml: bool,
}

pub fn run_cli(args: Args, data_dir_override: Option<&Path>) -> anyhow::Result<()> {
    let data_dir = resolve_data_dir(data_dir_override)?;
    let name = args.name.clone().unwrap_or_else(|| infer_model_name(&args.reactions));
    let opts = DraftOptions {
        model_id: name.clone(),
        biomass: args.biomass,
        high_evi_rxn_bs: args.high_evi_rxn_bs,
        min_bs_for_core: args.min_bs_for_core,
        gapsmith_version: Some(format!("gapsmith {}", env!("CARGO_PKG_VERSION"))),
        seqdb_version: None,
    };

    tracing::info!(
        reactions = %args.reactions.display(),
        transporter = %args.transporter.display(),
        data_dir = %data_dir.display(),
        model = %name,
        "gapsmith draft starting"
    );

    let report = run(&args.reactions, &args.transporter, &data_dir, &opts)?;

    std::fs::create_dir_all(&args.out_dir)?;
    let cbor_path = args.out_dir.join(format!("{name}-draft.gmod.cbor"));
    write_model_cbor(&report.model, &cbor_path)?;
    eprintln!(
        "wrote {} (mets={}, rxns={}, nnz={})",
        cbor_path.display(),
        report.model.met_count(),
        report.model.rxn_count(),
        report.model.s.nnz()
    );

    if !args.no_sbml {
        let xml_path = args.out_dir.join(format!("{name}-draft.xml"));
        let wopts = WriteOptions { pretty: true, ..WriteOptions::default() };
        write_sbml(&report.model, &xml_path, &wopts)?;
        eprintln!("wrote {}", xml_path.display());
    }

    let good_blast = report
        .candidates
        .rows
        .iter()
        .filter(|c| c.best_bitscore.unwrap_or(0.0) >= opts.high_evi_rxn_bs)
        .count();
    eprintln!(
        "  selected {} core SEED reactions; {} had high-evidence hits",
        report.selected_seed_ids.len(),
        good_blast
    );
    Ok(())
}

fn infer_model_name(path: &Path) -> String {
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("draft");
    // Strip `-<tag>-Reactions.tbl` suffix: ecore-PWY-6587-Reactions.tbl
    // → `ecore`.
    if let Some(idx) = stem.find("-Reactions") {
        let head = &stem[..idx];
        if let Some(dash) = head.rfind('-') {
            return head[..dash].to_string();
        }
        return head.to_string();
    }
    stem.to_string()
}
