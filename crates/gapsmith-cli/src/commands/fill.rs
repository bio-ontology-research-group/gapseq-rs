//! `gapsmith fill` — iterative gap-filling via pFBA + KO essentiality.
//!
//! Takes a draft model (from `gapsmith draft`), a medium CSV, and the
//! associated `*-Reactions.tbl` used to build the draft. Emits a gap-filled
//! CBOR + SBML plus an added-reactions TSV.
//!
//! Runs the 4-phase suite (Step 1 + 2 + 2b by default). Add
//! `--full-suite` for Steps 3 (energy-source screen) and 4 (fermentation
//! product screen) — those iterate ~200 exchange compounds each calling
//! pFBA and take 10-30 minutes on a whole bacterium.
//!
//! # Example
//!
//! ```bash
//! # Quick fill on MM_glu
//! gapsmith fill draft.gmod.cbor -n MM_glu.csv -r Reactions.tbl -k 0.01
//!
//! # Full 4-phase suite (slower, richer model)
//! gapsmith fill draft.gmod.cbor -n MM_glu.csv -r Reactions.tbl --full-suite
//!
//! # Step 1 only — useful for debugging
//! gapsmith fill draft.gmod.cbor -n MM_glu.csv -r Reactions.tbl --step1-only
//!
//! # With the opt-in futile-cycle prune (40+ min on large pools)
//! gapsmith fill ... --prune-futile
//! ```
//!
//! See [`crate::commands::adapt`] for post-fill model edits and
//! [`crate::commands::fba`] for solving FBA / pFBA on an existing model
//! (useful for validating a draft before fill).

use clap::Parser;
use gapsmith_db::load_seed_reactions;
use gapsmith_fill::{
    read_medium, read_weights_from_reactions_tbl, run_suite, RxnWeights, SuiteOptions,
    DEFAULT_BCORE, DEFAULT_DUMMY_WEIGHT, DEFAULT_HIGH_EVI,
};
use gapsmith_io::{read_model_cbor, read_model_json, resolve_data_dir, write_model_cbor, ModelFormat};
use gapsmith_sbml::{write_sbml, WriteOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
pub struct Args {
    /// Draft model (`*-draft.gmod.cbor` or `.json`) from `gapsmith draft`.
    pub model: PathBuf,

    /// Medium CSV. Columns: `compounds,name,maxFlux`.
    #[arg(long, short = 'n')]
    pub media: PathBuf,

    /// `*-Reactions.tbl` used to build the draft. Needed to derive per-
    /// reaction gap-fill weights from homology bitscores.
    #[arg(long, short = 'r')]
    pub reactions: Option<PathBuf>,

    /// Target metabolite. Default: `cpd11416` (Biomass).
    #[arg(long, short = 't', default_value = "cpd11416")]
    pub target: String,

    /// Minimum growth rate to accept from gap-filling.
    #[arg(long, short = 'k', default_value_t = 0.01)]
    pub min_growth: f64,

    /// Minimum bitscore for a reaction to count as "core" during the
    /// pFBA-weight assignment.
    #[arg(long, short = 'b', default_value_t = DEFAULT_BCORE)]
    pub bcore: f64,

    /// Bitscore at which a hit hits the cheapest weight floor (`0.005`).
    #[arg(long, default_value_t = DEFAULT_HIGH_EVI)]
    pub high_evi: f64,

    /// Weight for reactions without any blast hit.
    #[arg(long, default_value_t = DEFAULT_DUMMY_WEIGHT)]
    pub dummy_weight: f64,

    /// Output directory (default: current dir).
    #[arg(long, short = 'o', default_value = ".")]
    pub out_dir: PathBuf,

    /// Skip SBML output.
    #[arg(long = "no-sbml")]
    pub no_sbml: bool,

    /// Run only Step 1 (skip Steps 2/2b). Useful for debugging.
    #[arg(long)]
    pub step1_only: bool,

    /// Run the full 4-phase suite (Steps 1 → 4). Default is the quick
    /// suite (Steps 1 + 2 + 2b only). Steps 3 and 4 iterate ~100 carbon
    /// sources / fermentation products each calling pFBA, so the full
    /// suite typically takes 30+ minutes on a whole-bacterium genome.
    #[arg(long = "full-suite")]
    pub full_suite: bool,

    /// Enable the futile-cycle prune before gap-filling. Drops candidate
    /// reactions that saturate at 0.99·max_flux with all boundaries closed.
    /// Expensive on large candidate pools (≈8k rxns → tens of minutes);
    /// off by default.
    #[arg(long = "prune-futile")]
    pub prune_futile: bool,
}

pub fn run_cli(args: Args, data_dir_override: Option<&Path>) -> anyhow::Result<()> {
    let data_dir = resolve_data_dir(data_dir_override)?;

    tracing::info!(
        model = %args.model.display(),
        media = %args.media.display(),
        data_dir = %data_dir.display(),
        "gapsmith fill starting"
    );

    let draft = match ModelFormat::from_path(&args.model) {
        ModelFormat::Cbor => read_model_cbor(&args.model)?,
        ModelFormat::Json => read_model_json(&args.model)?,
    };
    draft.check_shape()?;

    let medium = read_medium(&args.media)?;
    eprintln!("loaded medium: {} compounds", medium.len());

    let weights = if let Some(rxns_path) = &args.reactions {
        let w = read_weights_from_reactions_tbl(
            rxns_path,
            args.bcore,
            args.high_evi,
            args.dummy_weight,
        )?;
        eprintln!(
            "loaded {} SEED rxn weights from {}",
            w.by_seed.len(),
            rxns_path.display()
        );
        w
    } else {
        eprintln!("no --reactions provided — all candidates weighted at dummy_weight");
        RxnWeights {
            by_seed: Default::default(),
            bcore: args.bcore,
            high_evi: args.high_evi,
            dummy_weight: args.dummy_weight,
        }
    };

    let seed_rxns_path = data_dir.join("seed_reactions_corrected.tsv");
    let seed_rxns = load_seed_reactions(&seed_rxns_path)?;
    eprintln!("loaded {} SEED reactions", seed_rxns.len());

    let suite_opts = SuiteOptions {
        target_cpd: args.target.clone(),
        min_growth: args.min_growth,
        bcore: args.bcore,
        high_evi: args.high_evi,
        dummy_weight: args.dummy_weight,
        quick: !args.full_suite,
        step1_only: args.step1_only,
        prune_futile: args.prune_futile,
    };

    let (filled, report) = run_suite(&draft, &medium, &weights, &seed_rxns, &data_dir, &suite_opts)?;

    eprintln!("\nGapfill summary");
    eprintln!("  Step 1  added: {:>3}  growth={:.4}", report.step1_added.len(), report.step1_growth);
    if !args.step1_only {
        eprintln!("  Step 2  added: {:>3}  growth={:.4}", report.step2_added.len(), report.step2_growth);
        eprintln!("  Step 2b added: {:>3}  growth={:.4}", report.step2b_added.len(), report.step2b_growth);
    }
    if args.full_suite && !args.step1_only {
        eprintln!("  Step 3  added: {:>3}  carbon sources tested: {}",
                  report.step3_added.len(), report.carbon_sources.len());
        eprintln!("  Step 4  added: {:>3}  ferm products tested: {}",
                  report.step4_added.len(), report.ferm_products.len());
    }
    eprintln!("  Total added:   {:>3}", report.total_added());
    eprintln!("  Final growth:  {:.6}", report.final_growth);

    // Emit outputs.
    std::fs::create_dir_all(&args.out_dir)?;
    let base = infer_model_name(&args.model);
    let cbor_path = args.out_dir.join(format!("{base}-filled.gmod.cbor"));
    write_model_cbor(&filled, &cbor_path)?;
    eprintln!(
        "\nwrote {} (mets={}, rxns={}, nnz={})",
        cbor_path.display(),
        filled.met_count(),
        filled.rxn_count(),
        filled.s.nnz()
    );

    if !args.no_sbml {
        let xml_path = args.out_dir.join(format!("{base}-filled.xml"));
        let wopts = WriteOptions { pretty: true, ..WriteOptions::default() };
        write_sbml(&filled, &xml_path, &wopts)?;
        eprintln!("wrote {}", xml_path.display());
    }

    // TSV summary of what was added.
    let tsv_path = args.out_dir.join(format!("{base}-filled-added.tsv"));
    let mut tf = std::fs::File::create(&tsv_path)?;
    writeln!(tf, "rxn_id\tstep")?;
    for id in &report.step1_added {
        writeln!(tf, "{id}\t1")?;
    }
    for id in &report.step2_added {
        writeln!(tf, "{id}\t2")?;
    }
    for id in &report.step2b_added {
        writeln!(tf, "{id}\t2b")?;
    }
    for id in &report.step3_added {
        writeln!(tf, "{id}\t3")?;
    }
    for id in &report.step4_added {
        writeln!(tf, "{id}\t4")?;
    }
    eprintln!("wrote {}", tsv_path.display());

    Ok(())
}

fn infer_model_name(path: &Path) -> String {
    let stem = path.file_name().and_then(|s| s.to_str()).unwrap_or("model");
    let cleaned = stem
        .trim_end_matches(".gmod.cbor")
        .trim_end_matches(".cbor")
        .trim_end_matches(".json")
        .trim_end_matches("-draft");
    cleaned.to_string()
}
