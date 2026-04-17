//! `gapsmith medium` — rule-based medium inference.
//!
//! Given a draft model + the companion `*-Pathways.tbl`, evaluate
//! `dat/medium_prediction_rules.tsv` and emit a medium CSV suitable for
//! `gapsmith fill -n`.
//!
//! # Example
//!
//! ```bash
//! gapsmith medium \
//!     -m drafts/ecoli-draft.gmod.cbor \
//!     -p ecoli-all-Pathways.tbl \
//!     -o ecoli-medium.csv
//!
//! # Manual overrides: force oxygen to 0 + disable CO2
//! gapsmith medium ... -c "cpd00007:0;cpd00011:0"
//! ```
//!
//! # Algorithm
//!
//! Each row of `dat/medium_prediction_rules.tsv` has a `rule` column
//! containing a boolean expression over pathway ids (in double quotes),
//! SEED reaction ids (`rxnNNNNN`), and SEED compound ids (`cpdNNNNN`).
//! Supported operators: `| & !` plus `< > <= >= == +` for "count-of-TRUE"
//! rules like `rxn02213 + rxn01255 < 3`. For every rule that evaluates
//! true, its `cpd.id` and `maxFlux` enter the medium. Duplicate cpds
//! across multiple triggered rules are averaged; Saccharides / Organic
//! acids categories keep only the first rule. Finally the proton
//! balancer adds H+ if the net charge of contributing entries is
//! negative.
//!
//! **Parity:** byte-identical to `toy/ecoli-medium.csv` on the upstream
//! toy ecoli example (21 compounds).

use clap::Parser;
use gapsmith_db::load_seed_metabolites;
use gapsmith_io::{read_model_cbor, read_model_json, resolve_data_dir, ModelFormat};
use gapsmith_medium::{load_rules, parse_manual_flux, predict_medium};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
pub struct Args {
    /// Draft model (`.gmod.cbor`, `.cbor`, or `.json`).
    #[arg(long, short = 'm')]
    pub model: PathBuf,

    /// `*-Pathways.tbl` emitted by `gapsmith find`.
    #[arg(long, short = 'p')]
    pub pathways: PathBuf,

    /// Optional manual flux overrides: `cpd00007:18.5;cpd00011:0`.
    #[arg(long, short = 'c')]
    pub manual_flux: Option<String>,

    /// Output CSV path. Default: replace the input stem's `-draft` suffix
    /// with `-medium.csv` in the `--out-dir`.
    #[arg(long, short = 'o')]
    pub output: Option<PathBuf>,

    /// Output directory (default: current dir).
    #[arg(long, short = 'f', default_value = ".")]
    pub out_dir: PathBuf,
}

pub fn run_cli(args: Args, data_dir_override: Option<&Path>) -> anyhow::Result<()> {
    let data_dir = resolve_data_dir(data_dir_override)?;

    let model = match ModelFormat::from_path(&args.model) {
        ModelFormat::Cbor => read_model_cbor(&args.model)?,
        ModelFormat::Json => read_model_json(&args.model)?,
    };
    model.check_shape()?;

    let predicted = read_predicted_pathways(&args.pathways)?;
    eprintln!("predicted pathways: {}", predicted.len());

    let rules_path = data_dir.join("medium_prediction_rules.tsv");
    let rules = load_rules(&rules_path)?;
    eprintln!("rules loaded: {}", rules.len());

    let seed_cpds = load_seed_metabolites(data_dir.join("seed_metabolites_edited.tsv"))?;

    let manual = match args.manual_flux.as_deref() {
        Some(s) => parse_manual_flux(s)?,
        None => Vec::new(),
    };

    let medium = predict_medium(&model, &predicted, &rules, &manual, &seed_cpds)?;

    std::fs::create_dir_all(&args.out_dir)?;
    let out_path = match args.output {
        Some(p) => p,
        None => {
            let stem = args
                .model
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("model");
            let cleaned = stem
                .trim_end_matches(".gmod.cbor")
                .trim_end_matches(".cbor")
                .trim_end_matches(".json")
                .trim_end_matches("-draft");
            args.out_dir.join(format!("{cleaned}-medium.csv"))
        }
    };

    let mut f = std::fs::File::create(&out_path)?;
    medium.write_csv(&mut f)?;
    eprintln!("wrote {} ({} compounds)", out_path.display(), medium.compounds.len());
    Ok(())
}

/// Parse the `Prediction` column from a `*-Pathways.tbl` and return the
/// set of MetaCyc-style pathway ids whose prediction is `TRUE`.
fn read_predicted_pathways(path: &Path) -> anyhow::Result<HashSet<String>> {
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(path)?;
    let r = BufReader::new(f);
    let mut iter = r.lines();
    // Skip comment rows (begin with `#`). The first non-comment line is
    // the header.
    let header = loop {
        match iter.next() {
            Some(Ok(line)) => {
                if line.starts_with('#') || line.trim().is_empty() {
                    continue;
                }
                break line;
            }
            Some(Err(e)) => return Err(e.into()),
            None => anyhow::bail!("Pathways.tbl is empty"),
        }
    };
    let cols: Vec<&str> = header.split('\t').collect();
    let i_id = cols
        .iter()
        .position(|c| *c == "ID")
        .ok_or_else(|| anyhow::anyhow!("missing ID column"))?;
    let i_pred = cols
        .iter()
        .position(|c| *c == "Prediction")
        .ok_or_else(|| anyhow::anyhow!("missing Prediction column"))?;
    let mut out = HashSet::new();
    for line in iter.map_while(Result::ok) {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() <= i_id.max(i_pred) {
            continue;
        }
        if matches!(parts[i_pred], "TRUE" | "True" | "true") {
            out.insert(parts[i_id].to_string());
        }
    }
    Ok(out)
}
