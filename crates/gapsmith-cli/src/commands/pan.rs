//! `gapsmith pan` — pan-draft model from N drafts.
//!
//! # Example
//!
//! ```bash
//! # Union from a directory of drafts, default 6% frequency threshold
//! gapsmith pan -m drafts/ -f pan_out/
//!
//! # Explicit list with a stricter threshold
//! gapsmith pan -m a-draft.gmod.cbor,b-draft.gmod.cbor,c-draft.gmod.cbor -t 0.5
//!
//! # Just the binary presence/absence table — no pan-draft model
//! gapsmith pan -m drafts/ -b -f pan_out/
//! ```
//!
//! Outputs:
//! - `pan-draft.gmod.cbor` + `pan-draft.xml` — the union model.
//! - `pan-rxn-presence.tsv` — rxn × model × frequency table.
//!
//! Port of `src/pan-draft.R` (plain-union variant). For every SEED
//! reaction found in any input draft, count how many drafts contain it;
//! retain those whose presence frequency ≥ `--min-freq` (default 0.06, the
//! gapseq default). The retained reactions — plus all metabolites they
//! touch, all exchanges in any draft, and bio1 from the first draft —
//! form the pan-draft.
//!
//! Not yet ported:
//!
//! - Per-reaction weight medianing (`custom_median` in
//!   `pan-draft_functions.R`) — we emit the pan-model without merged
//!   rxnWeights metadata, which `gapsmith fill` would need to re-derive
//!   from the source `*-Reactions.tbl` files.

use clap::Parser;
use gapsmith_core::{Model, StoichMatrix};
use gapsmith_io::{read_model_cbor, read_model_json, write_model_cbor, ModelFormat};
use gapsmith_sbml::{write_sbml, WriteOptions};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    /// Comma-separated or directory of draft model files
    /// (`*.gmod.cbor` or `*.json`). Glob-like patterns are expanded.
    #[arg(long, short = 'm')]
    pub models: String,

    /// Minimum reaction frequency across drafts (0.0–1.0) to include a
    /// reaction in the pan-draft. Default `0.06`.
    #[arg(long, short = 't', default_value_t = 0.06)]
    pub min_freq: f64,

    /// Only emit the binary presence/absence table, skip the pan-draft model.
    #[arg(long, short = 'b')]
    pub only_binary: bool,

    /// Output directory.
    #[arg(long, short = 'f', default_value = ".")]
    pub out_dir: PathBuf,

    /// Skip SBML output.
    #[arg(long = "no-sbml")]
    pub no_sbml: bool,
}

pub fn run_cli(args: Args) -> anyhow::Result<()> {
    let paths = expand_paths(&args.models)?;
    anyhow::ensure!(!paths.is_empty(), "no model files found for {}", args.models);
    eprintln!("loading {} models", paths.len());

    let models: Vec<(String, Model)> = paths
        .into_iter()
        .map(|p| -> anyhow::Result<(String, Model)> {
            let id = p
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("model")
                .trim_end_matches(".gmod.cbor")
                .trim_end_matches(".cbor")
                .trim_end_matches(".json")
                .trim_end_matches("-draft")
                .to_string();
            let m = match ModelFormat::from_path(&p) {
                ModelFormat::Cbor => read_model_cbor(&p)?,
                ModelFormat::Json => read_model_json(&p)?,
            };
            Ok((id, m))
        })
        .collect::<Result<_, _>>()?;

    // Frequency table: rxn_id → count.
    let n = models.len();
    let mut counts: HashMap<String, usize> = HashMap::new();
    for (_, m) in &models {
        let seen_here: HashSet<&str> = m.rxns.iter().map(|r| r.id.as_str()).collect();
        for rxn_id in seen_here {
            *counts.entry(rxn_id.to_string()).or_insert(0) += 1;
        }
    }

    std::fs::create_dir_all(&args.out_dir)?;

    // -- Binary table --
    let mut all_rxns: Vec<String> = counts.keys().cloned().collect();
    all_rxns.sort();
    let table_path = args.out_dir.join("pan-rxn-presence.tsv");
    let mut f = std::fs::File::create(&table_path)?;
    write!(f, "rxn_id\tfreq")?;
    for (id, _) in &models {
        write!(f, "\t{id}")?;
    }
    writeln!(f)?;
    for rxn_id in &all_rxns {
        let c = counts.get(rxn_id).copied().unwrap_or(0);
        let freq = c as f64 / n as f64;
        write!(f, "{rxn_id}\t{freq:.3}")?;
        for (_, m) in &models {
            let present = m.rxns.iter().any(|r| r.id.as_str() == rxn_id);
            write!(f, "\t{}", if present { 1 } else { 0 })?;
        }
        writeln!(f)?;
    }
    eprintln!("wrote {} ({} rxns × {} models)", table_path.display(), all_rxns.len(), n);

    if args.only_binary {
        return Ok(());
    }

    // -- Pan-draft model --
    let retained: HashSet<String> = counts
        .iter()
        .filter(|(_, c)| (**c as f64 / n as f64) >= args.min_freq)
        .map(|(k, _)| k.clone())
        .collect();
    eprintln!("retained {} reactions above {:.2} frequency", retained.len(), args.min_freq);

    let pan = build_pan_model(&models, &retained)?;
    let cbor_path = args.out_dir.join("pan-draft.gmod.cbor");
    write_model_cbor(&pan, &cbor_path)?;
    eprintln!(
        "wrote {} (mets={}, rxns={}, nnz={})",
        cbor_path.display(),
        pan.met_count(),
        pan.rxn_count(),
        pan.s.nnz()
    );
    if !args.no_sbml {
        let xml_path = args.out_dir.join("pan-draft.xml");
        let wopts = WriteOptions { pretty: true, ..WriteOptions::default() };
        write_sbml(&pan, &xml_path, &wopts)?;
        eprintln!("wrote {}", xml_path.display());
    }

    Ok(())
}

/// Expand a comma-separated list, directory path, or single file into a
/// flat `Vec<PathBuf>` of model files.
fn expand_paths(spec: &str) -> anyhow::Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    for piece in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let p = PathBuf::from(piece);
        if p.is_dir() {
            for entry in std::fs::read_dir(&p)? {
                let entry = entry?;
                let p = entry.path();
                if let Some(ext) = p.file_name().and_then(|s| s.to_str()) {
                    if ext.ends_with(".gmod.cbor") || ext.ends_with(".cbor") || ext.ends_with(".json") {
                        out.push(p);
                    }
                }
            }
        } else if p.is_file() {
            out.push(p);
        } else {
            anyhow::bail!("path not found: {}", p.display());
        }
    }
    out.sort();
    Ok(out)
}

/// Assemble a pan-draft from `models`, keeping only reactions listed in
/// `retained`. bio1 (if present in any model) is taken from the first.
fn build_pan_model(
    models: &[(String, Model)],
    retained: &HashSet<String>,
) -> anyhow::Result<Model> {
    let mut pan = Model::new("pan_draft");
    if let Some((_, first)) = models.first() {
        pan.compartments = first.compartments.clone();
        pan.annot = first.annot.clone();
    }
    pan.annot.id = "pan_draft".into();

    // Collect mets lazily as rxns are added.
    let mut met_ix: HashMap<String, usize> = HashMap::new();

    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();

    // 1. Retained reactions (plus any EX / bio1 from every model).
    let mut added: HashSet<String> = HashSet::new();
    let mut add_reaction = |pan: &mut Model,
                            met_ix: &mut HashMap<String, usize>,
                            triplets: &mut Vec<(usize, usize, f64)>,
                            model: &Model,
                            rxn_idx: usize| {
        let r = &model.rxns[rxn_idx];
        if added.contains(r.id.as_str()) {
            return;
        }
        let new_col = pan.rxns.len();
        pan.rxns.push(r.clone());
        for (row, coef) in model.s.column(rxn_idx) {
            let met = &model.mets[row];
            let met_id = met.id.as_str();
            let idx = match met_ix.get(met_id) {
                Some(&i) => i,
                None => {
                    pan.mets.push(met.clone());
                    met_ix.insert(met_id.to_string(), pan.mets.len() - 1);
                    pan.mets.len() - 1
                }
            };
            triplets.push((idx, new_col, coef));
        }
        added.insert(r.id.as_str().to_string());
    };

    for (_, m) in models {
        for (i, r) in m.rxns.iter().enumerate() {
            let keep = retained.contains(r.id.as_str())
                || r.is_exchange
                || r.id.as_str().starts_with("EX_")
                || r.id.as_str().starts_with("DM_")
                || r.is_biomass
                || r.id.as_str() == "bio1";
            if keep {
                add_reaction(&mut pan, &mut met_ix, &mut triplets, m, i);
            }
        }
    }

    let n_mets = pan.mets.len();
    let n_rxns = pan.rxns.len();
    pan.s = StoichMatrix::from_triplets(n_mets, n_rxns, triplets);
    Ok(pan)
}
