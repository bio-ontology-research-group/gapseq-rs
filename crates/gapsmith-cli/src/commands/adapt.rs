//! `gapsmith adapt` — add or remove reactions / pathways on an existing model.
//!
//! # Example
//!
//! ```bash
//! # Add specific SEED reactions
//! gapsmith adapt -m model.gmod.cbor -a rxn00011,rxn00102
//!
//! # Add every reaction in a MetaCyc pathway
//! gapsmith adapt -m model.gmod.cbor -a "PWY-6587"
//!
//! # Remove reactions
//! gapsmith adapt -m model.gmod.cbor -r rxn00011
//!
//! # Force growth on glucose — runs gap-fill
//! gapsmith adapt -m model.gmod.cbor -w cpd00027:TRUE -b Reactions.tbl
//!
//! # Disable growth on glucose — closes the EX uptake bound
//! gapsmith adapt -m model.gmod.cbor -w cpd00027:FALSE
//! ```
//!
//! Scope (subset of R `adapt.R`):
//!
//! - `-a <ids>` — add reactions. Each id is resolved to SEED rxn id(s) via:
//!   (1) direct SEED `rxn.....` match, (2) MetaCyc pathway id (add every
//!   reaction in the pathway's `reaId` column).
//! - `-r <ids>` — remove reactions. Resolved the same way.
//! - `-w <cpd>:TRUE|FALSE` — force growth / no-growth on a compound via
//!   gap-fill or KO. Uses the existing `gapfill4` heuristic.
//!
//! Not yet ported:
//!
//! - EC-number / KEGG / enzyme-name resolution (handled in R by
//!   `getDBhit`; our dbhit index is in `gapsmith-find` but plumbing it into
//!   adapt is deferred).

use clap::Parser;
use gapsmith_core::{CompartmentId, Metabolite, Model, Reaction, StoichMatrix};
use gapsmith_db::{load_seed_reactions, pathway::PwySource, PathwayTable, SeedRxnRow};
use gapsmith_fill::{
    apply_medium_to_full, build_full_model, fba, gapfill4, read_medium,
    read_weights_from_reactions_tbl, FbaOptions, GapfillOptions, RxnWeights, SolveStatus,
    DEFAULT_BCORE, DEFAULT_DUMMY_WEIGHT, DEFAULT_HIGH_EVI,
};
use gapsmith_io::{read_model_cbor, read_model_json, resolve_data_dir, write_model_cbor, ModelFormat};
use gapsmith_sbml::{write_sbml, WriteOptions};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
pub struct Args {
    /// Model file (`.gmod.cbor` or `.json`).
    #[arg(long, short = 'm')]
    pub model: PathBuf,

    /// Comma-separated reactions / pathway ids to add.
    #[arg(long, short = 'a')]
    pub add: Option<String>,

    /// Comma-separated reactions / pathway ids to remove.
    #[arg(long, short = 'r')]
    pub remove: Option<String>,

    /// Growth adjustments, comma-separated `cpdNNNNN:TRUE|FALSE`.
    #[arg(long, short = 'w')]
    pub growth: Option<String>,

    /// `*-Reactions.tbl` needed for growth-adjustment gap-filling (to
    /// derive pFBA weights).
    #[arg(long, short = 'b')]
    pub reactions: Option<PathBuf>,

    /// Minimum growth rate target for growth-adjustment fills.
    #[arg(long, short = 'k', default_value_t = 0.01)]
    pub min_growth: f64,

    /// Output directory.
    #[arg(long, short = 'f', default_value = ".")]
    pub out_dir: PathBuf,

    /// Skip SBML output.
    #[arg(long = "no-sbml")]
    pub no_sbml: bool,
}

pub fn run_cli(args: Args, data_dir_override: Option<&Path>) -> anyhow::Result<()> {
    if args.add.is_none() && args.remove.is_none() && args.growth.is_none() {
        anyhow::bail!("pass at least one of -a (add), -r (remove), -w (growth)");
    }

    let data_dir = resolve_data_dir(data_dir_override)?;
    let mut model = match ModelFormat::from_path(&args.model) {
        ModelFormat::Cbor => read_model_cbor(&args.model)?,
        ModelFormat::Json => read_model_json(&args.model)?,
    };
    model.check_shape()?;

    let seed_rxns = load_seed_reactions(data_dir.join("seed_reactions_corrected.tsv"))?;

    // Pathway resolver — MetaCyc + custom merge.
    let mut pwys = PathwayTable::load(&data_dir.join("meta_pwy.tbl"), PwySource::MetaCyc)
        .map(|t| t.rows)
        .unwrap_or_default();
    if let Ok(t) = PathwayTable::load(&data_dir.join("custom_pwy.tbl"), PwySource::Custom) {
        for r in t.rows {
            pwys.retain(|p| p.id != r.id);
            pwys.push(r);
        }
    }
    let resolver = IdResolver { pathways: &pwys };

    if let Some(ids) = &args.add {
        let seeds = resolver.resolve_many(ids);
        eprintln!("resolving -a: {} → {} SEED rxns", ids, seeds.len());
        let added = add_reactions(&mut model, &seeds, &seed_rxns);
        eprintln!("added: {}", added.join(", "));
    }

    if let Some(ids) = &args.remove {
        let seeds = resolver.resolve_many(ids);
        eprintln!("resolving -r: {} → {} SEED rxns", ids, seeds.len());
        let removed = remove_reactions(&mut model, &seeds);
        eprintln!("removed: {}", removed.join(", "));
    }

    if let Some(growth_spec) = &args.growth {
        let entries = parse_growth_spec(growth_spec)?;
        let weights = match &args.reactions {
            Some(rp) => read_weights_from_reactions_tbl(
                rp,
                DEFAULT_BCORE,
                DEFAULT_HIGH_EVI,
                DEFAULT_DUMMY_WEIGHT,
            )?,
            None => {
                eprintln!("warning: no --reactions provided; all candidates weighted at dummy_weight");
                RxnWeights::new()
            }
        };
        for (cpd, want_growth) in entries {
            adapt_growth(
                &mut model,
                &cpd,
                want_growth,
                &weights,
                &seed_rxns,
                args.min_growth,
            )?;
        }
    }

    // -- Write outputs --
    std::fs::create_dir_all(&args.out_dir)?;
    let stem = args
        .model
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("model")
        .trim_end_matches(".gmod.cbor")
        .trim_end_matches(".cbor")
        .trim_end_matches(".json")
        .trim_end_matches("-draft")
        .to_string();
    let cbor_path = args.out_dir.join(format!("{stem}-adapt.gmod.cbor"));
    write_model_cbor(&model, &cbor_path)?;
    eprintln!(
        "wrote {} (mets={}, rxns={}, nnz={})",
        cbor_path.display(),
        model.met_count(),
        model.rxn_count(),
        model.s.nnz()
    );
    if !args.no_sbml {
        let xml_path = args.out_dir.join(format!("{stem}-adapt.xml"));
        let wopts = WriteOptions { pretty: true, ..WriteOptions::default() };
        write_sbml(&model, &xml_path, &wopts)?;
        eprintln!("wrote {}", xml_path.display());
    }

    Ok(())
}

struct IdResolver<'a> {
    pathways: &'a [gapsmith_db::PathwayRow],
}

impl IdResolver<'_> {
    fn resolve_many(&self, ids_csv: &str) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        let mut seen: HashSet<String> = HashSet::new();
        for id in ids_csv.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            for seed in self.resolve_one(id) {
                if seen.insert(seed.clone()) {
                    out.push(seed);
                }
            }
        }
        out
    }

    fn resolve_one(&self, id: &str) -> Vec<String> {
        // 1. Direct SEED rxn id.
        if id.starts_with("rxn") && id[3..].chars().all(|c| c.is_ascii_digit()) {
            return vec![id.to_string()];
        }
        // 2. MetaCyc pathway id (strip leading/trailing `|`).
        let stripped = id.trim_matches('|');
        if let Some(p) = self
            .pathways
            .iter()
            .find(|p| p.id.trim_matches('|') == stripped)
        {
            // `reaId` is a comma-separated list. We match MetaCyc rxn
            // names against SEED db; for simplicity, treat entries that
            // look like rxnNNNNN directly and defer the rest.
            let mut out = Vec::new();
            for r in p.rea_id.split(',').map(str::trim) {
                if r.starts_with("rxn") && r.len() >= 4 {
                    out.push(r.to_string());
                }
            }
            return out;
        }
        Vec::new()
    }
}

fn add_reactions(model: &mut Model, seeds: &[String], seed_rxns: &[SeedRxnRow]) -> Vec<String> {
    let existing: HashSet<String> = model
        .rxns
        .iter()
        .map(|r| r.id.as_str().trim_end_matches("_c0").to_string())
        .collect();
    let mut added = Vec::new();
    for seed in seeds {
        if existing.contains(seed) {
            continue;
        }
        if let Some(row) = seed_rxns.iter().find(|r| r.id.as_str() == seed) {
            gapsmith_draft::builder::add_seed_reaction(model, row, Some(10));
            added.push(format!("{seed}_c0"));
        }
    }
    gapsmith_draft::builder::rebuild_s_matrix(model);
    added
}

fn remove_reactions(model: &mut Model, seeds: &[String]) -> Vec<String> {
    let removal_ids: HashSet<String> =
        seeds.iter().map(|s| format!("{s}_c0")).collect();
    let present: Vec<String> = model
        .rxns
        .iter()
        .filter(|r| removal_ids.contains(r.id.as_str()))
        .map(|r| r.id.as_str().to_string())
        .collect();
    gapsmith_fill::drop_reactions(model, &present.iter().cloned().collect());
    present
}

fn parse_growth_spec(s: &str) -> anyhow::Result<Vec<(String, bool)>> {
    let mut out = Vec::new();
    for term in s.split(',').map(str::trim).filter(|t| !t.is_empty()) {
        let (cpd, flag) = term
            .split_once(':')
            .ok_or_else(|| anyhow::anyhow!("bad growth spec `{term}` (expected cpdNNNNN:TRUE|FALSE)"))?;
        let cpd = cpd.trim();
        let flag = flag.trim();
        if !cpd.starts_with("cpd") {
            anyhow::bail!("growth spec must start with cpd: `{term}`");
        }
        let want = match flag.to_ascii_uppercase().as_str() {
            "TRUE" | "T" | "1" => true,
            "FALSE" | "F" | "0" => false,
            _ => anyhow::bail!("bad growth flag `{flag}` in `{term}`"),
        };
        out.push((cpd.to_string(), want));
    }
    Ok(out)
}

/// Force `want_growth` on `cpd`: if TRUE, add an EX reaction for the cpd
/// at `[-10, 1000]`, ensure a transporter exists, and gap-fill. If FALSE,
/// close the EX reaction's lower bound.
fn adapt_growth(
    model: &mut Model,
    cpd: &str,
    want_growth: bool,
    weights: &RxnWeights,
    seed_rxns: &[SeedRxnRow],
    min_growth: f64,
) -> anyhow::Result<()> {
    let ex_id = format!("EX_{cpd}_e0");
    let met_e = format!("{cpd}_e0");

    if !want_growth {
        if let Some(r) = model.rxns.iter_mut().find(|r| r.id.as_str() == ex_id) {
            r.lb = 0.0;
            eprintln!("closed uptake on {ex_id}");
        } else {
            eprintln!("model already has no uptake on {ex_id}");
        }
        return Ok(());
    }

    // Ensure the external metabolite exists.
    let met_idx = match model.mets.iter().position(|m| m.id.as_str() == met_e) {
        Some(i) => i,
        None => {
            let m = Metabolite::new(met_e.as_str(), cpd, CompartmentId::EXTRACELLULAR);
            model.mets.push(m);
            model.mets.len() - 1
        }
    };

    if !model.rxns.iter().any(|r| r.id.as_str() == ex_id) {
        let mut r = Reaction::new(ex_id.as_str(), format!("{cpd} uptake"), -10.0, 1000.0);
        r.is_exchange = true;
        r.gs_origin = Some(7);
        model.rxns.push(r);
        let n_mets = model.mets.len();
        let n_rxns = model.rxns.len();
        let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(model.s.nnz() + 1);
        let old_cols = model.s.cols();
        for c in 0..old_cols.min(n_rxns - 1) {
            for (row, v) in model.s.column(c) {
                triplets.push((row, c, v));
            }
        }
        triplets.push((met_idx, n_rxns - 1, -1.0));
        model.s = StoichMatrix::from_triplets(n_mets, n_rxns, triplets);
    } else if let Some(r) = model.rxns.iter_mut().find(|r| r.id.as_str() == ex_id) {
        r.lb = -10.0;
    }

    // Probe: does the model already grow on this compound?
    let sol = fba(
        model,
        &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
    )?;
    if matches!(sol.status, SolveStatus::Optimal) && sol.objective >= min_growth {
        eprintln!("model already grows on {cpd} (bio1={:.4})", sol.objective);
        return Ok(());
    }

    // Read a minimal medium from the current open exchanges; then gap-fill.
    let medium: Vec<gapsmith_fill::MediumEntry> = model
        .rxns
        .iter()
        .filter(|r| (r.is_exchange || r.id.as_str().starts_with("EX_")) && r.lb < 0.0)
        .map(|r| gapsmith_fill::MediumEntry {
            compound: r
                .id
                .as_str()
                .trim_start_matches("EX_")
                .trim_end_matches("_e0")
                .to_string(),
            name: r.name.clone(),
            max_flux: -r.lb,
        })
        .collect();
    let _ = read_medium; // keep the import alive for downstream doc

    let (mut full, _) = build_full_model(model, seed_rxns, weights)
        .map_err(|e| anyhow::anyhow!(e))?;
    apply_medium_to_full(&mut full, &medium);
    let opts = GapfillOptions::new(min_growth, full.rxn_count());
    let report = gapfill4(model, &full, weights, seed_rxns, &opts)?;
    if matches!(report.status, SolveStatus::Optimal) {
        *model = report.model;
        eprintln!(
            "grew on {cpd}: +{} rxns, bio1={:.4}",
            report.rxns_added.len(),
            report.growth_rate
        );
    } else {
        eprintln!("failed to gap-fill growth on {cpd}");
    }
    Ok(())
}
