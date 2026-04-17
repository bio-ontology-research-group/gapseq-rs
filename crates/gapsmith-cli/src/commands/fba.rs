//! `gapsmith fba` — solve FBA or pFBA on an existing CBOR/JSON model.
//!
//! Prints the solver status, objective value, biomass flux, and the top-N
//! absolute fluxes. Useful for sanity-checking a draft model before /
//! after gap-filling, and for debugging objective-scaling issues.

use clap::Parser;
use gapsmith_core::Model;
use gapsmith_fill::{fba, pfba, FbaOptions, PfbaOptions, SolveStatus};
use gapsmith_io::{read_model_cbor, read_model_json, ModelFormat};
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    /// Model file (`.gmod.cbor`, `.cbor`, or `.json`).
    pub model: PathBuf,

    /// Objective reaction id. When absent, reactions with `obj_coef != 0`
    /// in the model are used (normally `bio1_c0` for gapseq drafts).
    #[arg(short = 'r', long)]
    pub objective: Option<String>,

    /// Run parsimonious FBA (`Σ w_r · |v_r|` with `w_r = 1` for non-obj)
    /// instead of plain FBA.
    #[arg(long)]
    pub pfba: bool,

    /// pFBA coefficient (trade-off between total flux and biomass).
    #[arg(long, default_value_t = 1e-3)]
    pub pfba_coef: f64,

    /// Minimum required growth (biomass flux floor) when `--pfba` is set.
    #[arg(long, default_value_t = 0.0)]
    pub min_growth: f64,

    /// Print this many top-|flux| reactions (default: 20).
    #[arg(long, default_value_t = 20)]
    pub top: usize,

    /// Minimise instead of maximise (rarely useful).
    #[arg(long)]
    pub minimise: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let fmt = ModelFormat::from_path(&args.model);
    let mut model = match fmt {
        ModelFormat::Cbor => read_model_cbor(&args.model)?,
        ModelFormat::Json => read_model_json(&args.model)?,
    };
    model.check_shape()?;

    if let Some(obj_id) = &args.objective {
        set_objective(&mut model, obj_id)?;
    }

    if args.pfba {
        let weights: Vec<f64> = model
            .rxns
            .iter()
            .map(|r| if r.obj_coef.abs() >= f64::EPSILON { 0.0 } else { 1.0 })
            .collect();
        let sol = pfba(
            &model,
            &PfbaOptions {
                weights,
                pfba_coef: args.pfba_coef,
                min_growth: args.min_growth,
                objective: None,
            },
        )?;
        print_summary(&model, sol.status, sol.growth, &sol.fluxes, args.top, "pFBA");
    } else {
        let sol = fba(
            &model,
            &FbaOptions {
                objective: None,
                maximise: !args.minimise,
                max_flux: 1000.0,
            },
        )?;
        print_summary(&model, sol.status, sol.objective, &sol.fluxes, args.top, "FBA");
    }

    Ok(())
}

fn set_objective(model: &mut Model, obj_id: &str) -> anyhow::Result<()> {
    let mut found = false;
    for r in &mut model.rxns {
        if r.id.as_str() == obj_id {
            r.obj_coef = 1.0;
            found = true;
        } else {
            r.obj_coef = 0.0;
        }
    }
    if !found {
        anyhow::bail!("objective reaction `{obj_id}` not found in model");
    }
    Ok(())
}

fn print_summary(
    model: &Model,
    status: SolveStatus,
    objective: f64,
    fluxes: &[f64],
    top: usize,
    label: &str,
) {
    println!("{label} solver status: {status:?}");
    println!("{label} objective      : {objective:.6}");
    let bio_flux: Vec<(String, f64)> = model
        .rxns
        .iter()
        .zip(fluxes)
        .filter(|(r, _)| r.is_biomass || r.obj_coef.abs() >= f64::EPSILON)
        .map(|(r, v)| (r.id.to_string(), *v))
        .collect();
    for (id, v) in &bio_flux {
        println!("  biomass flux    : {id} = {v:.6}");
    }

    if top > 0 {
        let mut ranked: Vec<(usize, f64)> = fluxes
            .iter()
            .enumerate()
            .map(|(i, v)| (i, *v))
            .collect();
        ranked.sort_by(|a, b| b.1.abs().partial_cmp(&a.1.abs()).unwrap_or(std::cmp::Ordering::Equal));
        println!("\ntop {} reactions by |flux|:", top.min(ranked.len()));
        for (i, v) in ranked.iter().take(top) {
            if v.abs() < 1e-9 {
                break;
            }
            println!("  {:16}  {:+.6}", model.rxns[*i].id.as_str(), v);
        }
    }
}
