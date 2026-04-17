//! Exchange / sink / diffusion reaction helpers. Port of
//! `src/add_missing_exRxns.R`.

use gapsmith_core::{CompartmentId, Metabolite, Model, Reaction, StoichMatrix};
use std::collections::HashSet;

/// For every extracellular metabolite not yet covered by an
/// `EX_<cpd>_e0` reaction, add one (`<cpd>[e0] → ∅`) with bounds
/// `[0, ub]`. Mirrors `add_missing_exchanges(mod)`.
pub fn add_missing_exchanges(model: &mut Model, ub: f64) {
    let existing: HashSet<String> = model
        .rxns
        .iter()
        .filter_map(|r| {
            r.id.as_str()
                .strip_prefix("EX_")
                .and_then(|s| s.strip_suffix("_e0"))
                .map(|s| s.to_string())
        })
        .collect();

    // Collect extracellular mets.
    let extracellular: Vec<(usize, Metabolite)> = model
        .mets
        .iter()
        .enumerate()
        .filter(|(_, m)| m.compartment == CompartmentId::EXTRACELLULAR)
        .map(|(i, m)| (i, m.clone()))
        .collect();

    for (met_idx, met) in extracellular {
        // Met id is `<cpd>_e0` in the M7 convention; strip the suffix.
        let base = met.id.as_str().trim_end_matches("_e0").to_string();
        if existing.contains(&base) {
            continue;
        }
        let rxn_id = format!("EX_{base}_e0");
        // Append reaction + column to S.
        add_boundary_reaction(
            model,
            &rxn_id,
            &format!("{} Exchange", met.name),
            met_idx,
            0.0,
            ub,
            Some(7),
        );
    }
}

/// Add a diffusion reaction from a pre-loaded table (`(met_id, rxn_id)`
/// pairs read from `dat/diffusion_mets.tsv`). For each `rxn_id`, look it
/// up in the SEED reactions DB and append it if not already present;
/// gs_origin = 8.
pub fn add_missing_diffusion(
    model: &mut Model,
    diffusion_rxn_ids: &[String],
    seed_reactions: &[gapsmith_db::SeedRxnRow],
) {
    let existing: HashSet<String> = model.rxns.iter().map(|r| r.id.as_str().to_string()).collect();
    for rxn_id in diffusion_rxn_ids {
        let full_id = format!("{rxn_id}_c0");
        if existing.contains(&full_id) {
            continue;
        }
        if let Some(row) = seed_reactions.iter().find(|r| r.id.as_str() == rxn_id) {
            crate::builder::add_seed_reaction(model, row, Some(8));
        } else {
            tracing::debug!(rxn_id, "diffusion reaction not in SEED DB; skipping");
        }
    }
    add_missing_exchanges(model, 1000.0);
}

/// Load `dat/diffusion_mets.tsv` — `met\tdiffrxn\tcomment` — and return
/// the list of `diffrxn` ids.
pub fn load_diffusion_rxns(path: &std::path::Path) -> std::io::Result<Vec<String>> {
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(path)?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    let mut header = true;
    for line in r.lines() {
        let line = line?;
        if header {
            header = false;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 || cols[1].is_empty() {
            continue;
        }
        out.push(cols[1].to_string());
    }
    Ok(out)
}


/// Create a single-metabolite boundary reaction (exchange / demand /
/// sink) and append it to the model, including extending the S matrix.
pub fn add_boundary_reaction(
    model: &mut Model,
    rxn_id: &str,
    name: &str,
    met_idx: usize,
    lb: f64,
    ub: f64,
    gs_origin: Option<i8>,
) {
    let mut r = Reaction::new(rxn_id, name, lb, ub);
    r.is_exchange = rxn_id.starts_with("EX_");
    r.gs_origin = gs_origin;
    model.rxns.push(r);

    let n_mets = model.mets.len();
    let n_rxns = model.rxns.len();
    // Append the new column to the S matrix. Simplest (and fine for
    // draft models up to ~5k reactions): rebuild the CSC from triplets.
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    for (c, _) in model.rxns.iter().enumerate().take(n_rxns - 1) {
        for (row, coef) in model.s.column(c) {
            triplets.push((row, c, coef));
        }
    }
    triplets.push((met_idx, n_rxns - 1, -1.0));
    model.s = StoichMatrix::from_triplets(n_mets, n_rxns, triplets);
}
