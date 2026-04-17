//! 4-phase gap-fill driver. Port of `src/gf.suite.R`.
//!
//! Five phases, all sharing the same underlying [`crate::gapfill4`] call:
//!
//! - **Step 1** — user medium + biomass target. Bulk of the first-pass fills.
//! - **Step 2** — per-biomass-component gap-fill on minimal medium + available
//!   carbon sources (gf.suite.R:285-372). Catches cases where the Step 1
//!   pFBA didn't explore a component because multiple pathways could
//!   satisfy biomass simultaneously.
//! - **Step 2b** — aerobic/anaerobic variant of Step 2. If user medium
//!   excludes O2, do the per-component fill under anaerobic conditions.
//! - **Step 3** — energy-source screen with ESP1–5 pseudo-reactions
//!   (gf.suite.R:480-581). Not yet implemented — skipped when requested.
//! - **Step 4** — fermentation product potential (gf.suite.R:585-683).
//!   Not yet implemented — skipped when requested.

use crate::error::FillError;
use crate::fba::{fba, FbaOptions};
use crate::futile::{detect_futile_cycles, FutileOptions};
use crate::gapfill::{drop_reactions, gapfill4, GapfillOptions};
use crate::medium::{apply_medium, read_medium, MediumEntry};
use crate::pfba::PfbaHeuristicOptions;
use crate::pool::{apply_medium_to_full, build_full_model, RxnWeights};
use crate::SolveStatus;
use gapsmith_core::{CompartmentId, Metabolite, Model, Reaction, StoichMatrix};
use gapsmith_db::SeedRxnRow;
use std::path::Path;

/// Per-phase options shared by Steps 1, 2, 2b.
#[derive(Debug, Clone)]
pub struct SuiteOptions {
    pub target_cpd: String,
    pub min_growth: f64,
    pub bcore: f64,
    pub high_evi: f64,
    pub dummy_weight: f64,
    /// Skip Steps 3 and 4. Matches gapseq's `--quick.gf`.
    pub quick: bool,
    /// Skip Steps 2/2b too (fastest path — Step 1 only).
    pub step1_only: bool,
    /// Drop thermodynamically-infeasible candidate reactions before
    /// gap-filling. Recent upstream feature (`cccbb6f0`). Default `true`.
    pub prune_futile: bool,
}

impl Default for SuiteOptions {
    fn default() -> Self {
        Self {
            target_cpd: "cpd11416".into(),
            min_growth: 0.01,
            bcore: crate::DEFAULT_BCORE,
            high_evi: crate::DEFAULT_HIGH_EVI,
            dummy_weight: crate::DEFAULT_DUMMY_WEIGHT,
            quick: true,
            step1_only: false,
            // Off by default: on large candidate pools (~8k+) the pairwise
            // LP probe takes tens of minutes. Opt in when you know the
            // pool is narrow or when running on a small genome.
            prune_futile: false,
        }
    }
}

/// Structured report summarising every phase of the fill.
#[derive(Debug, Clone, Default)]
pub struct SuiteReport {
    /// Reaction ids added by each step (in order).
    pub step1_added: Vec<String>,
    pub step2_added: Vec<String>,
    pub step2b_added: Vec<String>,
    pub step3_added: Vec<String>,
    pub step4_added: Vec<String>,
    pub step1_growth: f64,
    pub step2_growth: f64,
    pub step2b_growth: f64,
    pub step3_growth: f64,
    pub step4_growth: f64,
    pub final_growth: f64,
    /// Carbon sources accepted in Step 3.
    pub carbon_sources: Vec<(String, String, bool)>,
    /// Fermentation products tested in Step 4.
    pub ferm_products: Vec<(String, String, bool)>,
}

impl SuiteReport {
    pub fn total_added(&self) -> usize {
        self.step1_added.len()
            + self.step2_added.len()
            + self.step2b_added.len()
            + self.step3_added.len()
            + self.step4_added.len()
    }
}

/// Run Step 1 (user-medium biomass fill) followed by Steps 2 / 2b
/// (per-biomass-component fill on minimal medium + carbon sources). Steps 3
/// and 4 are skipped in this milestone.
///
/// `data_dir` must point at gapseq's `dat/`. The driver reads
/// `dat/media/MM_glu.csv` and `dat/subex.tbl` for Step 2.
pub fn run_suite(
    draft: &Model,
    user_medium: &[MediumEntry],
    weights: &RxnWeights,
    seed_rxns: &[SeedRxnRow],
    data_dir: &Path,
    opts: &SuiteOptions,
) -> Result<(Model, SuiteReport), FillError> {
    let mut report = SuiteReport::default();

    // -- Step 1: user medium + biomass target.
    let mut current = draft.clone();
    apply_medium(&mut current, user_medium, 1.0, 1000.0);
    add_target_sink_obj(&mut current, &opts.target_cpd);

    let (mut full, cands) = build_full_model(&current, seed_rxns, weights)
        .map_err(FillError::Solver)?;
    apply_medium_to_full(&mut full, user_medium);

    if opts.prune_futile {
        let candidate_rxn_ids: Vec<String> =
            cands.iter().map(|s| format!("{s}_c0")).collect();
        match detect_futile_cycles(&full, &candidate_rxn_ids, &FutileOptions::default()) {
            Ok(bad) => {
                tracing::info!(n = bad.len(), "futile-cycle prune: dropping reactions");
                drop_reactions(&mut full, &bad);
            }
            Err(e) => {
                tracing::warn!(err = %e, "futile-cycle prune failed; continuing without");
            }
        }
    }

    let step_opts = build_gapfill_options(opts.min_growth, full.rxn_count());
    let s1 = gapfill4(&current, &full, weights, seed_rxns, &step_opts)?;
    if s1.rxns_added.is_empty() {
        tracing::info!(growth = s1.growth_rate, "step 1: nothing added");
    }
    report.step1_added = s1.rxns_added.clone();
    report.step1_growth = s1.growth_rate;
    current = s1.model;

    if opts.step1_only {
        report.final_growth = current_growth(&current)?;
        return Ok((current, report));
    }

    // -- Step 2: per-biomass-component fill on MM_glu + carbon sources.
    let mm_glu_path = data_dir.join("media/MM_glu.csv");
    let mut mm_glu = match read_medium(&mm_glu_path) {
        Ok(m) => m,
        Err(e) => {
            tracing::warn!(err = %e, "Step 2 skipped: couldn't read MM_glu.csv");
            report.final_growth = current_growth(&current)?;
            return Ok((current, report));
        }
    };
    augment_with_carbon_sources(&mut mm_glu, &current, data_dir);

    // Aerobic Step 2.
    let (aerobic_model, added) = step2(
        &current,
        &mm_glu,
        weights,
        seed_rxns,
        opts,
    )?;
    report.step2_added = added;
    report.step2_growth = current_growth(&aerobic_model)?;
    current = aerobic_model;
    // Restore user-medium constraints after Step 2.
    apply_medium(&mut current, user_medium, 1.0, 1000.0);

    // -- Step 2b: aerobic/anaerobic variant depending on user medium.
    let mut mm_2b = match read_medium(&mm_glu_path) {
        Ok(m) => m,
        Err(_) => {
            report.final_growth = current_growth(&current)?;
            return Ok((current, report));
        }
    };
    augment_with_carbon_sources(&mut mm_2b, &current, data_dir);
    let anaerobic = user_is_anaerobic(user_medium);
    if anaerobic {
        mm_2b.retain(|e| e.compound != "cpd00007"); // drop O2
        tracing::info!("step 2b: anaerobic (O2 removed)");
    } else {
        tracing::info!("step 2b: aerobic (re-run on MM_glu)");
    }
    let (model_2b, added_2b) = step2(
        &current,
        &mm_2b,
        weights,
        seed_rxns,
        opts,
    )?;
    report.step2b_added = added_2b;
    report.step2b_growth = current_growth(&model_2b)?;
    current = model_2b;
    apply_medium(&mut current, user_medium, 1.0, 1000.0);

    if opts.quick {
        report.final_growth = current_growth(&current)?;
        return Ok((current, report));
    }

    // -- Step 3: energy-source screen. Iterate each exchange compound,
    // test as carbon source on minimal medium + ESP1..5 redox obj.
    let mut mm_step3 = match read_medium(&mm_glu_path) {
        Ok(m) => m,
        Err(_) => {
            report.final_growth = current_growth(&current)?;
            return Ok((current, report));
        }
    };
    // Drop glucose — we'll substitute each test compound.
    mm_step3.retain(|e| e.compound != "cpd00027");

    let (model_s3, added_s3, cs_table) =
        step3(&current, &mm_step3, weights, seed_rxns, opts)?;
    report.step3_added = added_s3;
    report.carbon_sources = cs_table;
    report.step3_growth = current_growth(&model_s3)?;
    current = model_s3;
    apply_medium(&mut current, user_medium, 1.0, 1000.0);

    // -- Step 4: fermentation products on user medium.
    let (model_s4, added_s4, ferm_table) =
        step4(&current, user_medium, weights, seed_rxns, opts)?;
    report.step4_added = added_s4;
    report.ferm_products = ferm_table;
    report.step4_growth = current_growth(&model_s4)?;
    current = model_s4;
    apply_medium(&mut current, user_medium, 1.0, 1000.0);

    report.final_growth = current_growth(&current)?;
    Ok((current, report))
}

/// Does the user medium disable oxygen? (Either cpd00007 not in medium, or
/// explicitly at 0.)
fn user_is_anaerobic(user_medium: &[MediumEntry]) -> bool {
    !user_medium
        .iter()
        .any(|e| e.compound == "cpd00007" && e.max_flux > 0.0)
}

/// Step 2: iterate every biomass substrate, test for production, gap-fill
/// if needed. Returns the resulting model plus the set of added reaction ids.
fn step2(
    model: &Model,
    medium: &[MediumEntry],
    weights: &RxnWeights,
    seed_rxns: &[SeedRxnRow],
    opts: &SuiteOptions,
) -> Result<(Model, Vec<String>), FillError> {
    let mut current = model.clone();
    apply_medium(&mut current, medium, 1.0, 1000.0);

    // Collect biomass substrates (rxn id == "bio1" column with negative
    // coefficients). Snapshot the set up-front so mutations to the model
    // inside the loop don't re-expose substrates we already handled.
    let bio1_col = match current.rxns.iter().position(|r| r.id.as_str() == "bio1") {
        Some(i) => i,
        None => {
            tracing::warn!("step 2 skipped: bio1 reaction not found");
            return Ok((current, Vec::new()));
        }
    };
    let bio_substrates: Vec<String> = current
        .s
        .column(bio1_col)
        .into_iter()
        .filter(|(_, v)| *v < 0.0)
        .map(|(r, _)| current.mets[r].id.as_str().to_string())
        .collect();
    tracing::info!(n = bio_substrates.len(), "step 2: iterating biomass substrates");

    let (mut full, _cands) = build_full_model(&current, seed_rxns, weights)
        .map_err(FillError::Solver)?;
    apply_medium_to_full(&mut full, medium);

    let mut added_accum: Vec<String> = Vec::new();
    for (i, cpd) in bio_substrates.iter().enumerate() {
        // Turn OFF every existing objective (bio1, etc.), add a sink on the
        // target metabolite, set obj=1 on the sink.
        let target = strip_c0(cpd);
        let (sink_added, sink_id) = add_sink_if_missing(&mut current, target);

        for r in &mut current.rxns { r.obj_coef = 0.0; }
        if let Some(r) = current.rxns.iter_mut().find(|r| r.id.as_str() == sink_id) {
            r.obj_coef = 1.0;
        }

        // Also add the sink to the full model so obj extraction lines up.
        let (full_sink_added, _) = add_sink_if_missing(&mut full, target);
        for r in &mut full.rxns { r.obj_coef = 0.0; }
        if let Some(r) = full.rxns.iter_mut().find(|r| r.id.as_str() == sink_id) {
            r.obj_coef = 1.0;
        }

        // Probe feasibility.
        let probe = fba(
            &current,
            &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
        )?;
        if matches!(probe.status, SolveStatus::Optimal) && probe.objective >= 1e-6 {
            if sink_added {
                remove_reaction(&mut current, &sink_id);
            }
            if full_sink_added {
                remove_reaction(&mut full, &sink_id);
            }
            continue;
        }

        // Run gapfill4 targeting this sink.
        let step_opts = build_gapfill_options(opts.min_growth, full.rxn_count());
        let rep = match gapfill4(&current, &full, weights, seed_rxns, &step_opts) {
            Ok(r) => r,
            Err(e) => {
                tracing::debug!(idx = i, cpd = %cpd, err = %e, "step 2: gapfill4 failed");
                if sink_added { remove_reaction(&mut current, &sink_id); }
                if full_sink_added { remove_reaction(&mut full, &sink_id); }
                continue;
            }
        };
        if !rep.rxns_added.is_empty() && matches!(rep.status, SolveStatus::Optimal) {
            tracing::debug!(
                idx = i,
                cpd = %cpd,
                added = rep.rxns_added.len(),
                growth = rep.growth_rate,
                "step 2: filled"
            );
            current = rep.model;
            // Rebuild the full model to reflect the newly added rxns.
            let (new_full, _) = build_full_model(&current, seed_rxns, weights)
                .map_err(FillError::Solver)?;
            full = new_full;
            apply_medium_to_full(&mut full, medium);
            for id in rep.rxns_added {
                if !added_accum.contains(&id) {
                    added_accum.push(id);
                }
            }
        }

        if sink_added {
            remove_reaction(&mut current, &sink_id);
        }
        if full_sink_added {
            remove_reaction(&mut full, &sink_id);
        }
    }

    // Restore bio1 as the objective.
    for r in &mut current.rxns { r.obj_coef = 0.0; }
    let target_sink = format!("EX_{}_c0", opts.target_cpd);
    if let Some(r) = current.rxns.iter_mut().find(|r| r.id.as_str() == target_sink) {
        r.obj_coef = 1.0;
    } else if let Some(r) = current.rxns.iter_mut().find(|r| r.id.as_str() == "bio1") {
        r.obj_coef = 1.0;
    }

    Ok((current, added_accum))
}

/// Compounds (and exchange rxns) ignored by Steps 3 & 4 — metals, gases,
/// cofactors that are always available in MM_glu. Matches gapseq's
/// `ignore` vector at `gf.suite.R:499`.
const STEP34_IGNORE: &[&str] = &[
    "EX_cpd17041_e0",
    "EX_cpd17042_e0",
    "EX_cpd17043_e0",
    "EX_cpd11416_e0",
    "rxn13782_c0",
    "rxn13783_c0",
    "EX_cpd00001_e0", // H2O
    "EX_cpd00007_e0", // O2
    "EX_cpd00009_e0", // Phosphate
    "EX_cpd00011_e0", // CO2
    "EX_cpd00012_e0", // PPi
    "EX_cpd00030_e0", // Mn
    "EX_cpd00034_e0", // Zn
    "EX_cpd00058_e0", // Cu
    "EX_cpd00063_e0", // Ca
    "EX_cpd00067_e0", // H+
    "EX_cpd00075_e0", // NO2
    "EX_cpd00099_e0", // Cl-
    "EX_cpd00149_e0", // Cobalt
    "EX_cpd00205_e0", // K+
    "EX_cpd00254_e0", // Mg
    "EX_cpd10515_e0", // Fe2+
    "EX_cpd00971_e0", // Sodium
    "EX_cpd01012_e0",
    "EX_cpd10516_e0", // Fe3+
    "EX_cpd11574_e0", // Mo
];

/// Step 3: energy-source screen. Iterate each exchange compound, test as
/// carbon source on `medium` + 5 ESP pseudo-reactions as the objective.
/// `medium` is typically MM_glu with glucose already removed.
fn step3(
    model: &Model,
    medium: &[MediumEntry],
    weights: &RxnWeights,
    seed_rxns: &[SeedRxnRow],
    opts: &SuiteOptions,
) -> Result<(Model, Vec<String>, Vec<(String, String, bool)>), FillError> {
    let mut current = model.clone();
    // Install ESP reactions on the draft.
    for esp in esp_reactions() {
        add_pseudo_reaction(&mut current, &esp);
    }

    // Build the full model (with ESP in it) once; re-use across carbon-source iterations.
    let (mut full, _cands) = build_full_model(&current, seed_rxns, weights)
        .map_err(FillError::Solver)?;
    for esp in esp_reactions() {
        add_pseudo_reaction(&mut full, &esp);
    }

    // Snapshot: exchange compounds present in the draft.
    let ignore: std::collections::HashSet<String> =
        STEP34_IGNORE.iter().map(|s| (*s).to_string()).collect();
    let candidates: Vec<(String, String)> = current
        .rxns
        .iter()
        .filter(|r| r.is_exchange || r.id.as_str().starts_with("EX_"))
        .filter(|r| !ignore.contains(r.id.as_str()))
        .filter_map(|r| {
            let cpd = r.id.as_str().trim_start_matches("EX_").trim_end_matches("_e0");
            if cpd.is_empty() { return None; }
            Some((cpd.to_string(), r.name.clone()))
        })
        .collect();
    tracing::info!(n = candidates.len(), "step 3: iterating carbon sources");

    let mut added_accum: Vec<String> = Vec::new();
    let mut cs_table: Vec<(String, String, bool)> = Vec::new();

    // Objective = sum of ESP1..5. We zero all other obj coefs.
    set_esp_objective(&mut current);
    set_esp_objective(&mut full);

    for (cpd, name) in &candidates {
        let mut med = medium.to_vec();
        med.push(MediumEntry { compound: cpd.clone(), name: name.clone(), max_flux: 100.0 });

        apply_medium(&mut current, &med, 1.0, 1000.0);
        apply_medium_to_full(&mut full, &med);

        let probe = fba(
            &current,
            &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
        )?;
        if matches!(probe.status, SolveStatus::Optimal) && probe.objective >= 1e-7 {
            cs_table.push((cpd.clone(), name.clone(), true));
            continue;
        }

        let step_opts = build_gapfill_options(opts.min_growth, full.rxn_count());
        match gapfill4(&current, &full, weights, seed_rxns, &step_opts) {
            Ok(rep) if matches!(rep.status, SolveStatus::Optimal) && !rep.rxns_added.is_empty() => {
                current = rep.model;
                let (new_full, _) = build_full_model(&current, seed_rxns, weights)
                    .map_err(FillError::Solver)?;
                full = new_full;
                for esp in esp_reactions() {
                    add_pseudo_reaction(&mut full, &esp);
                }
                set_esp_objective(&mut current);
                set_esp_objective(&mut full);
                for id in rep.rxns_added {
                    if !added_accum.contains(&id) {
                        added_accum.push(id);
                    }
                }
                cs_table.push((cpd.clone(), name.clone(), true));
            }
            _ => cs_table.push((cpd.clone(), name.clone(), false)),
        }
    }

    // Restore user-facing state: remove ESP pseudo-reactions, put back
    // the biomass-sink objective.
    for esp_id in ["ESP1", "ESP2", "ESP3", "ESP4", "ESP5"] {
        remove_reaction(&mut current, esp_id);
    }
    for r in &mut current.rxns { r.obj_coef = 0.0; }
    let target_sink = format!("EX_{}_c0", opts.target_cpd);
    if let Some(r) = current.rxns.iter_mut().find(|r| r.id.as_str() == target_sink) {
        r.obj_coef = 1.0;
    } else if let Some(r) = current.rxns.iter_mut().find(|r| r.id.as_str() == "bio1") {
        r.obj_coef = 1.0;
    }

    Ok((current, added_accum, cs_table))
}

/// Step 4: fermentation-product screen. For each exchange compound on the
/// user medium, set the exchange itself as the objective (maximise
/// secretion) and gap-fill if infeasible.
fn step4(
    model: &Model,
    medium: &[MediumEntry],
    weights: &RxnWeights,
    seed_rxns: &[SeedRxnRow],
    opts: &SuiteOptions,
) -> Result<(Model, Vec<String>, Vec<(String, String, bool)>), FillError> {
    let mut current = model.clone();
    apply_medium(&mut current, medium, 1.0, 1000.0);

    let (mut full, _cands) = build_full_model(&current, seed_rxns, weights)
        .map_err(FillError::Solver)?;
    apply_medium_to_full(&mut full, medium);

    let ignore: std::collections::HashSet<String> =
        STEP34_IGNORE.iter().map(|s| (*s).to_string()).collect();
    let candidates: Vec<(String, String, String)> = current
        .rxns
        .iter()
        .filter(|r| r.is_exchange || r.id.as_str().starts_with("EX_"))
        .filter(|r| !ignore.contains(r.id.as_str()))
        .map(|r| {
            let cpd = r.id.as_str().trim_start_matches("EX_").trim_end_matches("_e0").to_string();
            (r.id.as_str().to_string(), cpd, r.name.clone())
        })
        .filter(|(_, cpd, _)| !cpd.is_empty())
        .collect();
    tracing::info!(n = candidates.len(), "step 4: iterating fermentation products");

    let mut added_accum: Vec<String> = Vec::new();
    let mut ferm_table: Vec<(String, String, bool)> = Vec::new();

    for (rxn_id, cpd, name) in &candidates {
        // Zero everything and set obj=1 on this exchange (secretion).
        for r in &mut current.rxns { r.obj_coef = 0.0; }
        for r in &mut full.rxns { r.obj_coef = 0.0; }
        if let Some(r) = current.rxns.iter_mut().find(|r| r.id.as_str() == rxn_id) {
            r.obj_coef = 1.0;
        }
        if let Some(r) = full.rxns.iter_mut().find(|r| r.id.as_str() == rxn_id) {
            r.obj_coef = 1.0;
        }

        let probe = fba(
            &current,
            &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
        )?;
        if matches!(probe.status, SolveStatus::Optimal) && probe.objective >= 1e-7 {
            ferm_table.push((cpd.clone(), name.clone(), true));
            continue;
        }

        let step_opts = build_gapfill_options(opts.min_growth, full.rxn_count());
        match gapfill4(&current, &full, weights, seed_rxns, &step_opts) {
            Ok(rep) if matches!(rep.status, SolveStatus::Optimal) && !rep.rxns_added.is_empty() => {
                current = rep.model;
                let (new_full, _) = build_full_model(&current, seed_rxns, weights)
                    .map_err(FillError::Solver)?;
                full = new_full;
                apply_medium_to_full(&mut full, medium);
                for id in rep.rxns_added {
                    if !added_accum.contains(&id) {
                        added_accum.push(id);
                    }
                }
                ferm_table.push((cpd.clone(), name.clone(), true));
            }
            _ => ferm_table.push((cpd.clone(), name.clone(), false)),
        }
    }

    // Restore biomass-sink objective.
    for r in &mut current.rxns { r.obj_coef = 0.0; }
    let target_sink = format!("EX_{}_c0", opts.target_cpd);
    if let Some(r) = current.rxns.iter_mut().find(|r| r.id.as_str() == target_sink) {
        r.obj_coef = 1.0;
    } else if let Some(r) = current.rxns.iter_mut().find(|r| r.id.as_str() == "bio1") {
        r.obj_coef = 1.0;
    }

    Ok((current, added_accum, ferm_table))
}

/// Define the 5 ESP pseudo-reactions (gf.suite.R:506-516).
/// Each ESP drains a reduced/oxidized redox couple, so setting them as
/// objective lets the solver use whatever electron-acceptor chemistry the
/// carbon source enables.
struct EspSpec {
    id: &'static str,
    substrates: &'static [(&'static str, f64)],
    products: &'static [(&'static str, f64)],
}

fn esp_reactions() -> [EspSpec; 5] {
    [
        EspSpec { id: "ESP1", substrates: &[("cpd15499_c0", 1.0)],
                  products: &[("cpd00067_c0", 2.0), ("cpd15500_c0", 1.0)] }, // menaquinone
        EspSpec { id: "ESP2", substrates: &[("cpd15561_c0", 1.0)],
                  products: &[("cpd00067_c0", 2.0), ("cpd15560_c0", 1.0)] }, // ubiquinone
        EspSpec { id: "ESP3", substrates: &[("cpd00004_c0", 1.0)],
                  products: &[("cpd00067_c0", 1.0), ("cpd00003_c0", 1.0)] }, // NADH -> NAD
        EspSpec { id: "ESP4", substrates: &[("cpd11620_c0", 1.0)],
                  products: &[("cpd11621_c0", 1.0)] }, // ferredoxin
        EspSpec { id: "ESP5", substrates: &[("cpd27796_c0", 1.0)],
                  products: &[("cpd00067_c0", 2.0), ("cpd27797_c0", 1.0)] }, // plastoquinone
    ]
}

fn add_pseudo_reaction(model: &mut Model, spec: &EspSpec) {
    // Skip if already present.
    if model.rxns.iter().any(|r| r.id.as_str() == spec.id) {
        return;
    }
    use std::collections::HashMap;
    let mut col: HashMap<usize, f64> = HashMap::new();
    for (cpd, coef) in spec.substrates {
        if let Some(idx) = ensure_met(model, cpd) {
            *col.entry(idx).or_default() -= *coef;
        }
    }
    for (cpd, coef) in spec.products {
        if let Some(idx) = ensure_met(model, cpd) {
            *col.entry(idx).or_default() += *coef;
        }
    }
    let mut r = Reaction::new(spec.id, spec.id, 0.0, 1000.0);
    r.gs_origin = Some(10); // synthetic ESP
    model.rxns.push(r);
    let n_mets = model.mets.len();
    let n_rxns = model.rxns.len();
    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(model.s.nnz() + col.len());
    let old_cols = model.s.cols();
    for c in 0..old_cols.min(n_rxns - 1) {
        for (row, v) in model.s.column(c) {
            triplets.push((row, c, v));
        }
    }
    for (row, v) in col {
        if v != 0.0 {
            triplets.push((row, n_rxns - 1, v));
        }
    }
    model.s = StoichMatrix::from_triplets(n_mets, n_rxns, triplets);
}

/// Look up (or create) a metabolite by id. Returns its index.
fn ensure_met(model: &mut Model, met_id: &str) -> Option<usize> {
    if let Some(i) = model.mets.iter().position(|m| m.id.as_str() == met_id) {
        return Some(i);
    }
    let comp = if met_id.ends_with("_c0") {
        CompartmentId::CYTOSOL
    } else if met_id.ends_with("_e0") {
        CompartmentId::EXTRACELLULAR
    } else {
        CompartmentId::PERIPLASM
    };
    let name = met_id.trim_end_matches("_c0").trim_end_matches("_e0").trim_end_matches("_p0");
    model.mets.push(Metabolite::new(met_id, name, comp));
    Some(model.mets.len() - 1)
}

fn set_esp_objective(model: &mut Model) {
    for r in &mut model.rxns {
        r.obj_coef = 0.0;
    }
    for esp_id in ["ESP1", "ESP2", "ESP3", "ESP4", "ESP5"] {
        if let Some(r) = model.rxns.iter_mut().find(|r| r.id.as_str() == esp_id) {
            r.obj_coef = 1.0;
        }
    }
}

/// Read `subex.tbl` and for every carbon-source row whose SEED id points at an
/// EX reaction already in the model (meaning the user's draft has this
/// compound available), add it to the medium at 100 u/g DW.
fn augment_with_carbon_sources(medium: &mut Vec<MediumEntry>, model: &Model, data_dir: &Path) {
    let subex_path = data_dir.join("subex.tbl");
    let text = match std::fs::read_to_string(&subex_path) {
        Ok(t) => t,
        Err(_) => return,
    };
    let mut lines = text.lines();
    let header = match lines.next() {
        Some(h) => h,
        None => return,
    };
    let cols: Vec<&str> = header.split('\t').collect();
    let ix = |name: &str| cols.iter().position(|c| *c == name);
    let (Some(i_name), Some(i_seed), Some(i_group)) = (ix("name"), ix("seed"), ix("group"))
    else {
        return;
    };

    let rxn_ids: std::collections::HashSet<String> =
        model.rxns.iter().map(|r| r.id.as_str().to_string()).collect();
    let already_in_medium: std::collections::HashSet<String> =
        medium.iter().map(|e| e.compound.clone()).collect();

    // gapseq's Step 2 only considers Carbohydrates / Polymers /
    // Carboxylic acids / Amino acids as fill carbon sources.
    let wanted_groups = [
        "Carbohydrates",
        "Polymers",
        "Carboxylic acids",
        "Amino acids",
    ];

    for line in lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() <= i_name.max(i_seed).max(i_group) {
            continue;
        }
        let seed = parts[i_seed];
        if seed.is_empty() || !rxn_ids.contains(seed) {
            continue;
        }
        if !wanted_groups.contains(&parts[i_group]) {
            continue;
        }
        // seed is `EX_cpdXXXXX_e0`; extract cpd.
        let cpd = seed
            .trim_start_matches("EX_")
            .trim_end_matches("_e0");
        if already_in_medium.contains(cpd) {
            continue;
        }
        medium.push(MediumEntry {
            compound: cpd.to_string(),
            name: parts[i_name].to_string(),
            max_flux: 100.0,
        });
    }
}

fn strip_c0(met_id: &str) -> &str {
    met_id
        .strip_suffix("_c0")
        .or_else(|| met_id.strip_suffix("_e0"))
        .or_else(|| met_id.strip_suffix("_p0"))
        .unwrap_or(met_id)
}

/// Add a `EX_<cpd>_c0` sink if not already present. Returns `(added, id)`.
fn add_sink_if_missing(model: &mut Model, cpd: &str) -> (bool, String) {
    let rxn_id = format!("EX_{cpd}_c0");
    if model.rxns.iter().any(|r| r.id.as_str() == rxn_id) {
        return (false, rxn_id);
    }
    let met_id = format!("{cpd}_c0");
    let met_idx = match model.mets.iter().position(|m| m.id.as_str() == met_id) {
        Some(i) => i,
        None => {
            let met = Metabolite::new(met_id.as_str(), cpd, CompartmentId::CYTOSOL);
            model.mets.push(met);
            model.mets.len() - 1
        }
    };

    let mut r = Reaction::new(rxn_id.as_str(), format!("Sink: {cpd}"), 0.0, 1000.0);
    r.is_exchange = true;
    r.gs_origin = Some(7);
    model.rxns.push(r);

    append_single_column(model, met_idx, -1.0);
    (true, rxn_id)
}

/// Remove a single reaction by id (brute-force rebuild of S).
fn remove_reaction(model: &mut Model, rxn_id: &str) {
    let mut removed = std::collections::HashSet::new();
    removed.insert(rxn_id.to_string());
    crate::gapfill::drop_reactions(model, &removed);
}

fn append_single_column(model: &mut Model, met_idx: usize, coef: f64) {
    let n_mets = model.mets.len();
    let n_rxns = model.rxns.len();
    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(model.s.nnz() + 1);
    let old_cols = model.s.cols();
    for c in 0..old_cols.min(n_rxns - 1) {
        for (row, v) in model.s.column(c) {
            triplets.push((row, c, v));
        }
    }
    triplets.push((met_idx, n_rxns - 1, coef));
    model.s = StoichMatrix::from_triplets(n_mets, n_rxns, triplets);
}

fn add_target_sink_obj(model: &mut Model, cpd: &str) {
    let (_, id) = add_sink_if_missing(model, cpd);
    for r in &mut model.rxns {
        r.obj_coef = 0.0;
    }
    if let Some(r) = model.rxns.iter_mut().find(|r| r.id.as_str() == id) {
        r.obj_coef = 1.0;
    }
}

fn build_gapfill_options(min_growth: f64, n_rxns_full: usize) -> GapfillOptions {
    let mut o = GapfillOptions::new(min_growth, n_rxns_full);
    o.pfba_heuristic = PfbaHeuristicOptions::new(vec![1.0; n_rxns_full], min_growth);
    o
}

/// Plain FBA growth, robust to solver hiccups (returns 0 on failure).
fn current_growth(model: &Model) -> Result<f64, FillError> {
    let sol = fba(
        model,
        &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
    )?;
    Ok(if matches!(sol.status, SolveStatus::Optimal) {
        sol.objective
    } else {
        0.0
    })
}
