//! Integration smoke-test for the gap-fill pipeline against the real
//! `dat/` and a pre-built ecoli draft. Skips when either is missing.
//!
//! Exercises the full stack: `build_full_model` → medium application →
//! `gapfill4` → KO loop. Confirms the ecoli draft grows on `MM_glu.csv`
//! after gap-filling.

use gapsmith_db::load_seed_reactions;
use gapsmith_fill::{
    apply_medium, apply_medium_to_full, build_full_model, fba, gapfill4, read_medium,
    read_weights_from_reactions_tbl, FbaOptions, GapfillOptions, PfbaHeuristicOptions,
    RxnWeights, SolveStatus, DEFAULT_BCORE, DEFAULT_DUMMY_WEIGHT, DEFAULT_HIGH_EVI,
};
use gapsmith_io::read_model_cbor;
use std::path::PathBuf;

fn gapseq_data() -> PathBuf {
    PathBuf::from(
        std::env::var("GAPSEQ_ROOT")
            .unwrap_or_else(|_| "/opt/gapseq".into()),
    )
    .join("dat")
}

fn draft_path() -> PathBuf {
    PathBuf::from(
        std::env::var("GAPSMITH_TEST_DIR")
            .unwrap_or_else(|_| "/tmp/gapsmith-test".into()),
    )
    .join("ecoli-draft.gmod.cbor")
}

fn reactions_tbl_path() -> PathBuf {
    PathBuf::from(
        std::env::var("GAPSMITH_TEST_DIR")
            .unwrap_or_else(|_| "/tmp/gapsmith-test".into()),
    )
    .join("ecoli-all-Reactions.tbl")
}

#[test]
fn ecoli_end_to_end_fill() {
    if !draft_path().exists() || !reactions_tbl_path().exists() || !gapseq_data().is_dir() {
        eprintln!("SKIP: ecoli draft or dat dir not available");
        return;
    }

    let mut draft = read_model_cbor(&draft_path()).unwrap();
    let medium =
        read_medium(&gapseq_data().join("media/MM_glu.csv")).expect("read MM_glu");
    let weights = read_weights_from_reactions_tbl(
        &reactions_tbl_path(),
        DEFAULT_BCORE,
        DEFAULT_HIGH_EVI,
        DEFAULT_DUMMY_WEIGHT,
    )
    .unwrap();
    let seed_rxns =
        load_seed_reactions(gapseq_data().join("seed_reactions_corrected.tsv")).unwrap();

    apply_medium(&mut draft, &medium, 1.0, 1000.0);
    add_target_sink(&mut draft, "cpd11416");

    let (mut full, _) = build_full_model(&draft, &seed_rxns, &weights).unwrap();
    apply_medium_to_full(&mut full, &medium);

    // Sanity: the full model should be feasible under the biomass
    // objective. Obj > 0 after the target-met sink is attached.
    let probe = fba(
        &full,
        &FbaOptions { objective: None, maximise: true, max_flux: 1000.0 },
    )
    .unwrap();
    assert_eq!(probe.status, SolveStatus::Optimal);
    assert!(probe.objective > 0.01, "full+MM_glu growth={}", probe.objective);

    // Gap-fill.
    let mut opts = GapfillOptions::new(0.01, full.rxn_count());
    opts.pfba_heuristic = PfbaHeuristicOptions::new(vec![1.0; full.rxn_count()], 0.01);
    let report = gapfill4(&draft, &full, &weights, &seed_rxns, &opts).unwrap();

    assert_eq!(report.status, SolveStatus::Optimal);
    assert!(
        report.growth_rate >= 0.01,
        "post-fill growth={} < min_growth",
        report.growth_rate
    );
    assert!(!report.rxns_added.is_empty(), "no reactions added");
    println!(
        "Filled ecoli: +{} rxns ({} core), growth={:.4}",
        report.rxns_added.len(),
        report.rxns_added_core.len(),
        report.growth_rate,
    );
}

/// Mirror of the CLI's `add_met_sink_objective` — attaches obj=1 to a
/// freshly-created `EX_<cpd>_c0` sink reaction.
fn add_target_sink(model: &mut gapsmith_core::Model, cpd: &str) {
    use gapsmith_core::{CompartmentId, Metabolite, Reaction, StoichMatrix};
    let met_id = format!("{cpd}_c0");
    let met_idx = match model.mets.iter().position(|m| m.id.as_str() == met_id) {
        Some(i) => i,
        None => {
            let met = Metabolite::new(met_id.as_str(), cpd, CompartmentId::CYTOSOL);
            model.mets.push(met);
            model.mets.len() - 1
        }
    };
    let rxn_id = format!("EX_{cpd}_c0");
    if !model.rxns.iter().any(|r| r.id.as_str() == rxn_id) {
        let mut r = Reaction::new(rxn_id.as_str(), format!("Sink: {cpd}"), 0.0, 1000.0);
        r.obj_coef = 1.0;
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
    }
    for r in &mut model.rxns {
        if r.id.as_str() != rxn_id {
            r.obj_coef = 0.0;
        }
    }
}
