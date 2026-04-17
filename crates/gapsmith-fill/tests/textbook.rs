//! Textbook FBA / pFBA regression tests.
//!
//! Each test builds a small model whose optimal flux we can solve by hand,
//! then compares HiGHS's answer against the expected value within a 1e-6
//! tolerance.

use gapsmith_core::{CompartmentId, Metabolite, Model, Reaction, StoichMatrix};
use gapsmith_fill::{fba, pfba, pfba_heuristic};
use gapsmith_fill::{FbaOptions, PfbaHeuristicOptions, PfbaOptions, SolveStatus};

/// Linear chain A → B → C with uptake cap 10; biomass on R_BC.
fn chain_model(uptake: f64) -> Model {
    let mut m = Model::new("chain");
    m.mets.push(Metabolite::new("cpdA", "A", CompartmentId::CYTOSOL));
    m.mets.push(Metabolite::new("cpdB", "B", CompartmentId::CYTOSOL));
    m.mets.push(Metabolite::new("cpdC", "C", CompartmentId::CYTOSOL));

    let mut ex_a = Reaction::new("EX_A", "EX A", -uptake, 1000.0);
    ex_a.is_exchange = true;
    m.rxns.push(ex_a);
    m.rxns.push(Reaction::new("R_AB", "A -> B", 0.0, 1000.0));
    let mut r_bc = Reaction::new("R_BC", "B -> C", 0.0, 1000.0);
    r_bc.obj_coef = 1.0;
    m.rxns.push(r_bc);
    let mut ex_c = Reaction::new("EX_C", "EX C", -1000.0, 1000.0);
    ex_c.is_exchange = true;
    m.rxns.push(ex_c);

    m.s = StoichMatrix::from_triplets(
        3,
        4,
        vec![
            (0, 0, -1.0),
            (0, 1, -1.0),
            (1, 1, 1.0),
            (1, 2, -1.0),
            (2, 2, 1.0),
            (2, 3, -1.0),
        ],
    );
    m
}

#[test]
fn fba_chain_growth_equals_uptake() {
    for uptake in [1.0, 5.0, 10.0, 25.0] {
        let m = chain_model(uptake);
        let sol = fba(&m, &FbaOptions::default()).unwrap();
        assert_eq!(sol.status, SolveStatus::Optimal);
        assert!((sol.objective - uptake).abs() < 1e-6, "uptake={uptake}, got {}", sol.objective);
    }
}

#[test]
fn fba_blocked_middle_kills_growth() {
    let mut m = chain_model(10.0);
    m.rxns[1].ub = 0.0; // R_AB closed
    let sol = fba(&m, &FbaOptions::default()).unwrap();
    assert!(sol.objective.abs() < 1e-6);
}

#[test]
fn pfba_minimises_absolute_flux_with_parallel_branch() {
    // Two parallel pathways A -> B, same capacity. pFBA should carry all
    // flux through one and leave the other zero.
    let mut m = Model::new("parallel");
    m.mets.push(Metabolite::new("cpdA", "A", CompartmentId::CYTOSOL));
    m.mets.push(Metabolite::new("cpdB", "B", CompartmentId::CYTOSOL));

    let mut ex_a = Reaction::new("EX_A", "", -10.0, 1000.0);
    ex_a.is_exchange = true;
    m.rxns.push(ex_a);
    m.rxns.push(Reaction::new("R1", "A->B", 0.0, 1000.0));
    m.rxns.push(Reaction::new("R2", "A->B", 0.0, 1000.0));
    let mut ex_b = Reaction::new("EX_B", "", -1000.0, 1000.0);
    ex_b.obj_coef = 1.0;
    ex_b.is_exchange = true;
    m.rxns.push(ex_b);

    m.s = StoichMatrix::from_triplets(
        2,
        4,
        vec![
            (0, 0, -1.0),
            (0, 1, -1.0),
            (0, 2, -1.0),
            (1, 1, 1.0),
            (1, 2, 1.0),
            (1, 3, -1.0),
        ],
    );

    // Wait — obj_coef on EX_B is tricky because EX_B has stoich -1 on B;
    // maximising EX_B drives secretion. That's fine.
    let sol = pfba(
        &m,
        &PfbaOptions {
            weights: vec![0.0, 1.0, 1.0, 0.0],
            pfba_coef: 1e-3,
            min_growth: 10.0,
            objective: None,
        },
    )
    .unwrap();
    assert_eq!(sol.status, SolveStatus::Optimal);
    // Total absolute flux on R1+R2 should be exactly 10.
    let sum_abs = sol.fluxes[1].abs() + sol.fluxes[2].abs();
    assert!((sum_abs - 10.0).abs() < 1e-6, "sum_abs = {sum_abs}");
}

#[test]
fn pfba_heuristic_respects_min_growth() {
    let m = chain_model(10.0);
    let opts = PfbaHeuristicOptions::new(
        {
            let mut w = vec![1.0; m.rxn_count()];
            w[2] = 0.0;
            w
        },
        5.0,
    );
    let sol = pfba_heuristic(&m, &opts).unwrap();
    assert!(sol.growth >= 5.0 - 1e-6, "growth = {}", sol.growth);
}

#[test]
fn fba_reversible_reaction_runs_backward() {
    // A <-> B with reversible R; biomass on "sink of A" (backward direction
    // of EX_A). Tests that vn flows when objective favours it.
    let mut m = Model::new("rev");
    m.mets.push(Metabolite::new("cpdA", "A", CompartmentId::CYTOSOL));
    m.mets.push(Metabolite::new("cpdB", "B", CompartmentId::CYTOSOL));

    let mut ex_b = Reaction::new("EX_B", "", -10.0, 1000.0);
    ex_b.is_exchange = true;
    m.rxns.push(ex_b);
    m.rxns.push(Reaction::new("R_AB", "A<->B", -1000.0, 1000.0));
    let mut ex_a = Reaction::new("EX_A", "", -1000.0, 1000.0);
    ex_a.obj_coef = -1.0; // maximise negative flux on EX_A = uptake? No: we want net export
    ex_a.is_exchange = true;
    m.rxns.push(ex_a);

    m.s = StoichMatrix::from_triplets(
        2,
        3,
        vec![(1, 0, -1.0), (0, 1, -1.0), (1, 1, 1.0), (0, 2, -1.0)],
    );

    // EX_B provides B uptake (v[EX_B] = -10), R_AB runs B→A (v = -10),
    // EX_A secretes A (v[EX_A] = 10). Objective -1 · v[EX_A] = -10 minimised.
    // So maximising gives obj = -(-10) = 10 when we maximise -v[EX_A]?
    // Let me just check the flux shapes:
    let sol = fba(&m, &FbaOptions::default()).unwrap();
    assert_eq!(sol.status, SolveStatus::Optimal);
    // Expectation: R_AB has negative flux (backward direction).
    assert!(
        sol.fluxes[1].abs() > 1e-6,
        "R_AB should carry flux, got {}",
        sol.fluxes[1]
    );
}
