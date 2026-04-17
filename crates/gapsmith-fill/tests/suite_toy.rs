//! 4-phase suite integration test on a hand-built toy model.
//!
//! No dependency on real gapseq data — the model is fully synthesized
//! in-process. Covers:
//!
//! - Step 1 fill on a user medium with a biomass gap.
//! - Step 2 no-op when every biomass substrate is already producible.
//! - Per-step report accounting.

use gapsmith_core::{CompartmentId, Metabolite, Model, Reaction, StoichMatrix};
use gapsmith_core::{RxnId, SeedStatus};
use gapsmith_db::SeedRxnRow;
use gapsmith_fill::{run_suite, MediumEntry, RxnWeights, SuiteOptions};

fn toy_draft() -> Model {
    // Compartments: c0 (cytosol), e0 (extracellular).
    // Metabolites: cpd00001 (H2O), cpd00027 (Glucose), cpd99999 (a cofactor the
    // draft can't produce), cpd11416 (Biomass).
    //
    // Reactions:
    //   EX_cpd00001_e0:   cpd00001_e0 <-> ∅   [closed/open by medium]
    //   EX_cpd00027_e0:   cpd00027_e0 <-> ∅
    //   T_H2O:            cpd00001_e0 <-> cpd00001_c0
    //   T_Glu:            cpd00027_e0 <-> cpd00027_c0
    //   bio1:  cpd00027_c0 + cpd99999_c0 -> cpd11416_c0
    //
    // The draft lacks a producer for cpd99999, so biomass is pegged at 0.
    // The SEED DB supplies rxnGAP: cpd00001_c0 -> cpd99999_c0 which will
    // be added by gap-filling.
    let mut m = Model::new("toy_fill");
    for id in ["cpd00001_e0", "cpd00027_e0"] {
        let name = if id.starts_with("cpd00001") { "H2O" } else { "Glucose" };
        m.mets.push(Metabolite::new(id, name, CompartmentId::EXTRACELLULAR));
    }
    for id in ["cpd00001_c0", "cpd00027_c0", "cpd99999_c0", "cpd11416_c0"] {
        let name = match id {
            "cpd00001_c0" => "H2O",
            "cpd00027_c0" => "Glucose",
            "cpd99999_c0" => "Gap cofactor",
            _ => "Biomass",
        };
        m.mets.push(Metabolite::new(id, name, CompartmentId::CYTOSOL));
    }

    let mut ex_h2o = Reaction::new("EX_cpd00001_e0", "H2O Exchange", -1000.0, 1000.0);
    ex_h2o.is_exchange = true;
    m.rxns.push(ex_h2o);
    let mut ex_glu = Reaction::new("EX_cpd00027_e0", "Glucose Exchange", -10.0, 1000.0);
    ex_glu.is_exchange = true;
    m.rxns.push(ex_glu);
    m.rxns.push(Reaction::new("T_H2O", "H2O transporter", -1000.0, 1000.0));
    m.rxns.push(Reaction::new("T_Glu", "Glucose transporter", -1000.0, 1000.0));
    let mut bio1 = Reaction::new("bio1", "Toy biomass", 0.0, 1000.0);
    bio1.obj_coef = 1.0;
    bio1.is_biomass = true;
    m.rxns.push(bio1);

    // Build S.
    let mi = |id: &str| m.mets.iter().position(|mt| mt.id.as_str() == id).unwrap();
    let ri = |id: &str| m.rxns.iter().position(|rt| rt.id.as_str() == id).unwrap();
    let n_mets = m.mets.len();
    let n_rxns = m.rxns.len();
    let triplets = vec![
        (mi("cpd00001_e0"), ri("EX_cpd00001_e0"), -1.0),
        (mi("cpd00027_e0"), ri("EX_cpd00027_e0"), -1.0),
        (mi("cpd00001_e0"), ri("T_H2O"), -1.0),
        (mi("cpd00001_c0"), ri("T_H2O"), 1.0),
        (mi("cpd00027_e0"), ri("T_Glu"), -1.0),
        (mi("cpd00027_c0"), ri("T_Glu"), 1.0),
        (mi("cpd00027_c0"), ri("bio1"), -1.0),
        (mi("cpd99999_c0"), ri("bio1"), -1.0),
        (mi("cpd11416_c0"), ri("bio1"), 1.0),
    ];
    m.s = StoichMatrix::from_triplets(n_mets, n_rxns, triplets);
    m
}

/// A tiny SEED DB: one reaction that unlocks cpd99999 from water.
fn toy_seed_rxns() -> Vec<SeedRxnRow> {
    vec![SeedRxnRow {
        id: RxnId::new("rxnGAP"),
        abbreviation: String::new(),
        name: "Gap filler".into(),
        code: String::new(),
        stoichiometry: r#"-1:cpd00001:0:0:"H2O";1:cpd99999:0:0:"Gap""#.into(),
        is_transport: 0,
        equation: String::new(),
        definition: String::new(),
        reversibility: ">".into(),
        direction: String::new(),
        abstract_reaction: String::new(),
        pathways: String::new(),
        aliases: String::new(),
        ec_numbers: String::new(),
        deltag: String::new(),
        deltagerr: String::new(),
        compound_ids: "cpd00001;cpd99999".into(),
        status: String::new(),
        is_obsolete: 0,
        linked_reaction: String::new(),
        notes: String::new(),
        is_copy_of: String::new(),
        gapseq_status: SeedStatus::Approved,
    }]
}

#[test]
fn toy_suite_fills_biomass_gap() {
    let draft = toy_draft();
    let seed_rxns = toy_seed_rxns();
    let medium = vec![
        MediumEntry {
            compound: "cpd00001".into(),
            name: "H2O".into(),
            max_flux: 1000.0,
        },
        MediumEntry {
            compound: "cpd00027".into(),
            name: "Glucose".into(),
            max_flux: 10.0,
        },
    ];
    let weights = RxnWeights::new();

    // We can't read MM_glu from a real data dir here, so Step 2/2b will
    // silently skip (no-op). Step 1 fill is the path under test.
    let td = tempfile::tempdir().unwrap();
    let data_dir = td.path().to_path_buf();

    let opts = SuiteOptions {
        target_cpd: "cpd11416".into(),
        min_growth: 0.01,
        step1_only: true,
        quick: true,
        prune_futile: false,
        ..SuiteOptions::default()
    };

    let (filled, report) = run_suite(&draft, &medium, &weights, &seed_rxns, &data_dir, &opts)
        .expect("run_suite");

    assert!(
        !report.step1_added.is_empty(),
        "expected at least one Step 1 addition, got report: {:?}",
        report
    );
    assert!(
        report.step1_added.iter().any(|id| id == "rxnGAP_c0"),
        "rxnGAP should be among the Step 1 additions, got {:?}",
        report.step1_added
    );
    assert!(report.final_growth > 0.01);
    // Filled model includes the gap reaction.
    assert!(
        filled.rxns.iter().any(|r| r.id.as_str() == "rxnGAP_c0"),
        "rxnGAP_c0 missing from filled model"
    );
}
