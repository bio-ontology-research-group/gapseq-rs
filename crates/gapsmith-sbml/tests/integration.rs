//! End-to-end validation of the emitted SBML against structural invariants.
//!
//! Without libsbml or cobrapy available in CI we cannot run the canonical
//! validator. Instead we assert a tight set of invariants a COBRA consumer
//! relies on:
//!
//! - root `<sbml>` with the three required namespaces and `fbc:strict`.
//! - every reaction has `fbc:lowerFluxBound` and `fbc:upperFluxBound`
//!   pointing at an emitted `<parameter>`.
//! - every `speciesReference` points at an emitted `<species>`.
//! - every `fbc:fluxObjective` points at an emitted reaction.
//! - every `fbc:geneProductRef` points at an emitted `<fbc:geneProduct>`.
//! - every `groups:member/@idRef` points at an emitted reaction.
//! - element counts match the input model.

use gapsmith_core::{CompartmentId, Metabolite, Model, Reaction, StoichMatrix};
use gapsmith_sbml::{write_to, WriteOptions};
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use std::collections::HashSet;

fn build_bigger_model() -> Model {
    // 5 mets, 4 rxns — covers: irreversible rxn, reversible exchange, custom
    // bound (0.01 lower), biomass-ish objective, multi-gene GPR, subsystems.
    let mut m = Model::new("bigger");
    m.annot.gapsmith_version = Some("0.1.0".into());
    m.annot.tax_domain = Some("Bacteria".into());
    for (i, (cpd, name)) in [
        ("cpd00001", "H2O"),
        ("cpd00002", "ATP"),
        ("cpd00007", "O2"),
        ("cpd00011", "CO2"),
        ("cpd11416", "biomass"),
    ]
    .iter()
    .enumerate()
    {
        let comp = if i == 2 { CompartmentId::EXTRACELLULAR } else { CompartmentId::CYTOSOL };
        let mut met = Metabolite::new(*cpd, *name, comp);
        met.charge = 0;
        m.mets.push(met);
    }

    let mut r1 = Reaction::new("rxn00001", "ATPase", 0.0, 1000.0);
    r1.gpr_raw = Some("b0001 and b0002".into());
    r1.subsystem = Some("Core".into());
    m.rxns.push(r1);

    let mut r2 = Reaction::new("rxn00099", "biomass", 0.01, 1000.0);
    r2.obj_coef = 1.0;
    r2.is_biomass = true;
    r2.subsystem = Some("Biomass".into());
    m.rxns.push(r2);

    let mut r3 = Reaction::new("EX_cpd00007_e0", "O2 exchange", -1000.0, 1000.0);
    r3.is_exchange = true;
    m.rxns.push(r3);

    let mut r4 = Reaction::new("rxn00007", "branched", -50.0, 50.0);
    r4.gpr_raw = Some("(b0003 or b0004) and b0005".into());
    r4.subsystem = Some("Core".into());
    m.rxns.push(r4);

    // Stoichiometry: meaningless but valid for emission.
    m.s = StoichMatrix::from_triplets(
        5,
        4,
        vec![
            (0, 0, 1.0),
            (1, 0, -1.0),
            (4, 1, 1.0),
            (2, 2, -1.0),
            (3, 3, -1.0),
            (1, 3, 1.0),
        ],
    );
    m
}

#[test]
fn structural_invariants_hold() {
    let m = build_bigger_model();
    let mut buf = Vec::new();
    write_to(&m, &mut buf, &WriteOptions::default()).unwrap();
    let xml = String::from_utf8(buf).unwrap();

    // Collect IDs by element type.
    let mut rdr = Reader::from_str(&xml);
    let mut species: HashSet<String> = HashSet::new();
    let mut reactions: HashSet<String> = HashSet::new();
    let mut parameters: HashSet<String> = HashSet::new();
    let mut gene_products: HashSet<String> = HashSet::new();
    let mut species_refs: Vec<String> = Vec::new();
    let mut rxn_lb_refs: Vec<String> = Vec::new();
    let mut rxn_ub_refs: Vec<String> = Vec::new();
    let mut flux_obj_rxns: Vec<String> = Vec::new();
    let mut gpr_refs: Vec<String> = Vec::new();
    let mut group_member_refs: Vec<String> = Vec::new();
    let mut root_attrs: Vec<(String, String)> = Vec::new();
    let mut model_attrs: Vec<(String, String)> = Vec::new();
    let mut buf2 = Vec::new();

    loop {
        match rdr.read_event_into(&mut buf2).unwrap() {
            Event::Start(e) | Event::Empty(e) => {
                let tag = String::from_utf8_lossy(e.name().as_ref()).to_string();
                let attrs: Vec<(String, String)> = e
                    .attributes()
                    .filter_map(Result::ok)
                    .map(|a| {
                        (
                            String::from_utf8_lossy(a.key.as_ref()).to_string(),
                            String::from_utf8_lossy(&a.value).to_string(),
                        )
                    })
                    .collect();
                let get = |k: &str| attrs.iter().find(|(ak, _)| ak == k).map(|(_, v)| v.clone());

                match tag.as_str() {
                    "sbml" => root_attrs = attrs.clone(),
                    "model" => model_attrs = attrs.clone(),
                    "species" => {
                        if let Some(id) = get("id") {
                            species.insert(id);
                        }
                    }
                    "reaction" => {
                        if let Some(id) = get("id") {
                            reactions.insert(id);
                        }
                        if let Some(v) = get("fbc:lowerFluxBound") {
                            rxn_lb_refs.push(v);
                        }
                        if let Some(v) = get("fbc:upperFluxBound") {
                            rxn_ub_refs.push(v);
                        }
                    }
                    "parameter" => {
                        if let Some(id) = get("id") {
                            parameters.insert(id);
                        }
                    }
                    "speciesReference" => {
                        if let Some(v) = get("species") {
                            species_refs.push(v);
                        }
                    }
                    "fbc:fluxObjective" => {
                        if let Some(v) = get("fbc:reaction") {
                            flux_obj_rxns.push(v);
                        }
                    }
                    "fbc:geneProduct" => {
                        if let Some(id) = get("fbc:id") {
                            gene_products.insert(id);
                        }
                    }
                    "fbc:geneProductRef" => {
                        if let Some(v) = get("fbc:geneProduct") {
                            gpr_refs.push(v);
                        }
                    }
                    "groups:member" => {
                        if let Some(v) = get("groups:idRef") {
                            group_member_refs.push(v);
                        }
                    }
                    _ => {}
                }
            }
            Event::Eof => break,
            _ => {}
        }
        buf2.clear();
    }

    // Root namespaces and fbc:strict.
    let root_has = |k: &str, v: &str| {
        root_attrs.iter().any(|(a, b)| a == k && b == v)
    };
    assert!(root_has("xmlns", "http://www.sbml.org/sbml/level3/version1/core"));
    assert!(root_has("xmlns:fbc", "http://www.sbml.org/sbml/level3/version1/fbc/version2"));
    assert!(root_has("xmlns:groups", "http://www.sbml.org/sbml/level3/version1/groups/version1"));
    assert!(model_attrs.iter().any(|(k, v)| k == "fbc:strict" && v == "true"));

    // Counts.
    assert_eq!(species.len(), 5);
    assert_eq!(reactions.len(), 4);

    // Every speciesReference resolves.
    for sr in &species_refs {
        assert!(species.contains(sr), "speciesRef `{sr}` missing from <listOfSpecies>");
    }
    // Every reaction flux-bound resolves to a parameter.
    for p in rxn_lb_refs.iter().chain(rxn_ub_refs.iter()) {
        assert!(parameters.contains(p), "flux bound `{p}` missing from <listOfParameters>");
    }
    // Custom bound was emitted.
    assert!(parameters.contains("R_rxn00099_lower_bound"));
    assert!(parameters.contains("R_rxn00007_lower_bound"));
    assert!(parameters.contains("R_rxn00007_upper_bound"));

    // Objective points at an emitted reaction.
    assert_eq!(flux_obj_rxns, vec!["R_rxn00099"]);

    // GPR references resolve.
    for g in &gpr_refs {
        assert!(gene_products.contains(g), "gene product `{g}` missing");
    }
    // All five genes present.
    assert_eq!(gene_products.len(), 5);

    // Group members resolve to reactions.
    for m in &group_member_refs {
        assert!(reactions.contains(m), "group member `{m}` missing");
    }
    // Exactly 3 reactions carry subsystems (rxn00001, rxn00099, rxn00007).
    assert_eq!(group_member_refs.len(), 3);
}
