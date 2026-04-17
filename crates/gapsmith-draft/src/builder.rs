//! Core model builder. Given a selected set of SEED reactions + a
//! biomass spec, assemble a [`gapsmith_core::Model`] with compartments,
//! metabolites, the sparse S matrix, and a `bio1` reaction.

use crate::biomass::BiomassSpec;
use crate::gpr::{build_gpr_string, GeneAssignment};
use gapsmith_core::{
    Compartment, CompartmentId, Metabolite, Model, Reaction, StoichMatrix,
};
use gapsmith_db::SeedRxnRow;
use std::collections::HashMap;

pub struct BuilderOptions {
    pub model_id: String,
    pub gapsmith_version: Option<String>,
    pub seqdb_version: Option<String>,
    pub tax_domain: Option<String>,
    pub gram: Option<String>,
}

pub fn build_model(
    opts: &BuilderOptions,
    selected_reactions: &[&SeedRxnRow],
    biomass: Option<&BiomassSpec>,
    gene_assignments_per_seed: &HashMap<String, Vec<(String, String)>>, // seed -> [(complex, gene)]
) -> Model {
    let mut m = Model::new(opts.model_id.clone());
    m.annot.gapsmith_version = opts.gapsmith_version.clone();
    m.annot.seqdb_version = opts.seqdb_version.clone();
    m.annot.tax_domain = opts.tax_domain.clone();
    m.annot.gram = opts.gram.clone();
    m.compartments = Compartment::default_three();

    for r in selected_reactions {
        let gene_rows = gene_assignments_per_seed.get(r.id.as_str()).cloned();
        add_seed_reaction_with_gpr(&mut m, r, None, gene_rows.as_deref());
    }

    if let Some(bm) = biomass {
        add_biomass(&mut m, bm);
    }
    rebuild_s_matrix(&mut m);
    m
}

/// Add one SEED reaction to the model. Creates any missing metabolites;
/// the S matrix column is NOT added here — call [`rebuild_s_matrix`]
/// once after all additions.
pub fn add_seed_reaction(model: &mut Model, row: &SeedRxnRow, gs_origin: Option<i8>) {
    add_seed_reaction_with_gpr(model, row, gs_origin, None);
}

fn add_seed_reaction_with_gpr(
    model: &mut Model,
    row: &SeedRxnRow,
    gs_origin: Option<i8>,
    gene_assignments: Option<&[(String, String)]>,
) {
    // If this id is already present, skip (caller dedup preferred).
    let rxn_id = format!("{}_c0", row.id.as_str());
    if model.rxns.iter().any(|r| r.id.as_str() == rxn_id) {
        return;
    }
    let terms = match row.parse_stoich() {
        Ok(t) if !t.is_empty() => t,
        _ => return,
    };

    // Ensure mets exist; collect stoich coefficients keyed by met_idx.
    // Metabolite IDs are compartment-unique (`<cpd>_<comp>`) so a single
    // compound living in multiple compartments gets distinct mets. The
    // SBML writer re-derives the `M_<cpd>_<comp>` species id from this.
    let mut col: Vec<(usize, f64)> = Vec::new();
    for t in &terms {
        let comp = match t.compartment {
            0 => CompartmentId::CYTOSOL,
            1 => CompartmentId::EXTRACELLULAR,
            _ => CompartmentId::PERIPLASM,
        };
        let comp_suffix = comp_suffix(comp);
        let unique_id = format!("{}_{comp_suffix}", t.cpd.as_str());
        let met_idx = match model.mets.iter().position(|m| m.id.as_str() == unique_id) {
            Some(i) => i,
            None => {
                let met = Metabolite::new(unique_id.as_str(), t.name.clone(), comp);
                model.mets.push(met);
                model.mets.len() - 1
            }
        };
        col.push((met_idx, t.coef));
    }

    // Decide bounds from reversibility.
    let (lb, ub) = match row.reversibility.as_str() {
        "=" => (-1000.0, 1000.0),
        "<" => (-1000.0, 0.0),
        _ => (0.0, 1000.0), // ">" or empty
    };

    // Build reaction object.
    let mut rxn = Reaction::new(rxn_id.as_str(), row.name.clone(), lb, ub);
    rxn.ec = row
        .ec_list()
        .into_iter()
        .map(|s| s.to_string())
        .collect();
    rxn.seed_status = row.gapseq_status;
    rxn.gs_origin = gs_origin;

    // Compose GPR if gene assignments were supplied.
    if let Some(assigns) = gene_assignments {
        if !assigns.is_empty() {
            let refs: Vec<GeneAssignment> = assigns
                .iter()
                .map(|(c, g)| GeneAssignment {
                    complex: if c.is_empty() { None } else { Some(c.as_str()) },
                    gene: g.as_str(),
                })
                .collect();
            let gpr = build_gpr_string(&refs);
            if !gpr.is_empty() {
                rxn.gpr_raw = Some(gpr);
            }
        }
    }

    // Stash the column on the model as a "pending" scratch (we rebuild
    // S at the end to avoid N² appends during draft construction). Use
    // reaction's `subsystem` field as the serialized stash — no: use a
    // side map returned from the builder's caller. Simpler: store
    // pending columns in a per-model vector that lives only at
    // construction time. We use a thread-local to avoid leaking this
    // detail through the public API.
    pending::push_column(model.rxns.len(), col);
    model.rxns.push(rxn);
}

fn comp_suffix(c: CompartmentId) -> &'static str {
    match c {
        CompartmentId::CYTOSOL => "c0",
        CompartmentId::EXTRACELLULAR => "e0",
        _ => "p0",
    }
}

fn add_biomass(model: &mut Model, bm: &BiomassSpec) {
    // Ensure every biomass metabolite exists.
    let mut col: Vec<(usize, f64)> = Vec::new();
    for e in &bm.entries {
        let comp = match e.compartment {
            'c' => CompartmentId::CYTOSOL,
            'e' => CompartmentId::EXTRACELLULAR,
            _ => CompartmentId::PERIPLASM,
        };
        let comp_suffix = comp_suffix(comp);
        let unique_id = format!("{}_{comp_suffix}", e.cpd);
        let met_idx = match model.mets.iter().position(|m| m.id.as_str() == unique_id) {
            Some(i) => i,
            None => {
                let mut met = Metabolite::new(unique_id.as_str(), e.name.clone(), comp);
                if !e.formula.is_empty() {
                    met.formula = Some(e.formula.clone());
                }
                model.mets.push(met);
                model.mets.len() - 1
            }
        };
        col.push((met_idx, e.scoef));
    }

    let mut r = Reaction::new("bio1", bm.name.clone(), 0.0, 1000.0);
    r.obj_coef = 1.0;
    r.is_biomass = true;
    r.gs_origin = Some(6);
    pending::push_column(model.rxns.len(), col);
    model.rxns.push(r);
}

/// Finalize the model's S matrix. Reads every non-zero currently in
/// `model.s` (if any) and merges any pending columns produced by
/// subsequent `add_seed_reaction` calls.
pub fn rebuild_s_matrix(model: &mut Model) {
    let n_mets = model.mets.len();
    let n_rxns = model.rxns.len();
    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(n_rxns * 4);

    // Preserve non-zeros already present in `model.s` (from earlier
    // flushes). `StoichMatrix::column` iterates the CSC column.
    let existing_cols = model.s.cols();
    for c in 0..existing_cols.min(n_rxns) {
        for (row, coef) in model.s.column(c) {
            if row < n_mets {
                triplets.push((row, c, coef));
            }
        }
    }

    // Apply pending updates — these may be for existing columns (which
    // should REPLACE any pre-existing entries) or for newly appended
    // columns. We handle that by collecting pending into a HashMap
    // keyed by (col, row), overlay on triplets.
    use std::collections::HashMap;
    let cols = pending::take();
    let mut overrides: HashMap<(usize, usize), f64> = HashMap::new();
    let mut overridden_cols: std::collections::HashSet<usize> =
        std::collections::HashSet::new();
    for (col, entries) in cols {
        overridden_cols.insert(col);
        for (met_idx, coef) in entries {
            overrides.insert((met_idx, col), coef);
        }
    }
    // Drop any earlier triplet in an overridden column — pending data
    // takes precedence.
    if !overridden_cols.is_empty() {
        triplets.retain(|(_, c, _)| !overridden_cols.contains(c));
        for ((row, col), coef) in overrides {
            triplets.push((row, col, coef));
        }
    }

    model.s = StoichMatrix::from_triplets(n_mets, n_rxns, triplets);
}

mod pending {
    //! Thread-local scratch area used while building a Model. Cleared
    //! at the end of each build via [`take`].

    use std::cell::RefCell;

    thread_local! {
        static PENDING: RefCell<Vec<(usize, Vec<(usize, f64)>)>> = const { RefCell::new(Vec::new()) };
    }

    pub fn push_column(rxn_idx: usize, entries: Vec<(usize, f64)>) {
        PENDING.with(|p| p.borrow_mut().push((rxn_idx, entries)));
    }

    pub fn take() -> Vec<(usize, Vec<(usize, f64)>)> {
        PENDING.with(|p| std::mem::take(&mut *p.borrow_mut()))
    }
}
