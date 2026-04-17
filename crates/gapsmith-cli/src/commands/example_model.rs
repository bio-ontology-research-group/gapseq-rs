//! `gapsmith example-model` — emit a tiny hand-built model as CBOR or JSON.
//!
//! Purpose: smoke-test serialization end-to-end without requiring the full
//! reference-data stack. Also handy as a minimal fixture for downstream
//! tooling that wants to load a [`Model`].

use clap::Parser;
use gapsmith_core::{CompartmentId, Metabolite, Model, Reaction, StoichMatrix};
use gapsmith_io::{write_model_cbor, write_model_json, ModelFormat};
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    /// Output path. Format inferred from extension (`.json` → JSON, else CBOR).
    pub output: PathBuf,

    /// Pretty-print if writing JSON.
    #[arg(long)]
    pub pretty: bool,

    /// Emit a larger sample with GPR trees, subsystems, custom bounds, and
    /// a biomass reaction. Useful for end-to-end SBML validation.
    #[arg(long)]
    pub complex: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let model = if args.complex { build_complex() } else { build_toy() };
    let fmt = ModelFormat::from_path(&args.output);
    match fmt {
        ModelFormat::Cbor => write_model_cbor(&model, &args.output)?,
        ModelFormat::Json => write_model_json(&model, &args.output, args.pretty)?,
    }
    eprintln!(
        "wrote example model to `{}` ({} mets, {} rxns, {} nnz)",
        args.output.display(),
        model.met_count(),
        model.rxn_count(),
        model.s.nnz()
    );
    Ok(())
}

/// Build the standard tiny demo model: three metabolites in cytosol, one
/// ATP-hydrolysis reaction, one O2 exchange reaction. Just enough to exercise
/// every field of the `Model` type on serialization.
pub fn build_toy() -> Model {
    let mut m = Model::new("toy_ecoli");
    m.annot.name = Some("toy E. coli sanity model".into());
    m.annot.gapsmith_version = Some(env!("CARGO_PKG_VERSION").into());
    m.annot.tax_domain = Some("Bacteria".into());
    m.annot.gram = Some("neg".into());

    m.mets.push(Metabolite::new("cpd00001", "H2O", CompartmentId::CYTOSOL));
    m.mets.push(Metabolite::new("cpd00002", "ATP", CompartmentId::CYTOSOL));
    m.mets.push(Metabolite::new("cpd00007", "O2", CompartmentId::EXTRACELLULAR));

    let mut r1 = Reaction::new("rxn00001", "ATP hydrolysis", 0.0, 1000.0);
    r1.obj_coef = 1.0;
    r1.bitscore = Some(312.5);
    r1.gs_origin = Some(0);
    m.rxns.push(r1);

    let mut ex = Reaction::new("EX_cpd00007_e0", "O2 exchange", -1000.0, 1000.0);
    ex.is_exchange = true;
    ex.gs_origin = Some(7);
    m.rxns.push(ex);

    // 3 mets × 2 rxns. Column 0: H2O produced, ATP consumed. Column 1: O2 boundary.
    m.s = StoichMatrix::from_triplets(
        3,
        2,
        vec![(0, 0, 1.0), (1, 0, -1.0), (2, 1, -1.0)],
    );
    m
}

/// Larger sample covering features absent from the bare toy: GPR trees,
/// subsystems (→ `groups:group`), custom per-reaction bounds, and a biomass
/// objective.
pub fn build_complex() -> Model {
    let mut m = Model::new("complex_demo");
    m.annot.name = Some("Complex demo model".into());
    m.annot.gapsmith_version = Some(env!("CARGO_PKG_VERSION").into());
    m.annot.tax_domain = Some("Bacteria".into());
    m.annot.gram = Some("neg".into());

    for (cpd, name, comp) in [
        ("cpd00001", "H2O", CompartmentId::CYTOSOL),
        ("cpd00002", "ATP", CompartmentId::CYTOSOL),
        ("cpd00007", "O2", CompartmentId::EXTRACELLULAR),
        ("cpd00011", "CO2", CompartmentId::CYTOSOL),
        ("cpd11416", "biomass", CompartmentId::CYTOSOL),
    ] {
        m.mets.push(Metabolite::new(cpd, name, comp));
    }

    let mut r1 = Reaction::new("rxn00001", "ATPase", 0.0, 1000.0);
    r1.gpr_raw = Some("b0001 and b0002".into());
    r1.subsystem = Some("Core metabolism".into());
    m.rxns.push(r1);

    let mut bio = Reaction::new("bio1", "biomass", 0.01, 1000.0);
    bio.obj_coef = 1.0;
    bio.is_biomass = true;
    bio.subsystem = Some("Biomass".into());
    m.rxns.push(bio);

    let mut ex = Reaction::new("EX_cpd00007_e0", "O2 exchange", -1000.0, 1000.0);
    ex.is_exchange = true;
    m.rxns.push(ex);

    let mut rb = Reaction::new("rxn00007", "branched reaction", -50.0, 50.0);
    rb.gpr_raw = Some("(b0003 or b0004) and b0005".into());
    rb.subsystem = Some("Core metabolism".into());
    m.rxns.push(rb);

    m.s = gapsmith_core::StoichMatrix::from_triplets(
        5,
        4,
        vec![
            (0, 0, 1.0),
            (1, 0, -1.0),
            (4, 1, 1.0),
            (1, 1, -1.0),
            (2, 2, -1.0),
            (3, 3, -1.0),
            (1, 3, 1.0),
        ],
    );
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toy_is_self_consistent() {
        let m = build_toy();
        m.check_shape().unwrap();
    }

    #[test]
    fn complex_is_self_consistent() {
        let m = build_complex();
        m.check_shape().unwrap();
    }
}
