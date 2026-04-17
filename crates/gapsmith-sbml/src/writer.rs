//! Streaming SBML writer.
//!
//! Design: a single pass over the model, driving a `quick_xml::Writer`. We
//! build no intermediate DOM — every element is emitted in order. That keeps
//! memory flat even for very large genome-scale models (~20k reactions) and
//! makes the code easy to follow top-to-bottom.

use crate::namespaces::*;
use crate::DEFAULT_BOUND;
use gapsmith_core::{Gpr, Model};
use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
use quick_xml::writer::Writer;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum SbmlError {
    #[error("i/o error on `{path}`: {source}")]
    Io {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("write error: {0}")]
    Write(#[from] std::io::Error),
    #[error("XML error: {0}")]
    Xml(#[from] quick_xml::Error),
    #[error("GPR parse error on reaction `{rxn}`: {source}")]
    Gpr {
        rxn: String,
        #[source]
        source: gapsmith_core::GprParseError,
    },
    #[error("model shape mismatch: {0}")]
    ShapeMismatch(#[from] gapsmith_core::ModelError),
}

/// Options for [`write_sbml`]. All fields are optional; defaults produce a
/// COBRApy-compatible document.
#[derive(Debug, Clone)]
pub struct WriteOptions {
    /// Whether to indent / pretty-print. Default: true.
    pub pretty: bool,
    /// Objective identifier. Default: `"obj"`.
    pub objective_id: String,
    /// Sense of the objective. Default: `"maximize"`.
    pub objective_sense: ObjectiveSense,
}

#[derive(Debug, Clone, Copy)]
pub enum ObjectiveSense {
    Maximize,
    Minimize,
}

impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            pretty: true,
            objective_id: "obj".into(),
            objective_sense: ObjectiveSense::Maximize,
        }
    }
}

/// Entry point: write `model` as SBML to `path`.
pub fn write_sbml(
    model: &Model,
    path: impl AsRef<Path>,
    opts: &WriteOptions,
) -> Result<(), SbmlError> {
    model.check_shape()?;
    let path = path.as_ref();
    let f = File::create(path).map_err(|e| SbmlError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    let mut sink = BufWriter::new(f);
    write_to(model, &mut sink, opts)?;
    sink.flush().map_err(|e| SbmlError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    Ok(())
}

/// Write the SBML document to any [`Write`] sink. Used both by [`write_sbml`]
/// (after opening the file) and by the test suite (to a memory buffer).
pub fn write_to<W: Write>(
    model: &Model,
    sink: &mut W,
    opts: &WriteOptions,
) -> Result<(), SbmlError> {
    let mut w: Writer<&mut W> = if opts.pretty {
        Writer::new_with_indent(sink, b' ', 2)
    } else {
        Writer::new(sink)
    };

    // XML declaration.
    w.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None)))?;

    // <sbml ...>
    let mut sbml = BytesStart::new("sbml");
    sbml.push_attribute(("xmlns", CORE));
    sbml.push_attribute((ATTR_XMLNS_FBC, FBC));
    sbml.push_attribute((ATTR_XMLNS_GROUPS, GROUPS));
    sbml.push_attribute(("level", "3"));
    sbml.push_attribute(("version", "1"));
    sbml.push_attribute((ATTR_FBC_REQUIRED, "false"));
    sbml.push_attribute((ATTR_GROUPS_REQUIRED, "false"));
    w.write_event(Event::Start(sbml))?;

    write_model(&mut w, model, opts)?;

    w.write_event(Event::End(BytesEnd::new("sbml")))?;
    Ok(())
}

fn write_model<W: Write>(
    w: &mut Writer<&mut W>,
    model: &Model,
    opts: &WriteOptions,
) -> Result<(), SbmlError> {
    let mut el = BytesStart::new("model");
    // SBML SIds require `[A-Za-z_][A-Za-z0-9_]*` — sanitize to avoid
    // COBRApy warnings when the model id contains `-`, `.`, etc. We
    // keep the original as the `name` so the provenance survives.
    let sanitized = sanitize_sid(model.annot.id.as_str());
    el.push_attribute(("id", sanitized.as_str()));
    if let Some(ref name) = model.annot.name {
        el.push_attribute(("name", name.as_str()));
    } else if sanitized != model.annot.id {
        el.push_attribute(("name", model.annot.id.as_str()));
    }
    // Model-level default units. libSBML emits a warning for every compartment
    // or species without declared units unless the enclosing model carries
    // defaults — so we set them once here.
    el.push_attribute(("substanceUnits", "substance"));
    el.push_attribute(("timeUnits", "time"));
    el.push_attribute(("volumeUnits", "volume"));
    el.push_attribute(("extentUnits", "substance"));
    el.push_attribute((ATTR_FBC_STRICT, "true"));
    w.write_event(Event::Start(el))?;

    write_notes(w, model)?;
    write_unit_definitions(w)?;
    write_compartments(w, model)?;
    write_species(w, model)?;
    let bounds = resolve_bound_parameters(model);
    write_parameters(w, &bounds)?;
    write_reactions(w, model, &bounds)?;
    write_objectives(w, model, opts)?;
    write_gene_products(w, model)?;
    write_groups(w, model)?;

    w.write_event(Event::End(BytesEnd::new("model")))?;
    Ok(())
}

fn write_notes<W: Write>(w: &mut Writer<&mut W>, model: &Model) -> Result<(), SbmlError> {
    if model.annot.gapsmith_version.is_none()
        && model.annot.seqdb_version.is_none()
        && model.annot.tax_domain.is_none()
        && model.annot.gram.is_none()
        && model.annot.notes.is_empty()
    {
        return Ok(());
    }
    w.write_event(Event::Start(BytesStart::new("notes")))?;

    let mut body = BytesStart::new("body");
    body.push_attribute(("xmlns", XHTML));
    w.write_event(Event::Start(body))?;

    if let Some(v) = &model.annot.gapsmith_version {
        write_note(w, "gapsmith version", v)?;
    }
    if let Some(v) = &model.annot.seqdb_version {
        write_note(w, "sequence DB version", v)?;
    }
    if let Some(v) = &model.annot.tax_domain {
        write_note(w, "tax_domain", v)?;
    }
    if let Some(v) = &model.annot.gram {
        write_note(w, "gram", v)?;
    }
    for n in &model.annot.notes {
        w.write_event(Event::Start(BytesStart::new("p")))?;
        w.write_event(Event::Text(BytesText::new(n)))?;
        w.write_event(Event::End(BytesEnd::new("p")))?;
    }

    w.write_event(Event::End(BytesEnd::new("body")))?;
    w.write_event(Event::End(BytesEnd::new("notes")))?;
    Ok(())
}

fn write_note<W: Write>(
    w: &mut Writer<&mut W>,
    key: &str,
    value: &str,
) -> Result<(), SbmlError> {
    w.write_event(Event::Start(BytesStart::new("p")))?;
    w.write_event(Event::Text(BytesText::new(&format!("{key}: {value}"))))?;
    w.write_event(Event::End(BytesEnd::new("p")))?;
    Ok(())
}

fn write_unit_definitions<W: Write>(w: &mut Writer<&mut W>) -> Result<(), SbmlError> {
    w.write_event(Event::Start(BytesStart::new("listOfUnitDefinitions")))?;

    // Flux unit used on every reaction bound parameter.
    emit_unit_def(w, "mmol_per_gDW_per_hr", &[
        ("mole", 1, -3, 1.0),
        ("gram", -1, 0, 1.0),
        ("second", -1, 0, 3600.0),
    ])?;

    // Default substance unit referenced by `model/@substanceUnits` —
    // millimoles (matches cobrar convention).
    emit_unit_def(w, "substance", &[("mole", 1, -3, 1.0)])?;

    // Default time unit — hours.
    emit_unit_def(w, "time", &[("second", 1, 0, 3600.0)])?;

    // Default volume unit — litres (placeholder; gapseq models don't use
    // volumes meaningfully but libSBML wants a declaration).
    emit_unit_def(w, "volume", &[("litre", 1, 0, 1.0)])?;

    w.write_event(Event::End(BytesEnd::new("listOfUnitDefinitions")))?;
    Ok(())
}

fn emit_unit_def<W: Write>(
    w: &mut Writer<&mut W>,
    id: &str,
    units: &[(&str, i32, i32, f64)],
) -> Result<(), SbmlError> {
    let mut ud = BytesStart::new("unitDefinition");
    ud.push_attribute(("id", id));
    w.write_event(Event::Start(ud))?;
    w.write_event(Event::Start(BytesStart::new("listOfUnits")))?;
    for (kind, exp, scale, mult) in units {
        write_unit(w, kind, *exp, *scale, *mult)?;
    }
    w.write_event(Event::End(BytesEnd::new("listOfUnits")))?;
    w.write_event(Event::End(BytesEnd::new("unitDefinition")))?;
    Ok(())
}

fn write_unit<W: Write>(
    w: &mut Writer<&mut W>,
    kind: &str,
    exponent: i32,
    scale: i32,
    multiplier: f64,
) -> Result<(), SbmlError> {
    let mut el = BytesStart::new("unit");
    el.push_attribute(("kind", kind));
    let e = exponent.to_string();
    let s = scale.to_string();
    let m = fmt_f64(multiplier);
    el.push_attribute(("exponent", e.as_str()));
    el.push_attribute(("scale", s.as_str()));
    el.push_attribute(("multiplier", m.as_str()));
    w.write_event(Event::Empty(el))?;
    Ok(())
}

fn write_compartments<W: Write>(
    w: &mut Writer<&mut W>,
    model: &Model,
) -> Result<(), SbmlError> {
    w.write_event(Event::Start(BytesStart::new("listOfCompartments")))?;
    for c in &model.compartments {
        let mut el = BytesStart::new("compartment");
        el.push_attribute(("id", c.id.as_str()));
        if !c.name.is_empty() {
            el.push_attribute(("name", c.name.as_str()));
        }
        el.push_attribute(("constant", "true"));
        el.push_attribute(("spatialDimensions", "3"));
        el.push_attribute(("size", "1"));
        el.push_attribute(("units", "volume"));
        w.write_event(Event::Empty(el))?;
    }
    w.write_event(Event::End(BytesEnd::new("listOfCompartments")))?;
    Ok(())
}

fn write_species<W: Write>(w: &mut Writer<&mut W>, model: &Model) -> Result<(), SbmlError> {
    w.write_event(Event::Start(BytesStart::new("listOfSpecies")))?;
    for m in &model.mets {
        let comp_id = model
            .compartments
            .get(m.compartment.0 as usize)
            .map(|c| c.id.as_str())
            .unwrap_or("c0");
        let sid = species_id(m.id.as_str(), comp_id);

        let mut el = BytesStart::new("species");
        el.push_attribute(("id", sid.as_str()));
        el.push_attribute(("name", m.name.as_str()));
        el.push_attribute(("compartment", comp_id));
        el.push_attribute(("hasOnlySubstanceUnits", "true"));
        el.push_attribute(("boundaryCondition", "false"));
        el.push_attribute(("constant", "false"));
        let charge = m.charge.to_string();
        el.push_attribute((ATTR_FBC_CHARGE, charge.as_str()));
        if let Some(f) = &m.formula {
            if !f.is_empty() && f != "null" {
                el.push_attribute((ATTR_FBC_FORMULA, f.as_str()));
            }
        }
        w.write_event(Event::Empty(el))?;
    }
    w.write_event(Event::End(BytesEnd::new("listOfSpecies")))?;
    Ok(())
}

/// Resolved flux-bound parameters for every reaction.
///
/// Each reaction has `(lb_param_id, ub_param_id)`. Shared defaults
/// (`cobra_default_lb`, `cobra_default_ub`, `cobra_0_bound`) are reused; any
/// remaining numeric bounds get their own `R_<id>_<lower|upper>_bound`
/// parameter, emitted inside `<listOfParameters>`.
struct BoundResolution {
    per_rxn: HashMap<usize, (String, String)>,
    /// Deterministic (id, value) pairs for custom parameters.
    custom: Vec<(String, f64)>,
}

fn resolve_bound_parameters(model: &Model) -> BoundResolution {
    let mut per_rxn = HashMap::new();
    let mut customs: BTreeMap<String, f64> = BTreeMap::new();
    for (i, r) in model.rxns.iter().enumerate() {
        let lb_id = bound_param_id(r.lb, r.id.as_str(), "lower_bound");
        let ub_id = bound_param_id(r.ub, r.id.as_str(), "upper_bound");
        if is_custom(&lb_id) {
            customs.insert(lb_id.clone(), r.lb);
        }
        if is_custom(&ub_id) {
            customs.insert(ub_id.clone(), r.ub);
        }
        per_rxn.insert(i, (lb_id, ub_id));
    }
    BoundResolution { per_rxn, custom: customs.into_iter().collect() }
}

fn bound_param_id(v: f64, rxn_id: &str, suffix: &str) -> String {
    if v == -DEFAULT_BOUND {
        "cobra_default_lb".into()
    } else if v == DEFAULT_BOUND {
        "cobra_default_ub".into()
    } else if v == 0.0 {
        "cobra_0_bound".into()
    } else {
        format!("R_{rxn_id}_{suffix}")
    }
}

fn is_custom(id: &str) -> bool {
    !matches!(id, "cobra_default_lb" | "cobra_default_ub" | "cobra_0_bound")
}

fn write_parameters<W: Write>(
    w: &mut Writer<&mut W>,
    bounds: &BoundResolution,
) -> Result<(), SbmlError> {
    w.write_event(Event::Start(BytesStart::new("listOfParameters")))?;
    write_shared_parameter(w, "cobra_default_lb", -DEFAULT_BOUND)?;
    write_shared_parameter(w, "cobra_default_ub", DEFAULT_BOUND)?;
    write_shared_parameter(w, "cobra_0_bound", 0.0)?;
    for (id, v) in &bounds.custom {
        write_shared_parameter(w, id, *v)?;
    }
    w.write_event(Event::End(BytesEnd::new("listOfParameters")))?;
    Ok(())
}

fn write_shared_parameter<W: Write>(
    w: &mut Writer<&mut W>,
    id: &str,
    value: f64,
) -> Result<(), SbmlError> {
    let mut el = BytesStart::new("parameter");
    el.push_attribute(("id", id));
    let v = fmt_f64(value);
    el.push_attribute(("value", v.as_str()));
    el.push_attribute(("constant", "true"));
    el.push_attribute(("sboTerm", SBO_FLUX_BOUND));
    el.push_attribute(("units", "mmol_per_gDW_per_hr"));
    w.write_event(Event::Empty(el))?;
    Ok(())
}

fn write_reactions<W: Write>(
    w: &mut Writer<&mut W>,
    model: &Model,
    bounds: &BoundResolution,
) -> Result<(), SbmlError> {
    w.write_event(Event::Start(BytesStart::new("listOfReactions")))?;
    for (i, r) in model.rxns.iter().enumerate() {
        let (lb_id, ub_id) = &bounds.per_rxn[&i];
        let mut el = BytesStart::new("reaction");
        let rid = reaction_id(r.id.as_str());
        el.push_attribute(("id", rid.as_str()));
        if !r.name.is_empty() {
            el.push_attribute(("name", r.name.as_str()));
        }
        let reversible = matches!(r.reversibility(), gapsmith_core::Reversibility::Reversible);
        el.push_attribute(("reversible", if reversible { "true" } else { "false" }));
        el.push_attribute(("fast", "false"));
        el.push_attribute((ATTR_FBC_LOWER, lb_id.as_str()));
        el.push_attribute((ATTR_FBC_UPPER, ub_id.as_str()));
        w.write_event(Event::Start(el))?;

        write_stoich_refs(w, model, i)?;
        if let Some(raw) = &r.gpr_raw {
            if !raw.trim().is_empty() {
                let gpr: Gpr = raw.parse().map_err(|source| SbmlError::Gpr {
                    rxn: r.id.as_str().to_string(),
                    source,
                })?;
                write_gpr(w, &gpr)?;
            }
        }

        w.write_event(Event::End(BytesEnd::new("reaction")))?;
    }
    w.write_event(Event::End(BytesEnd::new("listOfReactions")))?;
    Ok(())
}

fn write_stoich_refs<W: Write>(
    w: &mut Writer<&mut W>,
    model: &Model,
    rxn_idx: usize,
) -> Result<(), SbmlError> {
    let col = model.s.column(rxn_idx);
    let (reactants, products): (Vec<_>, Vec<_>) = col.into_iter().partition(|(_, v)| *v < 0.0);

    emit_refs(w, model, &reactants, "listOfReactants")?;
    emit_refs(w, model, &products, "listOfProducts")?;
    Ok(())
}

fn emit_refs<W: Write>(
    w: &mut Writer<&mut W>,
    model: &Model,
    list: &[(usize, f64)],
    tag: &str,
) -> Result<(), SbmlError> {
    if list.is_empty() {
        return Ok(());
    }
    w.write_event(Event::Start(BytesStart::new(tag)))?;
    for &(row, v) in list {
        let met = &model.mets[row];
        let comp_id = model
            .compartments
            .get(met.compartment.0 as usize)
            .map(|c| c.id.as_str())
            .unwrap_or("c0");
        let sid = species_id(met.id.as_str(), comp_id);
        let mut s = BytesStart::new("speciesReference");
        s.push_attribute(("species", sid.as_str()));
        let stoich = fmt_f64(v.abs());
        s.push_attribute(("stoichiometry", stoich.as_str()));
        s.push_attribute(("constant", "true"));
        w.write_event(Event::Empty(s))?;
    }
    w.write_event(Event::End(BytesEnd::new(tag)))?;
    Ok(())
}

fn write_gpr<W: Write>(w: &mut Writer<&mut W>, gpr: &Gpr) -> Result<(), SbmlError> {
    w.write_event(Event::Start(BytesStart::new(TAG_FBC_GPA)))?;
    write_gpr_node(w, gpr)?;
    w.write_event(Event::End(BytesEnd::new(TAG_FBC_GPA)))?;
    Ok(())
}

fn write_gpr_node<W: Write>(w: &mut Writer<&mut W>, gpr: &Gpr) -> Result<(), SbmlError> {
    match gpr {
        Gpr::Gene { id } => {
            let mut el = BytesStart::new(TAG_FBC_GENEPRODUCTREF);
            let gid = gene_id(id.as_str());
            el.push_attribute((ATTR_FBC_GENEPRODUCT, gid.as_str()));
            w.write_event(Event::Empty(el))?;
        }
        Gpr::And { operands } => {
            w.write_event(Event::Start(BytesStart::new(TAG_FBC_AND)))?;
            for op in operands {
                write_gpr_node(w, op)?;
            }
            w.write_event(Event::End(BytesEnd::new(TAG_FBC_AND)))?;
        }
        Gpr::Or { operands } => {
            w.write_event(Event::Start(BytesStart::new(TAG_FBC_OR)))?;
            for op in operands {
                write_gpr_node(w, op)?;
            }
            w.write_event(Event::End(BytesEnd::new(TAG_FBC_OR)))?;
        }
    }
    Ok(())
}

fn write_objectives<W: Write>(
    w: &mut Writer<&mut W>,
    model: &Model,
    opts: &WriteOptions,
) -> Result<(), SbmlError> {
    let mut el = BytesStart::new(TAG_FBC_OBJECTIVES);
    el.push_attribute((ATTR_FBC_ACTIVE_OBJECTIVE, opts.objective_id.as_str()));
    w.write_event(Event::Start(el))?;

    let sense = match opts.objective_sense {
        ObjectiveSense::Maximize => "maximize",
        ObjectiveSense::Minimize => "minimize",
    };
    let mut obj = BytesStart::new(TAG_FBC_OBJECTIVE);
    obj.push_attribute((ATTR_FBC_ID, opts.objective_id.as_str()));
    obj.push_attribute((ATTR_FBC_TYPE, sense));
    w.write_event(Event::Start(obj))?;

    w.write_event(Event::Start(BytesStart::new(TAG_FBC_FLUX_OBJECTIVES)))?;
    for r in &model.rxns {
        if r.obj_coef != 0.0 {
            let mut fo = BytesStart::new(TAG_FBC_FLUX_OBJECTIVE);
            let rid = reaction_id(r.id.as_str());
            fo.push_attribute((ATTR_FBC_REACTION, rid.as_str()));
            let c = fmt_f64(r.obj_coef);
            fo.push_attribute((ATTR_FBC_COEFFICIENT, c.as_str()));
            w.write_event(Event::Empty(fo))?;
        }
    }
    w.write_event(Event::End(BytesEnd::new(TAG_FBC_FLUX_OBJECTIVES)))?;

    w.write_event(Event::End(BytesEnd::new(TAG_FBC_OBJECTIVE)))?;
    w.write_event(Event::End(BytesEnd::new(TAG_FBC_OBJECTIVES)))?;
    Ok(())
}

fn write_gene_products<W: Write>(
    w: &mut Writer<&mut W>,
    model: &Model,
) -> Result<(), SbmlError> {
    // Collect unique genes from reaction GPRs and from model.genes.
    let mut seen: std::collections::BTreeSet<String> = Default::default();
    for g in &model.genes {
        seen.insert(g.as_str().to_string());
    }
    for r in &model.rxns {
        if let Some(raw) = &r.gpr_raw {
            if let Ok(gpr) = raw.parse::<Gpr>() {
                let mut acc = Vec::new();
                gpr.collect_genes(&mut acc);
                for g in acc {
                    seen.insert(g.as_str().to_string());
                }
            }
        }
    }
    if seen.is_empty() {
        return Ok(());
    }
    w.write_event(Event::Start(BytesStart::new(TAG_FBC_GENEPRODUCTS)))?;
    for g in seen {
        let mut el = BytesStart::new(TAG_FBC_GENEPRODUCT);
        let gid = gene_id(&g);
        el.push_attribute((ATTR_FBC_ID, gid.as_str()));
        el.push_attribute((ATTR_FBC_LABEL, g.as_str()));
        w.write_event(Event::Empty(el))?;
    }
    w.write_event(Event::End(BytesEnd::new(TAG_FBC_GENEPRODUCTS)))?;
    Ok(())
}

fn write_groups<W: Write>(w: &mut Writer<&mut W>, model: &Model) -> Result<(), SbmlError> {
    let mut by_subsys: BTreeMap<&str, Vec<usize>> = Default::default();
    for (i, r) in model.rxns.iter().enumerate() {
        if let Some(ss) = &r.subsystem {
            if !ss.is_empty() {
                by_subsys.entry(ss.as_str()).or_default().push(i);
            }
        }
    }
    if by_subsys.is_empty() {
        return Ok(());
    }
    w.write_event(Event::Start(BytesStart::new(TAG_GROUPS_GROUPS)))?;
    for (idx, (subsys, rxn_idxs)) in by_subsys.iter().enumerate() {
        let gid = format!("g_{}", idx + 1);
        let mut gel = BytesStart::new(TAG_GROUPS_GROUP);
        gel.push_attribute((ATTR_GROUPS_ID, gid.as_str()));
        gel.push_attribute((ATTR_GROUPS_NAME, *subsys));
        gel.push_attribute((ATTR_GROUPS_KIND, "partonomy"));
        gel.push_attribute(("sboTerm", "SBO:0000633"));
        w.write_event(Event::Start(gel))?;

        w.write_event(Event::Start(BytesStart::new(TAG_GROUPS_MEMBERS)))?;
        for &ri in rxn_idxs {
            let rid = reaction_id(model.rxns[ri].id.as_str());
            let mut m = BytesStart::new(TAG_GROUPS_MEMBER);
            m.push_attribute((ATTR_GROUPS_IDREF, rid.as_str()));
            w.write_event(Event::Empty(m))?;
        }
        w.write_event(Event::End(BytesEnd::new(TAG_GROUPS_MEMBERS)))?;

        w.write_event(Event::End(BytesEnd::new(TAG_GROUPS_GROUP)))?;
    }
    w.write_event(Event::End(BytesEnd::new(TAG_GROUPS_GROUPS)))?;
    Ok(())
}

// -- Helpers --

/// Canonicalize a metabolite id + compartment into a COBRApy-compatible
/// `M_<cpd>_<comp>` species id. Idempotent in three ways:
///
/// - Keeps an already-`M_`-prefixed id verbatim.
/// - Does not append `_<comp>` when the cpd id already ends with it
///   (gapseq's draft stores metabolites as `cpd00001_c0`, so the
///   compartment is baked into the id itself).
/// - Trims trailing separator if the cpd ends with `_`.
fn species_id(cpd: &str, comp: &str) -> String {
    if cpd.starts_with("M_") {
        return cpd.to_string();
    }
    let suffix = format!("_{comp}");
    if cpd.ends_with(&suffix) {
        format!("M_{cpd}")
    } else {
        format!("M_{cpd}{suffix}")
    }
}

fn reaction_id(rxn: &str) -> String {
    if rxn.starts_with("R_") {
        rxn.to_string()
    } else {
        format!("R_{rxn}")
    }
}

/// Sanitize a string into a valid SBML SId: `[A-Za-z_][A-Za-z0-9_]*`.
/// Non-conforming characters (`-`, `.`, `:`, space, …) are replaced
/// with `_`. Leading digits get an underscore prefix.
fn sanitize_sid(s: &str) -> String {
    if s.is_empty() {
        return "_".to_string();
    }
    let mut out = String::with_capacity(s.len() + 1);
    let mut chars = s.chars();
    let first = chars.next().unwrap();
    if first.is_ascii_alphabetic() || first == '_' {
        out.push(first);
    } else if first.is_ascii_digit() {
        out.push('_');
        out.push(first);
    } else {
        out.push('_');
    }
    for c in chars {
        if c.is_ascii_alphanumeric() || c == '_' {
            out.push(c);
        } else {
            out.push('_');
        }
    }
    out
}

fn gene_id(gene: &str) -> String {
    if gene.starts_with("G_") {
        gene.to_string()
    } else {
        format!("G_{gene}")
    }
}

/// Format an `f64` without a trailing `.0` if it's an integer value.
/// Mirrors cobrar / COBRApy output style.
fn fmt_f64(v: f64) -> String {
    if v.fract() == 0.0 && v.is_finite() && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        format!("{v}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gapsmith_core::{CompartmentId, Metabolite, Reaction, StoichMatrix};

    fn toy() -> Model {
        let mut m = Model::new("toy_model");
        m.annot.name = Some("Test model".into());
        m.annot.gapsmith_version = Some("0.1.0".into());
        m.mets.push(Metabolite::new("cpd00001", "H2O", CompartmentId::CYTOSOL));
        m.mets.push(Metabolite::new("cpd00002", "ATP", CompartmentId::CYTOSOL));
        m.mets.push({
            let mut x = Metabolite::new("cpd00007", "O2", CompartmentId::EXTRACELLULAR);
            x.formula = Some("O2".into());
            x.charge = 0;
            x
        });
        let mut r1 = Reaction::new("rxn00001", "ATPase", 0.0, 1000.0);
        r1.obj_coef = 1.0;
        r1.subsystem = Some("Central metabolism".into());
        r1.gpr_raw = Some("(b0001 and b0002) or b0003".into());
        m.rxns.push(r1);
        let mut ex = Reaction::new("EX_cpd00007_e0", "O2 exchange", -1000.0, 1000.0);
        ex.is_exchange = true;
        m.rxns.push(ex);
        m.s = StoichMatrix::from_triplets(
            3,
            2,
            vec![(0, 0, 1.0), (1, 0, -1.0), (2, 1, -1.0)],
        );
        m
    }

    #[test]
    fn toy_writes_valid_sbml() {
        let m = toy();
        let mut buf = Vec::new();
        write_to(&m, &mut buf, &WriteOptions::default()).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("<?xml"), "output: {s}");
        assert!(s.contains("<sbml"));
        assert!(s.contains("xmlns:fbc=\""));
        assert!(s.contains("xmlns:groups=\""));
        assert!(s.contains("<model id=\"toy_model\""));
        assert!(s.contains("fbc:strict=\"true\""));
        assert!(s.contains("<species id=\"M_cpd00001_c0\""));
        assert!(s.contains("<reaction id=\"R_rxn00001\""));
        assert!(s.contains("fbc:lowerFluxBound=\"cobra_0_bound\""));
        assert!(s.contains("fbc:upperFluxBound=\"cobra_default_ub\""));
        assert!(s.contains("<fbc:geneProductAssociation>"));
        assert!(s.contains("<fbc:or>"));
        assert!(s.contains("<fbc:and>"));
        assert!(s.contains("G_b0001"));
        assert!(s.contains("fbc:fluxObjective fbc:reaction=\"R_rxn00001\""));
        assert!(s.contains("groups:group "));
        assert!(s.contains("Central metabolism"));
    }

    #[test]
    fn parse_back_with_quick_xml() {
        use quick_xml::events::Event as E;
        use quick_xml::reader::Reader;
        let m = toy();
        let mut buf = Vec::new();
        write_to(&m, &mut buf, &WriteOptions::default()).unwrap();
        let mut rdr = Reader::from_reader(std::io::Cursor::new(buf));
        let mut stack: Vec<String> = Vec::new();
        let mut species = 0usize;
        let mut rxns = 0usize;
        let mut compartments = 0usize;
        let mut bufx = Vec::new();
        loop {
            match rdr.read_event_into(&mut bufx).unwrap() {
                E::Start(e) => {
                    let name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                    stack.push(name);
                }
                E::Empty(e) => {
                    let name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                    if name == "species" {
                        species += 1;
                    } else if name == "compartment" {
                        compartments += 1;
                    }
                }
                E::End(_) => {
                    if let Some(n) = stack.pop() {
                        if n == "reaction" {
                            rxns += 1;
                        }
                    }
                }
                E::Eof => break,
                _ => {}
            }
            bufx.clear();
        }
        assert_eq!(species, 3);
        assert_eq!(rxns, 2);
        // Three compartments always emitted (c0, e0, p0).
        assert_eq!(compartments, 3);
    }
}
