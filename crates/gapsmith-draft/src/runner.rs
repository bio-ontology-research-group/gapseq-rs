//! End-to-end driver for `gapsmith draft`. Loads every input, selects
//! reactions, builds the `Model`, emits diagnostic info.

use crate::biomass::{parse_biomass_json, BiomassSpec};
use crate::builder::{add_seed_reaction, build_model, BuilderOptions};
use crate::candidate::{build_candidates, CandidateOptions, CandidateTable};
use crate::exchanges::{add_missing_diffusion, load_diffusion_rxns};
use crate::reactions_tbl::{read_reactions_tbl, read_transporter_tbl};
use crate::stoich_hash::rxn_stoich_hash;
use gapsmith_core::Model;
use gapsmith_db::{load_seed_metabolites, load_seed_reactions, SeedRxnRow};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

#[derive(Debug, thiserror::Error)]
pub enum DraftError {
    #[error("i/o error on `{path}`: {source}")]
    Io { path: PathBuf, #[source] source: std::io::Error },
    #[error("DB error: {0}")]
    Db(#[from] gapsmith_db::DbError),
    #[error("biomass error: {0}")]
    Biomass(#[from] crate::biomass::BiomassError),
    #[error("stoich hash error: {0}")]
    Stoich(String),
    #[error("bad config: {0}")]
    BadConfig(String),
}

#[derive(Debug, Clone)]
pub struct DraftOptions {
    pub model_id: String,
    /// Biomass choice: `pos`, `neg`, `archaea`, or path to a custom
    /// JSON template. `auto` resolves to gram from the header of the
    /// Reactions.tbl.
    pub biomass: String,
    pub high_evi_rxn_bs: f32,
    pub min_bs_for_core: f32,
    pub gapsmith_version: Option<String>,
    pub seqdb_version: Option<String>,
}

impl Default for DraftOptions {
    fn default() -> Self {
        Self {
            model_id: "draft".into(),
            biomass: "auto".into(),
            high_evi_rxn_bs: 200.0,
            min_bs_for_core: 50.0,
            gapsmith_version: None,
            seqdb_version: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DraftReport {
    pub model: Model,
    pub candidates: CandidateTable,
    pub selected_seed_ids: Vec<String>,
    pub biomass_spec: Option<BiomassSpec>,
}

pub fn run(
    reactions_tbl: &Path,
    transporter_tbl: &Path,
    data_dir: &Path,
    opts: &DraftOptions,
) -> Result<DraftReport, DraftError> {
    // -- 1. Load find/find-transport output.
    let reactions = read_reactions_tbl(reactions_tbl).map_err(|e| DraftError::Io {
        path: reactions_tbl.to_path_buf(),
        source: e,
    })?;
    let transporters = read_transporter_tbl(transporter_tbl).map_err(|e| DraftError::Io {
        path: transporter_tbl.to_path_buf(),
        source: e,
    })?;
    tracing::info!(
        reactions = reactions.len(),
        transporters = transporters.len(),
        "find-tables loaded"
    );

    // -- 2. Extract per-reaction gene assignments + subunit complex info.
    let mut gene_assignments: HashMap<String, Vec<(String, String)>> = HashMap::new();
    for r in &reactions {
        if r.stitle.is_empty() || r.bitscore.is_none() {
            continue;
        }
        if r.bitscore.unwrap_or(0.0) < opts.high_evi_rxn_bs {
            continue;
        }
        for seed in r.dbhit.split_whitespace().filter(|s| !s.is_empty()) {
            let gene = first_token(&r.stitle);
            gene_assignments
                .entry(seed.to_string())
                .or_default()
                .push((r.complex.clone(), gene.to_string()));
        }
    }
    // Dedup (gene, complex) per seed.
    for v in gene_assignments.values_mut() {
        let mut seen = HashSet::new();
        v.retain(|(c, g)| seen.insert((c.clone(), g.clone())));
    }

    // -- 3. Candidate reaction table.
    let cand_opts = CandidateOptions {
        high_evi_rxn_bs: opts.high_evi_rxn_bs,
        min_bs_for_core: opts.min_bs_for_core,
    };
    let candidates = build_candidates(&reactions, &transporters, &cand_opts);

    // Selection: keep reactions with bitscore above threshold OR predicted
    // pathway present (core reactions).
    let selected_seed_ids: Vec<String> = candidates
        .rows
        .iter()
        .filter(|c| {
            let bs = c.best_bitscore.unwrap_or(0.0);
            bs >= opts.high_evi_rxn_bs
                || matches!(
                    c.pathway_status.as_str(),
                    "full" | "threshold" | "keyenzyme"
                )
        })
        .filter(|c| {
            // Drop incomplete complexes unless pathway gives topology support.
            if c.is_complex && c.complex_status.is_none() {
                return matches!(
                    c.pathway_status.as_str(),
                    "full" | "threshold" | "keyenzyme"
                );
            }
            true
        })
        .map(|c| c.seed.clone())
        .collect();
    tracing::info!(selected = selected_seed_ids.len(), "reactions selected for draft");

    // -- 4. Load SEED DB + metabolite map.
    let seed_rxns_all = load_seed_reactions(data_dir.join("seed_reactions_corrected.tsv"))?;
    let seed_cpds = load_seed_metabolites(data_dir.join("seed_metabolites_edited.tsv"))?;

    // Filter to approved/corrected.
    let seed_usable: Vec<&SeedRxnRow> = seed_rxns_all
        .iter()
        .filter(|r| r.gapseq_status.is_usable())
        .collect();

    // Intersect with selected ids.
    let selected_set: HashSet<&str> =
        selected_seed_ids.iter().map(|s| s.as_str()).collect();
    let mut draft_rxns: Vec<&SeedRxnRow> = seed_usable
        .iter()
        .copied()
        .filter(|r| selected_set.contains(r.id.as_str()))
        .collect();

    // -- 5. Dedup by stoichiometric hash.
    let mut by_hash: HashMap<String, &SeedRxnRow> = HashMap::new();
    for r in &draft_rxns {
        let rev = r.reversibility.as_str();
        let key = rxn_stoich_hash(&r.stoichiometry, rev)
            .map_err(DraftError::Stoich)?;
        by_hash.entry(key).or_insert(*r);
    }
    draft_rxns = by_hash.into_values().collect();
    draft_rxns.sort_by(|a, b| a.id.as_str().cmp(b.id.as_str()));
    tracing::info!(after_dedup = draft_rxns.len(), "reactions after stoich dedup");

    // -- 6. Biomass.
    let biomass_kind = resolve_biomass(&opts.biomass, reactions_tbl)?;
    let bm_path = match biomass_kind.as_str() {
        "neg" | "Gram_neg" => data_dir.join("biomass/biomass_Gram_neg.json"),
        "pos" | "Gram_pos" => data_dir.join("biomass/biomass_Gram_pos.json"),
        "archaea" | "Archaea" => data_dir.join("biomass/biomass_archaea.json"),
        other if Path::new(other).exists() => PathBuf::from(other),
        other => {
            return Err(DraftError::BadConfig(format!(
                "unknown biomass `{other}` and not a file path"
            )))
        }
    };
    let biomass_spec = parse_biomass_json(&bm_path, &seed_cpds)?;

    // -- 7. Build Model.
    let tax_domain = biomass_spec.domain.clone();
    let gram = match biomass_kind.as_str() {
        "neg" | "Gram_neg" => Some("neg".to_string()),
        "pos" | "Gram_pos" => Some("pos".to_string()),
        _ => None,
    };
    let bopts = BuilderOptions {
        model_id: opts.model_id.clone(),
        gapsmith_version: opts.gapsmith_version.clone(),
        seqdb_version: opts.seqdb_version.clone(),
        tax_domain: Some(tax_domain.clone()),
        gram,
    };
    let mut model =
        build_model(&bopts, &draft_rxns, Some(&biomass_spec), &gene_assignments);

    // -- 8. Add diffusion + missing exchanges.
    let diff_path = data_dir.join("diffusion_mets.tsv");
    let diffs = load_diffusion_rxns(&diff_path).map_err(|e| DraftError::Io {
        path: diff_path,
        source: e,
    })?;
    add_missing_diffusion(&mut model, &diffs, &seed_rxns_all);

    // A few conditional transporters — gapseq adds these when specific
    // source reactions are present in the model (origin = 5).
    add_conditional_transporters(&mut model, &seed_rxns_all);

    // Flush pending S-matrix updates accumulated during exchange /
    // diffusion / conditional-transporter passes.
    crate::builder::rebuild_s_matrix(&mut model);

    Ok(DraftReport { model, candidates, selected_seed_ids, biomass_spec: Some(biomass_spec) })
}

/// Peek at the Reactions.tbl header line to find `gram=pos|neg` when the
/// user requested `auto`.
fn resolve_biomass(biomass: &str, reactions_tbl: &Path) -> Result<String, DraftError> {
    match biomass.to_ascii_lowercase().as_str() {
        "auto" => {
            let meta = read_gram_from_header(reactions_tbl);
            Ok(meta.unwrap_or_else(|| "neg".into()))
        }
        "pos" | "gram_pos" => Ok("pos".into()),
        "neg" | "gram_neg" => Ok("neg".into()),
        "archaea" | "bacteria" => Ok(biomass.to_ascii_lowercase()),
        other => Ok(other.to_string()),
    }
}

fn read_gram_from_header(path: &Path) -> Option<String> {
    use std::io::BufRead;
    let f = std::fs::File::open(path).ok()?;
    let r = std::io::BufReader::new(f);
    for line in r.lines().take(3).flatten() {
        if let Some(idx) = line.find("gram=") {
            let rest = &line[idx + 5..];
            let end = rest.find(';').unwrap_or(rest.len());
            return Some(rest[..end].trim().to_string());
        }
    }
    None
}

/// Conditional transporter additions — rxn05683 (butyrate), rxn90116
/// (IPA), rxn09182 (PPA), rxn90117 (Phloretate). Each depends on a set
/// of source reactions already being in the model.
fn add_conditional_transporters(
    model: &mut Model,
    seed_rxns: &[gapsmith_db::SeedRxnRow],
) {
    let present: HashSet<String> = model
        .rxns
        .iter()
        .filter_map(|r| r.id.as_str().strip_suffix("_c0").map(|s| s.to_string()))
        .collect();
    let has = |s: &str| present.contains(s);

    let rules: &[(&[&str], &str)] = &[
        (&["rxn90001"], "rxn05683"),
        (
            &["rxn43343", "rxn45361", "rxn00483", "rxn01447"],
            "rxn90116",
        ),
        (
            &["rxn00493", "rxn00997", "rxn07603", "rxn40746"],
            "rxn09182",
        ),
        (
            &["rxn00527", "rxn02393", "rxn46948", "rxn46031"],
            "rxn90117",
        ),
    ];

    for (sources, target) in rules {
        if sources.iter().all(|s| has(s)) && !has(target) {
            if let Some(row) = seed_rxns.iter().find(|r| r.id.as_str() == *target) {
                add_seed_reaction(model, row, Some(5));
            }
        }
    }
}

fn first_token(s: &str) -> &str {
    match s.find(char::is_whitespace) {
        Some(i) => &s[..i],
        None => s,
    }
}
