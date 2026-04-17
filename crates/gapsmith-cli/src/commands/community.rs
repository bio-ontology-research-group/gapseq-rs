//! `gapsmith community` — community-level optimisation for metagenomes or
//! multi-MAG reconstructions.
//!
//! Two modes, picked via subcommand:
//!
//! - **`per-mag`** (scalable): each MAG keeps its own model and its own
//!   FBA; we only build a shared medium (union of per-MAG media with the
//!   max flux per compound) and report weighted community growth.
//!   Complexity is linear in `N` with no solver-side changes — the mode
//!   to reach for at 100+ MAGs.
//!
//! - **`cfba`** (accurate, opt-in): compose every draft into one
//!   block-diagonal community model with a shared `_e0` pool and a
//!   weighted-sum biomass objective. Optionally enforce balanced growth
//!   (`v(bio_k) == v(bio_community)` for all k). The single community LP
//!   then captures obligate cross-feeding. Solver cost grows roughly
//!   `N × (rxns + mets)` so this is the mode for small synthetic
//!   communities or curated co-cultures, not 1000-MAG metagenomes.
//!
//! Both modes accept either an explicit list of draft CBORs or a gspa run
//! directory (for the latter, abundance is taken from `genomes.tsv`'s
//! optional third column).

use clap::{Parser, Subcommand};
use gapsmith_align::GspaManifest;
use gapsmith_core::Model;
use gapsmith_fill::{
    apply_medium, compose_models, add_community_biomass, fba, per_mag_weights, read_medium,
    union_medium, weighted_growth, FbaOptions, MediumEntry, Organism, SolveStatus,
};
use gapsmith_io::{read_model_cbor, write_model_cbor};
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
pub struct Args {
    #[command(subcommand)]
    pub mode: Mode,
}

#[derive(Debug, Subcommand)]
pub enum Mode {
    /// Per-MAG FBA with a shared (union) medium. Scales to 1000+ MAGs.
    #[command(name = "per-mag")]
    PerMag(PerMagArgs),
    /// Full community FBA — one LP across a composed community model.
    Cfba(CfbaArgs),
}

#[derive(Debug, Parser)]
pub struct PerMagArgs {
    /// Directory of per-MAG draft CBORs (`<genome_id>-draft.gmod.cbor`
    /// or `<genome_id>-filled.gmod.cbor`). Mutually exclusive with
    /// `--drafts-list` and `--gspa-run`.
    #[arg(long, conflicts_with_all = ["drafts_list", "gspa_run"])]
    pub drafts_dir: Option<PathBuf>,

    /// TSV of `<genome_id>\t<cbor_path>[\t<medium_csv>]`. The 3rd column
    /// is optional — when absent we look for `<out_dir>/<genome_id>-medium.csv`
    /// next to each draft (the name `gapsmith medium` / `doall` writes).
    #[arg(long, conflicts_with_all = ["drafts_dir", "gspa_run"])]
    pub drafts_list: Option<PathBuf>,

    /// Gspa manifest — pulls the genome list (and abundances) from there.
    /// Per-MAG drafts are expected under `<drafts_root>/<genome_id>/` with
    /// the standard `doall` naming.
    #[arg(long, conflicts_with_all = ["drafts_dir", "drafts_list"])]
    pub gspa_run: Option<PathBuf>,

    /// Root directory where per-genome doall outputs live when resolving
    /// paths from a gspa manifest. Defaults to the current directory.
    #[arg(long, default_value = ".")]
    pub drafts_root: PathBuf,

    /// Output directory. Writes `community-medium.csv` + `per-mag-growth.tsv`.
    #[arg(long, short = 'o', default_value = ".")]
    pub out_dir: PathBuf,

    /// Skip inferring a shared medium — read it from `--medium` instead.
    #[arg(long, short = 'm')]
    pub medium: Option<PathBuf>,
}

#[derive(Debug, Parser)]
pub struct CfbaArgs {
    /// Directory of per-MAG draft CBORs, or a TSV list, or a gspa-run.
    /// Exactly one of the three must be supplied; same semantics as `per-mag`.
    #[arg(long, conflicts_with_all = ["drafts_list", "gspa_run"])]
    pub drafts_dir: Option<PathBuf>,

    #[arg(long, conflicts_with_all = ["drafts_dir", "gspa_run"])]
    pub drafts_list: Option<PathBuf>,

    #[arg(long, conflicts_with_all = ["drafts_dir", "drafts_list"])]
    pub gspa_run: Option<PathBuf>,

    #[arg(long, default_value = ".")]
    pub drafts_root: PathBuf,

    /// Abundance TSV: `<genome_id>\t<abundance>` (header optional). Applied
    /// after gspa-run abundance resolution — this file overrides. Missing
    /// genomes fall back to uniform weights.
    #[arg(long)]
    pub abundance: Option<PathBuf>,

    /// Reaction id used as each organism's biomass (default `bio1`).
    #[arg(long, default_value = "bio1")]
    pub biomass_rxn: String,

    /// Enforce balanced growth: every organism's biomass flux equals
    /// the community biomass flux. Matches classical cFBA semantics.
    #[arg(long)]
    pub balanced_growth: bool,

    /// Shared medium CSV. When omitted the composed model runs with each
    /// organism's own `EX_*_e0` bounds as merged by `compose_models`.
    #[arg(long, short = 'm')]
    pub medium: Option<PathBuf>,

    /// Output directory. Writes `community.gmod.cbor` + `community-fba.tsv`.
    #[arg(long, short = 'o', default_value = ".")]
    pub out_dir: PathBuf,
}

// --------------------------------------------------------------------------
// Entry points
// --------------------------------------------------------------------------

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.mode {
        Mode::PerMag(a) => run_per_mag(a),
        Mode::Cfba(a) => run_cfba(a),
    }
}

fn run_per_mag(args: PerMagArgs) -> anyhow::Result<()> {
    let entries = resolve_entries(
        args.drafts_dir.as_deref(),
        args.drafts_list.as_deref(),
        args.gspa_run.as_deref(),
        &args.drafts_root,
    )?;
    if entries.is_empty() {
        anyhow::bail!("no drafts resolved");
    }
    std::fs::create_dir_all(&args.out_dir)?;

    // 1. Resolve per-MAG media.
    let mut per_mag_media: Vec<(String, Vec<MediumEntry>)> = Vec::with_capacity(entries.len());
    for e in &entries {
        let med_path = match &e.medium {
            Some(p) => p.clone(),
            None => guess_medium_path(&e.draft_path, &e.id),
        };
        match read_medium(&med_path) {
            Ok(m) => per_mag_media.push((e.id.clone(), m)),
            Err(err) => {
                tracing::warn!(genome = %e.id, path = %med_path.display(), err = %err, "no medium found — skipping");
            }
        }
    }

    // 2. Compute shared medium.
    let shared = if let Some(path) = &args.medium {
        read_medium(path)?
    } else {
        let slices: Vec<&[MediumEntry]> = per_mag_media.iter().map(|(_, m)| m.as_slice()).collect();
        union_medium(&slices)
    };
    let medium_out = args.out_dir.join("community-medium.csv");
    write_medium_csv(&medium_out, &shared)?;
    eprintln!(
        "community medium: {} compounds → {}",
        shared.len(),
        medium_out.display()
    );

    // 3. Per-MAG FBA under the shared medium.
    let growth_out = args.out_dir.join("per-mag-growth.tsv");
    let mut w = std::io::BufWriter::new(std::fs::File::create(&growth_out)?);
    writeln!(w, "genome_id\tabundance\tgrowth_rate\tstatus")?;
    let mut per_organism: Vec<(String, f64)> = Vec::with_capacity(entries.len());
    let mut abundances: Vec<(String, Option<f64>)> = Vec::with_capacity(entries.len());
    for e in &entries {
        abundances.push((e.id.clone(), e.abundance));
        let mut model: Model = read_model_cbor(&e.draft_path)?;
        apply_medium(&mut model, &shared, 1.0, 1000.0);
        let sol = fba(&model, &FbaOptions::default())?;
        let g = if matches!(sol.status, SolveStatus::Optimal) { sol.objective } else { 0.0 };
        let a_str = match e.abundance {
            Some(v) => format!("{v}"),
            None => String::from("NA"),
        };
        writeln!(
            w,
            "{}\t{}\t{:.6}\t{:?}",
            e.id, a_str, g, sol.status
        )?;
        per_organism.push((e.id.clone(), g));
    }
    w.flush()?;

    // 4. Weighted summary.
    let weights = per_mag_weights(&abundances);
    let (mean_growth, _) = weighted_growth(&per_organism, &weights);
    eprintln!("weighted community growth: {mean_growth:.6}");
    eprintln!("per-MAG growth → {}", growth_out.display());

    Ok(())
}

fn run_cfba(args: CfbaArgs) -> anyhow::Result<()> {
    let mut entries = resolve_entries(
        args.drafts_dir.as_deref(),
        args.drafts_list.as_deref(),
        args.gspa_run.as_deref(),
        &args.drafts_root,
    )?;
    if entries.is_empty() {
        anyhow::bail!("no drafts resolved");
    }
    std::fs::create_dir_all(&args.out_dir)?;

    // Overlay abundance TSV (overrides anything from gspa-run).
    if let Some(p) = &args.abundance {
        let overrides = read_abundance_tsv(p)?;
        for e in &mut entries {
            if let Some(a) = overrides.get(&e.id) {
                e.abundance = Some(*a);
            }
        }
    }

    // Build normalised weights for each organism.
    let abundances: Vec<(String, Option<f64>)> =
        entries.iter().map(|e| (e.id.clone(), e.abundance)).collect();
    let weights = per_mag_weights(&abundances);

    // Load drafts into Organism records.
    let mut organisms: Vec<Organism> = Vec::with_capacity(entries.len());
    for (e, (_, w)) in entries.iter().zip(weights.iter()) {
        let model: Model = read_model_cbor(&e.draft_path)?;
        organisms.push(Organism {
            id: e.id.clone(),
            model,
            biomass_rxn: args.biomass_rxn.clone(),
            weight: *w,
        });
    }

    eprintln!(
        "cFBA: composing {} organisms (balanced_growth = {})",
        organisms.len(),
        args.balanced_growth
    );
    let mut comm = compose_models(&organisms)?;
    add_community_biomass(&mut comm, args.balanced_growth)?;

    // Optional shared medium.
    if let Some(p) = &args.medium {
        let m = read_medium(p)?;
        apply_medium(&mut comm.model, &m, 1.0, 1000.0);
    }

    let sol = fba(&comm.model, &FbaOptions::default())?;
    eprintln!(
        "community growth: {:.6} ({:?})",
        sol.objective, sol.status
    );

    // Output: community model + per-organism growth rates.
    let model_out = args.out_dir.join("community.gmod.cbor");
    write_model_cbor(&comm.model, &model_out)?;
    let report_out = args.out_dir.join("community-fba.tsv");
    let mut w = std::io::BufWriter::new(std::fs::File::create(&report_out)?);
    writeln!(w, "genome_id\tweight\tbiomass_flux")?;
    let rxn_idx = comm.model.rxn_index();
    for (org, bio_id) in &comm.biomass_rxn_ids {
        let i = rxn_idx[&gapsmith_core::RxnId::new(bio_id.as_str())];
        let weight = comm
            .weights
            .iter()
            .find(|(g, _)| g == org)
            .map(|(_, w)| *w)
            .unwrap_or(0.0);
        writeln!(w, "{}\t{:.6}\t{:.6}", org, weight, sol.fluxes[i])?;
    }
    writeln!(w, "#community_biomass\t{:.6}\t{:.6}", 1.0, sol.objective)?;
    w.flush()?;

    eprintln!(
        "  community model → {}\n  per-organism table → {}",
        model_out.display(),
        report_out.display()
    );
    Ok(())
}

// --------------------------------------------------------------------------
// Entry resolution shared by both modes.
// --------------------------------------------------------------------------

struct DraftEntry {
    id: String,
    draft_path: PathBuf,
    abundance: Option<f64>,
    medium: Option<PathBuf>,
}

fn resolve_entries(
    dir: Option<&Path>,
    list: Option<&Path>,
    gspa_run: Option<&Path>,
    drafts_root: &Path,
) -> anyhow::Result<Vec<DraftEntry>> {
    if let Some(d) = dir {
        let mut out = Vec::new();
        for entry in std::fs::read_dir(d)? {
            let e = entry?;
            let p = e.path();
            if !p.is_file() {
                continue;
            }
            let name = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
            if !name.ends_with(".gmod.cbor") {
                continue;
            }
            // `<id>-draft.gmod.cbor` or `<id>-filled.gmod.cbor` → id = stem up to the first `-`.
            let stem = name.trim_end_matches(".gmod.cbor");
            let id = stem
                .trim_end_matches("-filled")
                .trim_end_matches("-draft")
                .to_string();
            out.push(DraftEntry { id, draft_path: p, abundance: None, medium: None });
        }
        out.sort_by(|a, b| a.id.cmp(&b.id));
        return Ok(out);
    }
    if let Some(l) = list {
        let text = std::fs::read_to_string(l)?;
        let mut out = Vec::new();
        for (ln, line) in text.lines().enumerate() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 2 {
                anyhow::bail!("drafts-list line {} needs `<id>\\t<cbor>[\\t<medium>]`: `{line}`", ln + 1);
            }
            out.push(DraftEntry {
                id: cols[0].trim().into(),
                draft_path: PathBuf::from(cols[1].trim()),
                abundance: None,
                medium: cols.get(2).and_then(|s| {
                    let t = s.trim();
                    if t.is_empty() { None } else { Some(PathBuf::from(t)) }
                }),
            });
        }
        return Ok(out);
    }
    if let Some(gr) = gspa_run {
        let m = GspaManifest::load(gr)?;
        let mut out = Vec::with_capacity(m.genomes.len());
        for g in m.genomes {
            let draft = drafts_root.join(&g.id).join(format!("{}-draft.gmod.cbor", g.id));
            let filled = drafts_root.join(&g.id).join(format!("{}-filled.gmod.cbor", g.id));
            let path = if filled.exists() {
                filled
            } else {
                draft
            };
            out.push(DraftEntry {
                id: g.id,
                draft_path: path,
                abundance: g.abundance,
                medium: None,
            });
        }
        return Ok(out);
    }
    anyhow::bail!("one of --drafts-dir / --drafts-list / --gspa-run must be provided")
}

fn guess_medium_path(draft_path: &Path, id: &str) -> PathBuf {
    draft_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join(format!("{id}-medium.csv"))
}

fn read_abundance_tsv(path: &Path) -> anyhow::Result<std::collections::HashMap<String, f64>> {
    let text = std::fs::read_to_string(path)?;
    let mut out = std::collections::HashMap::new();
    for (ln, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let (id, val) = line.split_once('\t').ok_or_else(|| {
            anyhow::anyhow!(
                "abundance TSV line {} needs `<id>\\t<value>`: `{line}`",
                ln + 1
            )
        })?;
        if id.trim() == "genome_id" || val.trim() == "abundance" {
            continue;
        }
        let v: f64 = val.trim().parse().map_err(|_| {
            anyhow::anyhow!("abundance TSV line {} value `{val}` is not a float", ln + 1)
        })?;
        out.insert(id.trim().to_string(), v);
    }
    Ok(out)
}

fn write_medium_csv(path: &Path, medium: &[MediumEntry]) -> anyhow::Result<()> {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "compounds,name,maxFlux")?;
    for e in medium {
        writeln!(
            f,
            "{},{},{}",
            e.compound,
            // Names may contain commas; replace cheaply to avoid quoting.
            e.name.replace(',', ";"),
            e.max_flux
        )?;
    }
    f.flush()?;
    Ok(())
}
