//! `gapsmith find` — pathway/reaction detection pipeline.
//!
//! Thin CLI wrapper around [`gapsmith_find::run_find`]. Mirrors the short
//! flags of `src/gapseq_find.sh` where reasonable.

use clap::{Parser, ValueEnum};
use gapsmith_align::{
    AlignOpts, Aligner, BlastpAligner, DiamondAligner, Mmseqs2Aligner, PrecomputedTsvAligner,
};
use gapsmith_db::{exception, pathway::PwySource, ComplexSubunitTable, PathwayTable};
use gapsmith_find::{
    dbhit::DbhitIndex, pathways::MatchMode, run_find, taxonomy::valid_tax_ids_for,
    write_pathways_tbl, write_reactions_tbl, FindOptions, SeqfileOptions,
};
use gapsmith_io::{resolve_data_dir, resolve_seq_dir};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
pub struct Args {
    /// Path to the input genome's protein FASTA.
    pub genome: PathBuf,

    /// Keyword filter on pathway id / name / hierarchy. `all` matches
    /// every active pathway.
    #[arg(long, short = 'p', default_value = "all")]
    pub pathways: String,

    /// Pathway database(s) to query. Comma-separated subset of
    /// `metacyc,kegg,seed,custom,all`. Default matches gapseq:
    /// `metacyc,custom` (custom shadows metacyc on ID collision).
    #[arg(long, short = 'l', default_value = "metacyc,custom")]
    pub pathway_db: String,

    /// Taxonomy — maps to the seq-dir subfolder (`Bacteria`, `Archaea`).
    #[arg(long, short = 't', default_value = "Bacteria")]
    pub taxonomy: String,

    /// Aligner backend.
    #[arg(long, short = 'A', value_enum, default_value_t = AlignerArg::Diamond)]
    pub aligner: AlignerArg,

    /// Optional pre-computed alignment TSV (skips the aligner run).
    #[arg(long, short = 'P')]
    pub precomputed: Option<PathBuf>,

    /// Bitscore cutoff.
    #[arg(long, short = 'b', default_value_t = 200.0)]
    pub bitcutoff: f32,

    /// Identity cutoff (0–100).
    #[arg(long, short = 'i', default_value_t = 0.0)]
    pub identcutoff: f32,

    /// Identity cutoff for exception-table ECs (false-friend enzymes).
    #[arg(long, default_value_t = 70.0)]
    pub ident_exception: f32,

    /// Coverage cutoff (0–100).
    #[arg(long, short = 'c', default_value_t = 75)]
    pub coverage: u32,

    /// Main pathway-completeness cutoff (percent).
    #[arg(long, short = 'a', default_value_t = 80.0)]
    pub completeness_main: f32,

    /// Lowered pathway-completeness cutoff applied when all key reactions
    /// are present (percent).
    #[arg(long, short = 'k', default_value_t = 66.0)]
    pub completeness_hints: f32,

    /// Strict mode — disable key-reaction / completeness heuristics.
    #[arg(long, short = 's')]
    pub strict: bool,

    /// Include superpathways. Default (`false`) matches gapseq's
    /// `noSuperpathways=true`.
    #[arg(long, short = 'n')]
    pub include_superpathways: bool,

    /// Output directory. Creates `<stem>-Reactions.tbl` and
    /// `<stem>-Pathways.tbl` inside.
    #[arg(long, short = 'o', default_value = ".")]
    pub out_dir: PathBuf,

    /// Optional suffix used in output file names (`<stem>-<suffix>-Reactions.tbl`).
    #[arg(long, short = 'u')]
    pub suffix: Option<String>,

    /// Thread count for the aligner.
    #[arg(long, short = 'K')]
    pub threads: Option<usize>,
}


#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AlignerArg {
    Blastp,
    Diamond,
    Mmseqs2,
    Precomputed,
}

pub fn run(
    args: Args,
    data_dir_override: Option<&Path>,
    seq_dir_override: Option<&Path>,
) -> anyhow::Result<()> {
    let data_dir = resolve_data_dir(data_dir_override)?;
    let seq_dir = resolve_seq_dir(seq_dir_override, &data_dir)?;
    let tax_root = seq_dir.join(&args.taxonomy);

    tracing::info!(
        data_dir = %data_dir.display(),
        seq_dir = %seq_dir.display(),
        taxonomy = %args.taxonomy,
        pathway_db = ?args.pathway_db,
        keyword = %args.pathways,
        "gapsmith find starting"
    );

    // 1. Load the requested pathway database(s). We follow gapseq's merge
    //    semantics exactly (`gapseq_find.sh:520-532`): concatenate every
    //    selected table, then for any `id` that appears more than once
    //    keep only the `custom` row.
    let dbs_lc = args.pathway_db.to_ascii_lowercase();
    let include = |name: &str| dbs_lc == "all" || dbs_lc.split(',').any(|s| s.trim() == name);
    let mut table = PathwayTable::default();
    table.source = Some(PwySource::MetaCyc);
    let mut custom_rows: Vec<gapsmith_db::PathwayRow> = Vec::new();
    if include("metacyc") {
        let t = PathwayTable::load(&data_dir.join("meta_pwy.tbl"), PwySource::MetaCyc)?;
        table.rows.extend(t.rows);
    }
    if include("kegg") {
        let t = PathwayTable::load(&data_dir.join("kegg_pwy.tbl"), PwySource::Kegg)?;
        table.rows.extend(t.rows);
    }
    if include("seed") {
        let t = PathwayTable::load(&data_dir.join("seed_pwy.tbl"), PwySource::Seed)?;
        table.rows.extend(t.rows);
    }
    if include("custom") {
        let t = PathwayTable::load(&data_dir.join("custom_pwy.tbl"), PwySource::Custom)?;
        custom_rows = t.rows;
    }
    // Merge-with-custom-wins: if an id appears in both, drop every earlier
    // copy and splice the custom row in its place.
    if !custom_rows.is_empty() {
        let custom_ids: std::collections::HashSet<String> =
            custom_rows.iter().map(|r| r.id.clone()).collect();
        table.rows.retain(|r| !custom_ids.contains(&r.id));
        table.rows.extend(custom_rows);
    }
    tracing::info!(pathways = table.rows.len(), "pathway table loaded");

    // 2. Load exception table + complex subunit dictionary + dbhit index.
    let exception_rows = exception::load(data_dir.join("exception.tbl"))?;
    let subunit_dict = ComplexSubunitTable::load(data_dir.join("complex_subunit_dict.tsv"))?;
    let dbhit_index = DbhitIndex::load(&data_dir)?;
    let valid_tax_ids = valid_tax_ids_for(&data_dir.join("taxonomy.tbl"), &args.taxonomy)
        .unwrap_or_default();

    // 3. Prepare aligner.
    let aligner: Box<dyn Aligner> = match args.aligner {
        AlignerArg::Blastp => Box::new(BlastpAligner::new()),
        AlignerArg::Diamond => Box::new(DiamondAligner::new()),
        AlignerArg::Mmseqs2 => Box::new(Mmseqs2Aligner::new()),
        AlignerArg::Precomputed => {
            let p = args.precomputed.clone().ok_or_else(|| {
                anyhow::anyhow!("--precomputed required when --aligner precomputed")
            })?;
            Box::new(PrecomputedTsvAligner::new_percentage(p))
        }
    };

    let align_opts = AlignOpts {
        threads: args.threads.unwrap_or_else(|| {
            std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
        }),
        coverage_pct: args.coverage,
        evalue: None,
        extra_args: Vec::new(),
        quiet: false,
    };

    // 4. Run find pipeline.
    // user/ lookups must always hit the built-in gapseq data tree, not the
    // user's `--seq-dir` override — matches `prepare_batch_alignments.R:217`.
    let builtin_tax_root = data_dir.join("seq").join(&args.taxonomy);
    let seq_opts = SeqfileOptions {
        tax_root: tax_root.clone(),
        user_root: if builtin_tax_root.exists() { builtin_tax_root } else { tax_root.clone() },
        seq_src: gapsmith_find::seqfile::SeqSrc::PreferRev,
    };
    // For predefined shorthands (`amino`, `nucl`, …) gapseq's shell uses
    // `pwyKeyCol=hierarchy`; for `all` and the `custom` group it uses
    // hierarchy too; everything else falls through to the regex path. We
    // mirror that dispatch here.
    let (match_mode, exclude_superpathways) = match args.pathways.as_str() {
        "all" | "amino" | "nucl" | "cofactor" | "carbo" | "carbo-deg" | "polyamine"
        | "fatty" | "energy" | "terpenoid" | "degradation" | "kegg" => {
            (MatchMode::Hierarchy, !args.include_superpathways)
        }
        "core" | "min" => (MatchMode::Regex, false), // gapseq sets noSuperpathways=false for these
        _ => (MatchMode::Regex, !args.include_superpathways),
    };
    let find_opts = FindOptions {
        keyword: &args.pathways,
        match_mode,
        exclude_superpathways,
        only_active: true,
        bitcutoff: args.bitcutoff,
        identcutoff: args.identcutoff,
        ident_exception: args.ident_exception,
        coverage_pct: args.coverage,
        vague_cutoff: 0.3,
        completeness_hint_off: args.completeness_main / 100.0,
        completeness_hint_on: args.completeness_hints / 100.0,
        strict_candidates: args.strict,
        subunit_cutoff: 0.5,
        valid_tax_ids: &valid_tax_ids,
    };

    let workdir = tempfile::tempdir()?;
    let report = run_find(
        &table,
        &exception_rows,
        &subunit_dict,
        &dbhit_index,
        &seq_opts,
        &args.genome,
        aligner.as_ref(),
        &align_opts,
        &find_opts,
        workdir.path(),
    )?;

    // 5. Stamp pathway names by joining with pathway table.
    let mut by_id: std::collections::HashMap<&str, &str> =
        table.rows.iter().map(|r| (r.id.as_str(), r.name.as_str())).collect();
    let mut pathways = report.pathways;
    for p in &mut pathways {
        if let Some(n) = by_id.remove(p.id.as_str()) {
            p.name = n.to_string();
        }
    }

    // 6. Write output tables.
    std::fs::create_dir_all(&args.out_dir)?;
    let stem = args
        .genome
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("genome")
        .trim_end_matches(".faa")
        .trim_end_matches(".fasta")
        .to_string();
    let tag = args.suffix.unwrap_or_else(|| args.pathways.clone());
    let rxn_path = args.out_dir.join(format!("{stem}-{tag}-Reactions.tbl"));
    let pwy_path = args.out_dir.join(format!("{stem}-{tag}-Pathways.tbl"));
    write_reactions_tbl(&report.reactions, &rxn_path)?;
    write_pathways_tbl(&pathways, &pwy_path)?;

    let n_found_pwy = pathways.iter().filter(|p| p.prediction).count();
    let n_rxn_good = report
        .reactions
        .iter()
        .filter(|r| r.status == gapsmith_find::HitStatus::GoodBlast)
        .count();
    eprintln!(
        "wrote {} reactions to {}",
        report.reactions.len(),
        rxn_path.display()
    );
    eprintln!("wrote {} pathways to {}", pathways.len(), pwy_path.display());
    eprintln!("  good-blast reactions: {n_rxn_good}");
    eprintln!("  predicted-present pathways: {n_found_pwy}");
    Ok(())
}
