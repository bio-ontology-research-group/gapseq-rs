//! `gapsmith find-transport` — transporter detection.
//!
//! Port of `src/transporter.sh`. Loads subex + seed metabolite name map,
//! runs the transporter pipeline from `gapsmith-transport`, writes a
//! `<stem>-Transporter.tbl`.

use clap::{Parser, ValueEnum};
use gapsmith_align::{
    AlignOpts, Aligner, BlastpAligner, DiamondAligner, GspaRunAligner, Mmseqs2Aligner,
    PrecomputedTsvAligner,
};
use gapsmith_db::{load_seed_metabolites, subex};
use gapsmith_io::{resolve_data_dir, resolve_seq_dir};
use gapsmith_transport::{run, write_transporter_tbl, TransportOptions};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
pub struct Args {
    pub genome: PathBuf,

    #[arg(long, short = 'A', value_enum, default_value_t = AlignerArg::Blastp)]
    pub aligner: AlignerArg,

    #[arg(long, short = 'P')]
    pub precomputed: Option<PathBuf>,

    /// Optional gspa run directory (see `docs/multi-genome.md`). Skips
    /// alignment; expands rep-level hits onto this genome's members.
    #[arg(long, conflicts_with = "precomputed")]
    pub gspa_run: Option<PathBuf>,

    /// Genome id matching a row in the gspa manifest's `genomes.tsv`.
    /// Defaults to the FASTA's file stem.
    #[arg(long)]
    pub gspa_genome_id: Option<String>,

    /// Treat the gspa alignment TSV's coverage column as a 0–1 fraction
    /// (mmseqs2 native). Default is 0–100 (diamond/blastp).
    #[arg(long)]
    pub gspa_coverage_fraction: bool,

    #[arg(long, short = 'b', default_value_t = 50.0)]
    pub bitcutoff: f32,

    #[arg(long, short = 'i', default_value_t = 0.0)]
    pub identcutoff: f32,

    #[arg(long, short = 'c', default_value_t = 75)]
    pub coverage: u32,

    #[arg(long, short = 'a')]
    pub nouse_alternatives: bool,

    /// Only check this substrate keyword.
    #[arg(long, short = 'm')]
    pub only_met: Option<String>,

    #[arg(long, short = 'o', default_value = ".")]
    pub out_dir: PathBuf,

    #[arg(long, short = 'u')]
    pub suffix: Option<String>,

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

pub fn run_cli(
    args: Args,
    data_dir_override: Option<&Path>,
    seq_dir_override: Option<&Path>,
) -> anyhow::Result<()> {
    let data_dir = resolve_data_dir(data_dir_override)?;
    let seq_dir = resolve_seq_dir(seq_dir_override, &data_dir)?;

    tracing::info!(
        data_dir = %data_dir.display(),
        seq_dir = %seq_dir.display(),
        "gapsmith find-transport starting"
    );

    // Reference FASTAs: tcdb.fasta + transporter.fasta (both in `seq_dir`).
    let tcdb = seq_dir.join("tcdb.fasta");
    let transporter_fa = seq_dir.join("transporter.fasta");
    let reference: Vec<&Path> =
        [&tcdb, &transporter_fa].iter().map(|p| p.as_path()).collect();

    // Subex + metabolite name map.
    let subex_rows = subex::load(data_dir.join("subex.tbl"))?;
    let seed_cpds = load_seed_metabolites(data_dir.join("seed_metabolites_edited.tsv"))?;
    let metabolite_name_by_id: HashMap<String, String> =
        seed_cpds.iter().map(|c| (c.id.as_str().to_string(), c.name.clone())).collect();

    // Aligner.
    let aligner: Box<dyn Aligner> = if let Some(gspa_root) = &args.gspa_run {
        let gid = match &args.gspa_genome_id {
            Some(s) if !s.is_empty() => s.clone(),
            _ => args
                .genome
                .file_stem()
                .and_then(|s| s.to_str())
                .map(|s| {
                    s.trim_end_matches(".faa")
                        .trim_end_matches(".fasta")
                        .to_string()
                })
                .unwrap_or_default(),
        };
        if gid.is_empty() {
            anyhow::bail!("--gspa-run requires --gspa-genome-id (or a genome FASTA with a usable file stem)");
        }
        Box::new(GspaRunAligner::new(gspa_root, gid, args.gspa_coverage_fraction)?)
    } else {
        match args.aligner {
            AlignerArg::Blastp => Box::new(BlastpAligner::new()),
            AlignerArg::Diamond => Box::new(DiamondAligner::new()),
            AlignerArg::Mmseqs2 => Box::new(Mmseqs2Aligner::new()),
            AlignerArg::Precomputed => {
                let p = args.precomputed.clone().ok_or_else(|| {
                    anyhow::anyhow!("--precomputed required when --aligner precomputed")
                })?;
                Box::new(PrecomputedTsvAligner::new_percentage(p))
            }
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
    let opts = TransportOptions {
        bitcutoff: args.bitcutoff,
        identcutoff: args.identcutoff,
        coverage_pct: args.coverage,
        nouse_alternatives: args.nouse_alternatives,
        only_met: args.only_met.as_deref(),
    };

    let workdir = tempfile::tempdir()?;
    let report = run(
        &reference,
        &data_dir.join("seed_transporter.tbl"),
        &data_dir.join("seed_transporter_custom.tbl"),
        &data_dir.join("tcdb_substrates.tbl"),
        &data_dir.join("tcdb_custom.tbl"),
        &subex_rows,
        &metabolite_name_by_id,
        &args.genome,
        aligner.as_ref(),
        &align_opts,
        &opts,
        workdir.path(),
    )?;

    std::fs::create_dir_all(&args.out_dir)?;
    let stem = args
        .genome
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("genome")
        .trim_end_matches(".faa")
        .trim_end_matches(".fasta")
        .to_string();
    let tag = args.suffix.unwrap_or_default();
    let out_name = if tag.is_empty() {
        format!("{stem}-Transporter.tbl")
    } else {
        format!("{stem}-{tag}-Transporter.tbl")
    };
    let out_path = args.out_dir.join(out_name);
    write_transporter_tbl(&report.rows, &out_path)?;
    eprintln!("wrote {} rows to {}", report.rows.len(), out_path.display());
    Ok(())
}
