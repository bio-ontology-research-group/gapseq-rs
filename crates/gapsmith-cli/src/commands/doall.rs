//! `gapsmith doall` — chain find → find-transport → draft → medium → fill.
//!
//! Port of `src/doall.sh`. Accepts a genome FASTA + common options and
//! drives the full reconstruction pipeline in one invocation. By default
//! the medium is inferred via `gapsmith medium` (pass `-m <file>` to skip).
//!
//! # Example
//!
//! ```bash
//! # Minimal reconstruction, fully automatic medium:
//! gapsmith doall genome.faa.gz -f out/
//!
//! # Diamond + custom medium:
//! gapsmith doall genome.faa.gz -f out/ -A diamond -m MM_glu.csv
//!
//! # Include fill Steps 3 + 4 (slow, but richer model):
//! gapsmith doall genome.faa.gz -f out/ --full-suite
//! ```
//!
//! # Flow (each step is an in-process call into the matching subcommand)
//!
//! 1. `gapsmith find -p all` → `<stem>-all-{Reactions,Pathways}.tbl`.
//! 2. `gapsmith find-transport` → `<stem>-Transporter.tbl`.
//! 3. `gapsmith draft` → `<stem>-draft.{gmod.cbor,xml}`.
//! 4. `gapsmith medium` (unless `-m <csv>` supplied) → `<stem>-medium.csv`.
//! 5. `gapsmith fill` → `<stem>-filled.{gmod.cbor,xml,-added.tsv}`.
//!
//! Gzipped FASTA is auto-decompressed into a `tempfile::tempdir` before
//! step 1 — `makeblastdb` and `mmseqs createdb` can't read `.gz` directly.

use clap::Parser;
use std::path::{Path, PathBuf};

use crate::commands::{draft, fill, find, find_transport, medium};

#[derive(Debug, Parser)]
pub struct Args {
    /// Input genome (protein FASTA; .gz accepted).
    pub genome: PathBuf,

    /// Aligner backend for `find` and `find-transport`.
    #[arg(long, short = 'A', value_enum, default_value_t = find::AlignerArg::Blastp)]
    pub aligner: find::AlignerArg,

    /// Bitscore cutoff passed to `find` / `find-transport`.
    #[arg(long, short = 'b', default_value_t = 200.0)]
    pub bitcutoff: f32,

    /// Coverage cutoff (0–100).
    #[arg(long, short = 'c', default_value_t = 75)]
    pub coverage: u32,

    /// Core-reaction bitscore cutoff (passed to `draft` and `fill`).
    #[arg(long, short = 'l', default_value_t = 50.0)]
    pub min_bs_core: f64,

    /// Taxonomy: `Bacteria`, `Archaea`, or `auto`. `auto` currently maps
    /// to `Bacteria` since HMM-based detection isn't ported yet.
    #[arg(long, short = 't', default_value = "auto")]
    pub taxonomy: String,

    /// Growth medium: an explicit CSV or `auto` (infer via `gapsmith medium`).
    #[arg(long, short = 'm', default_value = "auto")]
    pub medium: String,

    /// Output directory.
    #[arg(long, short = 'f', default_value = ".")]
    pub out_dir: PathBuf,

    /// Thread count (default: all cores).
    #[arg(long, short = 'K')]
    pub threads: Option<usize>,

    /// Skip Steps 3 + 4 of the fill suite (default: skipped for speed).
    /// Pass `--full-suite` to enable them.
    #[arg(long = "full-suite")]
    pub full_suite: bool,

    /// Minimum growth rate for gap-filling.
    #[arg(long, short = 'k', default_value_t = 0.01)]
    pub min_growth: f64,
}

pub fn run_cli(
    args: Args,
    data_dir_override: Option<&Path>,
    seq_dir_override: Option<&Path>,
) -> anyhow::Result<()> {
    let out_dir = args.out_dir.clone();
    std::fs::create_dir_all(&out_dir)?;

    let stem = infer_stem(&args.genome);
    eprintln!("gapsmith doall — genome: {} (stem: {})", args.genome.display(), stem);

    // makeblastdb / mmseqs2 createdb don't read gzipped FASTA directly.
    // Decompress to a temp copy when the input ends in `.gz`.
    let _tmp_fasta; // keep alive for the duration of the call
    let genome: std::path::PathBuf = if args.genome
        .extension()
        .is_some_and(|e| e == "gz")
    {
        let td = tempfile::tempdir_in(std::env::temp_dir())?;
        let out = td.path().join(format!("{stem}.faa"));
        decompress_gz(&args.genome, &out)?;
        eprintln!("  decompressed {} → {}", args.genome.display(), out.display());
        _tmp_fasta = td;
        out
    } else {
        _tmp_fasta = tempfile::tempdir_in(std::env::temp_dir())?;
        args.genome.clone()
    };

    // Map `taxonomy = auto` to Bacteria (HMM domain prediction not yet ported).
    let taxonomy = if args.taxonomy == "auto" {
        "Bacteria".to_string()
    } else {
        args.taxonomy.clone()
    };

    // ---------- 1. find -p all ----------
    eprintln!("\n[1/5] gapsmith find -p all");
    let find_args = find::Args {
        genome: genome.clone(),
        pathways: "all".into(),
        pathway_db: "metacyc,custom".into(),
        taxonomy: taxonomy.clone(),
        aligner: args.aligner,
        precomputed: None,
        bitcutoff: args.bitcutoff,
        identcutoff: 0.0,
        ident_exception: 70.0,
        coverage: args.coverage,
        completeness_main: 80.0,
        completeness_hints: 66.0,
        strict: false,
        include_superpathways: false,
        out_dir: out_dir.clone(),
        suffix: Some("all".into()),
        threads: args.threads,
    };
    find::run(find_args, data_dir_override, seq_dir_override)?;

    let reactions_tbl = out_dir.join(format!("{stem}-all-Reactions.tbl"));
    let pathways_tbl = out_dir.join(format!("{stem}-all-Pathways.tbl"));

    // ---------- 2. find-transport ----------
    eprintln!("\n[2/5] gapsmith find-transport");
    let tr_aligner = match args.aligner {
        find::AlignerArg::Blastp => find_transport::AlignerArg::Blastp,
        find::AlignerArg::Diamond => find_transport::AlignerArg::Diamond,
        find::AlignerArg::Mmseqs2 => find_transport::AlignerArg::Mmseqs2,
        find::AlignerArg::Precomputed => find_transport::AlignerArg::Precomputed,
    };
    let tr_args = find_transport::Args {
        genome: genome.clone(),
        aligner: tr_aligner,
        precomputed: None,
        bitcutoff: 50.0,
        identcutoff: 0.0,
        coverage: args.coverage,
        nouse_alternatives: false,
        only_met: None,
        out_dir: out_dir.clone(),
        suffix: None,
        threads: args.threads,
    };
    find_transport::run_cli(tr_args, data_dir_override, seq_dir_override)?;
    let transporter_tbl = out_dir.join(format!("{stem}-Transporter.tbl"));

    // ---------- 3. draft ----------
    eprintln!("\n[3/5] gapsmith draft");
    let draft_args = draft::Args {
        reactions: reactions_tbl.clone(),
        transporter: transporter_tbl.clone(),
        biomass: "auto".into(),
        name: Some(stem.clone()),
        high_evi_rxn_bs: args.bitcutoff,
        min_bs_for_core: args.min_bs_core as f32,
        out_dir: out_dir.clone(),
        no_sbml: false,
    };
    draft::run_cli(draft_args, data_dir_override)?;
    let draft_cbor = out_dir.join(format!("{stem}-draft.gmod.cbor"));

    // ---------- 4. medium (auto) or pass through --------
    let medium_path = if args.medium == "auto" {
        eprintln!("\n[4/5] gapsmith medium (auto)");
        let medium_args = medium::Args {
            model: draft_cbor.clone(),
            pathways: pathways_tbl.clone(),
            manual_flux: None,
            output: None,
            out_dir: out_dir.clone(),
        };
        medium::run_cli(medium_args, data_dir_override)?;
        out_dir.join(format!("{stem}-medium.csv"))
    } else {
        eprintln!("\n[4/5] using user-supplied medium {}", args.medium);
        PathBuf::from(&args.medium)
    };

    // ---------- 5. fill ----------
    eprintln!("\n[5/5] gapsmith fill");
    let fill_args = fill::Args {
        model: draft_cbor.clone(),
        media: medium_path.clone(),
        reactions: Some(reactions_tbl.clone()),
        target: "cpd11416".into(),
        min_growth: args.min_growth,
        bcore: args.min_bs_core,
        high_evi: args.bitcutoff as f64,
        dummy_weight: gapsmith_fill::DEFAULT_DUMMY_WEIGHT,
        out_dir: out_dir.clone(),
        no_sbml: false,
        step1_only: false,
        full_suite: args.full_suite,
        prune_futile: false,
    };
    fill::run_cli(fill_args, data_dir_override)?;

    eprintln!("\ngapseq doall complete. Outputs in {}", out_dir.display());
    Ok(())
}

/// Extract a stem id from the genome filename, stripping `.gz` / `.faa` /
/// `.fa` extensions. Mirrors `src/doall.sh:`:
/// ```sh
/// base=$(basename "$file"); id="${base%.*}"; [[ $file == *.gz ]] && id="${id%.*}"
/// ```
fn infer_stem(path: &Path) -> String {
    let base = path.file_name().and_then(|s| s.to_str()).unwrap_or("model");
    let gz_stripped = base.strip_suffix(".gz").unwrap_or(base);
    let dot = gz_stripped.rfind('.').unwrap_or(gz_stripped.len());
    gz_stripped[..dot].to_string()
}

fn decompress_gz(src: &Path, dst: &Path) -> std::io::Result<()> {
    use std::io::copy;
    let f = std::fs::File::open(src)?;
    let mut gz = flate2::read::GzDecoder::new(f);
    let mut out = std::fs::File::create(dst)?;
    copy(&mut gz, &mut out)?;
    Ok(())
}
