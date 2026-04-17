//! `gapsmith align` — standalone alignment utility.
//!
//! Thin wrapper around the `gapsmith-align` crate. Useful for exercising any
//! of the backends (blastp, tblastn, diamond, mmseqs2, precomputed) without
//! the rest of the find/transport pipeline. Intended as a debugging and
//! plumbing-sanity tool; production flows drive the aligner via `gapseq
//! find` / `gapsmith find-transport`.

use clap::{Parser, ValueEnum};
use gapsmith_align::{
    AlignOpts, Aligner, BlastpAligner, DiamondAligner, Mmseqs2Aligner, PrecomputedTsvAligner,
    TblastnAligner,
};
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    /// Which aligner backend to invoke.
    #[arg(long, short = 'A', value_enum, default_value_t = Backend::Diamond)]
    pub aligner: Backend,

    /// Query FASTA (protein). Required unless `--aligner precomputed`.
    #[arg(long, short)]
    pub query: Option<PathBuf>,

    /// Target FASTA (protein for blastp/diamond/mmseqs2, nucleotide for
    /// tblastn). Required unless `--aligner precomputed`.
    #[arg(long, short)]
    pub target: Option<PathBuf>,

    /// Pre-computed alignment TSV (only with `--aligner precomputed`).
    #[arg(long, short = 'P')]
    pub precomputed: Option<PathBuf>,

    /// Treat precomputed `qcov` as 0–1 fraction (mmseqs native). Otherwise
    /// interpreted as 0–100.
    #[arg(long)]
    pub fraction_coverage: bool,

    /// Output TSV path. `-` writes to stdout. Default: stdout.
    #[arg(long, short, default_value = "-")]
    pub out: PathBuf,

    /// Thread count for the aligner.
    #[arg(long, short = 'K')]
    pub threads: Option<usize>,

    /// Query-coverage cutoff (0–100). Default 75.
    #[arg(long, short = 'c', default_value_t = 75)]
    pub coverage: u32,

    /// E-value cutoff. Default: tool's default.
    #[arg(long, short = 'e')]
    pub evalue: Option<f64>,

    /// Additional arguments forwarded to the underlying tool.
    #[arg(long = "extra", num_args = 0..)]
    pub extra: Vec<String>,
}

/// Format an e-value the way blast/diamond CLIs do: exact `0`, exact integers,
/// otherwise scientific notation. Avoids f64's default Display emitting
/// multi-line subnormals like `0.000…0477`.
fn fmt_evalue(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else if v.abs() < 1e-3 || v.abs() >= 1e5 {
        format!("{v:.3e}")
    } else {
        format!("{v}")
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Backend {
    Blastp,
    Tblastn,
    Diamond,
    Mmseqs2,
    Precomputed,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let opts = AlignOpts {
        threads: args.threads.unwrap_or_else(|| {
            std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
        }),
        coverage_pct: args.coverage,
        evalue: args.evalue,
        extra_args: args.extra,
        quiet: false,
    };

    let (aligner, query, target): (Box<dyn Aligner>, PathBuf, PathBuf) = match args.aligner {
        Backend::Blastp => (
            Box::new(BlastpAligner::new()),
            args.query.ok_or_else(|| anyhow::anyhow!("--query required"))?,
            args.target.ok_or_else(|| anyhow::anyhow!("--target required"))?,
        ),
        Backend::Tblastn => (
            Box::new(TblastnAligner::new()),
            args.query.ok_or_else(|| anyhow::anyhow!("--query required"))?,
            args.target.ok_or_else(|| anyhow::anyhow!("--target required"))?,
        ),
        Backend::Diamond => (
            Box::new(DiamondAligner::new()),
            args.query.ok_or_else(|| anyhow::anyhow!("--query required"))?,
            args.target.ok_or_else(|| anyhow::anyhow!("--target required"))?,
        ),
        Backend::Mmseqs2 => (
            Box::new(Mmseqs2Aligner::new()),
            args.query.ok_or_else(|| anyhow::anyhow!("--query required"))?,
            args.target.ok_or_else(|| anyhow::anyhow!("--target required"))?,
        ),
        Backend::Precomputed => {
            let p = args.precomputed.ok_or_else(|| {
                anyhow::anyhow!("--precomputed required when --aligner precomputed")
            })?;
            let a: Box<dyn Aligner> = if args.fraction_coverage {
                Box::new(PrecomputedTsvAligner::new_fraction(&p))
            } else {
                Box::new(PrecomputedTsvAligner::new_percentage(&p))
            };
            (a, PathBuf::from("_"), PathBuf::from("_"))
        }
    };

    tracing::info!(
        aligner = aligner.name(),
        query = %query.display(),
        target = %target.display(),
        "running alignment"
    );
    let hits = aligner.align(&query, &target, &opts)?;

    let stdout = std::io::stdout();
    let mut sink: Box<dyn Write> = if args.out.as_os_str() == "-" {
        Box::new(stdout.lock())
    } else {
        Box::new(std::io::BufWriter::new(std::fs::File::create(&args.out)?))
    };
    for h in &hits {
        writeln!(
            sink,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            h.qseqid,
            h.pident,
            fmt_evalue(h.evalue),
            h.bitscore,
            h.qcov,
            h.stitle,
            h.sstart,
            h.send
        )?;
    }
    eprintln!("wrote {} hits", hits.len());
    Ok(())
}
