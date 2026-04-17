//! `gapsmith batch-align` — cluster N genomes, align once, expand per-genome.
//!
//! See `crates/gapsmith-align/src/batch.rs` for the algorithm. This
//! subcommand is the user-facing entry point; it accepts a directory of
//! genome FASTAs (or a list) and a query FASTA, and writes per-genome TSVs
//! into an output directory.

use clap::{Parser, ValueEnum};
use gapsmith_align::{
    AlignOpts, Aligner, BatchClusterAligner, DiamondAligner, GenomeInput, Mmseqs2Aligner,
};
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Parser)]
pub struct Args {
    /// Query FASTA (protein) — e.g. concatenated reference sequences.
    #[arg(long, short)]
    pub query: PathBuf,

    /// Directory of genome FASTAs. Every `*.faa`, `*.fasta`, or `*.fa`
    /// inside is treated as one genome; the file stem becomes the genome ID.
    #[arg(long, short = 'g', conflicts_with = "genome_list")]
    pub genomes_dir: Option<PathBuf>,

    /// Alternative: TSV with `<genome_id>\t<fasta_path>` rows.
    #[arg(long = "genomes-list", conflicts_with = "genomes_dir")]
    pub genome_list: Option<PathBuf>,

    /// Output directory. One `<genome_id>.tsv` is written per input genome.
    #[arg(long, short)]
    pub out: PathBuf,

    /// Persistent work directory for the concatenated FASTA, cluster files,
    /// and alignment result. Defaults to a temp dir cleaned up on exit.
    #[arg(long)]
    pub workdir: Option<PathBuf>,

    /// Inner aligner backend (query vs. cluster representatives).
    #[arg(long, value_enum, default_value_t = Inner::Diamond)]
    pub aligner: Inner,

    /// mmseqs `--min-seq-id` cutoff for clustering (0–1). Default 0.5.
    #[arg(long, default_value_t = 0.5)]
    pub cluster_identity: f32,

    /// mmseqs `-c` cutoff for clustering coverage (0–1). Default 0.8.
    #[arg(long, default_value_t = 0.8)]
    pub cluster_coverage: f32,

    /// Query coverage cutoff for the inner alignment (0–100). Default 75.
    #[arg(long, short = 'c', default_value_t = 75)]
    pub coverage: u32,

    /// Thread count.
    #[arg(long, short = 'K')]
    pub threads: Option<usize>,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Inner {
    Diamond,
    Mmseqs2,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let genomes = collect_genomes(&args)?;
    if genomes.is_empty() {
        anyhow::bail!("no genomes found under --genomes-dir / --genomes-list");
    }
    tracing::info!(n = genomes.len(), "genomes to batch-align");

    let inner: Box<dyn Aligner> = match args.aligner {
        Inner::Diamond => Box::new(DiamondAligner::new()),
        Inner::Mmseqs2 => Box::new(Mmseqs2Aligner::new()),
    };
    let batcher = BatchClusterAligner {
        inner,
        cluster_identity: args.cluster_identity,
        cluster_coverage: args.cluster_coverage,
    };

    let opts = AlignOpts {
        threads: args.threads.unwrap_or_else(|| {
            std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
        }),
        coverage_pct: args.coverage,
        evalue: None,
        extra_args: Vec::new(),
        quiet: false,
    };

    std::fs::create_dir_all(&args.out)?;

    let _tmp; // kept for lifetime if we allocate one
    let workdir_path: PathBuf = match args.workdir.clone() {
        Some(p) => {
            std::fs::create_dir_all(&p)?;
            p
        }
        None => {
            _tmp = tempfile::tempdir()?;
            _tmp.path().to_path_buf()
        }
    };

    let per_genome = batcher.align_genomes(&args.query, &genomes, &workdir_path, &opts)?;

    for gs in &per_genome {
        let out_path = args.out.join(format!("{}.tsv", gs.genome_id));
        let mut f = std::io::BufWriter::new(std::fs::File::create(&out_path)?);
        for h in &gs.hits {
            writeln!(
                f,
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
        eprintln!("  {}: {} hits -> {}", gs.genome_id, gs.hits.len(), out_path.display());
    }
    Ok(())
}

fn collect_genomes(args: &Args) -> anyhow::Result<Vec<GenomeInput>> {
    collect_genomes_from_paths(args.genomes_dir.as_deref(), args.genome_list.as_deref())
}

/// Shared helper used by `batch-align` and `doall-batch`.
///
/// Priority: if `dir` is `Some`, scan it for `*.fa` / `*.fasta` / `*.faa`
/// (and the `.gz` variants); if `list` is `Some`, parse a TSV of
/// `<id>\t<path>` rows. Exactly one must be supplied — both `None` is an
/// error, both `Some` is also an error.
pub fn collect_genomes_from_paths(
    dir: Option<&std::path::Path>,
    list: Option<&std::path::Path>,
) -> anyhow::Result<Vec<GenomeInput>> {
    match (dir, list) {
        (Some(_), Some(_)) => {
            anyhow::bail!("pass either --genomes-dir or --genomes-list, not both")
        }
        (None, None) => anyhow::bail!("one of --genomes-dir / --genomes-list is required"),
        (Some(d), None) => {
            let mut out = Vec::new();
            for entry in std::fs::read_dir(d)? {
                let e = entry?;
                let p = e.path();
                if !p.is_file() {
                    continue;
                }
                let name = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
                // Accept both `genome.faa` and `genome.faa.gz`; reject
                // everything else so we don't pick up sidecar files
                // (.dmnd, .phr, etc.) that commonly live beside FASTAs.
                let is_fasta = matches!(
                    name.rsplit_once('.').map(|(_, e)| e),
                    Some("faa") | Some("fasta") | Some("fa")
                ) || name.ends_with(".faa.gz")
                    || name.ends_with(".fasta.gz")
                    || name.ends_with(".fa.gz");
                if !is_fasta {
                    continue;
                }
                let id = strip_fasta_ext(name);
                if id.is_empty() {
                    continue;
                }
                out.push(GenomeInput { id: id.to_string(), fasta: p });
            }
            out.sort_by(|a, b| a.id.cmp(&b.id));
            Ok(out)
        }
        (None, Some(list_path)) => {
            let data = std::fs::read_to_string(list_path)?;
            let mut out = Vec::new();
            for line in data.lines() {
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }
                let mut parts = line.split('\t');
                let id = parts.next().unwrap_or("").trim();
                let path = parts.next().unwrap_or("").trim();
                if id.is_empty() || path.is_empty() {
                    anyhow::bail!("--genomes-list row needs `<id>\\t<path>`: `{line}`");
                }
                // Optional 3rd column (abundance) is accepted and ignored
                // here — doall-batch picks it up separately via --gspa-run.
                out.push(GenomeInput { id: id.into(), fasta: PathBuf::from(path) });
            }
            Ok(out)
        }
    }
}

fn strip_fasta_ext(name: &str) -> &str {
    for ext in [".faa.gz", ".fasta.gz", ".fa.gz", ".faa", ".fasta", ".fa"] {
        if let Some(stem) = name.strip_suffix(ext) {
            return stem;
        }
    }
    name
}

fn fmt_evalue(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else if v.abs() < 1e-3 || v.abs() >= 1e5 {
        format!("{v:.3e}")
    } else {
        format!("{v}")
    }
}
