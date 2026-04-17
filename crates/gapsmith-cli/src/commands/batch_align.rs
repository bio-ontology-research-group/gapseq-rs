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
    if let Some(dir) = &args.genomes_dir {
        let mut out = Vec::new();
        for entry in std::fs::read_dir(dir)? {
            let e = entry?;
            let p = e.path();
            if !p.is_file() {
                continue;
            }
            let ext = p.extension().and_then(|s| s.to_str()).unwrap_or("");
            if !matches!(ext, "faa" | "fasta" | "fa") {
                continue;
            }
            let id = p.file_stem().and_then(|s| s.to_str()).unwrap_or("").to_string();
            if id.is_empty() {
                continue;
            }
            out.push(GenomeInput { id, fasta: p });
        }
        out.sort_by(|a, b| a.id.cmp(&b.id));
        Ok(out)
    } else if let Some(list) = &args.genome_list {
        let data = std::fs::read_to_string(list)?;
        let mut out = Vec::new();
        for line in data.lines() {
            if line.trim().is_empty() {
                continue;
            }
            let (id, path) = line.split_once('\t').ok_or_else(|| {
                anyhow::anyhow!("--genomes-list row missing tab: `{line}`")
            })?;
            out.push(GenomeInput { id: id.trim().into(), fasta: PathBuf::from(path.trim()) });
        }
        Ok(out)
    } else {
        anyhow::bail!("one of --genomes-dir / --genomes-list is required");
    }
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
