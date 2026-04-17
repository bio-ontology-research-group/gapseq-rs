//! `gapsmith doall-batch` — run [`doall`] across N genomes with rayon and
//! SLURM-array-friendly sharding.
//!
//! # When to reach for this vs. `doall`
//!
//! - **One genome** → `gapsmith doall <fasta>`.
//! - **A few dozen genomes on one machine** → `gapsmith doall-batch --threads <cores>`.
//! - **Thousands of genomes on a local box, with a gspa pre-cluster to
//!   amortize the alignment cost** → `gapsmith doall-batch --gspa-run <dir>
//!   --threads <cores>`.
//! - **1 M genomes on a SLURM cluster** → launch a 1024-task array; every
//!   task runs the same `doall-batch` with a different
//!   `--shard i/1024`. See `docs/multi-genome.md` for the full recipe.
//!
//! # What the command does
//!
//! 1. Resolve the full genome list (directory or TSV), sort deterministically.
//! 2. Apply `--shard i/N` → keep only the genomes whose index satisfies
//!    `idx.rem_euclid(N) == i`.
//! 3. Build a rayon pool sized by `--threads` and run [`doall::run_cli`]
//!    once per genome.
//! 4. Collect a per-genome pass / fail report; summary printed at the end.
//!
//! Each genome writes to its own subdirectory under `--out-dir/<genome_id>/`,
//! so concurrent workers never collide.
//!
//! # Shard semantics
//!
//! `--shard i/N` (with `0 ≤ i < N`) selects every genome whose *sorted*
//! 0-based index `k` satisfies `k.rem_euclid(N) == i`. That keeps
//! workload balance across shards even when the input isn't
//! pre-shuffled; no shared state needs to be read by sibling SLURM
//! tasks, so scaling to thousands of tasks costs nothing.

use clap::Parser;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::commands::{batch_align::collect_genomes_from_paths, doall, find};

/// Bundle of per-shard configuration. Parsed from `"i/N"`.
#[derive(Debug, Clone, Copy)]
pub struct ShardSpec {
    pub index: usize,
    pub count: usize,
}

impl ShardSpec {
    pub fn parse(s: &str) -> anyhow::Result<Self> {
        let (i, n) = s.split_once('/').ok_or_else(|| {
            anyhow::anyhow!("--shard must be `i/N` (got `{s}`)")
        })?;
        let index: usize = i.parse().map_err(|_| anyhow::anyhow!("shard index `{i}` is not an integer"))?;
        let count: usize = n.parse().map_err(|_| anyhow::anyhow!("shard count `{n}` is not an integer"))?;
        if count == 0 {
            anyhow::bail!("--shard count must be ≥ 1");
        }
        if index >= count {
            anyhow::bail!("--shard index {index} must be < count {count}");
        }
        Ok(Self { index, count })
    }

    pub fn keep(&self, idx: usize) -> bool {
        idx.rem_euclid(self.count) == self.index
    }
}

#[derive(Debug, Parser)]
pub struct Args {
    /// Directory of genome FASTAs (any `*.faa*`, `*.fasta*`, `*.fa*`,
    /// with or without `.gz`). Mutually exclusive with `--genomes-list`.
    #[arg(long, short = 'g', conflicts_with = "genomes_list")]
    pub genomes_dir: Option<PathBuf>,

    /// TSV with 2–3 columns: `<genome_id>\t<fasta_path>[\t<abundance>]`.
    /// Mutually exclusive with `--genomes-dir`. Takes precedence over any
    /// `genomes.tsv` found inside `--gspa-run`.
    #[arg(long = "genomes-list", conflicts_with = "genomes_dir")]
    pub genomes_list: Option<PathBuf>,

    /// Output directory. One `<genome_id>/` subdirectory is written per
    /// genome; inside, the normal `doall` outputs appear.
    #[arg(long, short = 'f')]
    pub out_dir: PathBuf,

    /// Worker threads — size of the rayon pool. Default: all detected cores.
    /// Each worker runs one `doall` at a time.
    #[arg(long, short = 'j')]
    pub jobs: Option<usize>,

    /// Shard selector for SLURM array jobs, as `<index>/<count>`.
    /// Genomes are split deterministically across shards after sorting.
    /// Omit to process every genome in one invocation.
    #[arg(long)]
    pub shard: Option<String>,

    /// Continue on per-genome `doall` failures. Default behaviour is to
    /// abort on the first error (same as `doall`). Set this for large
    /// batches where losing the whole run because of one bad FASTA is
    /// unacceptable.
    #[arg(long)]
    pub continue_on_error: bool,

    /// Optional gspa run directory. When set, every genome's `find`
    /// stage reads precomputed rep-hits from the manifest instead of
    /// running DIAMOND/BLASTp. See `docs/multi-genome.md`. Mutually
    /// compatible with `--aligner` (the aligner flag is ignored when a
    /// gspa run is supplied).
    #[arg(long)]
    pub gspa_run: Option<PathBuf>,

    // -- Passthrough options to `doall` ---------------------------------
    /// Aligner backend when `--gspa-run` is NOT set.
    #[arg(long, short = 'A', value_enum, default_value_t = find::AlignerArg::Diamond)]
    pub aligner: find::AlignerArg,

    #[arg(long, short = 'b', default_value_t = 200.0)]
    pub bitcutoff: f32,

    #[arg(long, short = 'c', default_value_t = 75)]
    pub coverage: u32,

    #[arg(long, short = 'l', default_value_t = 50.0)]
    pub min_bs_core: f64,

    #[arg(long, short = 't', default_value = "auto")]
    pub taxonomy: String,

    #[arg(long, short = 'm', default_value = "auto")]
    pub medium: String,

    #[arg(long, short = 'k', default_value_t = 0.01)]
    pub min_growth: f64,

    #[arg(long = "full-suite")]
    pub full_suite: bool,
}

pub fn run(
    args: Args,
    data_dir_override: Option<&Path>,
    seq_dir_override: Option<&Path>,
) -> anyhow::Result<()> {
    let mut genomes = resolve_genome_list(&args)?;
    genomes.sort_by(|a, b| a.0.cmp(&b.0));

    if genomes.is_empty() {
        anyhow::bail!("no genomes resolved from --genomes-dir / --genomes-list / --gspa-run");
    }

    let shard = match args.shard.as_deref() {
        Some(s) => Some(ShardSpec::parse(s)?),
        None => None,
    };
    let selected: Vec<(String, PathBuf)> = genomes
        .into_iter()
        .enumerate()
        .filter(|(i, _)| shard.map_or(true, |s| s.keep(*i)))
        .map(|(_, g)| g)
        .collect();

    eprintln!(
        "gapsmith doall-batch — {} genome{} {}",
        selected.len(),
        if selected.len() == 1 { "" } else { "s" },
        match shard {
            Some(s) => format!("(shard {}/{})", s.index, s.count),
            None => String::new(),
        }
    );
    std::fs::create_dir_all(&args.out_dir)?;

    let pool = build_pool(args.jobs)?;

    let succeeded = AtomicUsize::new(0);
    let failed = AtomicUsize::new(0);
    let continue_on_error = args.continue_on_error;

    // Errors from *inside* the rayon closure can't be surfaced cleanly
    // when `continue_on_error` is true, so we collect them as strings.
    let errors: std::sync::Mutex<Vec<(String, String)>> = std::sync::Mutex::new(Vec::new());

    pool.install(|| {
        selected.par_iter().for_each(|(gid, fasta)| {
            let genome_out = args.out_dir.join(gid);
            if let Err(e) = std::fs::create_dir_all(&genome_out) {
                fail(&errors, gid, format!("mkdir {genome_out:?}: {e}"));
                failed.fetch_add(1, Ordering::Relaxed);
                return;
            }

            let doall_args = doall::Args {
                genome: fasta.clone(),
                aligner: args.aligner,
                bitcutoff: args.bitcutoff,
                coverage: args.coverage,
                min_bs_core: args.min_bs_core,
                taxonomy: args.taxonomy.clone(),
                medium: args.medium.clone(),
                out_dir: genome_out,
                threads: Some(1), // one thread per genome — rayon handles parallelism
                full_suite: args.full_suite,
                min_growth: args.min_growth,
                gspa_run: args.gspa_run.clone(),
                gspa_genome_id: Some(gid.clone()),
            };

            match doall::run_cli(doall_args, data_dir_override, seq_dir_override) {
                Ok(()) => {
                    succeeded.fetch_add(1, Ordering::Relaxed);
                }
                Err(e) => {
                    failed.fetch_add(1, Ordering::Relaxed);
                    fail(&errors, gid, format!("{e}"));
                }
            }
        });
    });

    let ok = succeeded.load(Ordering::Relaxed);
    let bad = failed.load(Ordering::Relaxed);
    eprintln!("\ndoall-batch complete: {ok} succeeded, {bad} failed");

    let errs = errors.into_inner().unwrap_or_default();
    if !errs.is_empty() {
        let report = args.out_dir.join("doall-batch-errors.tsv");
        if let Ok(mut f) = std::fs::File::create(&report) {
            use std::io::Write;
            for (g, msg) in &errs {
                let _ = writeln!(f, "{g}\t{}", msg.replace('\t', " ").replace('\n', " | "));
            }
            eprintln!("error log: {}", report.display());
        }
        if !continue_on_error {
            anyhow::bail!("{bad} genome(s) failed — see doall-batch-errors.tsv");
        }
    }
    Ok(())
}

fn fail(
    errors: &std::sync::Mutex<Vec<(String, String)>>,
    gid: &str,
    msg: String,
) {
    eprintln!("  [FAIL] {gid}: {msg}");
    if let Ok(mut guard) = errors.lock() {
        guard.push((gid.to_string(), msg));
    }
}

fn build_pool(jobs: Option<usize>) -> anyhow::Result<rayon::ThreadPool> {
    let n = jobs.unwrap_or_else(|| {
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    });
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .map_err(|e| anyhow::anyhow!("failed to build rayon pool ({n} threads): {e}"))
}

/// Resolve the `(genome_id, fasta)` list with the following priority:
///
/// 1. `--genomes-list` (explicit TSV).
/// 2. `--genomes-dir` (scan a directory for `*.fa*` / `*.faa*`).
/// 3. `--gspa-run`'s built-in `genomes.tsv` — saves the user an extra flag
///    when they've already set up a manifest.
fn resolve_genome_list(args: &Args) -> anyhow::Result<Vec<(String, PathBuf)>> {
    if args.genomes_list.is_some() || args.genomes_dir.is_some() {
        return collect_genomes_from_paths(args.genomes_dir.as_deref(), args.genomes_list.as_deref())
            .map(|v| v.into_iter().map(|g| (g.id, g.fasta)).collect());
    }
    if let Some(gr) = &args.gspa_run {
        let m = gapsmith_align::GspaManifest::load(gr)?;
        return Ok(m.genomes.into_iter().map(|g| (g.id, g.fasta)).collect());
    }
    anyhow::bail!("one of --genomes-dir / --genomes-list / --gspa-run must be provided");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shard_parses_and_keeps_balanced() {
        let s = ShardSpec::parse("3/10").unwrap();
        assert_eq!(s.index, 3);
        assert_eq!(s.count, 10);
        assert!(s.keep(3));
        assert!(s.keep(13));
        assert!(!s.keep(4));
    }

    #[test]
    fn shard_rejects_bad_input() {
        assert!(ShardSpec::parse("3").is_err());
        assert!(ShardSpec::parse("3/0").is_err());
        assert!(ShardSpec::parse("10/5").is_err());
    }
}
