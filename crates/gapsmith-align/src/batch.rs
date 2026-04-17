//! Batch-cluster alignment across N genomes.
//!
//! When annotating many genomes against the same gapseq reference seqdb,
//! the per-genome loop is dominated by redundant alignments — closely
//! related genomes share a lot of proteins. The batch-cluster path
//! amortizes that cost:
//!
//! 1. [`concat_genomes`] merges every genome's proteome into one FASTA,
//!    rewriting each header as `<GENOMEID>|<orig_header>` so we can trace
//!    any downstream hit back to its genome of origin.
//! 2. [`cluster_with_mmseqs`] calls `mmseqs easy-cluster` to produce
//!    (a) a representative FASTA and (b) a two-column `<rep>\t<member>`
//!    TSV.
//! 3. [`BatchClusterAligner::align_genomes`] then runs a single alignment
//!    of the gapseq query FASTA against the representatives and expands
//!    each rep-hit into its cluster members, bucketing the results by
//!    genome.
//!
//! The expansion is an approximation: a member inherits its
//! representative's bitscore/identity. Users who need per-member precision
//! can re-align the affected members with any of the standard aligners —
//! typically a tiny fraction of the original N-genome cost.

use crate::error::{io_err, AlignError};
use crate::hit::Hit;
use crate::{AlignOpts, Aligner};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Separator between the genome ID and the original FASTA header. A pipe
/// `|` matches the convention used in the R `prepare_batch_alignments.R`
/// pipeline; the same character also appears in UniProt-style `sp|P12345|NAME`
/// accessions, so parsers must split on the *first* `|` only.
pub const GENOME_SEP: char = '|';

/// A genome's protein FASTA plus a short identifier used to prefix every
/// header on concatenation.
#[derive(Debug, Clone)]
pub struct GenomeInput {
    pub id: String,
    pub fasta: PathBuf,
}

/// Hits for one genome after batch-cluster expansion. Emitted in original
/// input order so downstream code can pair them up with `GenomeInput` by
/// index without extra lookups.
#[derive(Debug, Clone)]
pub struct GenomeHitSet {
    pub genome_id: String,
    pub hits: Vec<Hit>,
}

/// Driver for the batch-cluster mode.
///
/// `inner` is any per-pair [`Aligner`] — typically [`crate::DiamondAligner`]
/// or [`crate::Mmseqs2Aligner`]. It aligns query-vs-representative exactly
/// once; the expansion layer then rebuilds per-genome hit sets.
pub struct BatchClusterAligner {
    pub inner: Box<dyn Aligner>,
    pub cluster_identity: f32,
    pub cluster_coverage: f32,
}

impl BatchClusterAligner {
    /// Sensible defaults: 0.5 identity / 0.8 coverage mirror mmseqs2's own
    /// default `easy-cluster` parameters.
    pub fn new(inner: Box<dyn Aligner>) -> Self {
        Self { inner, cluster_identity: 0.5, cluster_coverage: 0.8 }
    }

    /// Run the full batch-cluster pipeline. `workdir` is where intermediate
    /// files live (concatenated FASTA, cluster outputs, alignment TSV);
    /// pass a `tempfile::TempDir::path()` if you don't need to persist them.
    pub fn align_genomes(
        &self,
        query_fasta: &Path,
        genomes: &[GenomeInput],
        workdir: &Path,
        opts: &AlignOpts,
    ) -> Result<Vec<GenomeHitSet>, AlignError> {
        std::fs::create_dir_all(workdir).map_err(|e| io_err(workdir, e))?;

        // 1) Concatenate.
        let merged = workdir.join("all_genomes.faa");
        concat_genomes(genomes, &merged)?;

        // 2) Cluster.
        let cluster = cluster_with_mmseqs(
            &merged,
            workdir,
            self.cluster_identity,
            self.cluster_coverage,
            opts,
        )?;

        tracing::info!(
            genomes = genomes.len(),
            members = cluster.total_members,
            representatives = cluster.rep_to_members.len(),
            "clustering complete"
        );

        // 3) Single alignment vs representatives.
        let rep_hits = self.inner.align(query_fasta, &cluster.representatives_fasta, opts)?;

        // 4) Expand per-genome.
        let per_genome = expand_hits(&rep_hits, &cluster, genomes);
        Ok(per_genome)
    }
}

// -- Concatenation --

/// Produce a single FASTA at `out` whose headers are rewritten to
/// `>GENOMEID|ORIGHEADER`. Fails if any genome ID contains `|` (because
/// the downstream split would be ambiguous).
pub fn concat_genomes(genomes: &[GenomeInput], out: &Path) -> Result<(), AlignError> {
    for g in genomes {
        if g.id.contains(GENOME_SEP) || g.id.is_empty() {
            return Err(AlignError::BadArg(format!(
                "genome id `{}` must be non-empty and must not contain `{GENOME_SEP}`",
                g.id
            )));
        }
    }
    let f = File::create(out).map_err(|e| io_err(out, e))?;
    let mut w = BufWriter::new(f);
    for g in genomes {
        let r = File::open(&g.fasta).map_err(|e| io_err(&g.fasta, e))?;
        let rdr = BufReader::new(r);
        for line in rdr.lines() {
            let line = line.map_err(|e| io_err(&g.fasta, e))?;
            if let Some(rest) = line.strip_prefix('>') {
                writeln!(w, ">{}{}{}", g.id, GENOME_SEP, rest).map_err(|e| io_err(out, e))?;
            } else {
                writeln!(w, "{line}").map_err(|e| io_err(out, e))?;
            }
        }
    }
    w.flush().map_err(|e| io_err(out, e))?;
    Ok(())
}

/// Pull `GENOMEID` out of a concatenated header (first `|`-separated token).
/// Returns `None` on a malformed header.
pub fn split_genome_prefix(header: &str) -> Option<(&str, &str)> {
    header.split_once(GENOME_SEP)
}

// -- Clustering --

/// Result of an `mmseqs easy-cluster` run.
pub struct ClusterResult {
    pub representatives_fasta: PathBuf,
    /// Map: representative-header → list of member headers (including the
    /// representative itself, which mmseqs includes as one of the members).
    pub rep_to_members: HashMap<String, Vec<String>>,
    pub total_members: usize,
}

pub fn cluster_with_mmseqs(
    merged_fasta: &Path,
    workdir: &Path,
    min_identity: f32,
    min_coverage: f32,
    opts: &AlignOpts,
) -> Result<ClusterResult, AlignError> {
    require_mmseqs()?;

    let prefix = workdir.join("cluster");
    let scratch = workdir.join("cluster_tmp");
    std::fs::create_dir_all(&scratch).map_err(|e| io_err(&scratch, e))?;

    let verbosity = if opts.quiet { "1" } else { "3" };
    let mut cmd = Command::new("mmseqs");
    cmd.arg("easy-cluster")
        .arg(merged_fasta)
        .arg(&prefix)
        .arg(&scratch)
        .arg("--min-seq-id")
        .arg(format!("{min_identity:.2}"))
        .arg("-c")
        .arg(format!("{min_coverage:.2}"))
        .arg("--threads")
        .arg(opts.threads.to_string())
        .arg("-v")
        .arg(verbosity);

    tracing::debug!(?cmd, "running mmseqs easy-cluster");
    let out = cmd.output().map_err(|e| io_err(Path::new("mmseqs"), e))?;
    if !out.status.success() {
        return Err(AlignError::ToolFailed {
            tool: "mmseqs2",
            status: out.status,
            stderr: String::from_utf8_lossy(&out.stderr).to_string(),
        });
    }

    // easy-cluster writes:
    //   <prefix>_rep_seq.fasta       (representatives)
    //   <prefix>_cluster.tsv         (two cols: rep\tmember)
    //   <prefix>_all_seqs.fasta      (full clustered fasta, ignored here)
    let rep_fa: PathBuf = append_ext(&prefix, "_rep_seq.fasta");
    let cluster_tsv: PathBuf = append_ext(&prefix, "_cluster.tsv");

    let rep_to_members = parse_cluster_tsv(&cluster_tsv)?;
    let total_members: usize = rep_to_members.values().map(|v| v.len()).sum();

    Ok(ClusterResult { representatives_fasta: rep_fa, rep_to_members, total_members })
}

fn append_ext(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

pub fn parse_cluster_tsv(path: &Path) -> Result<HashMap<String, Vec<String>>, AlignError> {
    let f = File::open(path).map_err(|e| io_err(path, e))?;
    let rdr = BufReader::new(f);
    let mut out: HashMap<String, Vec<String>> = HashMap::new();
    for (i, line) in rdr.lines().enumerate() {
        let line = line.map_err(|e| io_err(path, e))?;
        if line.is_empty() {
            continue;
        }
        let (rep, member) = line.split_once('\t').ok_or_else(|| AlignError::TsvParse {
            line: (i + 1) as u64,
            msg: "cluster TSV row has no tab".into(),
        })?;
        out.entry(rep.to_string()).or_default().push(member.to_string());
    }
    Ok(out)
}

// -- Expansion --

fn expand_hits(
    rep_hits: &[Hit],
    cluster: &ClusterResult,
    genomes: &[GenomeInput],
) -> Vec<GenomeHitSet> {
    // Preserve genome order for stable output.
    let mut per_genome: BTreeMap<usize, Vec<Hit>> = BTreeMap::new();
    let index: HashMap<&str, usize> =
        genomes.iter().enumerate().map(|(i, g)| (g.id.as_str(), i)).collect();

    // mmseqs's cluster TSV keys are the representative's *first*
    // whitespace-delimited token (accession only); diamond/blastp's
    // `stitle` field carries the full header. Normalize both sides to
    // the accession form before the lookup.
    let to_accession = |h: &str| h.split_whitespace().next().unwrap_or("").to_string();
    let key_by_acc: HashMap<String, &Vec<String>> = cluster
        .rep_to_members
        .iter()
        .map(|(k, v)| (to_accession(k), v))
        .collect();

    for rep_hit in rep_hits {
        let rep_key = to_accession(&rep_hit.stitle);
        let members = match key_by_acc.get(rep_key.as_str()) {
            Some(m) => m.as_slice(),
            None => {
                tracing::warn!(
                    rep = %rep_key,
                    "rep-hit has no matching cluster — dropping"
                );
                continue;
            }
        };
        for m in members {
            let (genome_id, orig_header) = match split_genome_prefix(m) {
                Some((g, o)) => (g, o),
                None => continue,
            };
            let Some(&gidx) = index.get(genome_id) else { continue };
            let mut h = rep_hit.clone();
            h.stitle = orig_header.to_string();
            per_genome.entry(gidx).or_default().push(h);
        }
    }
    drop(key_by_acc);
    (0..genomes.len())
        .map(|i| GenomeHitSet {
            genome_id: genomes[i].id.clone(),
            hits: per_genome.remove(&i).unwrap_or_default(),
        })
        .collect()
}

fn require_mmseqs() -> Result<(), AlignError> {
    if crate::blast::which("mmseqs").is_some() {
        Ok(())
    } else {
        Err(AlignError::ToolMissing { tool: "mmseqs" })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn concat_rewrites_headers() {
        let dir = tempfile::tempdir().unwrap();
        let a = dir.path().join("A.faa");
        let b = dir.path().join("B.faa");
        std::fs::write(&a, ">seq1 descA\nMKLMA\n>seq2\nEEPS\n").unwrap();
        std::fs::write(&b, ">x\nAAA\n").unwrap();
        let out = dir.path().join("all.faa");
        concat_genomes(
            &[
                GenomeInput { id: "g1".into(), fasta: a },
                GenomeInput { id: "g2".into(), fasta: b },
            ],
            &out,
        )
        .unwrap();
        let mut s = String::new();
        File::open(&out).unwrap().read_to_string(&mut s).unwrap();
        assert!(s.contains(">g1|seq1 descA\n"));
        assert!(s.contains(">g1|seq2\n"));
        assert!(s.contains(">g2|x\n"));
        assert!(s.contains("MKLMA"));
    }

    #[test]
    fn concat_rejects_bad_genome_id() {
        let dir = tempfile::tempdir().unwrap();
        let a = dir.path().join("A.faa");
        std::fs::write(&a, ">s\nA\n").unwrap();
        let err = concat_genomes(
            &[GenomeInput { id: "g|1".into(), fasta: a }],
            &dir.path().join("out.faa"),
        )
        .unwrap_err();
        assert!(matches!(err, AlignError::BadArg(_)));
    }

    #[test]
    fn parses_cluster_tsv() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("c.tsv");
        std::fs::write(&p, "rep1\trep1\nrep1\tmemA\nrep2\trep2\n").unwrap();
        let m = parse_cluster_tsv(&p).unwrap();
        assert_eq!(m.get("rep1").unwrap(), &vec!["rep1".to_string(), "memA".to_string()]);
        assert_eq!(m.get("rep2").unwrap(), &vec!["rep2".to_string()]);
    }

    #[test]
    fn split_prefix_first_pipe_only() {
        assert_eq!(split_genome_prefix("g1|sp|P12345|X"), Some(("g1", "sp|P12345|X")));
        assert_eq!(split_genome_prefix("unpre"), None);
    }
}
