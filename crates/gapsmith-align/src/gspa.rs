//! Cluster-aware aligner backed by a gspa run directory.
//!
//! # Motivation
//!
//! When annotating thousands of genomes (or MAGs from metagenomes) the
//! dominant cost is redundant protein alignment — closely related genomes
//! share a lot of sequence. Upstream of gapsmith, the
//! [gspa](https://github.com/bio-ontology-research-group/gspa) pipeline
//! already clusters proteins across N genomes with mmseqs2/linclust and
//! runs every predictor exactly once against the cluster representatives.
//! gapsmith's `find` stage can then *consume* those precomputed alignments
//! instead of re-running DIAMOND / BLASTp.
//!
//! [`GspaRunAligner`] makes that wiring a one-liner: give it a gspa run
//! directory and a target genome id; it returns hits for *that* genome by
//! expanding every representative hit to the members that belong to the
//! genome in question.
//!
//! # Manifest layout (v1)
//!
//! ```text
//! <gspa-run>/
//!   gspa-manifest.toml      # optional; version marker + creation metadata
//!   clusters.tsv            # 3 columns: rep_id<TAB>member_id<TAB>genome_id
//!   alignment/reference.tsv # 8-column gapsmith-style TSV keyed on rep_id
//!   genomes.tsv             # 2–3 columns: genome_id<TAB>faa_path[<TAB>abundance]
//! ```
//!
//! The single `alignment/reference.tsv` is the simplest case; multiple
//! alignment files are concatenated transparently. The manifest is plain
//! text — easy to produce from any upstream pipeline (gspa-nf, a shell
//! script, a Snakemake rule, etc.).
//!
//! # Compatibility with gspa's native output
//!
//! gspa writes its cross-genome clusters via the same mmseqs2 `easy-cluster`
//! two-column TSV that [`crate::batch`] consumes. The third column
//! (`genome_id`) is added either by gspa-nf's join step or by the tiny
//! `prepare-gspa-run.sh` helper shipped in this repo (see
//! `docs/multi-genome.md`). Either way, the on-disk contract is stable and
//! parseable without extra dependencies.
//!
//! # What this file deliberately does NOT do
//!
//! - It does **not** run mmseqs2 itself; clustering is gspa's job.
//! - It does **not** merge evidence from multiple gspa predictors (diamond
//!   vs. mmseqs2 vs. foldseek) — that is gspa's integration layer. gapsmith
//!   only consumes one alignment TSV at a time, the same way BLASTp is one
//!   alignment. If you want to use multiple predictors, concatenate them
//!   upstream (they share the 8-column format).

use crate::error::{io_err, AlignError};
use crate::hit::Hit;
use crate::tsv::parse_tsv;
use crate::{AlignOpts, Aligner};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// Default file names inside a gspa run directory.
pub const CLUSTERS_FILE: &str = "clusters.tsv";
pub const GENOMES_FILE: &str = "genomes.tsv";
pub const ALIGNMENT_DIR: &str = "alignment";

/// Parsed gspa run manifest. Cheap to clone (all paths are owned);
/// the expensive per-rep hit fan-out happens inside [`GspaRunAligner`].
#[derive(Debug, Clone)]
pub struct GspaManifest {
    pub root: PathBuf,
    /// `rep_id` → list of `(member_id, genome_id)`. A single rep can map
    /// members in many genomes; per-genome expansion filters this at
    /// align-time.
    pub rep_to_members: HashMap<String, Vec<(String, String)>>,
    /// `genome_id` → (faa_path, optional abundance weight).
    /// Abundance is passed through to downstream cFBA / community code;
    /// the aligner itself doesn't use it.
    pub genomes: Vec<GenomeRow>,
    /// Resolved paths of every alignment TSV under `alignment/`.
    pub alignment_files: Vec<PathBuf>,
}

#[derive(Debug, Clone)]
pub struct GenomeRow {
    pub id: String,
    pub fasta: PathBuf,
    /// Optional relative abundance (counts / RPKM / TPM / whatever is
    /// meaningful to the caller). `None` when the column is absent.
    pub abundance: Option<f64>,
}

impl GspaManifest {
    /// Read `<root>/clusters.tsv`, `<root>/genomes.tsv`, and every `*.tsv`
    /// inside `<root>/alignment/`. The optional `gspa-manifest.toml` is
    /// *not* required — its absence is OK.
    pub fn load(root: impl Into<PathBuf>) -> Result<Self, AlignError> {
        let root: PathBuf = root.into();
        if !root.is_dir() {
            return Err(AlignError::BadArg(format!(
                "gspa-run `{}` is not a directory",
                root.display()
            )));
        }
        let rep_to_members = parse_clusters(&root.join(CLUSTERS_FILE))?;
        let genomes = parse_genomes(&root.join(GENOMES_FILE))?;
        let alignment_files = collect_alignment_files(&root.join(ALIGNMENT_DIR))?;
        Ok(Self { root, rep_to_members, genomes, alignment_files })
    }

    pub fn genome_ids(&self) -> Vec<&str> {
        self.genomes.iter().map(|g| g.id.as_str()).collect()
    }

    pub fn genome(&self, id: &str) -> Option<&GenomeRow> {
        self.genomes.iter().find(|g| g.id == id)
    }
}

/// Cluster-aware aligner. Per [`Aligner`] contract, `align()` ignores the
/// `query_fasta` / `target_fasta` arguments — the work was done upstream
/// by gspa. The caller selects *which* genome's hits to return at
/// construction time.
pub struct GspaRunAligner {
    manifest: GspaManifest,
    target_genome_id: String,
    coverage_is_fraction: bool,
    /// Parsed rep-level hits, built lazily on first call.
    cache: OnceLock<Vec<Hit>>,
}

impl GspaRunAligner {
    /// Construct from an already-parsed manifest. `coverage_is_fraction`
    /// follows the same convention as [`crate::PrecomputedTsvAligner`]:
    /// diamond/blastp emit 0–100, mmseqs2 emits 0–1.
    pub fn from_manifest(
        manifest: GspaManifest,
        target_genome_id: impl Into<String>,
        coverage_is_fraction: bool,
    ) -> Self {
        Self {
            manifest,
            target_genome_id: target_genome_id.into(),
            coverage_is_fraction,
            cache: OnceLock::new(),
        }
    }

    pub fn new(
        manifest_root: impl Into<PathBuf>,
        target_genome_id: impl Into<String>,
        coverage_is_fraction: bool,
    ) -> Result<Self, AlignError> {
        let m = GspaManifest::load(manifest_root)?;
        Ok(Self::from_manifest(m, target_genome_id, coverage_is_fraction))
    }

    fn load_rep_hits(&self) -> Result<&Vec<Hit>, AlignError> {
        if let Some(h) = self.cache.get() {
            return Ok(h);
        }
        let mut all = Vec::new();
        for path in &self.manifest.alignment_files {
            let f = File::open(path).map_err(|e| io_err(path, e))?;
            let rdr = BufReader::new(f);
            let hits = parse_tsv(rdr, self.coverage_is_fraction)?;
            all.extend(hits);
        }
        let _ = self.cache.set(all);
        Ok(self.cache.get().unwrap())
    }
}

impl Aligner for GspaRunAligner {
    fn name(&self) -> &'static str {
        "gspa-run"
    }

    fn align(
        &self,
        _query_fasta: &Path,
        _target_fasta: &Path,
        _opts: &AlignOpts,
    ) -> Result<Vec<Hit>, AlignError> {
        let rep_hits = self.load_rep_hits()?;
        let members = build_member_index(&self.manifest.rep_to_members, &self.target_genome_id);
        let mut out = Vec::with_capacity(rep_hits.len());
        for h in rep_hits {
            let Some(expanded_ids) = members.get(h.qseqid.as_str()) else {
                continue;
            };
            for mid in expanded_ids {
                let mut copy = h.clone();
                copy.qseqid = (*mid).clone();
                out.push(copy);
            }
        }
        tracing::debug!(
            genome = %self.target_genome_id,
            reps_with_hits = out_reps(rep_hits),
            expanded = out.len(),
            "gspa-run: fan-out complete"
        );
        Ok(out)
    }
}

fn out_reps(hits: &[Hit]) -> usize {
    let mut seen = std::collections::HashSet::new();
    for h in hits {
        seen.insert(h.qseqid.as_str());
    }
    seen.len()
}

type MemberIndex<'a> = HashMap<&'a str, Vec<&'a String>>;

fn build_member_index<'a>(
    rep_to_members: &'a HashMap<String, Vec<(String, String)>>,
    target_genome_id: &str,
) -> MemberIndex<'a> {
    let mut idx: MemberIndex<'a> = HashMap::with_capacity(rep_to_members.len());
    for (rep, members) in rep_to_members {
        let mut v: Vec<&String> = Vec::new();
        for (m, g) in members {
            if g == target_genome_id {
                v.push(m);
            }
        }
        if !v.is_empty() {
            idx.insert(rep.as_str(), v);
        }
    }
    idx
}

// --- Parsers ---

fn parse_clusters(
    path: &Path,
) -> Result<HashMap<String, Vec<(String, String)>>, AlignError> {
    let f = File::open(path).map_err(|e| io_err(path, e))?;
    let rdr = BufReader::new(f);
    let mut out: HashMap<String, Vec<(String, String)>> = HashMap::new();
    for (i, line) in rdr.lines().enumerate() {
        let line = line.map_err(|e| io_err(path, e))?;
        let l = line.trim_end_matches('\r');
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let mut it = l.splitn(3, '\t');
        let rep = it.next().unwrap_or("").trim().to_string();
        let member = it.next().unwrap_or("").trim().to_string();
        let genome = it.next().unwrap_or("").trim().to_string();
        if rep.is_empty() || member.is_empty() || genome.is_empty() {
            return Err(AlignError::TsvParse {
                line: (i + 1) as u64,
                msg: format!(
                    "clusters.tsv row needs 3 tab-separated columns (rep, member, genome); got `{l}`"
                ),
            });
        }
        out.entry(rep).or_default().push((member, genome));
    }
    Ok(out)
}

fn parse_genomes(path: &Path) -> Result<Vec<GenomeRow>, AlignError> {
    let f = File::open(path).map_err(|e| io_err(path, e))?;
    let rdr = BufReader::new(f);
    let mut out = Vec::new();
    for (i, line) in rdr.lines().enumerate() {
        let line = line.map_err(|e| io_err(path, e))?;
        let l = line.trim_end_matches('\r');
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        // Accept 2 or 3 tab-separated columns. An optional header row with
        // literal "genome_id" in column 0 is silently skipped.
        let cols: Vec<&str> = l.split('\t').collect();
        if cols.len() < 2 {
            return Err(AlignError::TsvParse {
                line: (i + 1) as u64,
                msg: format!("genomes.tsv row needs ≥2 columns; got `{l}`"),
            });
        }
        if cols[0].trim() == "genome_id" {
            continue;
        }
        let abundance = if cols.len() >= 3 && !cols[2].trim().is_empty() {
            Some(cols[2].trim().parse::<f64>().map_err(|_| AlignError::TsvParse {
                line: (i + 1) as u64,
                msg: format!("genomes.tsv abundance `{}` is not a float", cols[2]),
            })?)
        } else {
            None
        };
        out.push(GenomeRow {
            id: cols[0].trim().to_string(),
            fasta: PathBuf::from(cols[1].trim()),
            abundance,
        });
    }
    Ok(out)
}

fn collect_alignment_files(dir: &Path) -> Result<Vec<PathBuf>, AlignError> {
    if !dir.is_dir() {
        return Err(AlignError::BadArg(format!(
            "gspa-run alignment directory `{}` missing",
            dir.display()
        )));
    }
    let mut out = Vec::new();
    for entry in std::fs::read_dir(dir).map_err(|e| io_err(dir, e))? {
        let e = entry.map_err(|e| io_err(dir, e))?;
        let p = e.path();
        if !p.is_file() {
            continue;
        }
        let ext = p.extension().and_then(|s| s.to_str()).unwrap_or("");
        if matches!(ext, "tsv" | "tab" | "txt") {
            out.push(p);
        }
    }
    out.sort();
    if out.is_empty() {
        return Err(AlignError::BadArg(format!(
            "gspa-run alignment directory `{}` has no .tsv files",
            dir.display()
        )));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write(p: &Path, body: &str) {
        let mut f = std::fs::File::create(p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
    }

    fn mini_manifest(td: &Path) {
        std::fs::create_dir_all(td.join("alignment")).unwrap();
        write(
            &td.join("clusters.tsv"),
            "rep1\trep1\tgA\nrep1\tmemA2\tgA\nrep1\tmemB1\tgB\nrep2\trep2\tgB\n",
        );
        write(&td.join("genomes.tsv"), "gA\t/data/gA.faa\ngB\t/data/gB.faa\t0.3\n");
        write(
            &td.join("alignment/x.tsv"),
            "rep1\t90\t1e-40\t250\t80\trxn1_rep\t1\t200\n\
             rep2\t95\t1e-50\t300\t85\trxn2_rep\t1\t200\n",
        );
    }

    #[test]
    fn loads_manifest() {
        let td = tempfile::tempdir().unwrap();
        mini_manifest(td.path());
        let m = GspaManifest::load(td.path()).unwrap();
        assert_eq!(m.genomes.len(), 2);
        assert_eq!(m.genomes[1].abundance, Some(0.3));
        assert!(m.rep_to_members.contains_key("rep1"));
        assert_eq!(m.alignment_files.len(), 1);
    }

    #[test]
    fn expands_hits_per_genome() {
        let td = tempfile::tempdir().unwrap();
        mini_manifest(td.path());
        let a = GspaRunAligner::new(td.path(), "gA", false).unwrap();
        let hits = a
            .align(Path::new("unused"), Path::new("unused"), &AlignOpts::default())
            .unwrap();
        // gA has two members in rep1 (rep1 itself + memA2) and zero in rep2.
        assert_eq!(hits.len(), 2);
        let ids: Vec<_> = hits.iter().map(|h| h.qseqid.as_str()).collect();
        assert!(ids.contains(&"rep1"));
        assert!(ids.contains(&"memA2"));

        let b = GspaRunAligner::new(td.path(), "gB", false).unwrap();
        let hits_b = b
            .align(Path::new("unused"), Path::new("unused"), &AlignOpts::default())
            .unwrap();
        // gB has one member in rep1 (memB1) and one in rep2 (rep2 itself).
        assert_eq!(hits_b.len(), 2);
    }

    #[test]
    fn unknown_genome_yields_zero_hits() {
        let td = tempfile::tempdir().unwrap();
        mini_manifest(td.path());
        let a = GspaRunAligner::new(td.path(), "nonexistent", false).unwrap();
        let hits = a
            .align(Path::new("x"), Path::new("y"), &AlignOpts::default())
            .unwrap();
        assert!(hits.is_empty());
    }

    #[test]
    fn rejects_short_cluster_row() {
        let td = tempfile::tempdir().unwrap();
        std::fs::create_dir_all(td.path().join("alignment")).unwrap();
        write(&td.path().join("clusters.tsv"), "rep1\tmemA\n");
        write(&td.path().join("genomes.tsv"), "gA\t/x.faa\n");
        write(&td.path().join("alignment/x.tsv"), "");
        let err = GspaManifest::load(td.path()).unwrap_err();
        matches!(err, AlignError::TsvParse { .. });
    }
}
