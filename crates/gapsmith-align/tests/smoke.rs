//! Smoke tests that shell out to blast+ / diamond / mmseqs2 when present.
//!
//! Each test opens a temp dir, writes a two-sequence protein FASTA, and
//! aligns it against itself. Expected: ≥ 2 self-hits (one per sequence)
//! with pident ≈ 100. If the required binary is missing, the test is
//! skipped (`cargo test` reports it as passed but with a logged message).

use gapsmith_align::{
    AlignOpts, Aligner, BlastpAligner, DiamondAligner, Mmseqs2Aligner,
    PrecomputedTsvAligner,
};
use std::io::Write;
use std::path::PathBuf;

fn which(name: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path) {
        let c = dir.join(name);
        if c.is_file() {
            return Some(c);
        }
    }
    None
}

const FASTA: &str = "\
>seqA some protein
MSKVILVSGAAGYIGSHTCVELLEAGYDVVVLDNLCNSKRSVLPVIEKLGGKHPTFVEGD
IRNEALMTEIFAQHAIDTVIHFAGLKAVGESVAKPLEYYDNNVNGTLRLISAMR
>seqB another protein
MTIKVAIVGAGPAGLLLAQMLHKAGISHEIYERDSDAKAREYRKLPDRAHNISLVTNNVS
LVGQELFQIGFGYDKAIAETHKLFPEAKCYFNHKIKAVELRGKDTHLCIMKCG
";

fn write_fasta(dir: &std::path::Path) -> PathBuf {
    let p = dir.join("pair.faa");
    let mut f = std::fs::File::create(&p).unwrap();
    f.write_all(FASTA.as_bytes()).unwrap();
    p
}

#[test]
fn blastp_self_alignment() {
    if which("blastp").is_none() || which("makeblastdb").is_none() {
        eprintln!("SKIP: blastp/makeblastdb not on PATH");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let fa = write_fasta(dir.path());
    let hits = BlastpAligner::new()
        .align(
            &fa,
            &fa,
            &AlignOpts { threads: 1, coverage_pct: 50, quiet: true, ..Default::default() },
        )
        .unwrap();
    assert!(hits.len() >= 2, "expected ≥2 self-hits, got {}", hits.len());
    for h in &hits {
        assert!(h.pident > 90.0, "self-hit should be ~100%, got {}", h.pident);
    }
}

#[test]
fn diamond_self_alignment() {
    if which("diamond").is_none() {
        eprintln!("SKIP: diamond not on PATH");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let fa = write_fasta(dir.path());
    let hits = DiamondAligner::new()
        .align(
            &fa,
            &fa,
            &AlignOpts { threads: 1, coverage_pct: 50, quiet: true, ..Default::default() },
        )
        .unwrap();
    assert!(hits.len() >= 2, "expected ≥2 self-hits, got {}", hits.len());
    for h in &hits {
        assert!(h.pident > 90.0, "self-hit should be ~100%, got {}", h.pident);
        assert!(h.qcov > 0.0 && h.qcov <= 100.0, "qcov should be on 0-100 scale, got {}", h.qcov);
    }
}

#[test]
fn mmseqs2_self_alignment() {
    if which("mmseqs").is_none() {
        eprintln!("SKIP: mmseqs not on PATH");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let fa = write_fasta(dir.path());
    let hits = Mmseqs2Aligner::new()
        .align(
            &fa,
            &fa,
            &AlignOpts { threads: 1, coverage_pct: 50, quiet: true, ..Default::default() },
        )
        .unwrap();
    assert!(hits.len() >= 2, "expected ≥2 self-hits, got {}", hits.len());
    for h in &hits {
        assert!(h.pident > 90.0, "self-hit should be ~100%, got {}", h.pident);
        // We rescale fraction→percentage in tsv.rs; verify.
        assert!(h.qcov > 0.0 && h.qcov <= 100.0, "qcov should be on 0-100 scale, got {}", h.qcov);
    }
}

#[test]
fn precomputed_aligner_reads_tsv_we_generate() {
    // Generate a TSV with diamond if possible, else blastp — whichever
    // binary exists — then feed that TSV into the precomputed aligner and
    // confirm it reports the same number of hits.
    let dir = tempfile::tempdir().unwrap();
    let fa = write_fasta(dir.path());

    let reference_hits: Vec<gapsmith_align::Hit>;
    let aligner_name: &str;
    if which("diamond").is_some() {
        aligner_name = "diamond";
        reference_hits = DiamondAligner::new()
            .align(
                &fa,
                &fa,
                &AlignOpts { threads: 1, coverage_pct: 50, quiet: true, ..Default::default() },
            )
            .unwrap();
    } else if which("blastp").is_some() && which("makeblastdb").is_some() {
        aligner_name = "blastp";
        reference_hits = BlastpAligner::new()
            .align(
                &fa,
                &fa,
                &AlignOpts { threads: 1, coverage_pct: 50, quiet: true, ..Default::default() },
            )
            .unwrap();
    } else {
        eprintln!("SKIP: neither diamond nor blastp on PATH");
        return;
    }

    // Write reference_hits as an 8-column TSV and round-trip via
    // PrecomputedTsvAligner.
    let tsv = dir.path().join("precomputed.tsv");
    {
        let mut f = std::fs::File::create(&tsv).unwrap();
        for h in &reference_hits {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                h.qseqid, h.pident, h.evalue, h.bitscore, h.qcov, h.stitle, h.sstart, h.send
            )
            .unwrap();
        }
    }
    let a = PrecomputedTsvAligner::new_percentage(&tsv);
    let hits =
        a.align(std::path::Path::new("x"), std::path::Path::new("y"), &AlignOpts::default()).unwrap();
    assert_eq!(hits.len(), reference_hits.len());
    eprintln!("precomputed round-trip OK via {aligner_name}: {} hits", hits.len());
}
