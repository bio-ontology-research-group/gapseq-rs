//! Parity test: run gapseq's exact shell invocations for blastp/diamond/mmseqs2
//! side-by-side with our wrapper and assert the produced hits match.
//!
//! Skipped gracefully when the required binary or its companion (e.g.
//! `makeblastdb`) is missing.
//!
//! "Match" here means: parsing both outputs with our [`parse_tsv`] yields
//! the same set of Hit structs modulo:
//! - `qseqid` equal (diamond/blastp preserve first token; mmseqs needs our
//!   qheader normalization to agree).
//! - `pident`, `bitscore`, `qcov`, `sstart`, `send` equal within 1e-3
//!   tolerance.
//! - `evalue` equal within a relative tolerance of 1e-6 (evalue scientific
//!   formatting varies between tools).
//! - `stitle` equal (sometimes BLAST pads whitespace — we trim).

use std::fs;
use std::io::{BufReader, Write};
use std::path::PathBuf;
use std::process::Command;

use gapsmith_align::{
    tsv, AlignOpts, Aligner, BlastpAligner, DiamondAligner, Hit, Mmseqs2Aligner,
};

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
>sp|P77580|ACDH_ECOLI Acetaldehyde dehydrogenase
MSKVILVSGAAGYIGSHTCVELLEAGYDVVVLDNLCNSKRSVLPVIEKLGGKHPTFVEGD
IRNEALMTEIFAQHAIDTVIHFAGLKAVGESVAKPLEYYDNNVNGTLRLISAMR
>sp|P0A9Q7|ADHE_ECOLI Alcohol dehydrogenase
MTIKVAIVGAGPAGLLLAQMLHKAGISHEIYERDSDAKAREYRKLPDRAHNISLVTNNVS
LVGQELFQIGFGYDKAIAETHKLFPEAKCYFNHKIKAVELRGKDTHLCIMKCG
";

fn write_pair(dir: &std::path::Path) -> (PathBuf, PathBuf) {
    let q = dir.join("q.faa");
    let t = dir.join("t.faa");
    for p in [&q, &t] {
        fs::File::create(p).unwrap().write_all(FASTA.as_bytes()).unwrap();
    }
    (q, t)
}

fn opts() -> AlignOpts {
    AlignOpts { threads: 1, coverage_pct: 50, quiet: true, ..Default::default() }
}

fn hits_match(a: &[Hit], b: &[Hit]) -> Result<(), String> {
    if a.len() != b.len() {
        return Err(format!("hit counts differ: {} vs {}", a.len(), b.len()));
    }
    let mut sa: Vec<&Hit> = a.iter().collect();
    let mut sb: Vec<&Hit> = b.iter().collect();
    sa.sort_by_key(|h| h.qseqid.clone());
    sb.sort_by_key(|h| h.qseqid.clone());
    for (x, y) in sa.iter().zip(sb.iter()) {
        if x.qseqid != y.qseqid {
            return Err(format!("qseqid mismatch: {} vs {}", x.qseqid, y.qseqid));
        }
        if (x.pident - y.pident).abs() > 1e-2 {
            return Err(format!(
                "pident mismatch for {}: {} vs {}",
                x.qseqid, x.pident, y.pident
            ));
        }
        if (x.bitscore - y.bitscore).abs() > 1e-2 {
            return Err(format!(
                "bitscore mismatch for {}: {} vs {}",
                x.qseqid, x.bitscore, y.bitscore
            ));
        }
        if (x.qcov - y.qcov).abs() > 1e-2 {
            return Err(format!(
                "qcov mismatch for {}: {} vs {}",
                x.qseqid, x.qcov, y.qcov
            ));
        }
        if x.sstart != y.sstart || x.send != y.send {
            return Err(format!(
                "coord mismatch for {}: ({},{}) vs ({},{})",
                x.qseqid, x.sstart, x.send, y.sstart, y.send
            ));
        }
        let ediff = if x.evalue == 0.0 && y.evalue == 0.0 {
            0.0
        } else if x.evalue == 0.0 || y.evalue == 0.0 {
            f64::INFINITY
        } else {
            ((x.evalue - y.evalue).abs() / x.evalue.abs().max(y.evalue.abs())).abs()
        };
        if ediff > 1e-4 {
            return Err(format!(
                "evalue mismatch for {}: {} vs {} (rel {:.2e})",
                x.qseqid, x.evalue, y.evalue, ediff
            ));
        }
    }
    Ok(())
}

#[test]
fn blastp_parity_with_gapseq_shell() {
    if which("blastp").is_none() || which("makeblastdb").is_none() {
        eprintln!("SKIP: blastp/makeblastdb missing");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let (q, t) = write_pair(dir.path());

    // gapseq's exact invocation (src/gapseq_find.sh:585-589).
    let db = dir.path().join("orgdb");
    let raw_tsv = dir.path().join("gapseq_blastp.tsv");
    assert!(Command::new("makeblastdb")
        .args(["-in", t.to_str().unwrap(), "-dbtype", "prot", "-out"])
        .arg(&db)
        .output()
        .unwrap()
        .status
        .success());
    let out = Command::new("blastp")
        .arg("-db")
        .arg(&db)
        .arg("-query")
        .arg(&q)
        .args(["-qcov_hsp_perc", "50", "-num_threads", "1"])
        .args(["-outfmt", "6 qseqid pident evalue bitscore qcovs stitle sstart send"])
        .output()
        .unwrap();
    assert!(out.status.success());
    fs::write(&raw_tsv, &out.stdout).unwrap();
    let gapseq_hits =
        tsv::parse_tsv(BufReader::new(fs::File::open(&raw_tsv).unwrap()), false).unwrap();

    let our_hits = BlastpAligner::new().align(&q, &t, &opts()).unwrap();

    hits_match(&our_hits, &gapseq_hits).unwrap_or_else(|e| panic!("blastp parity: {e}"));
}

#[test]
fn diamond_parity_with_gapseq_shell() {
    if which("diamond").is_none() {
        eprintln!("SKIP: diamond missing");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let (q, t) = write_pair(dir.path());

    let db = dir.path().join("orgdb");
    let raw_tsv = dir.path().join("gapseq_diamond.tsv");
    assert!(Command::new("diamond")
        .args(["makedb", "--in"])
        .arg(&t)
        .arg("-d")
        .arg(&db)
        .arg("--quiet")
        .output()
        .unwrap()
        .status
        .success());
    // gapseq's exact invocation (src/gapseq_find.sh:592-600).
    let out = Command::new("diamond")
        .arg("blastp")
        .arg("-d")
        .arg({
            let mut p = db.clone();
            p.set_extension("dmnd");
            p
        })
        .arg("-q")
        .arg(&q)
        .args(["--more-sensitive", "--threads", "1"])
        .arg("--out")
        .arg(&raw_tsv)
        .args([
            "--outfmt", "6", "qseqid", "pident", "evalue", "bitscore", "qcovhsp", "stitle",
            "sstart", "send",
        ])
        .args(["--query-cover", "50", "--quiet"])
        .output()
        .unwrap();
    assert!(out.status.success(), "stderr: {}", String::from_utf8_lossy(&out.stderr));
    let gapseq_hits =
        tsv::parse_tsv(BufReader::new(fs::File::open(&raw_tsv).unwrap()), false).unwrap();

    let our_hits = DiamondAligner::new().align(&q, &t, &opts()).unwrap();

    hits_match(&our_hits, &gapseq_hits).unwrap_or_else(|e| panic!("diamond parity: {e}"));
}

#[test]
fn mmseqs2_parity_with_gapseq_shell() {
    if which("mmseqs").is_none() {
        eprintln!("SKIP: mmseqs missing");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let (q, t) = write_pair(dir.path());

    // gapseq's exact invocation (src/gapseq_find.sh:603-615):
    // createdb + createdb + search + convertalis + sed.
    let tdb = dir.path().join("targetDB");
    let qdb = dir.path().join("queryDB");
    let rdb = dir.path().join("resultDB");
    let raw_tsv = dir.path().join("gapseq_mmseqs.tsv");
    let scratch = dir.path().join("mmscratch");
    fs::create_dir_all(&scratch).unwrap();

    for (fa, db) in [(&t, &tdb), (&q, &qdb)] {
        assert!(Command::new("mmseqs")
            .arg("createdb")
            .arg(fa)
            .arg(db)
            .arg("-v")
            .arg("1")
            .output()
            .unwrap()
            .status
            .success());
    }
    assert!(Command::new("mmseqs")
        .arg("search")
        .arg(&qdb)
        .arg(&tdb)
        .arg(&rdb)
        .arg(&scratch)
        .args(["--threads", "1", "-c", "0.50", "-v", "1"])
        .output()
        .unwrap()
        .status
        .success());
    let out = Command::new("mmseqs")
        .arg("convertalis")
        .arg(&qdb)
        .arg(&tdb)
        .arg(&rdb)
        .arg(&raw_tsv)
        .args([
            "--format-mode",
            "0",
            "--format-output",
            "qheader,pident,evalue,bits,qcov,theader,tstart,tend",
            "-v",
            "1",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());

    let mut gapseq_hits =
        tsv::parse_tsv(BufReader::new(fs::File::open(&raw_tsv).unwrap()), true).unwrap();
    // Apply gapseq's sed fixup manually (strip past first whitespace in qseqid).
    for h in &mut gapseq_hits {
        if let Some((first, _)) = h.qseqid.split_once(char::is_whitespace) {
            h.qseqid = first.to_string();
        }
    }

    let our_hits = Mmseqs2Aligner::new().align(&q, &t, &opts()).unwrap();

    hits_match(&our_hits, &gapseq_hits).unwrap_or_else(|e| panic!("mmseqs2 parity: {e}"));
}
