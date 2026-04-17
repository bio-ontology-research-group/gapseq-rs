//! End-to-end parity test: run both real gapsmith find-transport and
//! gapsmith find-transport on `toy/ecore.faa` and compare the emitted
//! `*-Transporter.tbl` for matching distinct TC ids + row counts.
//!
//! Skipped when Rscript / blastp / gapseq checkout / gapsmith release
//! binary are missing.

use std::path::{Path, PathBuf};
use std::process::Command;

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

fn gapseq_root() -> PathBuf {
    std::env::var("GAPSEQ_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/opt/gapseq"))
}

fn gapseq_rs_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/release/gapsmith")
}

fn prerequisites_ok(root: &Path) -> bool {
    which("blastp").is_some()
        && which("makeblastdb").is_some()
        && which("Rscript").is_some()
        && root.join("gapseq").is_file()
        && gapseq_rs_bin().is_file()
}

fn decompress_ecore(root: &Path, out: &Path) -> std::io::Result<()> {
    let bytes = std::fs::read(root.join("toy/ecore.faa.gz"))?;
    use std::io::Read;
    let mut s = String::new();
    flate2::read::GzDecoder::new(&bytes[..]).read_to_string(&mut s)?;
    std::fs::write(out, s)?;
    Ok(())
}

#[test]
fn ecore_transporter_parity() {
    let root = gapseq_root();
    if !prerequisites_ok(&root) {
        eprintln!("SKIP: missing prerequisites (blastp, Rscript, gapseq, or release binary)");
        return;
    }
    let td = tempfile::tempdir().unwrap();
    let genome = td.path().join("ecore.faa");
    decompress_ecore(&root, &genome).unwrap();

    let real_out = td.path().join("real");
    let rs_out = td.path().join("rs");
    std::fs::create_dir_all(&real_out).unwrap();
    std::fs::create_dir_all(&rs_out).unwrap();

    // Real gapseq.
    let status = Command::new(root.join("gapseq"))
        .arg("find-transport")
        .args(["-b", "50", "-c", "50"])
        .args(["-f", real_out.to_str().unwrap()])
        .arg(&genome)
        .status()
        .unwrap();
    assert!(status.success(), "real gapsmith find-transport failed");

    // gapsmith.
    let status = Command::new(gapseq_rs_bin())
        .args(["--data-dir", root.join("dat").to_str().unwrap()])
        .arg("find-transport")
        .args(["-A", "blastp", "-b", "50", "-c", "50"])
        .args(["-o", rs_out.to_str().unwrap()])
        .args(["-u", "rs"])
        .arg(&genome)
        .status()
        .unwrap();
    assert!(status.success(), "gapsmith find-transport failed");

    let real_path = real_out.join("ecore-Transporter.tbl");
    let rs_path = rs_out.join("ecore-rs-Transporter.tbl");
    let real = std::fs::read_to_string(&real_path).unwrap();
    let rs = std::fs::read_to_string(&rs_path).unwrap();

    let real_rows: Vec<&str> = real.lines().skip_while(|l| l.starts_with('#')).skip(1).collect();
    let rs_rows: Vec<&str> = rs.lines().skip(1).collect();

    assert_eq!(
        real_rows.len(),
        rs_rows.len(),
        "row count differs: real={} vs rs={}",
        real_rows.len(),
        rs_rows.len()
    );

    let tc_of = |row: &&str| row.split('\t').nth(1).unwrap_or("").to_string();
    let real_tcs: std::collections::HashSet<String> = real_rows.iter().map(tc_of).collect();
    let rs_tcs: std::collections::HashSet<String> = rs_rows.iter().map(tc_of).collect();
    assert_eq!(
        real_tcs,
        rs_tcs,
        "distinct TC sets differ: real-only={:?} rs-only={:?}",
        real_tcs.difference(&rs_tcs).take(10).collect::<Vec<_>>(),
        rs_tcs.difference(&real_tcs).take(10).collect::<Vec<_>>()
    );
}
