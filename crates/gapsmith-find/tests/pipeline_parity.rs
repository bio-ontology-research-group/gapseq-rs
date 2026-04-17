//! End-to-end parity test: run both real gapseq and gapsmith on the
//! same input and compare `-Pathways.tbl` byte-for-byte.
//!
//! Requires: `Rscript` + gapseq checkout on `$GAPSEQ_ROOT` (default
//! `/opt/gapseq`), `blastp`/`makeblastdb` on
//! PATH. Skipped when any is missing.
//!
//! What the test covers:
//!
//! - Single-pathway (`-p PWY-6587`): pathway selection, reference-FASTA
//!   resolution (including `user/` fallback to the built-in
//!   `dat/seq/<tax>/user/` directory), hit classification, and
//!   completeness scoring on 3 reactions with 1 blast hit.
//! - Multi-pathway shorthand (`-p amino`): hierarchy-mode pathway
//!   selection, superpathway filter, `metacyc+custom` database merge
//!   with custom-row-wins semantics.

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
    which("Rscript").is_some()
        && which("blastp").is_some()
        && which("makeblastdb").is_some()
        && root.join("gapseq").is_file()
        && gapseq_rs_bin().is_file()
}

fn prepare_seqdb(root: &Path, dest: &Path) -> std::io::Result<()> {
    use std::io::Write;
    std::fs::create_dir_all(dest.join("Bacteria/rev"))?;
    std::fs::create_dir_all(dest.join("Bacteria/unrev"))?;
    std::fs::create_dir_all(dest.join("Bacteria/rxn"))?;
    std::fs::create_dir_all(dest.join("Bacteria/user"))?;

    let mut f = std::fs::File::create(dest.join("Bacteria/version_seqDB.json"))?;
    writeln!(
        f,
        r#"{{"version": ["test"], "date": ["2026-04-16"], "zenodoID": ["0"]}}"#
    )?;

    // Extract two proteins from toy/ecore into rev/EC.fasta files.
    let ecore = read_ecore(root)?;
    let p77580 = extract_fasta(&ecore, "sp|P77580|ACDH_ECOLI");
    let p0a9q7 = extract_fasta(&ecore, "sp|P0A9Q7|ADHE_ECOLI");

    let mut rev1 = std::fs::File::create(dest.join("Bacteria/rev/1.2.1.10.fasta"))?;
    rev1.write_all(p77580.as_bytes())?;
    rev1.write_all(p0a9q7.as_bytes())?;
    drop(rev1);

    let mut rev2 = std::fs::File::create(dest.join("Bacteria/rev/1.1.1.1.fasta"))?;
    rev2.write_all(p0a9q7.as_bytes())?;
    Ok(())
}

fn read_ecore(root: &Path) -> std::io::Result<String> {
    let gz = root.join("toy/ecore.faa.gz");
    let bytes = std::fs::read(&gz)?;
    use std::io::Read;
    let mut s = String::new();
    flate2::read::GzDecoder::new(&bytes[..]).read_to_string(&mut s)?;
    Ok(s)
}

fn extract_fasta(ecore: &str, header_marker: &str) -> String {
    let mut out = String::new();
    let mut capturing = false;
    for line in ecore.lines() {
        if line.starts_with('>') {
            capturing = line.contains(header_marker);
        }
        if capturing {
            out.push_str(line);
            out.push('\n');
        }
    }
    out
}

fn decompress_ecore_to(root: &Path, out: &Path) -> std::io::Result<()> {
    let ecore = read_ecore(root)?;
    std::fs::write(out, ecore)?;
    Ok(())
}

fn read_minus_metadata(path: &Path) -> std::io::Result<String> {
    let text = std::fs::read_to_string(path)?;
    Ok(text.lines().skip(3).collect::<Vec<_>>().join("\n"))
}

/// Build a case run under `workdir` and return the per-pathway pathway
/// file from each tool so the caller can diff them.
fn run_case(
    root: &Path,
    workdir: &Path,
    keyword: &str,
    suffix: &str,
) -> std::io::Result<(PathBuf, PathBuf)> {
    std::fs::create_dir_all(workdir)?;
    let seqdb = workdir.join("seqdb");
    prepare_seqdb(root, &seqdb)?;
    let genome = workdir.join("ecore.faa");
    decompress_ecore_to(root, &genome)?;

    let real_out = workdir.join("real_out");
    let rs_out = workdir.join("rs_out");
    std::fs::create_dir_all(&real_out)?;
    std::fs::create_dir_all(&rs_out)?;

    let real_status = Command::new(root.join("gapseq"))
        .arg("find")
        .args(["-p", keyword])
        .args(["-t", "Bacteria"])
        .arg("-O")
        .args(["-b", "200"])
        .args(["-c", "50"])
        .args(["-D", seqdb.to_str().unwrap()])
        .args(["-f", real_out.to_str().unwrap()])
        .arg(&genome)
        .status()?;
    assert!(real_status.success(), "real gapsmith find failed");

    let rs_status = Command::new(gapseq_rs_bin())
        .args(["--data-dir", root.join("dat").to_str().unwrap()])
        .args(["--seq-dir", seqdb.to_str().unwrap()])
        .arg("find")
        .args(["-p", keyword])
        .args(["-t", "Bacteria"])
        .args(["-A", "blastp"])
        .args(["-b", "200"])
        .args(["-c", "50"])
        .args(["-o", rs_out.to_str().unwrap()])
        .args(["-u", suffix])
        .arg(&genome)
        .status()?;
    assert!(rs_status.success(), "gapsmith find failed");

    let real_pwy = real_out.join(format!("ecore-{suffix}-Pathways.tbl"));
    let rs_pwy = rs_out.join(format!("ecore-{suffix}-Pathways.tbl"));
    Ok((real_pwy, rs_pwy))
}

#[test]
fn pwy_6587_pathways_byte_identical() {
    let root = gapseq_root();
    if !prerequisites_ok(&root) {
        eprintln!(
            "SKIP: missing one of Rscript, blastp, makeblastdb, gapseq checkout, or gapsmith release binary at {}",
            gapseq_rs_bin().display()
        );
        return;
    }
    let td = tempfile::tempdir().unwrap();
    let (real, rs) = run_case(&root, td.path(), "PWY-6587", "PWY-6587").unwrap();
    let real_body = read_minus_metadata(&real).unwrap();
    let rs_body = std::fs::read_to_string(&rs).unwrap();
    assert_eq!(real_body.trim_end(), rs_body.trim_end(), "PWY-6587 Pathways.tbl differs");
}

#[test]
fn amino_pathways_byte_identical() {
    let root = gapseq_root();
    if !prerequisites_ok(&root) {
        eprintln!("SKIP: missing prerequisites (see other parity test)");
        return;
    }
    let td = tempfile::tempdir().unwrap();
    let (real, rs) = run_case(&root, td.path(), "amino", "amino").unwrap();
    let real_body = read_minus_metadata(&real).unwrap();
    let rs_body = std::fs::read_to_string(&rs).unwrap();
    assert_eq!(real_body.trim_end(), rs_body.trim_end(), "amino Pathways.tbl differs");
}
