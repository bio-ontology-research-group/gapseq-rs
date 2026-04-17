//! Point-by-point parity test of `complex::detect_subunits` against the
//! real R `complex_detection.R` implementation.
//!
//! Skipped when `Rscript` is not on PATH. Otherwise: feed the same batch of
//! FASTA descriptors to both implementations, assert identical results
//! element-by-element.

use gapsmith_db::{ComplexSubunitTable, DbError};
use gapsmith_find::complex::detect_subunits;
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

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

fn load_dict() -> Result<ComplexSubunitTable, DbError> {
    let path = PathBuf::from(
        std::env::var("GAPSEQ_ROOT")
            .unwrap_or_else(|_| "/opt/gapseq".into()),
    )
    .join("dat/complex_subunit_dict.tsv");
    ComplexSubunitTable::load(&path)
}

fn run_r(batch: &[(String, Vec<String>)]) -> Vec<Vec<Option<String>>> {
    let mut tmpfile = tempfile::NamedTempFile::new().unwrap();
    for (rxn, descs) in batch {
        writeln!(tmpfile, "{rxn}\t{}", descs.join("|||")).unwrap();
    }
    tmpfile.flush().unwrap();
    let script = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tools/r_complex_detection.R");
    let mut cmd = Command::new("Rscript");
    cmd.arg(script).arg(tmpfile.path()).stdout(Stdio::piped()).stderr(Stdio::piped());
    let out = cmd.output().expect("Rscript failed to spawn");
    assert!(
        out.status.success(),
        "Rscript failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = String::from_utf8_lossy(&out.stdout).to_string();
    text.lines()
        .map(|line| {
            line.split("|||")
                .map(|s| if s.is_empty() { None } else { Some(s.to_string()) })
                .collect()
        })
        .collect()
}

fn cases() -> Vec<(String, Vec<String>)> {
    vec![
        (
            "RXN-A".into(),
            vec![
                "DNA polymerase subunit alpha OS=Ecoli".into(),
                "DNA polymerase subunit beta OS=Ecoli".into(),
                "DNA polymerase subunit gamma OS=Ecoli".into(),
                "Unrelated protein OS=X".into(),
            ],
        ),
        (
            "RXN-B".into(),
            vec![
                "ATP synthase alpha chain".into(),
                "ATP synthase beta chain".into(),
                "ATP synthase gamma chain".into(),
                "ATP synthase delta chain".into(),
                "ATP synthase epsilon chain".into(),
            ],
        ),
        (
            "RXN-C".into(),
            vec![
                "Ribulose bisphosphate carboxylase large chain".into(),
                "Ribulose bisphosphate carboxylase small chain".into(),
            ],
        ),
        (
            "RXN-D".into(),
            vec![
                "Cytochrome c oxidase subunit 1".into(),
                "Cytochrome c oxidase subunit 2".into(),
                "Cytochrome c oxidase subunit 3".into(),
            ],
        ),
        (
            "RXN-E".into(),
            vec![
                "DNA-directed RNA polymerase subunit alpha".into(),
                "DNA-directed RNA polymerase subunit beta".into(),
                "DNA-directed RNA polymerase subunit beta'".into(),
                "DNA-directed RNA polymerase subunit omega".into(),
                "DNA-directed RNA polymerase subunit sigma-70".into(),
            ],
        ),
        // Edge: ≤20% hit → everything blanked.
        (
            "RXN-F".into(),
            {
                let mut v = vec!["random protein".to_string(); 9];
                v.push("DNA polymerase subunit alpha".into());
                v
            },
        ),
        // Edge: 66%+ numbered subunits — non-numbered gets dropped only if
        // its subunit label is non-numeric *after* numeral mapping, which
        // would effectively never happen. Still good as a sanity case.
        (
            "RXN-G".into(),
            (1..=7)
                .map(|i| format!("DNA pol subunit {i}"))
                .chain(std::iter::once("Ribonuclease alpha chain".into()))
                .collect(),
        ),
        // Edge: latin numerals.
        (
            "RXN-H".into(),
            vec![
                "Complex V subunit I".into(),
                "Complex V subunit II".into(),
                "Complex V subunit III".into(),
                "Complex V subunit IV".into(),
                "Complex V subunit V".into(),
            ],
        ),
        // Edge: descriptor with quotes + OS= tail.
        (
            "RXN-I".into(),
            vec![
                "ATPase F1 alpha-subunit OS=Ecoli".into(),
                "ATPase F1 beta-subunit OS=Ecoli".into(),
            ],
        ),
    ]
}

#[test]
fn r_vs_rust_parity() {
    if which("Rscript").is_none() {
        eprintln!("SKIP: Rscript not on PATH");
        return;
    }
    let dict = match load_dict() {
        Ok(d) => d,
        Err(e) => {
            eprintln!("SKIP: complex_subunit_dict.tsv: {e}");
            return;
        }
    };
    let batch = cases();
    let r_out = run_r(&batch);
    assert_eq!(r_out.len(), batch.len());

    for (i, (rxn, descs)) in batch.iter().enumerate() {
        let refs: Vec<&str> = descs.iter().map(|s| s.as_str()).collect();
        let rs_out = detect_subunits(rxn, &refs, &dict);
        let r = &r_out[i];
        assert_eq!(
            rs_out.len(),
            r.len(),
            "length mismatch for {rxn}: rs={:?} r={:?}",
            rs_out,
            r
        );
        for (j, (rust_v, r_v)) in rs_out.iter().zip(r.iter()).enumerate() {
            assert_eq!(
                rust_v, r_v,
                "case {rxn} [{j}] descriptor `{}`: rust={:?} r={:?}",
                descs[j], rust_v, r_v
            );
        }
    }
}
