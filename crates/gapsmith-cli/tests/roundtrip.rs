//! End-to-end smoke test for the `gapsmith` binary:
//!
//! `example-model` → `convert` → `convert` and confirm CBOR output is
//! bit-identical across the round-trip. Guards against regressions in the
//! CLI dispatch layer, format auto-detection, and serde derives.

use std::process::Command;

/// Locate the test binary built by cargo. `CARGO_BIN_EXE_<name>` is set for
/// us as long as the binary target is named `gapsmith` — see Cargo.toml.
fn gapsmith_bin() -> &'static str {
    env!("CARGO_BIN_EXE_gapsmith")
}

#[test]
fn example_model_cbor_json_cbor_roundtrip() {
    let tmp = tempfile::tempdir().unwrap();
    let cbor1 = tmp.path().join("toy.gmod.cbor");
    let json = tmp.path().join("toy.json");
    let cbor2 = tmp.path().join("toy_roundtrip.gmod.cbor");

    let ok = Command::new(gapsmith_bin())
        .args(["example-model", cbor1.to_str().unwrap()])
        .status()
        .unwrap()
        .success();
    assert!(ok, "example-model failed");

    let ok = Command::new(gapsmith_bin())
        .args([
            "convert",
            cbor1.to_str().unwrap(),
            json.to_str().unwrap(),
            "--pretty",
        ])
        .status()
        .unwrap()
        .success();
    assert!(ok, "convert cbor->json failed");

    let ok = Command::new(gapsmith_bin())
        .args(["convert", json.to_str().unwrap(), cbor2.to_str().unwrap()])
        .status()
        .unwrap()
        .success();
    assert!(ok, "convert json->cbor failed");

    let a = std::fs::read(&cbor1).unwrap();
    let b = std::fs::read(&cbor2).unwrap();
    assert_eq!(a, b, "CBOR round-trip not bit-identical");

    // Also confirm the JSON parses and reports a sensible shape.
    let txt = std::fs::read_to_string(&json).unwrap();
    assert!(txt.contains("\"toy_ecoli\""));
    assert!(txt.contains("\"cpd00001\""));
    assert!(txt.contains("\"EX_cpd00007_e0\""));
}
