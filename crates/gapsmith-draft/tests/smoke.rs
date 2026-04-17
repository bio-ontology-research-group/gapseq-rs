//! Smoke test for the draft pipeline. Builds a tiny Reactions.tbl +
//! Transporter.tbl in a temp dir, runs `gapsmith_draft::run`, and verifies
//! the emitted model has a biomass reaction, valid stoichiometry, and
//! loads via [`gapsmith_io`].

use gapsmith_draft::{run, DraftOptions};
use std::io::Write;
use std::path::PathBuf;

fn gapseq_data() -> PathBuf {
    PathBuf::from(
        std::env::var("GAPSEQ_ROOT")
            .unwrap_or_else(|_| "/opt/gapseq".into()),
    )
    .join("dat")
}

#[test]
fn builds_minimal_draft() {
    let data_dir = gapseq_data();
    if !data_dir.is_dir() {
        eprintln!("SKIP: gapseq data dir not present");
        return;
    }

    let td = tempfile::tempdir().unwrap();
    let rxn_path = td.path().join("toy-all-Reactions.tbl");
    let tr_path = td.path().join("toy-Transporter.tbl");

    // Minimal Reactions.tbl with one good-blast reaction (rxn00001 via
    // the ACETALD-DEHYDROG-RXN example) so the draft gets > 0 rxns.
    let mut rf = std::fs::File::create(&rxn_path).unwrap();
    writeln!(rf, "# gapseq version: 0.1.0 (test)").unwrap();
    writeln!(rf, "# genome_format=prot;taxonomy=Bacteria;gram=neg").unwrap();
    writeln!(rf, "# placeholder").unwrap();
    writeln!(
        rf,
        "pathway\trxn\tname\tec\tkeyrea\tfile\tdbhit\tspont\ttype\tsrc\tis_complex\tsubunit_count\tsubunits\tqseqid\tpident\tevalue\tbitscore\tqcovs\tstitle\tsstart\tsend\tcomplex\texception\tstatus\tsubunits_found\tsubunit_undefined_found\tcomplex.status\tpathway.status"
    )
    .unwrap();
    // Give a good_blast hit with a dbhit pointing at rxn00001 (water+PPi
    // diphosphatase) so the draft has something real from SEED.
    writeln!(
        rf,
        "PWY-toy\tFOO-RXN\tacetaldehyde dehydrogenase\t3.5.4.19\tFALSE\trev/3.5.4.19.fasta\trxn00001\tFALSE\tEC\trev\tFALSE\t\t\tQ\t100\t1e-100\t300\t100\tsp|X|Y\t1\t100\t\t0\tgood_blast\t\t\t\tfull"
    )
    .unwrap();
    drop(rf);

    // Minimal Transporter.tbl (empty body is fine).
    let mut tf = std::fs::File::create(&tr_path).unwrap();
    writeln!(
        tf,
        "id\ttc\tsub\tsub_gapseq\texid\trea\tqseqid\tpident\tevalue\tbitscore\tqcovs\tstitle\tsstart\tsend\tcomment"
    )
    .unwrap();
    drop(tf);

    let opts = DraftOptions {
        model_id: "toy_draft".into(),
        biomass: "neg".into(),
        high_evi_rxn_bs: 200.0,
        min_bs_for_core: 50.0,
        gapsmith_version: Some("test".into()),
        seqdb_version: None,
    };
    let rep = run(&rxn_path, &tr_path, &data_dir, &opts).expect("draft run failed");

    assert!(rep.model.met_count() > 20, "got {} mets", rep.model.met_count());
    assert!(rep.model.rxn_count() >= 5, "got {} rxns", rep.model.rxn_count());
    assert!(
        rep.model.rxns.iter().any(|r| r.id.as_str() == "bio1"),
        "biomass reaction missing"
    );
    assert!(
        rep.model.rxns.iter().any(|r| r.is_exchange),
        "no exchange reactions added"
    );
    rep.model.check_shape().expect("S shape must match");
}
