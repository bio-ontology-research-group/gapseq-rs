//! Small helper binary: build the `small.fasta` we would use for
//! find-transport and dump its header list so we can diff against
//! gapseq's shell pipeline.

use gapsmith_db::{subex, SubexRow};
use gapsmith_transport::{build_small_fasta, load_tcdb_all};
use std::path::{Path, PathBuf};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let data_dir = PathBuf::from(args.get(1).cloned().unwrap_or_else(|| {
        std::env::var("GAPSMITH_DATA_DIR")
            .expect("pass data-dir as arg 1 or set GAPSMITH_DATA_DIR")
    }));
    let seq_dir = data_dir.join("seq");
    let out_dir = PathBuf::from(args.get(2).cloned().unwrap_or_else(|| {
        "/tmp/ftp_test/inspect".to_string()
    }));
    std::fs::create_dir_all(&out_dir).unwrap();

    let subex_rows: Vec<SubexRow> = subex::load(data_dir.join("subex.tbl")).unwrap();
    let tcdb_all = load_tcdb_all(
        &data_dir.join("tcdb_substrates.tbl"),
        &data_dir.join("tcdb_custom.tbl"),
    )
    .unwrap();
    let refs = [seq_dir.join("tcdb.fasta"), seq_dir.join("transporter.fasta")];
    let refs_paths: Vec<&Path> = refs.iter().map(|p| p.as_path()).collect();

    let result =
        build_small_fasta(&refs_paths, &subex_rows, &tcdb_all, None, &out_dir).unwrap();
    eprintln!(
        "wrote {} headers to {}",
        result.fasta_header_small.len(),
        result.small_fasta.display()
    );
    eprintln!("sub_keys={}", result.sub_key.len());
    eprintln!(
        "tcdb_all has 1.A.11.1.4 = {:?}",
        tcdb_all.get("1.A.11.1.4")
    );
    eprintln!(
        "sub_key has 'ammonia' = {}",
        result.sub_key.contains("ammonia")
    );
    eprintln!(
        "sub_key has 'maltodextrin' = {}",
        result.sub_key.contains("maltodextrin")
    );
    // Count selected headers containing each of several missing TCs.
    for tc in ["1.A.11.1.4", "3.A.1.1.16", "3.A.1.1.34", "2.A.56.1.10"] {
        let n = result
            .fasta_header_small
            .iter()
            .filter(|h| h.contains(tc))
            .count();
        eprintln!("headers with {tc}: {n}");
    }
}
