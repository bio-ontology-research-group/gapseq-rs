//! Minimal loader for `dat/taxonomy.tbl` — used to filter pathways by
//! their `taxrange` column.
//!
//! Format (TSV, with header `tax\tsuperkingdom\tgroup`):
//!
//! ```text
//! 2        Bacteria    Bacteria
//! 1224     Bacteria    Pseudomonadota
//! 4751     Eukaryota   Fungi
//! ```
//!
//! `validTaxIds(for: "Bacteria")` returns every `tax` column whose
//! `superkingdom` OR `group` contains the requested tag (case-insensitive,
//! matching `grep -i` in `gapseq_find.sh:555`).

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn valid_tax_ids_for(tax_table: &Path, taxonomy: &str) -> std::io::Result<Vec<String>> {
    let f = File::open(tax_table)?;
    let r = BufReader::new(f);
    let needle = taxonomy.to_ascii_lowercase();
    let mut out: Vec<String> = Vec::new();
    let mut header = true;
    for line in r.lines() {
        let line = line?;
        if header {
            header = false;
            continue;
        }
        let lower = line.to_ascii_lowercase();
        if lower.contains(&needle) {
            if let Some(id) = line.split('\t').next() {
                if !id.trim().is_empty() {
                    out.push(id.trim().to_string());
                }
            }
        }
    }
    Ok(out)
}
