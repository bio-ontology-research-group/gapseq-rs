//! Data loaders for the transport-specific tables.
//!
//! The SEED transporter table (`dat/seed_transporter.tbl` +
//! `dat/seed_transporter_custom.tbl`) has 5 columns:
//! `id, name, type, exmet, exmetnames`. We merge the two files and strip
//! the `[e0]`-style compartment suffix from `exmet` to produce a
//! `(type, metid) -> reaction ids` mapping.
//!
//! `tcdb_all` is the union of `dat/tcdb_substrates.tbl` and
//! `dat/tcdb_custom.tbl` — both have `<TC_id>\t<substrate list>` as bare
//! (header-less) rows. Substrates are pipe-separated pairs of
//! `CHEBI:<id>;<name>`. We keep only the semi-colon-joined human name
//! portion after cleaning the regex junk (charges, water/ion/etc.
//! filler words) that the R code strips in
//! `analyse_alignments_transport.R:24-34`.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeedTransporterRow {
    pub id: String,
    pub name: String,
    /// e.g. `"1.Channels and pores"`.
    pub transporter_type: String,
    /// Extracellular metabolite, `<cpd>[eN]` format.
    pub exmet: String,
    pub exmetnames: String,
}

/// Load + concatenate the SEED transporter tables.
pub fn load_seed_transporter(
    seed: &Path,
    custom: &Path,
) -> std::io::Result<Vec<SeedTransporterRow>> {
    let mut out = Vec::new();
    for path in [seed, custom] {
        let f = File::open(path)?;
        let r = BufReader::new(f);
        let mut header_skipped = false;
        for line in r.lines() {
            let line = line?;
            if !header_skipped {
                header_skipped = true;
                continue;
            }
            if line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 5 {
                continue;
            }
            out.push(SeedTransporterRow {
                id: cols[0].to_string(),
                name: cols[1].to_string(),
                transporter_type: cols[2].to_string(),
                exmet: cols[3].to_string(),
                exmetnames: cols[4].to_string(),
            });
        }
    }
    Ok(out)
}

/// Group SEED transporter rows by `(transporter_type, metid)` where
/// `metid` is the `[e0]`-stripped extracellular metabolite id. Reaction
/// ids for the same key are joined with `,` (matches R's
/// `paste(unique(id), collapse = ",")`).
pub fn group_seed_by_type_met(
    rows: &[SeedTransporterRow],
) -> HashMap<(String, String), Vec<String>> {
    let mut groups: HashMap<(String, String), Vec<String>> = HashMap::new();
    for r in rows {
        let metid = strip_compartment(&r.exmet);
        let key = (r.transporter_type.clone(), metid);
        let entry = groups.entry(key).or_default();
        if !entry.contains(&r.id) {
            entry.push(r.id.clone());
        }
    }
    groups
}

fn strip_compartment(exmet: &str) -> String {
    // R: `sub("\\[e.\\]$","",exmet)` — strip trailing `[eX]`.
    if let Some(idx) = exmet.rfind("[e") {
        if exmet.ends_with(']') {
            return exmet[..idx].to_string();
        }
    }
    exmet.to_string()
}

/// Load and clean `tcdb_substrates.tbl` + `tcdb_custom.tbl`. Returns a
/// map `TC_id -> ; -joined substrate names` ready for substring search.
pub fn load_tcdb_all(
    substrates: &Path,
    custom: &Path,
) -> std::io::Result<HashMap<String, String>> {
    let mut out: HashMap<String, String> = HashMap::new();
    for path in [substrates, custom] {
        let f = File::open(path)?;
        let r = BufReader::new(f);
        for line in r.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let mut cols = line.splitn(2, '\t');
            let tc = cols.next().unwrap_or("").trim().to_string();
            let raw = cols.next().unwrap_or("").to_string();
            if tc.is_empty() {
                continue;
            }
            let cleaned = clean_substrates(&raw);
            if cleaned.is_empty() {
                continue;
            }
            out.insert(tc, cleaned);
        }
    }
    Ok(out)
}

/// Port of `analyse_alignments_transport.R:24-34` — strip CHEBI ids, `|`
/// separators, charges, filler words, then trim and squish whitespace.
fn clean_substrates(raw: &str) -> String {
    use regex::Regex;
    use std::sync::OnceLock;

    static RE_CHEBI: OnceLock<Regex> = OnceLock::new();
    static RE_PIPE: OnceLock<Regex> = OnceLock::new();
    static RE_CHARGE: OnceLock<Regex> = OnceLock::new();
    static RE_STEREO: OnceLock<Regex> = OnceLock::new();
    static RE_FILLER: OnceLock<Regex> = OnceLock::new();
    static RE_DOUBLE_SEMI: OnceLock<Regex> = OnceLock::new();
    static RE_WS: OnceLock<Regex> = OnceLock::new();

    let chebi = RE_CHEBI.get_or_init(|| Regex::new(r"CHEBI:[0-9]+").unwrap());
    let pipe = RE_PIPE.get_or_init(|| Regex::new(r"\|").unwrap());
    let charge = RE_CHARGE.get_or_init(|| Regex::new(r"\([0-9](\+|-)\)").unwrap());
    let stereo = RE_STEREO.get_or_init(|| Regex::new(r"\((R|S|-)\)-").unwrap());
    let filler = RE_FILLER.get_or_init(|| {
        Regex::new(
            r"\bion\b|\bwater\b|\bmolecule\b|\bmetabolite\b|inorganic cation|metal cation|organic cation|\bcation\b|\banion\b|hydron|electron|proton",
        )
        .unwrap()
    });
    let dbl = RE_DOUBLE_SEMI.get_or_init(|| Regex::new(r";+").unwrap());
    let ws = RE_WS.get_or_init(|| Regex::new(r"\s+").unwrap());

    let mut s = chebi.replace_all(raw, "").to_string();
    s = pipe.replace_all(&s, "").to_string();
    s = charge.replace_all(&s, "").to_string();
    s = stereo.replace_all(&s, "").to_string();
    s = filler.replace_all(&s, "").to_string();
    s = s.replace("(ribonucleotide)nm", "");
    s = s.replace("(+)-", "");
    s = s.replace("(+)", "");
    s = dbl.replace_all(&s, ";").to_string();
    s = s.trim_matches(';').to_string();
    s = ws.replace_all(&s, " ").to_string();
    s.trim().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn cleans_substrate_junk() {
        let raw = "CHEBI:5651;ferroheme b|CHEBI:7435;ammonium|CHEBI:3283;carbon dioxide";
        let out = clean_substrates(raw);
        // After CHEBI removal and pipe removal we get sequences of ';'
        // separated names. The double-semi collapse turns `;;` into `;`.
        assert!(out.contains("ferroheme b"));
        assert!(out.contains("ammonium"));
        assert!(out.contains("carbon dioxide"));
        assert!(!out.contains("CHEBI"));
        assert!(!out.contains("|"));
    }

    #[test]
    fn groups_seed_by_type_met() {
        let rows = vec![
            SeedTransporterRow {
                id: "rxn1".into(),
                name: "X".into(),
                transporter_type: "1.Channels and pores".into(),
                exmet: "cpd00001[e0]".into(),
                exmetnames: "H2O".into(),
            },
            SeedTransporterRow {
                id: "rxn2".into(),
                name: "Y".into(),
                transporter_type: "1.Channels and pores".into(),
                exmet: "cpd00001[e0]".into(),
                exmetnames: "H2O".into(),
            },
        ];
        let g = group_seed_by_type_met(&rows);
        let k = ("1.Channels and pores".to_string(), "cpd00001".to_string());
        assert_eq!(g[&k], vec!["rxn1".to_string(), "rxn2".to_string()]);
    }

    #[test]
    fn strip_compartment_works() {
        assert_eq!(strip_compartment("cpd00001[e0]"), "cpd00001");
        assert_eq!(strip_compartment("cpd00001[e1]"), "cpd00001");
        assert_eq!(strip_compartment("cpd00001"), "cpd00001");
    }

    #[test]
    fn load_transporter_merges_both_files() {
        let d = tempfile::tempdir().unwrap();
        let a = d.path().join("a.tbl");
        let b = d.path().join("b.tbl");
        let mut fa = File::create(&a).unwrap();
        writeln!(fa, "id\tname\ttype\texmet\texmetnames").unwrap();
        writeln!(fa, "rxn1\tAlpha\t1.X\tcpd00001[e0]\tH2O").unwrap();
        drop(fa);
        let mut fb = File::create(&b).unwrap();
        writeln!(fb, "id\tname\ttype\texmet\texmetnames").unwrap();
        writeln!(fb, "rxn2\tBeta\t1.X\tcpd00002[e0]\tATP").unwrap();
        drop(fb);
        let rows = load_seed_transporter(&a, &b).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].id, "rxn1");
        assert_eq!(rows[1].id, "rxn2");
    }
}
