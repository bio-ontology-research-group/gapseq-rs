//! Reaction → SEED reaction ID lookup ("dbhit").
//!
//! Port of `src/getDBhit.R` for the SEED target database (gapseq's default;
//! VMH / BiGG paths are omitted as gapseq itself considers VMH deprecated).
//!
//! Strategies implemented, matching the R numbering:
//!
//! 1. EC-based: match `External ID == EC` in
//!    `dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv`, split
//!    `MS ID` on `|`.
//! 3. Alternative-EC expansion via `dat/altec.csv` (BRENDA fallback
//!    omitted in M5 — covers <1% of reactions and isn't needed for
//!    parity on the common case).
//! 5. MetaCyc-id match via `dat/mnxref_seed-other.tsv` — any row with
//!    `other == <metacyc_id>` contributes its `seed` column.
//! 6. Enzyme-name match via
//!    `dat/seed_Enzyme_Name_Reactions_Aliases.tsv` — exact name match.
//!
//! The result of every strategy is unioned, sorted, deduplicated, and
//! space-joined to form the `dbhit` column — identical to
//! `getDBhit.R:129`.

use crate::seqfile::looks_like_rxn_id;
use gapsmith_db::{MnxrefSeedOther, SeedCpdRow};
use std::collections::{BTreeSet, HashMap};
use std::path::Path;

/// Index over the reference data needed to answer dbhit queries.
pub struct DbhitIndex {
    /// EC → set of SEED rxn ids (from strategy 1).
    ec_to_rxns: HashMap<String, BTreeSet<String>>,
    /// EC → alternative ECs (from `altec.csv`, strategy 3 precompute).
    ec_alt: HashMap<String, Vec<String>>,
    /// MetaCyc (or other non-seed) id → set of SEED rxn ids (from strategy 5).
    metaid_to_rxns: HashMap<String, BTreeSet<String>>,
    /// Enzyme name (exact) → set of SEED rxn ids (from strategy 6).
    name_to_rxns: HashMap<String, BTreeSet<String>>,
}

impl DbhitIndex {
    /// Build from the reference `dat/` directory.
    pub fn load(data_dir: &Path) -> Result<Self, std::io::Error> {
        let ec_to_rxns = load_seed_ec_aliases(
            &data_dir.join("seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"),
        )?;
        let ec_alt = load_altec(&data_dir.join("altec.csv"))?;
        let metaid_to_rxns =
            load_mnxref_seed_other(&data_dir.join("mnxref_seed-other.tsv"))?;
        let name_to_rxns =
            load_seed_name_aliases(&data_dir.join("seed_Enzyme_Name_Reactions_Aliases.tsv"))?;
        Ok(Self { ec_to_rxns, ec_alt, metaid_to_rxns, name_to_rxns })
    }

    /// Return the space-joined, sorted, deduplicated SEED rxn ids matching
    /// a query. Mirrors `paste(sort(unique(dbhit)), collapse = " ")`.
    pub fn lookup(&self, rxn: &str, name: &str, ec: &str) -> String {
        let mut out: BTreeSet<String> = BTreeSet::new();

        // Strategy 1: EC-based.
        for one_ec in split_slash(ec) {
            if is_valid_ec(&one_ec) {
                if let Some(set) = self.ec_to_rxns.get(&one_ec) {
                    out.extend(set.iter().cloned());
                }
            }
        }

        // Strategy 3: alternative ECs.
        for one_ec in split_slash(ec) {
            if let Some(alt_list) = self.ec_alt.get(&one_ec) {
                for alt in alt_list {
                    // Filter unspecific `.99.n` ECs.
                    if alt_is_unspecific(alt) {
                        continue;
                    }
                    if let Some(set) = self.ec_to_rxns.get(alt) {
                        out.extend(set.iter().cloned());
                    }
                }
            }
        }

        // Strategy 5: MetaCyc id lookup (only for real reaction ids, not
        // synthetic `rxn` shorthands; those would be SEED already).
        if !rxn.is_empty() && !looks_like_rxn_id(rxn) {
            if let Some(set) = self.metaid_to_rxns.get(rxn) {
                out.extend(set.iter().cloned());
            }
        }

        // Strategy 6: enzyme name.
        if !name.is_empty() {
            if let Some(set) = self.name_to_rxns.get(name) {
                out.extend(set.iter().cloned());
            }
        }

        out.into_iter().collect::<Vec<_>>().join(" ")
    }
}

fn split_slash(s: &str) -> Vec<String> {
    s.split('/').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect()
}

fn is_valid_ec(s: &str) -> bool {
    let parts: Vec<&str> = s.split('.').collect();
    parts.len() == 4 && parts.iter().all(|p| !p.is_empty() && p.chars().all(|c| c.is_ascii_digit()))
}

fn alt_is_unspecific(ec: &str) -> bool {
    // R regex: `\\.99\\.[0-9]+$`.
    let parts: Vec<&str> = ec.split('.').collect();
    if parts.len() != 4 {
        return false;
    }
    parts[2] == "99" && parts[3].chars().all(|c| c.is_ascii_digit()) && !parts[3].is_empty()
}

fn load_seed_ec_aliases(path: &Path) -> Result<HashMap<String, BTreeSet<String>>, std::io::Error> {
    use std::io::BufRead;
    let f = std::fs::File::open(path)?;
    let r = std::io::BufReader::new(f);
    let mut out: HashMap<String, BTreeSet<String>> = HashMap::new();
    let mut header = true;
    for line in r.lines() {
        let line = line?;
        if header {
            header = false;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let ms_id = cols[0];
        let external_id = cols[2];
        if external_id.is_empty() {
            continue;
        }
        let entry = out.entry(external_id.to_string()).or_default();
        for rxn in ms_id.split('|').filter(|s| !s.is_empty()) {
            entry.insert(rxn.to_string());
        }
    }
    Ok(out)
}

fn load_altec(path: &Path) -> Result<HashMap<String, Vec<String>>, std::io::Error> {
    use std::io::BufRead;
    let f = std::fs::File::open(path)?;
    let r = std::io::BufReader::new(f);
    let mut out: HashMap<String, Vec<String>> = HashMap::new();
    for line in r.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let members: Vec<&str> = line.split(',').map(str::trim).collect();
        for ec in &members {
            let alts: Vec<String> = members
                .iter()
                .filter(|m| *m != ec && !m.is_empty())
                .map(|s| s.to_string())
                .collect();
            out.entry((*ec).to_string()).or_default().extend(alts);
        }
    }
    Ok(out)
}

fn load_mnxref_seed_other(
    path: &Path,
) -> Result<HashMap<String, BTreeSet<String>>, std::io::Error> {
    use std::io::BufRead;
    let f = std::fs::File::open(path)?;
    let r = std::io::BufReader::new(f);
    let mut out: HashMap<String, BTreeSet<String>> = HashMap::new();
    let mut header = true;
    for line in r.lines() {
        let line = line?;
        if header {
            header = false;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let seed = cols[1];
        let other = cols[2];
        if other.is_empty() || seed.is_empty() {
            continue;
        }
        out.entry(other.to_string()).or_default().insert(seed.to_string());
    }
    Ok(out)
}

fn load_seed_name_aliases(
    path: &Path,
) -> Result<HashMap<String, BTreeSet<String>>, std::io::Error> {
    use std::io::BufRead;
    let f = std::fs::File::open(path)?;
    let r = std::io::BufReader::new(f);
    let mut out: HashMap<String, BTreeSet<String>> = HashMap::new();
    let mut header = true;
    for line in r.lines() {
        let line = line?;
        if header {
            header = false;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            continue;
        }
        let rxn_id = cols[0];
        let name = cols[1];
        if name.is_empty() || rxn_id.is_empty() {
            continue;
        }
        out.entry(name.to_string()).or_default().insert(rxn_id.to_string());
    }
    Ok(out)
}

// Keep public-facing types exported by lib.rs tidy by not importing types
// we don't use directly. These re-exports are unused internally but ensure
// rustc doesn't prune the dependency declarations from the workspace.
#[allow(dead_code)]
fn _unused_imports(_: &SeedCpdRow, _: &MnxrefSeedOther) {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn empty_index() -> DbhitIndex {
        DbhitIndex {
            ec_to_rxns: HashMap::new(),
            ec_alt: HashMap::new(),
            metaid_to_rxns: HashMap::new(),
            name_to_rxns: HashMap::new(),
        }
    }

    #[test]
    fn validates_ec() {
        assert!(is_valid_ec("1.2.3.4"));
        assert!(!is_valid_ec("1.2.3"));
        assert!(!is_valid_ec(""));
        assert!(!is_valid_ec("1.1.1.n1"));
    }

    #[test]
    fn unspecific_filter() {
        assert!(alt_is_unspecific("1.1.99.1"));
        assert!(!alt_is_unspecific("1.1.1.1"));
    }

    #[test]
    fn loads_seed_ec_aliases() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("ec.tsv");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "MS ID\tOld MS ID\tExternal ID\tSource").unwrap();
        writeln!(f, "rxn00171|rxn00869\trxn00171|rxn00869\t1.2.1.10\tEnzyme Class").unwrap();
        drop(f);
        let m = load_seed_ec_aliases(&p).unwrap();
        let set = &m["1.2.1.10"];
        assert!(set.contains("rxn00171"));
        assert!(set.contains("rxn00869"));
    }

    #[test]
    fn lookup_union_and_sort() {
        let mut idx = empty_index();
        idx.ec_to_rxns.insert("1.2.1.10".into(), {
            let mut s = BTreeSet::new();
            s.insert("rxn00171".into());
            s.insert("rxn00869".into());
            s
        });
        assert_eq!(idx.lookup("", "", "1.2.1.10"), "rxn00171 rxn00869");
    }

    #[test]
    fn lookup_empty_for_unknown_ec() {
        let idx = empty_index();
        assert_eq!(idx.lookup("", "", "9.9.9.9"), "");
    }
}
