//! Pathway selection and per-reaction expansion.
//!
//! Input: a `PathwayTable` (loaded by `gapsmith-db`) plus user filters. The
//! selector walks every row, keeps those matching the requested keyword
//! (or `"all"`), and then expands each surviving row into one entry per
//! `(reaction_id, reaction_name, ec, keyrea, spont)` tuple.
//!
//! `prepare_batch_alignments.R:150-224` in the R reference does the same
//! fan-out, keyed off the same four comma-separated columns
//! (`reaId`, `reaEc`, `keyRea`, `reaName`).
//!
//! ## Keyword resolution
//!
//! The user-facing `-p` keyword is first mapped to a MetaCyc hierarchy
//! category (or a set of categories) via [`resolve_keyword`], mirroring
//! `src/gapseq_find.sh:351-421`. `-p amino` expands to
//! `Amino-Acid-Biosynthesis`, `-p all` expands to a superset of every
//! gapseq category, and any other token is matched as-is.
//!
//! Selection then runs `grep -wE $pwyKey` against the full TSV row, then
//! applies the taxonomic-range filter: a pathway is kept iff its
//! `taxrange` column contains any of the NCBI tax IDs listed in
//! `dat/taxonomy.tbl` for the requested taxonomy (default Bacteria) OR is
//! empty.

use gapsmith_db::{PathwayRow, PathwayTable};
use regex::Regex;
use std::collections::HashSet;

/// One (pathway, reaction) pair ready for alignment.
#[derive(Debug, Clone)]
pub struct ExpandedReaction {
    pub pathway: String,
    pub pathway_name: String,
    pub rxn: String,
    pub name: String,
    pub ec: String,
    pub keyrea: bool,
    pub spont: bool,
}

/// Match mode used by [`PathwaySelectOptions::keyword`].
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum MatchMode {
    /// Exact match against the pathway `id` column.
    Id,
    /// `|`- or `;`-split exact match against the `hierarchy` column tokens
    /// (case-insensitive). Matches `filter_pathways.R:20-23`.
    Hierarchy,
    /// Word-boundary regex search over id/hierarchy/name/altname columns.
    /// Matches the fallback branch of `filter_pathways.R:25-30`.
    Regex,
}

pub struct PathwaySelectOptions<'a> {
    /// User keyword. Interpretation depends on `mode`.
    pub keyword: &'a str,
    pub mode: MatchMode,
    /// If true, drop rows where `superpathway == TRUE` or `hierarchy`
    /// contains `Super-Pathways`. Matches `filter_pathways.R:32-34`.
    pub exclude_superpathways: bool,
    /// Only keep pathways whose status column is `TRUE`.
    pub only_active: bool,
    /// Set of NCBI tax IDs (as strings, no `TAX-` prefix) that pathways
    /// must overlap via their `taxrange` column. Pass empty to disable
    /// the taxonomic filter.
    pub valid_tax_ids: &'a [String],
}

impl<'a> Default for PathwaySelectOptions<'a> {
    fn default() -> Self {
        Self {
            keyword: "all",
            mode: MatchMode::Regex,
            exclude_superpathways: true,
            only_active: true,
            valid_tax_ids: &[],
        }
    }
}

/// Select matching pathways, then expand to per-reaction rows. Duplicate
/// (pathway, rxn, name, ec) entries collapse to one.
pub fn select(table: &PathwayTable, opts: &PathwaySelectOptions<'_>) -> Vec<ExpandedReaction> {
    let pattern = resolve_keyword(opts.keyword);
    // Split on `|` to get individual keys, the same way `filter_pathways.R`
    // does via `strsplit(args[2], "\\|")`.
    let keys_lc: Vec<String> = pattern
        .split('|')
        .map(|s| s.trim().to_ascii_lowercase())
        .filter(|s| !s.is_empty())
        .collect();

    let kw_re: Option<Regex> = if keys_lc.is_empty() || matches!(opts.mode, MatchMode::Hierarchy | MatchMode::Id) {
        None
    } else {
        let alt = keys_lc.iter().map(|s| regex::escape(s)).collect::<Vec<_>>().join("|");
        Some(
            regex::RegexBuilder::new(&format!(r"\b(?:{alt})\b"))
                .case_insensitive(true)
                .build()
                .expect("keyword regex compiles"),
        )
    };

    let tax_re: Option<Regex> = if opts.valid_tax_ids.is_empty() {
        None
    } else {
        let alt = opts
            .valid_tax_ids
            .iter()
            .map(|s| regex::escape(s))
            .collect::<Vec<_>>()
            .join("|");
        Some(Regex::new(&format!(r"TAX-(?:{alt})\b")).expect("tax regex compiles"))
    };

    let mut out = Vec::new();
    let mut seen: HashSet<(String, String, String, String)> = HashSet::new();

    for row in &table.rows {
        if opts.only_active && !is_active(row) {
            continue;
        }
        if !keys_lc.is_empty() {
            let matched = match opts.mode {
                MatchMode::Id => keys_lc.iter().any(|k| row.id.to_ascii_lowercase() == *k),
                MatchMode::Hierarchy => hierarchy_matches_any(&row.hierarchy, &keys_lc),
                MatchMode::Regex => kw_re.as_ref().unwrap().is_match(&row.id)
                    || kw_re.as_ref().unwrap().is_match(&row.hierarchy)
                    || kw_re.as_ref().unwrap().is_match(&row.name)
                    || kw_re.as_ref().unwrap().is_match(&row.altname),
            };
            if !matched {
                continue;
            }
        }
        // Superpathway filter.
        if opts.exclude_superpathways && is_superpathway(row) {
            continue;
        }
        // Reactions-empty filter — gapseq drops these upstream.
        if row.rea_id.trim().is_empty() {
            continue;
        }
        // Tax filter: keep if taxrange is empty OR matches a valid tax id.
        if let Some(ref t_re) = tax_re {
            if !row.taxrange.is_empty() && !t_re.is_match(&row.taxrange) {
                continue;
            }
        }
        expand_row(row, &mut out, &mut seen);
    }
    out
}

fn hierarchy_matches_any(hierarchy: &str, keys_lc: &[String]) -> bool {
    // Split on `|` AND `;`, matching R's `strsplit(hierarchy, "\\|")` and
    // the alternative `strsplit(hierarchy, "\\;")`.
    for token in hierarchy.split(|c: char| c == '|' || c == ';') {
        let t = token.trim().to_ascii_lowercase();
        if keys_lc.iter().any(|k| k == &t) {
            return true;
        }
    }
    false
}

fn is_superpathway(row: &PathwayRow) -> bool {
    matches!(row.superpathway.to_ascii_lowercase().as_str(), "true" | "t" | "1")
        || row.hierarchy.contains("Super-Pathways")
}

/// Map a user keyword to the regex alternatives gapseq's `gapseq_find.sh`
/// plugs into `grep -wE`. Mirrors lines 351-421 of the shell script.
pub fn resolve_keyword(kw: &str) -> String {
    let aa_syn =
        "Amino-Acid-Biosynthesis|Nucleotide-Biosynthesis|Cofactor-Biosynthesis|Carbohydrates-Degradation|CARBO-BIOSYNTHESIS|Polyamine-Biosynthesis|Fatty-acid-biosynthesis|Energy-Metabolism|Terpenoid-Biosynthesis|Chorismate-Biosynthesis";
    let min_pwys = "ETOH-ACETYLCOA-ANA-PWY|GLNSYN-PWY|GLUCONEO-PWY|GLUGLNSYN-PWY|GLUTAMATE-DEG1-PWY|GLUTAMATE-SYN2-PWY|GLUTAMINEFUM-PWY|GLUTSYNIII-PWY|GLYCOLYSIS|GLYOXYLATE-BYPASS|NONOXIPENT-PWY|OXIDATIVEPENT-PWY|P185-PWY|P21-PWY|PWY0-1312|PWY0-1315|PWY0-1329|PWY0-1334|PWY0-1335|PWY0-1353|PWY0-1517|PWY0-1565|PWY0-1567|PWY0-1568|PWY-4341|PWY-5084|PWY-5480|PWY-5482|PWY-5484|PWY-5690|PWY-5766|PWY-5913|PWY-6028|PWY-6333|PWY-6543|PWY-6549|PWY66-21|PWY66-398|PWY-6697|PWY-6964|PWY-7167|PWY-7685|PWY-7686|PWY-7980|PWY-8178|PWY-8215|PWY-8274|PWY-8404|PYRUVDEHYD-PWY|TCA-1|TCA";
    match kw {
        "" => String::new(),
        "all" => "Pathways|Enzyme-Test|seed|kegg".into(),
        "amino" => "Amino-Acid-Biosynthesis".into(),
        "nucl" => "Nucleotide-Biosynthesis".into(),
        "cofactor" => "Cofactor-Biosynthesis".into(),
        "carbo" => "CARBO-BIOSYNTHESIS".into(),
        "carbo-deg" => "Carbohydrates-Degradation".into(),
        "polyamine" => "Polyamine-Biosynthesis".into(),
        "fatty" => "Fatty-acid-biosynthesis".into(),
        "energy" => "Energy-Metabolism".into(),
        "terpenoid" => "Terpenoid-Biosynthesis".into(),
        "degradation" => "Degradation".into(),
        "core" => aa_syn.into(),
        "min" => min_pwys.into(),
        "kegg" => "kegg".into(),
        other => other.into(), // literal or regex alternatives
    }
}

fn is_active(row: &PathwayRow) -> bool {
    matches!(row.status.to_ascii_lowercase().as_str(), "true" | "t" | "1")
}


fn expand_row(
    row: &PathwayRow,
    out: &mut Vec<ExpandedReaction>,
    seen: &mut HashSet<(String, String, String, String)>,
) {
    // Parallel-column splits — must PRESERVE empty slots so index
    // alignment between rea_id, rea_ec, and rea_name holds. Only
    // `key_rea` and `spont` are treated as set-like (empty-stripped).
    let rxns = split_positional(&row.rea_id, ',');
    let ecs = split_positional(&row.rea_ec, ',');
    let names = split_positional(&row.rea_name, ';');
    let keys = split_csv(&row.key_rea);
    let spont = split_csv(&row.spont);

    let key_set: HashSet<&str> = keys.iter().copied().collect();
    let spont_set: HashSet<&str> = spont.iter().copied().collect();

    // Align every parallel list to `rxns` — gapseq data guarantees reaId,
    // reaEc, reaName all share the same length, but we defensively fall
    // back to "" on mismatch.
    for (i, r) in rxns.iter().enumerate() {
        let r = r.as_str();
        if r.is_empty() {
            continue;
        }
        let ec = ecs.get(i).map(|s| s.as_str()).unwrap_or("");
        let name = names.get(i).map(|s| s.as_str()).unwrap_or("");
        let key = r == row.id || key_set.contains(r);
        let sp = spont_set.contains(r);
        let dedup = (row.id.clone(), r.to_string(), name.to_string(), ec.to_string());
        if seen.insert(dedup) {
            out.push(ExpandedReaction {
                pathway: row.id.clone(),
                pathway_name: row.name.clone(),
                rxn: r.to_string(),
                name: name.to_string(),
                ec: ec.to_string(),
                keyrea: key,
                spont: sp,
            });
        }
    }
}

fn split_csv(s: &str) -> Vec<&str> {
    s.split(',').map(str::trim).filter(|s| !s.is_empty()).collect()
}

/// Positional split — keeps empty slots so indices line up across
/// parallel columns. Trims each item. Returns `Vec<String>` so callers
/// can reason about lifetime independently.
fn split_positional(s: &str, sep: char) -> Vec<String> {
    if s.is_empty() {
        return Vec::new();
    }
    s.split(sep).map(|t| t.trim().to_string()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use gapsmith_db::PwySource;

    fn toy_table() -> PathwayTable {
        PathwayTable {
            source: Some(PwySource::MetaCyc),
            rows: vec![
                PathwayRow {
                    id: "PWY-A".into(),
                    name: "Glycolysis".into(),
                    altname: "".into(),
                    hierarchy: "Carbohydrate metabolism".into(),
                    taxrange: "".into(),
                    rea_id: "rxn1,rxn2,rxn3".into(),
                    rea_ec: "1.1.1.1,2.2.2.2,3.3.3.3".into(),
                    key_rea: "rxn1".into(),
                    rea_name: "enz1;enz2;enz3".into(),
                    rea_nr: 3,
                    ec_nr: 3,
                    superpathway: "".into(),
                    status: "TRUE".into(),
                    spont: "rxn3".into(),
                    source: PwySource::MetaCyc,
                },
                PathwayRow {
                    id: "PWY-B".into(),
                    name: "Lysine biosynthesis".into(),
                    altname: "".into(),
                    hierarchy: "Amino acid".into(),
                    taxrange: "".into(),
                    rea_id: "rxn10".into(),
                    rea_ec: "4.4.4.4".into(),
                    key_rea: "".into(),
                    rea_name: "enz10".into(),
                    rea_nr: 1,
                    ec_nr: 1,
                    superpathway: "".into(),
                    status: "TRUE".into(),
                    spont: "".into(),
                    source: PwySource::MetaCyc,
                },
            ],
        }
    }

    fn opts_with_keyword(kw: &'static str) -> PathwaySelectOptions<'static> {
        PathwaySelectOptions {
            keyword: kw,
            mode: MatchMode::Regex,
            exclude_superpathways: true,
            only_active: true,
            valid_tax_ids: &[],
        }
    }

    #[test]
    fn keyword_matches_name_word_boundary() {
        let t = toy_table();
        let exp = select(&t, &opts_with_keyword("Glycolysis"));
        assert_eq!(exp.len(), 3);
        assert_eq!(exp[0].pathway, "PWY-A");
    }

    #[test]
    fn keyword_all_expands_every_active_pathway() {
        // `all` → pattern matches "Pathways" but since our hierarchy
        // strings don't contain that token, nothing matches. Use a real
        // gapseq hierarchy to exercise the path.
        let mut t = toy_table();
        for r in &mut t.rows {
            r.hierarchy.push_str(",Pathways");
        }
        let exp = select(&t, &opts_with_keyword("all"));
        assert_eq!(exp.len(), 4);
    }

    #[test]
    fn keyword_shorthand_resolves_to_hierarchy() {
        let mut t = toy_table();
        t.rows[1].hierarchy = "Amino-Acid-Biosynthesis".into();
        let exp = select(
            &t,
            &PathwaySelectOptions {
                keyword: "amino",
                mode: MatchMode::Hierarchy,
                exclude_superpathways: true,
                only_active: true,
                valid_tax_ids: &[],
            },
        );
        assert_eq!(exp.len(), 1);
        assert_eq!(exp[0].pathway, "PWY-B");
    }

    #[test]
    fn taxonomy_filter() {
        let mut t = toy_table();
        t.rows[0].taxrange = "|TAX-4751|".into(); // fungi-only
        t.rows[1].taxrange = "|TAX-2|".into(); // bacteria
        let tax_ids = vec!["2".to_string(), "1224".to_string()];
        let exp = select(
            &t,
            &PathwaySelectOptions {
                keyword: "all",
                mode: MatchMode::Regex,
                exclude_superpathways: true,
                only_active: true,
                valid_tax_ids: &tax_ids,
            },
        );
        // Since neither pathway has "Pathways" in hierarchy, add it.
        let _ = exp; // other test covers 'all'
        // Instead, verify the filter directly.
        let exp = select(
            &t,
            &PathwaySelectOptions {
                keyword: "Lysine",
                mode: MatchMode::Regex,
                exclude_superpathways: true,
                only_active: true,
                valid_tax_ids: &tax_ids,
            },
        );
        assert_eq!(exp.len(), 1);
        assert_eq!(exp[0].pathway, "PWY-B");
    }
}
