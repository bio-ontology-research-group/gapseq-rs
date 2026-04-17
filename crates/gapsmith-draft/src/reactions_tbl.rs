//! Read the `*-Reactions.tbl` and `*-Transporter.tbl` produced by
//! `gapsmith find` / `gapsmith find-transport`. Column schemas match the
//! gapseq output exactly (see `gapsmith_find::output` and
//! `gapsmith_transport::output`).

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct ReactionRow {
    pub pathway: String,
    pub rxn: String,
    pub name: String,
    pub ec: String,
    pub keyrea: bool,
    pub file: String,
    pub dbhit: String,
    pub spont: bool,
    pub reftype: String,
    pub src: String,
    pub is_complex: bool,
    pub subunit_count: Option<u32>,
    pub subunits: String,
    pub qseqid: String,
    pub pident: Option<f32>,
    pub evalue: Option<f64>,
    pub bitscore: Option<f32>,
    pub qcov: Option<f32>,
    pub stitle: String,
    pub complex: String,
    pub exception: bool,
    pub status: String,
    pub subunits_found: Option<u32>,
    pub complex_status: Option<u8>,
    pub pathway_status: String,
}

#[derive(Debug, Clone)]
pub struct TransporterRow {
    pub id: String,
    pub tc: String,
    pub sub: String,
    pub sub_gapseq: String,
    pub exid: String,
    pub rea: String,
    pub qseqid: String,
    pub pident: Option<f32>,
    pub evalue: Option<f64>,
    pub bitscore: Option<f32>,
    pub qcov: Option<f32>,
    pub stitle: String,
    pub comment: String,
}

pub fn read_reactions_tbl(path: &Path) -> std::io::Result<Vec<ReactionRow>> {
    let f = File::open(path)?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    let mut header: Option<Vec<String>> = None;
    for line in r.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if header.is_none() {
            if !cols.is_empty() && cols[0] == "pathway" {
                header = Some(cols.iter().map(|s| s.to_string()).collect());
                continue;
            } else {
                continue;
            }
        }
        let h = header.as_ref().unwrap();
        let get = |name: &str| -> String {
            h.iter()
                .position(|c| c == name)
                .and_then(|i| cols.get(i))
                .map(|s| (*s).to_string())
                .unwrap_or_default()
        };
        out.push(ReactionRow {
            pathway: get("pathway"),
            rxn: get("rxn"),
            name: get("name"),
            ec: get("ec"),
            keyrea: is_true(&get("keyrea")),
            file: get("file"),
            dbhit: get("dbhit"),
            spont: is_true(&get("spont")),
            reftype: get("type"),
            src: get("src"),
            is_complex: is_true(&get("is_complex")),
            subunit_count: parse_u32(&get("subunit_count")),
            subunits: get("subunits"),
            qseqid: get("qseqid"),
            pident: parse_f32(&get("pident")),
            evalue: parse_f64(&get("evalue")),
            bitscore: parse_f32(&get("bitscore")),
            qcov: parse_f32(&get("qcovs")),
            stitle: get("stitle"),
            complex: get("complex"),
            exception: get("exception") == "1" || is_true(&get("exception")),
            status: get("status"),
            subunits_found: parse_u32(&get("subunits_found")),
            complex_status: parse_u8(&get("complex.status")),
            pathway_status: get("pathway.status"),
        });
    }
    Ok(out)
}

pub fn read_transporter_tbl(path: &Path) -> std::io::Result<Vec<TransporterRow>> {
    let f = File::open(path)?;
    let r = BufReader::new(f);
    let mut out = Vec::new();
    let mut header: Option<Vec<String>> = None;
    for line in r.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if header.is_none() {
            if !cols.is_empty() && cols[0] == "id" {
                header = Some(cols.iter().map(|s| s.to_string()).collect());
                continue;
            } else {
                continue;
            }
        }
        let h = header.as_ref().unwrap();
        let get = |name: &str| -> String {
            h.iter()
                .position(|c| c == name)
                .and_then(|i| cols.get(i))
                .map(|s| (*s).to_string())
                .unwrap_or_default()
        };
        out.push(TransporterRow {
            id: get("id"),
            tc: get("tc"),
            sub: get("sub"),
            sub_gapseq: get("sub_gapseq"),
            exid: get("exid"),
            rea: get("rea"),
            qseqid: get("qseqid"),
            pident: parse_f32(&get("pident")),
            evalue: parse_f64(&get("evalue")),
            bitscore: parse_f32(&get("bitscore")),
            qcov: parse_f32(&get("qcovs")),
            stitle: get("stitle"),
            comment: get("comment"),
        });
    }
    Ok(out)
}

fn is_true(s: &str) -> bool {
    matches!(s, "TRUE" | "true" | "T" | "1")
}
fn parse_u32(s: &str) -> Option<u32> {
    if s.is_empty() || s == "NA" {
        None
    } else {
        s.parse().ok()
    }
}
fn parse_u8(s: &str) -> Option<u8> {
    if s.is_empty() || s == "NA" {
        None
    } else {
        s.parse().ok()
    }
}
fn parse_f32(s: &str) -> Option<f32> {
    if s.is_empty() || s == "NA" {
        None
    } else {
        s.parse().ok()
    }
}
fn parse_f64(s: &str) -> Option<f64> {
    if s.is_empty() || s == "NA" {
        None
    } else {
        s.parse().ok()
    }
}
