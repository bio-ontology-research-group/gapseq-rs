//! Emit `*-Reactions.tbl` and `*-Pathways.tbl` in gapseq's column order.
//!
//! Columns mirror `src/analyse_alignments.R` output exactly so existing
//! downstream consumers (generate_GSdraft.R, predict_medium.R) continue to
//! parse them.

use crate::types::{PathwayResult, ReactionHit};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Column order (28 fields) mirrors gapseq's actual output exactly:
///
/// `pathway rxn name ec keyrea file dbhit spont type src is_complex
///  subunit_count subunits qseqid pident evalue bitscore qcovs stitle
///  sstart send complex exception status subunits_found
///  subunit_undefined_found complex.status pathway.status`
pub fn write_reactions_tbl(
    rows: &[ReactionHit],
    path: &Path,
) -> Result<(), std::io::Error> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "pathway\trxn\tname\tec\tkeyrea\tfile\tdbhit\tspont\ttype\tsrc\tis_complex\tsubunit_count\tsubunits\tqseqid\tpident\tevalue\tbitscore\tqcovs\tstitle\tsstart\tsend\tcomplex\texception\tstatus\tsubunits_found\tsubunit_undefined_found\tcomplex.status\tpathway.status"
    )?;
    for r in rows {
        writeln!(
            w,
            "{pwy}\t{rxn}\t{name}\t{ec}\t{key}\t{file}\t{db}\t{spont}\t{rtype}\t{src}\t{cplx}\t{sc}\t{subs}\t{qseq}\t{pident}\t{evalue}\t{bit}\t{qcov}\t{stitle}\t{sstart}\t{send}\t{cpx_name}\t{exc}\t{status}\t{sf}\t{suf}\t{cpx_st}\t{pwy_st}",
            pwy = r.pathway,
            rxn = r.rxn,
            name = r.name,
            ec = r.ec,
            key = if r.keyrea { "TRUE" } else { "FALSE" },
            file = r.file.as_deref().unwrap_or(""),
            db = r.dbhit,
            spont = if r.spont { "TRUE" } else { "FALSE" },
            rtype = r.reftype,
            src = r.src,
            cplx = if r.is_complex { "TRUE" } else { "FALSE" },
            // gapseq emits `subunit_count` blank for non-complex reactions;
            // non-blank ≥ 2 for real complexes.
            sc = if r.subunit_count == 0 { String::new() } else { r.subunit_count.to_string() },
            // gapseq emits `subunits` as `NA` when complex detection ran
            // over the reference fastas but subunit_count came out 0
            // (i.e. non-complex). When no reference fasta existed, the
            // column stays blank. Runner clears the string in the latter
            // case; we only need to fill the `NA` here.
            subs = if r.is_complex {
                r.subunits.clone()
            } else if r.subunits.is_empty() && r.file.is_some() {
                "NA".to_string()
            } else {
                r.subunits.clone()
            },
            qseq = r.qseqid.as_deref().unwrap_or(""),
            pident = opt_f32(r.pident),
            evalue = opt_f64(r.evalue),
            bit = opt_f32(r.bitscore),
            qcov = opt_f32(r.qcov),
            stitle = r.stitle.as_deref().unwrap_or(""),
            sstart = r.sstart.map(|v| v.to_string()).unwrap_or_default(),
            send = r.send.map(|v| v.to_string()).unwrap_or_default(),
            cpx_name = r.complex.as_deref().unwrap_or(""),
            exc = if r.exception { "1" } else { "0" },
            status = r.status.as_str(),
            // subunits_found: gapseq writes the count for complexes and
            // leaves it blank for non-complex reactions.
            sf = r.subunits_found.map(|v| v.to_string()).unwrap_or_default(),
            suf = match r.subunit_undefined_found {
                Some(true) => "TRUE",
                Some(false) => "FALSE",
                None => "",
            },
            cpx_st = r.complex_status.map(|v| v.to_string()).unwrap_or_default(),
            pwy_st = r.pathway_status.map(|s| s.as_str()).unwrap_or(""),
        )?;
    }
    Ok(())
}

pub fn write_pathways_tbl(
    rows: &[PathwayResult],
    path: &Path,
) -> Result<(), std::io::Error> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "ID\tName\tPrediction\tCompleteness\tStatus\tNrReaction\tNrSpontaneous\tNrVague\tNrKeyReaction\tNrReactionFound\tNrKeyReactionFound\tReactionsFound\tSpontaneousReactions\tKeyReactions"
    )?;
    for r in rows {
        writeln!(
            w,
            "{id}\t{name}\t{pred}\t{comp}\t{status}\t{nr}\t{ns}\t{nv}\t{nk}\t{nf}\t{nkf}\t{found}\t{spont}\t{key}",
            id = r.id,
            name = r.name,
            pred = if r.prediction { "TRUE" } else { "FALSE" },
            comp = format_completeness(r.completeness),
            status = r.status.map(|s| s.as_str()).unwrap_or(""),
            nr = r.n_reaction,
            ns = r.n_spontaneous,
            nv = r.n_vague,
            nk = r.n_key_reaction,
            nf = r.n_reaction_found,
            nkf = r.n_key_reaction_found,
            found = r.reactions_found.join(" "),
            spont = r.spontaneous_reactions.join(" "),
            key = r.key_reactions.join(" "),
        )?;
    }
    Ok(())
}

fn opt_f32(v: Option<f32>) -> String {
    v.map(|x| format!("{x}")).unwrap_or_default()
}

/// Match R's default f64 → string conversion for completeness values like
/// `2/3 * 100 = 66.6666666666667` (R trims to 15 significant digits).
fn format_completeness(v: f64) -> String {
    if v.fract() == 0.0 {
        return format!("{}", v as i64);
    }
    // R's `format(x)` defaults to 15 significant digits, trimming trailing
    // zeros. We approximate with 13 decimal places and a zero trim — good
    // enough for values in the 0–100 range that completeness lives in.
    let precise = format!("{v:.13}");
    precise.trim_end_matches('0').trim_end_matches('.').to_string()
}

fn opt_f64(v: Option<f64>) -> String {
    v.map(|x| {
        if x == 0.0 {
            "0".to_string()
        } else if x.abs() < 1e-3 || x.abs() >= 1e5 {
            format!("{x:.3e}")
        } else {
            format!("{x}")
        }
    })
    .unwrap_or_default()
}
