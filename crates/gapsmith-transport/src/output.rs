//! Emit `*-Transporter.tbl` with the exact column order gapseq uses
//! (`analyse_alignments_transport.R:178-183`):
//!
//! `id tc sub sub_gapseq exid rea qseqid pident evalue bitscore qcovs stitle sstart send comment`

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransporterRow {
    pub id: String,
    pub tc: String,
    pub sub: String,
    pub sub_gapseq: String,
    pub exid: String,
    pub rea: String,
    pub qseqid: String,
    pub pident: f32,
    pub evalue: f64,
    pub bitscore: f32,
    pub qcov: f32,
    pub stitle: String,
    pub sstart: i32,
    pub send: i32,
    pub comment: Option<String>,
    /// Compartment-stripped metabolite id. Not written to the TSV, used
    /// internally.
    pub metid: String,
}

pub fn write_transporter_tbl(
    rows: &[TransporterRow],
    path: &Path,
) -> std::io::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "id\ttc\tsub\tsub_gapseq\texid\trea\tqseqid\tpident\tevalue\tbitscore\tqcovs\tstitle\tsstart\tsend\tcomment"
    )?;
    for r in rows {
        writeln!(
            w,
            "{id}\t{tc}\t{sub}\t{sg}\t{exid}\t{rea}\t{qseq}\t{pident}\t{evalue}\t{bit}\t{qcov}\t{stitle}\t{sstart}\t{send}\t{comment}",
            id = r.id,
            tc = r.tc,
            sub = r.sub,
            sg = r.sub_gapseq,
            exid = r.exid,
            rea = r.rea,
            qseq = r.qseqid,
            pident = r.pident,
            evalue = fmt_evalue(r.evalue),
            bit = r.bitscore,
            qcov = r.qcov,
            stitle = r.stitle,
            sstart = r.sstart,
            send = r.send,
            comment = r.comment.as_deref().unwrap_or("NA"),
        )?;
    }
    Ok(())
}

/// Format an e-value the way BLAST's `-outfmt 6` emits it: `0` for
/// exact zero, `<mantissa>e<exp>` with 2 significant digits below
/// 1e-3, plain decimal otherwise.
fn fmt_evalue(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let abs = v.abs();
    if abs < 1e-3 || abs >= 1e5 {
        // BLAST emits 2 digits after the decimal: "2.52e-41".
        let s = format!("{v:.2e}");
        // Rust emits "2.52e-41" the same as BLAST does in this range.
        return s;
    }
    // Plain decimal — trim to at most 3 significant digits after point.
    if abs >= 10.0 {
        return format!("{v:.1}");
    }
    if abs >= 1.0 {
        return format!("{v:.2}").trim_end_matches('0').trim_end_matches('.').to_string();
    }
    format!("{v:.3}").trim_end_matches('0').trim_end_matches('.').to_string()
}
