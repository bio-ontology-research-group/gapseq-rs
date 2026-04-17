//! TSV parser shared by every aligner backend.
//!
//! Input: lines of tab-separated values with the 8 columns documented on
//! the crate root (`qseqid pident evalue bitscore qcov stitle sstart send`).
//! Blank lines and `#`-prefixed comments are skipped.
//!
//! `coverage_is_fraction = true` indicates the source tool (e.g. mmseqs2)
//! reports `qcov` as a 0–1 fraction; we rescale to 0–100 for consistency.

use crate::error::AlignError;
use crate::hit::Hit;
use std::io::BufRead;

pub fn parse_tsv<R: BufRead>(
    rdr: R,
    coverage_is_fraction: bool,
) -> Result<Vec<Hit>, AlignError> {
    let mut out = Vec::new();
    for (i, line) in rdr.lines().enumerate() {
        let line = line.map_err(|e| AlignError::TsvParse {
            line: (i + 1) as u64,
            msg: format!("read error: {e}"),
        })?;
        let l = line.trim_end_matches('\r');
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let hit = parse_hit_line(l, coverage_is_fraction, (i + 1) as u64)?;
        out.push(hit);
    }
    Ok(out)
}

fn parse_hit_line(
    line: &str,
    coverage_is_fraction: bool,
    line_no: u64,
) -> Result<Hit, AlignError> {
    // 8 fields — but stitle can contain tabs(?) it shouldn't; BLAST's -outfmt 6
    // preserves spaces but replaces tabs in stitle. We split on tab directly.
    let cols: Vec<&str> = line.splitn(8, '\t').collect();
    if cols.len() < 8 {
        return Err(AlignError::TsvParse {
            line: line_no,
            msg: format!("expected 8 tab-separated fields, got {}", cols.len()),
        });
    }
    let parse_f32 = |s: &str, field: &str| -> Result<f32, AlignError> {
        s.trim().parse().map_err(|_| AlignError::TsvParse {
            line: line_no,
            msg: format!("{field} `{s}` is not a float"),
        })
    };
    let parse_f64 = |s: &str, field: &str| -> Result<f64, AlignError> {
        s.trim().parse().map_err(|_| AlignError::TsvParse {
            line: line_no,
            msg: format!("{field} `{s}` is not a float"),
        })
    };
    let parse_i = |s: &str, field: &str| -> Result<i32, AlignError> {
        s.trim().parse().map_err(|_| AlignError::TsvParse {
            line: line_no,
            msg: format!("{field} `{s}` is not an integer"),
        })
    };

    let mut qcov = parse_f32(cols[4], "qcov")?;
    if coverage_is_fraction {
        qcov *= 100.0;
    }
    Ok(Hit {
        qseqid: cols[0].to_string(),
        pident: parse_f32(cols[1], "pident")?,
        evalue: parse_f64(cols[2], "evalue")?,
        bitscore: parse_f32(cols[3], "bitscore")?,
        qcov,
        stitle: cols[5].to_string(),
        sstart: parse_i(cols[6], "sstart")?,
        send: parse_i(cols[7], "send")?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parses_blast_row() {
        let data = "q1\t95.5\t1e-100\t312.5\t90\ttarget1 some desc\t1\t200\n";
        let hits = parse_tsv(Cursor::new(data), false).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].qseqid, "q1");
        assert_eq!(hits[0].pident, 95.5);
        assert_eq!(hits[0].evalue, 1e-100);
        assert_eq!(hits[0].bitscore, 312.5);
        assert_eq!(hits[0].qcov, 90.0);
        assert_eq!(hits[0].stitle, "target1 some desc");
        assert_eq!(hits[0].sstart, 1);
        assert_eq!(hits[0].send, 200);
    }

    #[test]
    fn rescales_mmseqs_fraction_coverage() {
        let data = "q1\t95.5\t1e-100\t312.5\t0.9\ttarget1\t1\t200\n";
        let hits = parse_tsv(Cursor::new(data), true).unwrap();
        assert_eq!(hits[0].qcov, 90.0_f32);
    }

    #[test]
    fn skips_blank_and_comment() {
        let data = "\n#header\nq1\t95.5\t1e-100\t312.5\t90\tt1\t1\t200\n\n";
        let hits = parse_tsv(Cursor::new(data), false).unwrap();
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn error_reports_line_number() {
        let data = "q1\t95.5\t1e-100\t312.5\t90\tt1\t1\t200\nBAD_LINE\n";
        let err = parse_tsv(Cursor::new(data), false).unwrap_err();
        match err {
            AlignError::TsvParse { line, .. } => assert_eq!(line, 2),
            other => panic!("unexpected: {other}"),
        }
    }

    #[test]
    fn stitle_with_internal_spaces_preserved() {
        let data = "q1\t95.5\t1e-100\t312.5\t90\ta b c d\t1\t200\n";
        let h = parse_tsv(Cursor::new(data), false).unwrap();
        assert_eq!(h[0].stitle, "a b c d");
    }
}
