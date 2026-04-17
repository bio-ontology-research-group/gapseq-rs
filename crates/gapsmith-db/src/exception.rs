//! `dat/exception.tbl` loader — enzymes with known false-friends, requiring
//! stricter identity cutoffs (see `src/analyse_alignments.R:108–143`).
//!
//! Columns: `enzyme/reaction, comment`. Lines starting with `#` are comments.

use crate::common::{io_err, DbError};
use serde::{Deserialize, Serialize};
use std::io::BufRead;
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExceptionRow {
    pub id: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub comment: String,
}

pub fn load(path: impl AsRef<Path>) -> Result<Vec<ExceptionRow>, DbError> {
    let path = path.as_ref();
    let f = std::fs::File::open(path).map_err(|e| io_err(path, e))?;
    let rdr = std::io::BufReader::new(f);
    let mut out = Vec::new();
    let mut saw_header = false;
    for line in rdr.lines() {
        let line = line.map_err(|e| io_err(path, e))?;
        let trimmed = line.trim_start();
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }
        if !saw_header {
            saw_header = true;
            // First non-comment line is the header; skip.
            continue;
        }
        let mut cols = line.splitn(2, '\t');
        let id = cols.next().unwrap_or("").trim().to_string();
        let comment = cols.next().unwrap_or("").trim().to_string();
        if !id.is_empty() {
            out.push(ExceptionRow { id, comment });
        }
    }
    tracing::info!(path = %path.display(), rows = out.len(), "loaded exception table");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn skips_comment_and_header() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("e.tbl");
        let mut f = std::fs::File::create(&p).unwrap();
        writeln!(f, "# file with enzymes having false friends").unwrap();
        writeln!(f, "enzyme/reaction\tcomment").unwrap();
        writeln!(f, "7.1.1.9\tcytochrome-c oxidase").unwrap();
        writeln!(f, "1.11.1.6\tcatalase").unwrap();
        let r = load(&p).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].id, "7.1.1.9");
        assert_eq!(r[1].comment, "catalase");
    }
}
