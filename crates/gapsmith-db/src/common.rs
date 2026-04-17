//! Shared error type plus tiny helpers used across loaders.

use std::path::{Path, PathBuf};

#[derive(Debug, thiserror::Error)]
pub enum DbError {
    #[error("i/o error on `{path}`: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("CSV error on `{path}` at line {line:?}: {source}")]
    Csv {
        path: PathBuf,
        line: Option<u64>,
        #[source]
        source: csv::Error,
    },
    #[error("JSON error on `{path}`: {source}")]
    Json {
        path: PathBuf,
        #[source]
        source: serde_json::Error,
    },
    #[error("parse error on `{path}` line {line}: {msg}")]
    Parse { path: PathBuf, line: u64, msg: String },
    #[error("stoichiometry parse error: {0}")]
    Stoich(#[from] crate::stoich_parse::StoichParseError),
    #[error("biomass error: {0}")]
    Biomass(#[from] crate::biomass::BiomassError),
    #[error("file not found: `{0}`")]
    NotFound(PathBuf),
}

pub(crate) fn io_err(path: &Path, source: std::io::Error) -> DbError {
    DbError::Io { path: path.to_path_buf(), source }
}

pub(crate) fn csv_err(path: &Path, source: csv::Error) -> DbError {
    let line = source.position().map(|p| p.line());
    DbError::Csv { path: path.to_path_buf(), line, source }
}
