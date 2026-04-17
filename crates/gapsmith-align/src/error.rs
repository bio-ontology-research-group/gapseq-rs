//! Error type shared by every aligner backend.

use std::path::PathBuf;
use std::process::ExitStatus;

#[derive(Debug, thiserror::Error)]
pub enum AlignError {
    #[error("external tool `{tool}` not found on PATH")]
    ToolMissing { tool: &'static str },

    #[error("external tool `{tool}` exited with {status}: {stderr}")]
    ToolFailed {
        tool: &'static str,
        status: ExitStatus,
        stderr: String,
    },

    #[error("i/o error on `{path}`: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("TSV parse error on line {line}: {msg}")]
    TsvParse { line: u64, msg: String },

    #[error("unsupported argument: {0}")]
    BadArg(String),
}

pub(crate) fn io_err(path: &std::path::Path, source: std::io::Error) -> AlignError {
    AlignError::Io { path: path.to_path_buf(), source }
}
