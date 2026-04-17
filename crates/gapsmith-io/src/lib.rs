//! I/O for gapsmith.
//!
//! Currently provides:
//!
//! - [`cbor`]: CBOR (de)serialization of [`Model`](gapsmith_core::Model) via `ciborium`.
//!   Filename convention: `<id>.gmod.cbor`. Replaces R's `.RDS`.
//! - [`json`]: JSON (de)serialization (useful for inspection / debugging).
//! - [`paths`]: data-directory resolver (`--data-dir` flag → `GAPSMITH_DATA_DIR`
//!   env → XDG → exe-sibling → cwd).

pub mod cbor;
pub mod json;
pub mod paths;

pub use cbor::{read_model_cbor, write_model_cbor, CborError};
pub use json::{read_model_json, write_model_json, JsonError};
pub use paths::{resolve_data_dir, resolve_seq_dir, PathResolveError};

/// File-format detected from extension.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ModelFormat {
    Cbor,
    Json,
}

impl ModelFormat {
    /// Detect from path extension. `.gmod.cbor`, `.cbor` → CBOR; `.json` → JSON.
    /// Unknown extensions default to CBOR.
    pub fn from_path(path: &std::path::Path) -> Self {
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.ends_with(".json") {
            Self::Json
        } else {
            Self::Cbor
        }
    }
}
