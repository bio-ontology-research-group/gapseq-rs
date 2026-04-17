//! Data-directory resolution.
//!
//! Resolution order (first match wins):
//!
//! 1. Explicit `override` argument (`--data-dir` flag).
//! 2. `GAPSMITH_DATA_DIR` environment variable.
//! 3. `$XDG_DATA_HOME/gapsmith` (or `~/.local/share/gapsmith`).
//! 4. `<exe-dir>/../share/gapsmith/dat` (for system packaging).
//! 5. `./dat` (developer fallback in the checkout).
//!
//! A candidate is accepted if it is a directory and contains a sentinel file
//! (`seed_reactions_corrected.tsv`). Use [`resolve_data_dir`] for the main
//! reference-data root and [`resolve_seq_dir`] for the sequence subtree
//! (defaults to `<data-dir>/seq`, separately overridable because it is large
//! and often lives on a different disk).

use std::path::{Path, PathBuf};

#[derive(Debug, thiserror::Error)]
pub enum PathResolveError {
    #[error("could not locate gapsmith data directory. Tried: {candidates:?}")]
    NotFound { candidates: Vec<PathBuf> },
}

/// Sentinel file that identifies a valid `dat/` root.
const DATA_SENTINEL: &str = "seed_reactions_corrected.tsv";
/// Sentinel subdirectory for `dat/seq/`.
const SEQ_SENTINEL: &str = "Bacteria";

pub fn resolve_data_dir(override_path: Option<&Path>) -> Result<PathBuf, PathResolveError> {
    let mut tried: Vec<PathBuf> = Vec::new();
    let try_candidate = |c: PathBuf, tried: &mut Vec<PathBuf>| -> Option<PathBuf> {
        if c.join(DATA_SENTINEL).is_file() {
            Some(c)
        } else {
            tried.push(c);
            None
        }
    };

    if let Some(p) = override_path {
        if let Some(hit) = try_candidate(p.to_path_buf(), &mut tried) {
            return Ok(hit);
        }
    }
    if let Some(p) = std::env::var_os("GAPSMITH_DATA_DIR") {
        if let Some(hit) = try_candidate(PathBuf::from(p), &mut tried) {
            return Ok(hit);
        }
    }
    if let Some(xdg) = std::env::var_os("XDG_DATA_HOME") {
        let mut p = PathBuf::from(xdg);
        p.push("gapsmith");
        if let Some(hit) = try_candidate(p, &mut tried) {
            return Ok(hit);
        }
    } else if let Some(home) = std::env::var_os("HOME") {
        let mut p = PathBuf::from(home);
        p.push(".local/share/gapsmith");
        if let Some(hit) = try_candidate(p, &mut tried) {
            return Ok(hit);
        }
    }
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let mut p = dir.to_path_buf();
            p.push("../share/gapsmith/dat");
            if let Some(hit) = try_candidate(p, &mut tried) {
                return Ok(hit);
            }
        }
    }
    if let Ok(cwd) = std::env::current_dir() {
        let c = cwd.join("dat");
        if let Some(hit) = try_candidate(c, &mut tried) {
            return Ok(hit);
        }
        // Also try `../dat` (when running from inside `gapsmith/`).
        let c = cwd.join("../dat");
        if let Some(hit) = try_candidate(c, &mut tried) {
            return Ok(hit);
        }
    }
    Err(PathResolveError::NotFound { candidates: tried })
}

/// Resolve the sequence-database subtree. Order: override → `GAPSMITH_SEQ_DIR`
/// → `<data-dir>/seq`.
pub fn resolve_seq_dir(
    override_path: Option<&Path>,
    data_dir: &Path,
) -> Result<PathBuf, PathResolveError> {
    let mut tried: Vec<PathBuf> = Vec::new();
    let try_candidate = |c: PathBuf, tried: &mut Vec<PathBuf>| -> Option<PathBuf> {
        if c.join(SEQ_SENTINEL).is_dir() {
            Some(c)
        } else {
            tried.push(c);
            None
        }
    };
    if let Some(p) = override_path {
        if let Some(hit) = try_candidate(p.to_path_buf(), &mut tried) {
            return Ok(hit);
        }
    }
    if let Some(p) = std::env::var_os("GAPSMITH_SEQ_DIR") {
        if let Some(hit) = try_candidate(PathBuf::from(p), &mut tried) {
            return Ok(hit);
        }
    }
    let c = data_dir.join("seq");
    if let Some(hit) = try_candidate(c, &mut tried) {
        return Ok(hit);
    }
    Err(PathResolveError::NotFound { candidates: tried })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn explicit_override_wins() {
        let d = tempfile::tempdir().unwrap();
        fs::write(d.path().join(DATA_SENTINEL), "").unwrap();
        let p = resolve_data_dir(Some(d.path())).unwrap();
        assert_eq!(p, d.path());
    }

    #[test]
    fn missing_sentinel_fails() {
        let d = tempfile::tempdir().unwrap();
        let r = resolve_data_dir(Some(d.path()));
        assert!(r.is_err());
    }

    #[test]
    fn seq_dir_resolves_under_data_dir() {
        let d = tempfile::tempdir().unwrap();
        fs::write(d.path().join(DATA_SENTINEL), "").unwrap();
        fs::create_dir_all(d.path().join("seq").join(SEQ_SENTINEL)).unwrap();
        let seq = resolve_seq_dir(None, d.path()).unwrap();
        assert_eq!(seq, d.path().join("seq"));
    }
}
