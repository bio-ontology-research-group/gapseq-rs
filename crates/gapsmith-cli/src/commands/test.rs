//! `gapsmith test` — print resolved paths and sanity-check the environment.
//!
//! Bash equivalent: `src/test.sh`. The Rust version cannot check R packages,
//! but it reports which external binaries are on PATH and where reference
//! data resolved to.

use clap::Parser;
use gapsmith_io::{resolve_data_dir, resolve_seq_dir};
use std::path::Path;

#[derive(Debug, Parser)]
pub struct Args {}

pub fn run(_args: Args, data_dir_override: Option<&Path>) -> anyhow::Result<()> {
    println!("gapsmith {}", env!("CARGO_PKG_VERSION"));

    match resolve_data_dir(data_dir_override) {
        Ok(p) => {
            println!("data dir: {}", p.display());
            match resolve_seq_dir(None, &p) {
                Ok(s) => println!("seq dir:  {}", s.display()),
                Err(e) => println!("seq dir:  <unresolved: {e}>"),
            }
        }
        Err(e) => println!("data dir: <unresolved: {e}>"),
    }

    println!();
    println!("external tools on PATH:");
    for tool in [
        "blastp", "tblastn", "makeblastdb", "diamond", "mmseqs", "prodigal", "barrnap",
        "hmmsearch", "bedtools",
    ] {
        let found = which(tool);
        println!(
            "  {:<12} {}",
            tool,
            found.map(|p| p.display().to_string()).unwrap_or_else(|| "<missing>".into())
        );
    }

    Ok(())
}

/// Minimal `which` — splits `$PATH` and checks each entry for an executable
/// file with the given name. Avoids pulling in the `which` crate for a
/// one-page implementation.
fn which(name: &str) -> Option<std::path::PathBuf> {
    let path = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path) {
        let candidate = dir.join(name);
        if is_executable(&candidate) {
            return Some(candidate);
        }
    }
    None
}

#[cfg(unix)]
fn is_executable(p: &std::path::Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    p.metadata()
        .map(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
        .unwrap_or(false)
}

#[cfg(not(unix))]
fn is_executable(p: &std::path::Path) -> bool {
    p.is_file()
}
