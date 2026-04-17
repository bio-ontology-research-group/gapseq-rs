//! `gapsmith update-sequences` — sync the reference sequence database from Zenodo.
//!
//! # Example
//!
//! ```bash
//! # Check-only: what's installed vs upstream, no download
//! gapsmith update-sequences -c
//!
//! # Fetch the version this binary pins
//! gapsmith update-sequences -t Bacteria
//!
//! # Fetch the latest published version
//! gapsmith update-sequences -t Bacteria -Z latest
//!
//! # Install into a user-owned dir instead of the system path
//! gapsmith update-sequences -t Archaea -D ~/.gapseq/seq
//! ```
//!
//! Port of `src/update_sequences.sh`. The upstream layout on Zenodo is:
//!
//! ```text
//!   <record>/
//!     md5sums.txt                 — "<md5>  <taxonomy>/<path>" per line
//!     <Taxonomy>.tar.gz           — expands to <Taxonomy>/{rev,unrev,rxn}/sequences.tar.gz
//! ```
//!
//! We replicate the shell script:
//!
//! 1. Resolve `<record>` via the parent concept-DOI redirect (or use a
//!    user-supplied record id).
//! 2. Fetch `md5sums.txt`; diff against local files.
//! 3. If any file is missing or has the wrong checksum, download the
//!    taxonomy archive, extract it, then extract every nested
//!    `rev/unrev/rxn` `sequences.tar.gz`.
//! 4. Write / update `<seqdb>/<Taxonomy>/version_seqDB.json`.
//!
//! Flags:
//!
//! - `-c` check-only: report whether a newer version is available; no download.
//! - `-Z <id|latest>` pin a specific Zenodo record (default: the version
//!   pinned in this binary; `latest` resolves via the concept-DOI redirect).
//! - `-t Bacteria|Archaea` pick a taxonomy. Default `Bacteria`.
//! - `-D <dir>` override the seqdb location.

use clap::Parser;
use flate2::read::GzDecoder;
use md5::{Digest, Md5};
use serde::{Deserialize, Serialize};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

/// Parent concept DOI id. Individual versioned records are linked from
/// this concept record via the Zenodo API. Matches the default in
/// `src/update_sequences.sh`.
const ZENODO_CONCEPT_ID: &str = "10047603";
/// Record id that gapseq 2.0.0 ships against. Used when the user doesn't
/// pass `-Z`. Update when a new upstream release lands.
const ZENODO_PINNED_RECORD: &str = "16908828";

#[derive(Debug, Parser)]
pub struct Args {
    /// Taxonomy database to sync.
    #[arg(long, short = 't', default_value = "Bacteria")]
    pub taxonomy: String,

    /// Seqdb directory (default: `<data-dir>/seq` → `<gapseq>/dat/seq`
    /// via the global --seq-dir flag, if provided).
    #[arg(long, short = 'D')]
    pub seq_dir: Option<PathBuf>,

    /// Zenodo record id or the literal word `latest`. Default: a pinned
    /// record id that ships with this binary.
    #[arg(long, short = 'Z', default_value = ZENODO_PINNED_RECORD)]
    pub record: String,

    /// Check whether a newer database is available; do not download.
    #[arg(long, short = 'c')]
    pub check_only: bool,

    /// Suppress progress messages.
    #[arg(long, short = 'q')]
    pub quiet: bool,
}

pub fn run_cli(
    args: Args,
    data_dir_override: Option<&Path>,
    seq_dir_override: Option<&Path>,
) -> anyhow::Result<()> {
    let seq_dir = resolve_seq_dir(args.seq_dir.as_deref(), data_dir_override, seq_dir_override)?;
    let tax_dir = seq_dir.join(&args.taxonomy);
    std::fs::create_dir_all(tax_dir.join("rev"))?;
    std::fs::create_dir_all(tax_dir.join("unrev"))?;
    std::fs::create_dir_all(tax_dir.join("rxn"))?;

    let say = |msg: &str| {
        if !args.quiet {
            eprintln!("{msg}");
        }
    };

    say(&format!(
        "seqdb: {} (taxonomy: {})",
        seq_dir.display(),
        args.taxonomy
    ));

    if args.check_only {
        return check_for_update(&tax_dir, &args.taxonomy, &say);
    }

    // Resolve the record id (either explicit or via the concept DOI).
    let record_id = if args.record == "latest" {
        resolve_latest_record()?
    } else {
        args.record.clone()
    };
    say(&format!("target Zenodo record: {record_id}"));

    let meta = fetch_record_metadata(&record_id)?;
    verify_community(&meta)?;
    say(&format!(
        "Zenodo record id={} version={} date={}",
        meta.id,
        meta.metadata.version.as_deref().unwrap_or("(na)"),
        meta.created
    ));

    // Determine which files need download.
    let md5_path = format!("records/{record_id}/files/md5sums.txt/content");
    let md5_text = http_get_text(&format!("https://zenodo.org/api/{md5_path}"))?;
    let entries = parse_md5sums(&md5_text, &args.taxonomy);
    if entries.is_empty() {
        anyhow::bail!(
            "no md5sum entries for taxonomy `{}` on Zenodo record {record_id}",
            args.taxonomy
        );
    }
    say(&format!(
        "md5sums.txt: {} files for {}",
        entries.len(),
        args.taxonomy
    ));

    let mut needs_update = false;
    for e in &entries {
        let local = seq_dir.join(&e.relpath);
        if !local.exists() {
            needs_update = true;
            break;
        }
        let got = md5_file(&local)?;
        if got != e.md5 {
            needs_update = true;
            break;
        }
    }

    if !needs_update {
        say(&format!(
            "Reference sequences for {} are already at version {} (Zenodo {}).",
            args.taxonomy,
            meta.metadata.version.as_deref().unwrap_or("(na)"),
            meta.id
        ));
        write_version_stamp(&tax_dir, &meta)?;
        return Ok(());
    }

    say("Local sequences out of date — downloading archive …");
    let archive_url = format!(
        "https://zenodo.org/api/records/{record_id}/files/{}.tar.gz/content",
        args.taxonomy
    );
    let archive_path = seq_dir.join(format!("{}.tar.gz", args.taxonomy));
    download_to(&archive_url, &archive_path, &say)?;
    extract_tar_gz(&archive_path, &seq_dir, &say)?;
    std::fs::remove_file(&archive_path).ok();

    // The archive nests `rev/unrev/rxn` as further tar.gz files. We
    // leave the nested archives in place after extraction — upstream
    // gapseq's `update_sequences.sh` re-checks their md5sums on each
    // run against `md5sums.txt`, so deleting them would trigger a
    // spurious re-download next time.
    for sub in ["rev", "unrev", "rxn"] {
        let nested = tax_dir.join(sub).join("sequences.tar.gz");
        if nested.exists() {
            extract_tar_gz(&nested, &tax_dir.join(sub), &say)?;
        }
    }

    write_version_stamp(&tax_dir, &meta)?;
    say("Reference sequences updated.");
    Ok(())
}

fn resolve_seq_dir(
    explicit: Option<&Path>,
    data_dir_override: Option<&Path>,
    seq_dir_override: Option<&Path>,
) -> anyhow::Result<PathBuf> {
    if let Some(p) = explicit {
        return Ok(p.to_path_buf());
    }
    let data_dir = gapsmith_io::resolve_data_dir(data_dir_override)?;
    Ok(gapsmith_io::resolve_seq_dir(seq_dir_override, &data_dir)?)
}

fn check_for_update(
    tax_dir: &Path,
    taxonomy: &str,
    say: &impl Fn(&str),
) -> anyhow::Result<()> {
    let stamp = tax_dir.join("version_seqDB.json");
    let local = if stamp.exists() {
        let text = std::fs::read_to_string(&stamp)?;
        serde_json::from_str::<VersionStamp>(&text).ok()
    } else {
        None
    };
    match &local {
        Some(v) => say(&format!(
            "Local {taxonomy} sequence database: {} (zenodoID: {}, date: {})",
            v.version.as_deref().unwrap_or("(na)"),
            v.zenodo_id,
            v.date
        )),
        None => say(&format!(
            "No current {taxonomy} sequence database version found. Run `gapsmith update-sequences -t {taxonomy}`."
        )),
    };

    let latest_id = resolve_latest_record()?;
    let remote = fetch_record_metadata(&latest_id)?;
    verify_community(&remote)?;
    match local {
        Some(v) if v.zenodo_id == remote.id.to_string() => {
            say("Reference sequences are up-to-date.");
        }
        _ => say(&format!(
            "[NOTE] A newer {taxonomy} sequence database exists: {} (zenodoID: {}, date: {}).\n[NOTE] Run `gapsmith update-sequences -t {taxonomy} -Z latest` to update.",
            remote.metadata.version.as_deref().unwrap_or("(na)"),
            remote.id,
            remote.created,
        )),
    }
    Ok(())
}

/// Resolve `latest` via the concept-DOI redirect. We follow the
/// `https://zenodo.org/api/records/<concept>` redirect chain and extract
/// the final record id.
fn resolve_latest_record() -> anyhow::Result<String> {
    let url = format!("https://zenodo.org/api/records/{ZENODO_CONCEPT_ID}");
    let agent = http_agent();
    let resp = agent.get(&url).call().map_err(|e| anyhow::anyhow!("Zenodo concept fetch failed: {e}"))?;
    let final_url = resp.get_url().to_string();
    let id = final_url
        .rsplit('/')
        .find(|s| !s.is_empty())
        .ok_or_else(|| anyhow::anyhow!("couldn't parse record id from `{final_url}`"))?
        .to_string();
    Ok(id)
}

#[derive(Debug, Deserialize)]
struct ZenodoRecord {
    id: u64,
    created: String,
    metadata: ZenodoMeta,
}

#[derive(Debug, Deserialize)]
struct ZenodoMeta {
    #[serde(default)]
    version: Option<String>,
    #[serde(default)]
    communities: Vec<ZenodoCommunity>,
}

#[derive(Debug, Deserialize)]
struct ZenodoCommunity {
    #[serde(default)]
    id: String,
}

fn fetch_record_metadata(record_id: &str) -> anyhow::Result<ZenodoRecord> {
    let url = format!("https://zenodo.org/api/records/{record_id}");
    let text = http_get_text(&url)?;
    let rec: ZenodoRecord = serde_json::from_str(&text)?;
    Ok(rec)
}

fn verify_community(meta: &ZenodoRecord) -> anyhow::Result<()> {
    let ok = meta.metadata.communities.iter().any(|c| c.id == "gapseq");
    if !ok {
        anyhow::bail!("Zenodo record {} is not in the `gapseq` community", meta.id);
    }
    Ok(())
}

#[derive(Debug, Serialize, Deserialize)]
struct VersionStamp {
    #[serde(rename = "zenodoID")]
    zenodo_id: String,
    #[serde(default)]
    version: Option<String>,
    #[serde(default)]
    date: String,
}

fn write_version_stamp(tax_dir: &Path, meta: &ZenodoRecord) -> anyhow::Result<()> {
    let stamp = VersionStamp {
        zenodo_id: meta.id.to_string(),
        version: meta.metadata.version.clone(),
        date: meta.created.clone(),
    };
    let path = tax_dir.join("version_seqDB.json");
    let out = serde_json::to_string(&stamp)?;
    std::fs::write(&path, out)?;
    Ok(())
}

/// `md5sums.txt` is two columns: `<md5>  <taxonomy>/<path>`. We keep only
/// entries whose path starts with the requested taxonomy.
#[derive(Debug)]
struct Md5Entry {
    md5: String,
    relpath: String,
}

fn parse_md5sums(text: &str, taxonomy: &str) -> Vec<Md5Entry> {
    let prefix = format!("{taxonomy}/");
    text.lines()
        .filter_map(|line| {
            let line = line.trim();
            if line.is_empty() {
                return None;
            }
            let mut it = line.splitn(2, char::is_whitespace);
            let md5 = it.next()?.trim().to_string();
            let raw = it.next()?.trim();
            // Zenodo's md5sums.txt uses `./Bacteria/...`. Strip the leading
            // `./` so callers can resolve against the seqdb root.
            let path = raw.trim_start_matches("./").to_string();
            if !path.starts_with(&prefix) {
                return None;
            }
            Some(Md5Entry { md5, relpath: path })
        })
        .collect()
}

fn md5_file(path: &Path) -> std::io::Result<String> {
    let mut f = std::fs::File::open(path)?;
    let mut hasher = Md5::new();
    let mut buf = [0u8; 1 << 16];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    let digest = hasher.finalize();
    let mut s = String::with_capacity(32);
    for b in digest.iter() {
        s.push_str(&format!("{b:02x}"));
    }
    Ok(s)
}

fn http_agent() -> ureq::Agent {
    ureq::AgentBuilder::new()
        .redirects(16)
        .timeout_connect(std::time::Duration::from_secs(30))
        .timeout_read(std::time::Duration::from_secs(300))
        .build()
}

fn http_get_text(url: &str) -> anyhow::Result<String> {
    let agent = http_agent();
    let resp = agent
        .get(url)
        .call()
        .map_err(|e| anyhow::anyhow!("GET {url} failed: {e}"))?;
    let text = resp.into_string()?;
    Ok(text)
}

fn download_to(url: &str, dst: &Path, say: &impl Fn(&str)) -> anyhow::Result<()> {
    say(&format!("  fetching {url}"));
    let agent = http_agent();
    let resp = agent
        .get(url)
        .call()
        .map_err(|e| anyhow::anyhow!("GET {url} failed: {e}"))?;
    let len = resp
        .header("content-length")
        .and_then(|s| s.parse::<u64>().ok());
    say(&format!(
        "  writing {} ({})",
        dst.display(),
        len.map(|n| format!("{} bytes", n)).unwrap_or_else(|| "size unknown".into())
    ));
    let mut reader = resp.into_reader();
    let mut f = std::fs::File::create(dst)?;
    let mut buf = [0u8; 1 << 16];
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        f.write_all(&buf[..n])?;
    }
    Ok(())
}

fn extract_tar_gz(src: &Path, into: &Path, say: &impl Fn(&str)) -> anyhow::Result<()> {
    say(&format!("  extracting {} → {}", src.display(), into.display()));
    let f = std::fs::File::open(src)?;
    let dec = GzDecoder::new(f);
    let mut ar = tar::Archive::new(dec);
    ar.unpack(into)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_md5sums_filters_by_taxonomy() {
        let text = "\
abc123  Bacteria/rev/sequences.tar.gz
def456  Bacteria/rxn/sequences.tar.gz
ghi789  Archaea/rev/sequences.tar.gz
";
        let v = parse_md5sums(text, "Bacteria");
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].md5, "abc123");
        assert_eq!(v[0].relpath, "Bacteria/rev/sequences.tar.gz");
    }

    #[test]
    fn parse_md5sums_strips_dot_slash_prefix() {
        // Zenodo emits paths with `./` prefix; we need to canonicalise.
        let text = "\
ae174b52d48d9a82f0da1319f7ccc1e5  ./Bacteria/rev/sequences.tar.gz
d29bfc01e1a91cb6e8ae056a3e886ad7  ./Archaea/rxn/sequences.tar.gz
";
        let v = parse_md5sums(text, "Bacteria");
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].relpath, "Bacteria/rev/sequences.tar.gz");
    }

    #[test]
    fn md5_matches_reference() {
        use std::io::Write;
        let td = tempfile::tempdir().unwrap();
        let p = td.path().join("a.txt");
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(b"hello\n").unwrap();
        drop(f);
        // md5 of "hello\n"
        assert_eq!(md5_file(&p).unwrap(), "b1946ac92492d2347c6235b4d2611184");
    }
}
