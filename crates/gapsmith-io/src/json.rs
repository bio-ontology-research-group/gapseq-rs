//! JSON (de)serialization for [`Model`].
//!
//! Primary use: inspection and third-party tooling. CBOR remains the default
//! on-disk format.

use gapsmith_core::Model;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum JsonError {
    #[error("i/o error on `{path}`: {source}")]
    Io { path: String, #[source] source: std::io::Error },
    #[error("JSON (de)serialization error: {0}")]
    Serde(#[from] serde_json::Error),
}

pub fn write_model_json(
    model: &Model,
    path: impl AsRef<Path>,
    pretty: bool,
) -> Result<(), JsonError> {
    let p = path.as_ref();
    let f = File::create(p).map_err(|e| JsonError::Io {
        path: p.display().to_string(),
        source: e,
    })?;
    let mut w = BufWriter::new(f);
    if pretty {
        serde_json::to_writer_pretty(&mut w, model)?;
    } else {
        serde_json::to_writer(&mut w, model)?;
    }
    Ok(())
}

pub fn read_model_json(path: impl AsRef<Path>) -> Result<Model, JsonError> {
    let p = path.as_ref();
    let f = File::open(p).map_err(|e| JsonError::Io {
        path: p.display().to_string(),
        source: e,
    })?;
    let r = BufReader::new(f);
    let m: Model = serde_json::from_reader(r)?;
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gapsmith_core::{CompartmentId, Metabolite, Reaction, StoichMatrix};

    #[test]
    fn json_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("toy.json");
        let mut m = Model::new("toy");
        m.mets.push(Metabolite::new("cpdA", "A", CompartmentId::CYTOSOL));
        m.rxns.push(Reaction::new("r1", "drain", 0.0, 1000.0));
        m.s = StoichMatrix::from_triplets(1, 1, vec![(0, 0, -1.0)]);
        write_model_json(&m, &p, true).unwrap();
        let back = read_model_json(&p).unwrap();
        assert_eq!(back.met_count(), 1);
        assert_eq!(back.rxn_count(), 1);
        back.check_shape().unwrap();
    }
}
