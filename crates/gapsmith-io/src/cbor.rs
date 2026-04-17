//! CBOR (de)serialization for [`Model`].
//!
//! CBOR replaces R's RDS as the native on-disk representation. It is
//! self-describing, compact (~30% of equivalent JSON), and roundtrips cleanly
//! through serde — so `sprs::CsMat`, `Arc<str>`, and nested enums all work
//! without custom visitors.

use gapsmith_core::Model;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum CborError {
    #[error("i/o error on `{path}`: {source}")]
    Io { path: String, #[source] source: std::io::Error },
    #[error("CBOR encoding error: {0}")]
    Encode(#[from] ciborium::ser::Error<std::io::Error>),
    #[error("CBOR decoding error: {0}")]
    Decode(#[from] ciborium::de::Error<std::io::Error>),
}

pub fn write_model_cbor(model: &Model, path: impl AsRef<Path>) -> Result<(), CborError> {
    let p = path.as_ref();
    let f = File::create(p).map_err(|e| CborError::Io {
        path: p.display().to_string(),
        source: e,
    })?;
    let mut w = BufWriter::new(f);
    ciborium::ser::into_writer(model, &mut w)?;
    Ok(())
}

pub fn read_model_cbor(path: impl AsRef<Path>) -> Result<Model, CborError> {
    let p = path.as_ref();
    let f = File::open(p).map_err(|e| CborError::Io {
        path: p.display().to_string(),
        source: e,
    })?;
    let r = BufReader::new(f);
    let m: Model = ciborium::de::from_reader(r)?;
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gapsmith_core::{CompartmentId, Metabolite, Reaction, StoichMatrix};

    fn toy() -> Model {
        let mut m = Model::new("toy");
        m.mets.push(Metabolite::new("cpdA", "A", CompartmentId::CYTOSOL));
        m.mets.push(Metabolite::new("cpdB", "B", CompartmentId::CYTOSOL));
        m.rxns.push(Reaction::new("r1", "A -> B", 0.0, 1000.0));
        m.s = StoichMatrix::from_triplets(2, 1, vec![(0, 0, -1.0), (1, 0, 1.0)]);
        m
    }

    #[test]
    fn cbor_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("toy.gmod.cbor");
        let m = toy();
        write_model_cbor(&m, &p).unwrap();
        let back = read_model_cbor(&p).unwrap();
        assert_eq!(back.annot.id, "toy");
        assert_eq!(back.met_count(), 2);
        assert_eq!(back.rxn_count(), 1);
        assert_eq!(back.s.nnz(), 2);
        back.check_shape().unwrap();
    }
}
