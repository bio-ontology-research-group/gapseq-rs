//! Biomass JSON template loader.
//!
//! Input: `dat/biomass/biomass_{Gram_pos,Gram_neg,archaea}.json`.
//!
//! Schema (inferred from the real files — no schema document ships with
//! gapseq):
//!
//! ```json
//! {
//!   "id":        "Gram_neg",
//!   "name":      "Bacterial Gram-negative biomass reaction",
//!   "ref":       "derived from ...",
//!   "energy_GAM": 40,
//!   "domain":    "Bacteria",
//!   "met_groups": [
//!     {
//!       "group_name":      "DNA",
//!       "mass":            0.031,
//!       "unit_group":      "g",
//!       "unit_components": "MOLFRACTION",
//!       "components": [
//!         {"id":"cpd00115","name":"dATP","comp":"c","coef":0.246,
//!          "link":"cpd00012:-1"}
//!       ]
//!     }
//!   ]
//! }
//! ```
//!
//! Some fields carry gapseq-specific quirks (see `src/parse_BMjson.R:1–108`):
//!
//! - `link` is optional and empty-string means "no coupled product".
//!   It encodes a coupled metabolite: `"<cpd>:<coef>"`.
//! - `comp` is a single character (`c`, `e`, `p`) that maps to the standard
//!   compartments at model-build time.
//! - `unit_group` is the measurement unit for the *group* mass; gapseq
//!   assumes `"g"`.
//! - `unit_components` describes how component coefficients combine:
//!   `"MOLFRACTION"` → per-mole fraction, others exist but are rare.

use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum BiomassError {
    #[error("i/o error on `{path}`: {source}")]
    Io {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("JSON parse error on `{path}`: {source}")]
    Json {
        path: std::path::PathBuf,
        #[source]
        source: serde_json::Error,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiomassComponent {
    pub id: String,
    #[serde(default)]
    pub name: String,
    /// Compartment code: `c` (cytosol), `e` (extracellular), `p` (periplasm).
    pub comp: String,
    pub coef: f64,
    /// Optional coupled metabolite, encoded as `"<cpd>:<coef>"`. Empty when absent.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub link: String,
}

impl BiomassComponent {
    /// Decode the `link` field into (metabolite id, stoichiometric coefficient).
    ///
    /// Single link form. For entries like `"cpd01997:-1|cpd03422:-1"` with
    /// multiple coupled metabolites, use [`Self::links`] instead.
    pub fn link(&self) -> Option<(&str, f64)> {
        self.links().into_iter().next()
    }

    /// Decode the `link` field into every `(cpd, coef)` it encodes. The
    /// `|`-separator is how gapseq chains multiple couplings on one
    /// biomass component (e.g., Calomide → consume both cpd01997 and
    /// cpd03422).
    pub fn links(&self) -> Vec<(&str, f64)> {
        if self.link.is_empty() {
            return Vec::new();
        }
        self.link
            .split('|')
            .filter_map(|term| {
                let (cpd, coef) = term.split_once(':')?;
                coef.trim().parse::<f64>().ok().map(|c| (cpd.trim(), c))
            })
            .collect()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiomassGroup {
    pub group_name: String,
    pub mass: f64,
    #[serde(default)]
    pub unit_group: String,
    #[serde(default)]
    pub unit_components: String,
    #[serde(default)]
    pub components: Vec<BiomassComponent>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiomassTemplate {
    pub id: String,
    pub name: String,
    #[serde(default, rename = "ref", skip_serializing_if = "String::is_empty")]
    pub reference: String,
    #[serde(default, rename = "energy_GAM")]
    pub energy_gam: f64,
    #[serde(default)]
    pub domain: String,
    #[serde(default)]
    pub met_groups: Vec<BiomassGroup>,
}

impl BiomassTemplate {
    pub fn load(path: impl AsRef<Path>) -> Result<Self, BiomassError> {
        let path = path.as_ref();
        let f = std::fs::File::open(path).map_err(|e| BiomassError::Io {
            path: path.to_path_buf(),
            source: e,
        })?;
        let r = std::io::BufReader::new(f);
        let t: BiomassTemplate = serde_json::from_reader(r).map_err(|e| BiomassError::Json {
            path: path.to_path_buf(),
            source: e,
        })?;
        Ok(t)
    }

    /// Load returning `None` if the file doesn't exist. Useful when the
    /// user's `dat/` root only ships a subset of biomass templates.
    pub fn load_opt(path: impl AsRef<Path>) -> Result<Option<Self>, BiomassError> {
        let path = path.as_ref();
        if !path.exists() {
            tracing::warn!(path = %path.display(), "biomass template not found; skipping");
            return Ok(None);
        }
        Self::load(path).map(Some)
    }

    /// Return every component across all groups (convenience iterator).
    pub fn iter_components(&self) -> impl Iterator<Item = (&BiomassGroup, &BiomassComponent)> {
        self.met_groups
            .iter()
            .flat_map(|g| g.components.iter().map(move |c| (g, c)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const MINIMAL_JSON: &str = r#"
    { "id" : "Gram_neg",
      "name" : "Bacterial Gram-negative biomass reaction",
      "ref"  : "test",
      "energy_GAM" : 40,
      "domain" : "Bacteria",
      "met_groups" : [
        { "group_name" : "DNA",
          "mass" : 0.031,
          "unit_group" : "g",
          "unit_components" : "MOLFRACTION",
          "components" : [
            { "id":"cpd00115","name":"dATP","comp":"c","coef":0.246,"link":"cpd00012:-1" },
            { "id":"cpd00357","name":"dTTP","comp":"c","coef":0.246 }
          ]
        }
      ]
    }"#;

    #[test]
    fn parses_minimal_biomass() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("bm.json");
        std::fs::File::create(&p).unwrap().write_all(MINIMAL_JSON.as_bytes()).unwrap();
        let t = BiomassTemplate::load(&p).unwrap();
        assert_eq!(t.id, "Gram_neg");
        assert_eq!(t.met_groups.len(), 1);
        assert_eq!(t.met_groups[0].components.len(), 2);
        let (cpd, c) = t.met_groups[0].components[0].link().unwrap();
        assert_eq!(cpd, "cpd00012");
        assert_eq!(c, -1.0);
        assert!(t.met_groups[0].components[1].link.is_empty());
        assert!(t.met_groups[0].components[1].link().is_none());
    }

    #[test]
    fn load_opt_missing_file() {
        let d = tempfile::tempdir().unwrap();
        let t = BiomassTemplate::load_opt(d.path().join("missing.json")).unwrap();
        assert!(t.is_none());
    }
}
