//! Newtype wrappers for the three identifier namespaces used throughout gapseq.
//!
//! Keeping them distinct at the type level prevents passing a reaction id where
//! a compound id is expected. The underlying storage is `String` for now;
//! once the model-loading hot paths are profiled we can swap to interned
//! `Arc<str>` without touching call sites.

use serde::{Deserialize, Serialize};
use std::fmt;

macro_rules! id_newtype {
    ($name:ident, $doc:literal) => {
        #[doc = $doc]
        #[derive(
            Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize,
        )]
        #[serde(transparent)]
        pub struct $name(String);

        impl $name {
            pub fn new(s: impl Into<String>) -> Self {
                Self(s.into())
            }

            pub fn as_str(&self) -> &str {
                &self.0
            }

            pub fn into_string(self) -> String {
                self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }

        impl From<String> for $name {
            fn from(s: String) -> Self {
                Self(s)
            }
        }

        impl From<&str> for $name {
            fn from(s: &str) -> Self {
                Self(s.to_owned())
            }
        }

        impl AsRef<str> for $name {
            fn as_ref(&self) -> &str {
                &self.0
            }
        }
    };
}

id_newtype!(RxnId, "Reaction identifier (e.g. `rxn00148_c0`, `EX_cpd00007_e0`, `bio1`).");
id_newtype!(CpdId, "Compound identifier (e.g. `cpd00001`).");
id_newtype!(GeneId, "Gene identifier (e.g. `fig|83333.1.peg.3`, locus tag).");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distinct_types_do_not_cross() {
        let r = RxnId::new("rxn00001");
        let c = CpdId::new("cpd00001");
        assert_eq!(r.as_str(), "rxn00001");
        assert_eq!(c.as_str(), "cpd00001");
        // The following must NOT compile; kept as a doctest-style reminder:
        // let _: RxnId = c;
    }

    #[test]
    fn serde_transparent_json() {
        let r = RxnId::new("rxn00148");
        let j = serde_json::to_string(&r).unwrap();
        assert_eq!(j, "\"rxn00148\"");
        let back: RxnId = serde_json::from_str(&j).unwrap();
        assert_eq!(r, back);
    }
}
