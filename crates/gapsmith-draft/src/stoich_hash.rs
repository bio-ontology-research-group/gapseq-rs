//! Re-export of [`gapsmith_db::rxn_stoich_hash`].
//!
//! The function previously lived here. It now lives in `gapsmith-db` so
//! [`gapsmith_db::SeedRxnRow`] can cache its own hash at load time.
//! We keep the re-export here for API stability with existing callers.

pub use gapsmith_db::rxn_stoich_hash;
