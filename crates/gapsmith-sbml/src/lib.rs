//! SBML Level 3 Version 1 writer with FBC v2 and groups v1 extensions.
//!
//! Output is deliberately COBRApy-compatible: species get `M_` prefixes,
//! reactions `R_`, gene products `G_`. The `fbc:strict="true"` flag is set,
//! so every reaction must have `fbc:lowerFluxBound` and `fbc:upperFluxBound`
//! pointing at a `<parameter>`. Default-bound parameters are emitted once
//! (`cobra_default_lb`, `cobra_default_ub`, `cobra_0_bound`); reactions
//! whose bounds differ from defaults get dedicated per-reaction parameters.
//!
//! Only a writer is implemented in this milestone. Reading SBML is not
//! needed by gapsmith itself — every input comes from the CBOR/JSON
//! pipeline.

mod namespaces;
mod writer;

pub use writer::{write_sbml, write_to, ObjectiveSense, SbmlError, WriteOptions};

/// Default upper/lower bound magnitudes used when a reaction's own bounds
/// equal these and a shared parameter is preferred. Mirrors cobrar's choice.
pub const DEFAULT_BOUND: f64 = 1000.0;
