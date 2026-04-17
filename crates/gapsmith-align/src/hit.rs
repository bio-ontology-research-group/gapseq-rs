//! The canonical hit structure produced by every aligner backend.

use serde::{Deserialize, Serialize};

/// A single alignment hit, normalized across the three supported aligners
/// plus the pre-computed TSV path.
///
/// Coverage and identity are always on the 0–100 scale, regardless of
/// whether the source tool reported them as fractions or percentages.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Hit {
    pub qseqid: String,
    pub pident: f32,
    pub evalue: f64,
    pub bitscore: f32,
    pub qcov: f32,
    pub stitle: String,
    pub sstart: i32,
    pub send: i32,
}
