//! Hit classification — port of `src/analyse_alignments.R:108–143`.
//!
//! Inputs:
//! - per-reaction reference-FASTA presence (did we have sequences to search?)
//! - raw alignment hits for this reaction (possibly many, possibly none)
//! - bitscore / identity / coverage / exception-EC cutoffs
//!
//! Output: a [`HitStatus`] plus the best hit (largest bitscore) when present.

use crate::types::HitStatus;
use gapsmith_align::Hit;
use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct ClassifyOptions<'a> {
    pub bitcutoff: f32,
    pub identcutoff: f32,
    /// Identity cutoff applied only when the reaction's EC is listed in
    /// `dat/exception.tbl` (known "false-friend" enzymes).
    pub ident_exception: f32,
    /// Precomputed set of EC strings that require the exception cutoff.
    pub exception_ecs: &'a HashSet<String>,
}

/// Classify a single reaction's hits. Returns `(status, is_exception,
/// best_hit_index)` — the index points back into `hits` so the caller can
/// attach the winning hit's fields to the reaction row.
pub fn classify_reaction(
    hits: &[Hit],
    ec: &str,
    has_seq_data: bool,
    spont: bool,
    opts: &ClassifyOptions<'_>,
) -> (HitStatus, bool, Option<usize>) {
    let exception = ec_is_exception(ec, opts.exception_ecs);

    if hits.is_empty() {
        let status = if !has_seq_data && spont {
            HitStatus::Spontaneous
        } else if has_seq_data {
            HitStatus::NoBlast
        } else {
            HitStatus::NoSeqData
        };
        return (status, exception, None);
    }

    // Pick the top-bitscore hit. In the R pipeline ties would sort by
    // (-bitscore, stitle, complex) after the complex join; we keep a
    // single representative for M5's no-complex regime.
    let (idx, best) = hits
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.bitscore.partial_cmp(&b.1.bitscore).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap();

    let pass_main = best.bitscore >= opts.bitcutoff && best.pident >= opts.identcutoff;
    let pass_exc = if exception { best.pident >= opts.ident_exception } else { true };

    let status = if pass_main && pass_exc {
        HitStatus::GoodBlast
    } else {
        HitStatus::BadBlast
    };
    (status, exception, Some(idx))
}

/// Decide whether an EC belongs in the exception table. Handles both
/// exact-match EC codes and EC/EC-style entries (`enzyme/reaction` column
/// in `dat/exception.tbl`).
pub fn ec_is_exception(ec: &str, exception_set: &HashSet<String>) -> bool {
    for token in ec.split('/').map(str::trim).filter(|s| !s.is_empty()) {
        if exception_set.contains(token) {
            return true;
        }
    }
    false
}

/// Legacy one-shot classifier preserved for API callers that operate on a
/// slice of `(bitscore, pident)` pairs. Not used by the runner.
pub fn classify_hits(
    scores: &[(f32, f32)],
    ec: &str,
    has_seq_data: bool,
    spont: bool,
    opts: &ClassifyOptions<'_>,
) -> HitStatus {
    let exception = ec_is_exception(ec, opts.exception_ecs);
    if scores.is_empty() {
        if !has_seq_data && spont {
            return HitStatus::Spontaneous;
        }
        if has_seq_data {
            return HitStatus::NoBlast;
        }
        return HitStatus::NoSeqData;
    }
    let (b, p) = scores
        .iter()
        .cloned()
        .fold((f32::MIN, f32::MIN), |(ab, ap), (nb, np)| (ab.max(nb), ap.max(np)));
    if b >= opts.bitcutoff && p >= opts.identcutoff {
        if exception && p < opts.ident_exception {
            HitStatus::BadBlast
        } else {
            HitStatus::GoodBlast
        }
    } else {
        HitStatus::BadBlast
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(exc: &HashSet<String>) -> ClassifyOptions<'_> {
        ClassifyOptions {
            bitcutoff: 200.0,
            identcutoff: 30.0,
            ident_exception: 70.0,
            exception_ecs: exc,
        }
    }

    fn hit(bit: f32, pid: f32) -> Hit {
        Hit {
            qseqid: "q".into(),
            pident: pid,
            evalue: 0.0,
            bitscore: bit,
            qcov: 100.0,
            stitle: "t".into(),
            sstart: 1,
            send: 100,
        }
    }

    #[test]
    fn good_blast_above_cutoffs() {
        let e = HashSet::new();
        let (s, _, _) = classify_reaction(&[hit(250.0, 50.0)], "1.1.1.1", true, false, &opts(&e));
        assert_eq!(s, HitStatus::GoodBlast);
    }

    #[test]
    fn bad_blast_below_bitscore() {
        let e = HashSet::new();
        let (s, _, _) = classify_reaction(&[hit(100.0, 50.0)], "1.1.1.1", true, false, &opts(&e));
        assert_eq!(s, HitStatus::BadBlast);
    }

    #[test]
    fn exception_requires_higher_identity() {
        let mut e = HashSet::new();
        e.insert("7.1.1.9".into());
        let (s, exc, _) =
            classify_reaction(&[hit(250.0, 50.0)], "7.1.1.9", true, false, &opts(&e));
        assert_eq!(s, HitStatus::BadBlast);
        assert!(exc);
        let (s, _, _) = classify_reaction(&[hit(250.0, 80.0)], "7.1.1.9", true, false, &opts(&e));
        assert_eq!(s, HitStatus::GoodBlast);
    }

    #[test]
    fn no_blast_when_seq_data_but_no_hit() {
        let e = HashSet::new();
        let (s, _, _) = classify_reaction(&[], "1.1.1.1", true, false, &opts(&e));
        assert_eq!(s, HitStatus::NoBlast);
    }

    #[test]
    fn no_seq_data_when_nothing() {
        let e = HashSet::new();
        let (s, _, _) = classify_reaction(&[], "1.1.1.1", false, false, &opts(&e));
        assert_eq!(s, HitStatus::NoSeqData);
    }

    #[test]
    fn spontaneous_when_marked_and_no_data() {
        let e = HashSet::new();
        let (s, _, _) = classify_reaction(&[], "", false, true, &opts(&e));
        assert_eq!(s, HitStatus::Spontaneous);
    }
}
