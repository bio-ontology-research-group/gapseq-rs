//! Integration test for `BatchClusterAligner`.
//!
//! Builds three pseudo-genomes whose protein sets overlap, clusters them,
//! and asserts the expansion layer correctly attributes each hit to every
//! genome that carries the corresponding cluster member.

use gapsmith_align::{AlignOpts, BatchClusterAligner, DiamondAligner, GenomeInput, Hit};
use std::fs;
use std::io::Write;
use std::path::PathBuf;

fn which(name: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path) {
        let c = dir.join(name);
        if c.is_file() {
            return Some(c);
        }
    }
    None
}

const PROT_A: &str = "MSKVILVSGAAGYIGSHTCVELLEAGYDVVVLDNLCNSKRSVLPVIEKLGGKHPTFVEGD
IRNEALMTEIFAQHAIDTVIHFAGLKAVGESVAKPLEYYDNNVNGTLRLISAMR";
const PROT_B: &str = "MTIKVAIVGAGPAGLLLAQMLHKAGISHEIYERDSDAKAREYRKLPDRAHNISLVTNNVS
LVGQELFQIGFGYDKAIAETHKLFPEAKCYFNHKIKAVELRGKDTHLCIMKCG";
const PROT_C: &str = "MKQPVKGIVVLNKELKTILRKGDYVIEQEVRGKDLIRMEYSTTEEQAEELFRALDVEQDG
RVMLKIVYDPNGYKIDYSIAQGDDFAKLTLDTHMVTLIMGKEKPETVVTPIEMQPKE";

/// g1 carries A+B; g2 carries B+C; g3 carries A+C.
/// Every protein appears in exactly two genomes, so cluster members span
/// multiple genomes and the expansion must replicate hits accordingly.
fn write_three_genomes(dir: &std::path::Path) -> Vec<GenomeInput> {
    let make = |name: &str, contents: &[(&str, &str)]| {
        let p = dir.join(format!("{name}.faa"));
        let mut f = fs::File::create(&p).unwrap();
        for (id, seq) in contents {
            writeln!(f, ">{id}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
        p
    };
    vec![
        GenomeInput {
            id: "g1".into(),
            fasta: make("g1", &[("sp|A|PROT", PROT_A), ("sp|B|PROT", PROT_B)]),
        },
        GenomeInput {
            id: "g2".into(),
            fasta: make("g2", &[("sp|B|PROT", PROT_B), ("sp|C|PROT", PROT_C)]),
        },
        GenomeInput {
            id: "g3".into(),
            fasta: make("g3", &[("sp|A|PROT", PROT_A), ("sp|C|PROT", PROT_C)]),
        },
    ]
}

#[test]
fn batch_expansion_attributes_hits_to_every_genome_with_member() {
    if which("mmseqs").is_none() || which("diamond").is_none() {
        eprintln!("SKIP: mmseqs + diamond required");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let genomes = write_three_genomes(dir.path());

    // Query = A + B + C. Each should hit two of the three genomes.
    let query = dir.path().join("query.faa");
    let mut f = fs::File::create(&query).unwrap();
    writeln!(f, ">qA\n{PROT_A}\n>qB\n{PROT_B}\n>qC\n{PROT_C}").unwrap();
    drop(f);

    let batcher = BatchClusterAligner {
        inner: Box::new(DiamondAligner::new()),
        cluster_identity: 0.9,
        cluster_coverage: 0.8,
    };
    let opts = AlignOpts { threads: 1, coverage_pct: 50, quiet: true, ..Default::default() };

    let workdir = dir.path().join("work");
    let per_genome = batcher.align_genomes(&query, &genomes, &workdir, &opts).unwrap();

    assert_eq!(per_genome.len(), 3);
    fn count_by_qseqid(hits: &[Hit]) -> Vec<&str> {
        let mut q: Vec<&str> = hits.iter().map(|h| h.qseqid.as_str()).collect();
        q.sort();
        q.dedup();
        q
    }
    // Sort so we can index directly regardless of how BTreeMap expansion
    // came out.
    let by_id: std::collections::HashMap<&str, &gapsmith_align::GenomeHitSet> =
        per_genome.iter().map(|gs| (gs.genome_id.as_str(), gs)).collect();

    let g1_qs = count_by_qseqid(&by_id["g1"].hits);
    let g2_qs = count_by_qseqid(&by_id["g2"].hits);
    let g3_qs = count_by_qseqid(&by_id["g3"].hits);

    // g1 has PROT_A + PROT_B ⇒ should be hit by qA and qB.
    assert!(g1_qs.contains(&"qA"), "g1 missing qA: {g1_qs:?}");
    assert!(g1_qs.contains(&"qB"), "g1 missing qB: {g1_qs:?}");
    assert!(!g1_qs.contains(&"qC"), "g1 unexpectedly has qC: {g1_qs:?}");

    // g2: PROT_B + PROT_C ⇒ qB + qC.
    assert!(g2_qs.contains(&"qB"));
    assert!(g2_qs.contains(&"qC"));
    assert!(!g2_qs.contains(&"qA"));

    // g3: PROT_A + PROT_C ⇒ qA + qC.
    assert!(g3_qs.contains(&"qA"));
    assert!(g3_qs.contains(&"qC"));
    assert!(!g3_qs.contains(&"qB"));

    // Every hit stitle should no longer contain the genome prefix.
    for gs in &per_genome {
        for h in &gs.hits {
            assert!(
                !h.stitle.starts_with(&format!("{}|", gs.genome_id)),
                "stitle still prefixed: {}",
                h.stitle
            );
        }
    }
}
