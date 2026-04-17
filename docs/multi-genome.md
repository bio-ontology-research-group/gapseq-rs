# Multi-genome & metagenome workflows

Gapsmith's single-genome pipeline (`doall`) stays the default. This chapter
covers three **additional, opt-in** capabilities for running at scale:

1. **[`--gspa-run`](#1-consuming-a-gspa-manifest-gspa-run)** — consume
   precomputed protein clusters + alignments from the upstream
   [gspa](https://github.com/bio-ontology-research-group/gspa) pipeline
   instead of re-running DIAMOND/BLASTp per genome.
2. **[`doall-batch`](#2-batch-reconstruction-doall-batch)** — rayon- and
   SLURM-friendly driver that fans `doall` across N genomes (targets
   1 000 metagenomes; scales to ~1 M genomes via `--shard`).
3. **[`community`](#3-community-optimisation-community)** — community-level
   objectives for metagenomes: a cheap per-MAG mode with a shared union
   medium, and a full cFBA mode that composes N drafts into one model.

Nothing in here breaks single-genome usage. If you never pass `--gspa-run`
or use the new subcommands, behaviour is identical to earlier releases.

---

## 1. Consuming a gspa manifest (`--gspa-run`)

### Why

Annotating N genomes with a single mmseqs2/linclust cluster + one
alignment against cluster representatives is dramatically cheaper than
N independent BLASTp runs. gspa's Phase 11 already produces that
cluster; `--gspa-run` plugs the output into `gapsmith find` /
`find-transport` / `doall`.

### Manifest layout (v1)

gspa writes (or you assemble by hand) a directory with this shape:

```
<gspa-run>/
  clusters.tsv            # 3 cols: rep_id<TAB>member_id<TAB>genome_id
  alignment/<ref>.tsv     # 8-col gapsmith TSV keyed on rep_id
  genomes.tsv             # 2–3 cols: genome_id<TAB>faa_path[<TAB>abundance]
  gspa-manifest.toml      # optional — version/provenance marker
```

- `clusters.tsv`: every (rep, member, genome) triple emitted by the
  cross-genome mmseqs2 cluster. The rep appears on its own line as its
  own member.
- `alignment/*.tsv`: one or more 8-column files in gapsmith's canonical
  format (`qseqid pident evalue bitscore qcov stitle sstart send`). All
  `.tsv` / `.tab` / `.txt` files in that directory are concatenated
  transparently — useful when you run DIAMOND and mmseqs2 separately and
  want both sources of evidence.
- `genomes.tsv`: the authoritative list of organism ids + FASTA paths.
  The optional 3rd column carries a relative abundance (any positive
  float — RPKM / TPM / raw counts / …) that propagates into
  `community cfba`.

### Usage

```bash
# Single genome, alignment reused from an already-built gspa manifest:
gapsmith doall /path/to/genome.faa -f out/ \
    --gspa-run /path/to/gspa-run/ \
    --gspa-genome-id genome

# If the FASTA stem already matches the genome_id in genomes.tsv,
# --gspa-genome-id can be omitted.
```

`--gspa-coverage-fraction` re-scales mmseqs2-native coverage columns
(0–1) up to the gapsmith convention (0–100). Diamond/BLASTp TSVs are
already 0–100 and need no flag.

### What happens under the hood

`gapsmith_align::GspaRunAligner` implements the `Aligner` trait. On
`align()` it:

1. Reads `alignment/*.tsv` lazily (once, cached).
2. Looks up each rep hit in `clusters.tsv` filtered to
   `target_genome_id`.
3. Emits one `Hit` per member, rewriting `qseqid` to the member's
   original protein id so downstream find / find-transport see a normal
   per-genome hit table.

---

## 2. Batch reconstruction (`doall-batch`)

### When to reach for it

| Scale                   | Recommended tool                          |
|-------------------------|-------------------------------------------|
| 1 genome                | `gapsmith doall`                          |
| 10 – 200 genomes, 1 box | `gapsmith doall-batch --jobs <cores>`     |
| 1 k – 1 M genomes, cluster | `doall-batch --shard i/N` inside a SLURM array |

### Flags that matter

- `--genomes-dir <dir>` / `--genomes-list <tsv>` / `--gspa-run <dir>` —
  exactly one; the last option pulls the genome list straight out of
  the gspa manifest.
- `--jobs N` — rayon pool size. Each worker processes one genome at a
  time; we pin `doall`'s own `--threads` to 1 so rayon gets clean
  parallelism.
- `--shard i/N` — SLURM-array-friendly selector. After sorting the
  genome list alphabetically, only genomes whose 0-based index `k`
  satisfies `k mod N == i` are processed. Every shard is independent —
  no shared state, no cross-task coordination.
- `--continue-on-error` — log-and-skip instead of aborting the batch
  when a single genome's `doall` fails. Failures are written to
  `<out-dir>/doall-batch-errors.tsv`.

### Recipes

#### 100 genomes on one box, reusing a gspa cluster

```bash
gapsmith doall-batch \
    --gspa-run /data/gspa-run/ \
    -f /out/reconstructions/ \
    --jobs 32 \
    --continue-on-error
```

Each `<genome_id>/` subdirectory under `/out/reconstructions/` gets the
standard `doall` outputs.

#### 1 M genomes on a SLURM cluster

```bash
#!/bin/bash
#SBATCH --array=0-1023
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

gapsmith doall-batch \
    --gspa-run /scratch/gspa-run/ \
    -f /scratch/recon/ \
    --jobs $SLURM_CPUS_PER_TASK \
    --shard ${SLURM_ARRAY_TASK_ID}/1024 \
    --continue-on-error
```

- One gspa run on the full 1 M genome set covers the expensive
  alignment step globally.
- 1 024 SLURM tasks each process ~1 000 genomes in parallel; add tasks
  or add cores per task to change the balance.

---

## 3. Community optimisation (`community`)

Two subcommands; pick the one that matches your accuracy / scale
trade-off. You can run both on the same input and compare.

### `community per-mag` — scalable default

Each MAG keeps its own GEM and its own FBA. We (a) build a *shared*
medium that is the union of every per-MAG medium (largest `maxFlux` per
compound wins), (b) run FBA on every draft under that shared medium,
and (c) report the abundance-weighted mean growth rate.

```bash
gapsmith community per-mag \
    --gspa-run /data/gspa-run/ \
    --drafts-root /out/reconstructions/ \
    -o /out/community/
```

Outputs:

- `community-medium.csv` — the shared medium.
- `per-mag-growth.tsv` — `genome_id, abundance, growth_rate, status` per
  MAG.

Complexity: linear in `N`; no community-level LP. This is the mode for
≥ 100 MAGs where a full cFBA isn't tractable.

### `community cfba` — full community LP

Composes N drafts into one block-diagonal model. Private mets/rxns are
suffixed `__<genome_id>`; extracellular mets (`_e0`) and exchange
reactions (`EX_*_e0`) are **shared** so organisms cross-feed through a
single pool. Bounds on shared exchanges are widened to the most
permissive input.

A synthetic `bio_community` reaction is added whose flux equals
`Σ w_k · bio_k` (weights come from `--abundance` or `genomes.tsv`).

```bash
gapsmith community cfba \
    --gspa-run /data/gspa-run/ \
    --drafts-root /out/reconstructions/ \
    --abundance /data/abundances.tsv \
    --balanced-growth \
    -m /data/community-medium.csv \
    -o /out/community/
```

Outputs:

- `community.gmod.cbor` — composed community model.
- `community-fba.tsv` — per-organism biomass flux plus the
  community-level objective.

`--balanced-growth` adds one linear constraint per organism forcing
`v(bio_k) == v(bio_community)`, matching classical cFBA semantics. Drop
the flag for an abundance-weighted average without the equal-growth
constraint (faster, more realistic when abundances are very uneven).

Scaling: the community LP grows roughly `N × (rxns + mets)`. Expect
comfortable behaviour up to ~50 organisms; larger communities are
better served by `per-mag` first, with cFBA applied to selected
sub-communities.

---

## 4. Decision tree

```
                 "metagenome / multi-MAG run"
                             │
            ┌────────────────┴────────────────┐
            │                                 │
     small curated set                 thousands of MAGs
      (≤ ~50 organisms)                  or fast iteration
            │                                 │
    ┌───────┴───────┐                 ┌───────┴───────┐
    │ need obligate │                 │ need obligate │
    │ cross-feeding?│                 │ cross-feeding?│
    └───────┬───────┘                 └───────┬───────┘
      yes  │    │ no                    yes  │    │ no
           │    │                            │    │
  community cfba  community per-mag  community cfba on  community per-mag
  --balanced-growth                  selected cohorts
```

Rules of thumb:

- 1 — 50 organisms, synthetic or curated community: **`cfba`**.
- 50 — 1 000 organisms, real-world metagenome: **`per-mag`**.
- >1 000 organisms or iterating over many assemblies: always start with
  `doall-batch` + `per-mag`; drop to `cfba` only for sub-cohorts.

---

## 5. Producing a gspa run from outside gspa

Any tool that can emit the three files above works. A minimal
hand-rolled example:

```bash
# 1. Concatenate proteomes, tagging each header with its genome id.
mkdir -p run/
python tools/concat_with_prefix.py --in genomes/ --out run/all.faa

# 2. Cluster across genomes.
mmseqs easy-cluster run/all.faa run/cluster run/tmp \
    --min-seq-id 0.5 -c 0.8 --threads 32

# 3. Pivot the 2-col `run/cluster_cluster.tsv` into 3-col clusters.tsv
#    by peeling each member's `<genome_id>|` prefix.
awk -F$'\t' '{
  member=$2; rep=$1;
  n = split(member, a, "|"); gid=a[1];
  orig=substr(member, length(gid)+2);
  print rep "\t" orig "\t" gid
}' run/cluster_cluster.tsv > run/clusters.tsv

# 4. Run alignment ONCE against cluster reps.
mkdir -p run/alignment
diamond makedb --in reference.faa -d ref
diamond blastp -q run/cluster_rep_seq.fasta -d ref \
    --outfmt 6 qseqid pident evalue bitscore qcovhsp stitle sstart send \
    --more-sensitive -k 10 -e 1e-5 > run/alignment/reference.tsv

# 5. genomes.tsv (optional 3rd column = abundance).
for f in genomes/*.faa; do
  id=$(basename "$f" .faa)
  echo -e "$id\t$(realpath "$f")"
done > run/genomes.tsv
```

After that, `gapsmith doall-batch --gspa-run run/ --jobs 32 -f out/` is
all you need.
