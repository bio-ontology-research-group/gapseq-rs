# CLI reference

Every `gapseq` subcommand, exhaustive flag list. Generated from
`--help`, annotated per-section. Run any command with `-h` to see
these locally.

---

## Global options

These apply to every subcommand:

| Flag | Description |
|---|---|
| `--data-dir <PATH>` | Override gapseq reference-data directory. Resolution order when absent: `GAPSMITH_DATA_DIR` env → `$XDG_DATA_HOME/gapsmith` → `<exe-dir>/../dat` → `./dat`. |
| `--seq-dir <PATH>` | Override the sequence database directory. Defaults to `<data-dir>/seq`. |
| `-K, --threads <N>` | Worker threads. Default: all available cores. |
| `-v, -vv, -vvv` | Logging verbosity (info / debug / trace). |
| `-h, --help` | Print help. |

---

## `gapsmith doall`

End-to-end pipeline: find → find-transport → draft → medium → fill.

| Flag | Default | Description |
|---|---|---|
| `<GENOME>` | — | Required. Protein FASTA (`.faa` / `.faa.gz`). |
| `-A, --aligner` | `blastp` | `blastp`, `diamond`, `mmseqs2`, or `precomputed`. |
| `-b, --bitcutoff` | `200.0` | Bitscore cutoff for the find stages. |
| `-c, --coverage` | `75` | Coverage cutoff (0–100 %). |
| `-l, --min-bs-core` | `50.0` | Core-reaction bitscore cutoff. |
| `-t, --taxonomy` | `auto` | `Bacteria`, `Archaea`, or `auto` (currently maps to Bacteria). |
| `-m, --medium` | `auto` | Medium CSV, or `auto` to infer via `gapsmith medium`. |
| `-f, --out-dir` | `.` | Output directory. |
| `-K, --threads` | all cores | Thread count. |
| `--full-suite` | off | Include fill Steps 3 + 4 (slow). |
| `-k, --min-growth` | `0.01` | Minimum growth for gap-filling. |
| `--gspa-run <DIR>` | — | Read precomputed cluster-rep hits from a [gspa](multi-genome.md) manifest instead of running the aligner. |
| `--gspa-genome-id <ID>` | FASTA stem | Row id in the manifest's `genomes.tsv` matching this run. |

Gzipped inputs are auto-decompressed into a `tempfile::tempdir`.

---

## `gapsmith find`

Pathway and reaction detection.

| Flag | Default | Description |
|---|---|---|
| `<GENOME>` | — | Protein FASTA. |
| `-p, --pathways` | `all` | Pathway keyword or id (or `all`). |
| `-l, --pathway-db` | `metacyc,custom` | Comma-separated: `metacyc`, `kegg`, `seed`, `custom`, `all`. |
| `-t, --taxonomy` | `Bacteria` | Reference subfolder under `seq-dir`. |
| `-A, --aligner` | `diamond` | Aligner backend. |
| `-P, --precomputed` | — | Precomputed alignment TSV (skip aligner). |
| `--gspa-run <DIR>` | — | [gspa](multi-genome.md) manifest — expands rep-level hits onto this genome's cluster members. |
| `--gspa-genome-id <ID>` | FASTA stem | Genome id in the manifest matching this invocation. |
| `--gspa-coverage-fraction` | off | Treat the alignment TSV's qcov as 0–1 (mmseqs2 native). |
| `-b, --bitcutoff` | `200.0` | Bitscore cutoff. |
| `-i, --identcutoff` | `0.0` | Identity cutoff (%). |
| `-c, --coverage` | `75` | Coverage cutoff. |
| `-a, --completeness-main` | `80.0` | Pathway completeness cutoff. |
| `-k, --completeness-hints` | `66.0` | Relaxed cutoff when every key reaction is present. |
| `-s, --strict` | off | Disable key-reaction heuristic. |
| `-n, --include-superpathways` | off | Default is to filter them out (matches upstream). |
| `-o, --out-dir` | `.` | Output directory. |
| `-u, --suffix` | — | Filename suffix `<stem>-<suffix>-{Reactions,Pathways}.tbl`. |

---

## `gapsmith find-transport`

Transporter detection.

| Flag | Default | Description |
|---|---|---|
| `<GENOME>` | — | Protein FASTA. |
| `-A, --aligner` | `blastp` | Aligner backend. |
| `-P, --precomputed` | — | Precomputed TSV. |
| `-b, --bitcutoff` | `50.0` | Transport hits use a lower threshold than biosynthetic. |
| `-i, --identcutoff` | `0.0` | Identity cutoff. |
| `-c, --coverage` | `75` | Coverage cutoff. |
| `-a, --nouse-alternatives` | off | Disable the alt-transporter fallback. |
| `-m, --only-met` | — | Restrict to one substrate keyword. |
| `-o, --out-dir` | `.` | Output directory. |

---

## `gapsmith draft`

Build a draft metabolic model from `find` + `find-transport` output.

| Flag | Default | Description |
|---|---|---|
| `-r, --reactions` | required | `*-Reactions.tbl`. |
| `-t, --transporter` | required | `*-Transporter.tbl`. |
| `-b, --biomass` | `auto` | `auto` / `pos` / `neg` / `archaea` / `<path.json>`. |
| `-n, --name` | inferred | Model id. |
| `-u, --high-evi-rxn-bs` | `200.0` | Bitscore threshold for core reactions. |
| `-l, --min-bs-for-core` | `50.0` | Lower bound for candidate pool. |
| `-o, --out-dir` | `.` | Output directory. |
| `--no-sbml` | off | Skip SBML emission. |

---

## `gapsmith medium`

Rule-based medium inference.

| Flag | Default | Description |
|---|---|---|
| `-m, --model` | required | Draft CBOR / JSON. |
| `-p, --pathways` | required | `*-Pathways.tbl` from `find`. |
| `-c, --manual-flux` | — | `cpdXXXXX:val;cpdYYYYY:val` overrides. |
| `-o, --output` | inferred | Output CSV path. |
| `-f, --out-dir` | `.` | Output directory (when `-o` is absent). |

---

## `gapsmith fill`

Iterative gap-filling.

| Flag | Default | Description |
|---|---|---|
| `<MODEL>` | required | Draft CBOR / JSON. |
| `-n, --media` | required | Medium CSV. |
| `-r, --reactions` | — | `*-Reactions.tbl` for bitscore-weighted pFBA. |
| `-t, --target` | `cpd11416` | Biomass pseudo-metabolite id. |
| `-k, --min-growth` | `0.01` | Minimum growth rate floor. |
| `-b, --bcore` | `50.0` | Minimum bitscore for "core" classification. |
| `--high-evi` | `200.0` | Bitscore-to-weight upper calibration point. |
| `--dummy-weight` | `100.0` | Weight for reactions with no hit. |
| `-o, --out-dir` | `.` | Output directory. |
| `--no-sbml` | off | Skip SBML. |
| `--step1-only` | off | Run only Step 1 (debugging). |
| `--full-suite` | off | Include fill Steps 3 + 4. |
| `--prune-futile` | off | Thermodynamic futile-cycle prune (opt-in; slow). |

---

## `gapsmith fba`

FBA / pFBA on an existing model.

| Flag | Default | Description |
|---|---|---|
| `<MODEL>` | required | CBOR / JSON model. |
| `-r, --objective` | model's own | Objective reaction id. |
| `--pfba` | off | Parsimonious FBA. |
| `--pfba-coef` | `0.001` | pFBA biomass trade-off coefficient. |
| `--min-growth` | `0.0` | Biomass floor (pFBA only). |
| `--top` | `20` | Print the top-N highest-absolute-flux reactions. |
| `--minimise` | off | Minimise instead of maximise. |

---

## `gapsmith adapt`

Edit reactions or force growth on a compound.

| Flag | Default | Description |
|---|---|---|
| `-m, --model` | required | Model file. |
| `-a, --add` | — | Comma-separated `rxnNNNNN` or MetaCyc pathway ids. |
| `-r, --remove` | — | Same format; reactions are removed. |
| `-w, --growth` | — | `cpdNNNNN:TRUE` / `cpdNNNNN:FALSE`. |
| `-b, --reactions` | — | `*-Reactions.tbl` (needed for `-w ...:TRUE` gap-filling). |
| `-k, --min-growth` | `0.01` | Min growth floor for the `-w ...:TRUE` path. |
| `-f, --out-dir` | `.` | Output directory. |
| `--no-sbml` | off | Skip SBML. |

---

## `gapsmith pan`

Build a pan-draft from N drafts.

| Flag | Default | Description |
|---|---|---|
| `-m, --models` | required | Comma-separated / directory / single path. |
| `-t, --min-freq` | `0.06` | Minimum across-draft reaction frequency. |
| `-b, --only-binary` | off | Skip the pan-draft, emit only the rxn × model presence TSV. |
| `-f, --out-dir` | `.` | Output directory. |
| `--no-sbml` | off | Skip SBML. |

---

## `gapsmith update-sequences`

Zenodo seqdb sync.

| Flag | Default | Description |
|---|---|---|
| `-t, --taxonomy` | `Bacteria` | Which seqdb to sync. |
| `-D, --seq-dir` | `<data-dir>/seq` | Seqdb location. |
| `-Z, --record` | pinned | Zenodo record id, or `latest`. |
| `-c, --check-only` | off | Report version, no download. |
| `-q, --quiet` | off | Suppress progress messages. |

---

## `gapsmith convert`

CBOR ↔ JSON round-trip.

| Flag | Default | Description |
|---|---|---|
| `<INPUT>` | required | Path to source file. |
| `<OUTPUT>` | required | Destination path. |
| `--to` | from extension | `cbor` / `json` (overrides extension). |
| `--pretty` | off | Pretty-print JSON. |

---

## `gapsmith export-sbml`

Serialise a CBOR / JSON model as SBML L3V1 + FBC2 + groups.

| Flag | Default | Description |
|---|---|---|
| `<INPUT>` | required | CBOR / JSON model. |
| `<OUTPUT>` | required | Destination `.xml` path. |
| `--compact` | off | Omit pretty-print / nested indentation. |

---

## `gapsmith align`

Run a single aligner standalone. Useful for debugging reference-FASTA
issues or building a precomputed TSV to feed into `find --aligner precomputed`.

| Flag | Default | Description |
|---|---|---|
| `-A, --aligner` | required | `blastp` / `tblastn` / `diamond` / `mmseqs2` / `precomputed`. |
| `-q, --query` | required | Query FASTA. |
| `-t, --target` | required | Target FASTA. |
| `-P, --precomputed` | — | When `-A precomputed`, path to the TSV. |
| `-b, --bitcutoff` | `0.0` | Filter hits below this bitscore. |
| `-c, --coverage` | `0` | Filter below this coverage %. |
| `-e, --evalue` | `1e-5` | E-value cutoff. |
| `--extra` | — | Passthrough flags to the aligner binary. |
| `-o, --out` | stdout | Output TSV path. |

---

## `gapsmith batch-align`

Cluster N genomes + single alignment + per-genome TSV expansion.

| Flag | Default | Description |
|---|---|---|
| `-q, --query` | required | Reference query FASTA. |
| `-g, --genomes` | required | Directory containing `.faa(.gz)` files. |
| `-o, --out-dir` | required | Output `<genome>.tsv` directory. |
| `-A, --aligner` | `diamond` | Aligner backend (`blastp` / `diamond` / `mmseqs2`). |
| `--cluster-identity` | `0.5` | mmseqs2 cluster identity. |
| `--cluster-coverage` | `0.8` | mmseqs2 cluster coverage. |
| `-b, --bitcutoff` | `0.0` | Per-hit bitscore filter. |
| `-c, --coverage` | `0` | Per-hit coverage filter. |

---

## `gapsmith doall-batch`

Run `doall` across many genomes in parallel. See
[multi-genome.md](multi-genome.md) for the full recipe.

| Flag | Default | Description |
|---|---|---|
| `-g, --genomes-dir <DIR>` | — | Directory of protein FASTAs. |
| `--genomes-list <TSV>` | — | Explicit TSV list (`id<TAB>path[<TAB>abundance]`). |
| `--gspa-run <DIR>` | — | Pulls the genome list from a gspa manifest. |
| `-f, --out-dir <DIR>` | required | One `<genome_id>/` subdir is written per input. |
| `-j, --jobs <N>` | all cores | Rayon pool size. |
| `--shard <i/N>` | — | Select genomes where `index mod N == i`. |
| `--continue-on-error` | off | Log failed genomes to `doall-batch-errors.tsv` instead of aborting. |
| *(plus every passthrough flag from `doall`)* | | |

---

## `gapsmith community`

Community-level optimisation. Two subcommands.

### `gapsmith community per-mag`

Per-MAG FBA under a shared (union) medium. Linear in N.

| Flag | Default | Description |
|---|---|---|
| `--drafts-dir <DIR>` | — | Directory of `<id>-{draft,filled}.gmod.cbor`. |
| `--drafts-list <TSV>` | — | `<id><TAB><cbor>[<TAB><medium>]`. |
| `--gspa-run <DIR>` | — | Use manifest's genome list. |
| `--drafts-root <DIR>` | `.` | Where per-genome `doall` outputs live (paired with `--gspa-run`). |
| `-m, --medium <CSV>` | auto | Override the inferred shared medium. |
| `-o, --out-dir <DIR>` | `.` | Writes `community-medium.csv`, `per-mag-growth.tsv`. |

### `gapsmith community cfba`

Compose N drafts into one community model; solve weighted-sum biomass.

| Flag | Default | Description |
|---|---|---|
| `--drafts-dir <DIR>` / `--drafts-list <TSV>` / `--gspa-run <DIR>` | — | Exactly one (same semantics as `per-mag`). |
| `--drafts-root <DIR>` | `.` | Where `doall` outputs live. |
| `--abundance <TSV>` | — | `<id><TAB><abundance>` overrides. |
| `--biomass-rxn <ID>` | `bio1` | Per-organism biomass reaction id. |
| `--balanced-growth` | off | Add `v(bio_k) == v(bio_community)` constraints. |
| `-m, --medium <CSV>` | — | Optional shared medium CSV. |
| `-o, --out-dir <DIR>` | `.` | Writes `community.gmod.cbor`, `community-fba.tsv`. |

---

## `gapsmith db inspect`

Smoke-test the `--data-dir`. Loads every reference table and prints row
counts. No flags beyond `--data-dir`.

---

## `gapsmith test`

Print resolved paths + which external tool binaries are on `PATH`.

---

## `gapsmith example-model`

Emit a small hand-built toy model (3 metabolites / 2 reactions) as CBOR.
Useful for smoke-testing format round-trips and downstream tools.

| Flag | Default | Description |
|---|---|---|
| `<OUTPUT>` | required | Destination CBOR path. |
| `--complex` | off | Emit a slightly bigger "complex" variant. |
