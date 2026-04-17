#!/usr/bin/env bash
#
# Benchmark driver for gapsmith vs upstream R gapseq.
#
# Expects on PATH (or in CONDA_PREFIX/bin):
#   - gapseq             — the upstream bash dispatcher
#   - gapsmith binary at $GAPSMITH
#   - blastp / diamond / mmseqs
#
# Usage:
#   tools/bench/run_bench.sh <genome_dir> <out_dir>
#
# Produces out_dir/timings.tsv with columns:
#   genome  tool  subcommand  aligner  wall_s  rss_kb  exit
#
# Each run is `/usr/bin/time -v` so we capture wall + peak RSS.

set -u
set -o pipefail

GENOME_DIR=${1:-/mnt/data1/leechuck/gapseq-bench/genomes}
OUT_DIR=${2:-/mnt/data1/leechuck/gapseq-bench/results}
GAPSEQ_R=${GAPSEQ_R:-/mnt/data1/leechuck/gapseq-bench/gapseq/gapseq}
GAPSMITH=${GAPSMITH:-/mnt/data1/leechuck/gapseq-bench/gapsmith/target/release/gapsmith}
DATA_DIR=${DATA_DIR:-/mnt/data1/leechuck/gapseq-bench/gapseq/dat}
# Rust uses blastp/diamond/mmseqs2; R gapseq uses blast/diamond/mmseqs2.
# We normalise to a single input flag and map per-tool below.
ALIGNERS=${ALIGNERS:-blast}

mkdir -p "$OUT_DIR"
RESULTS="$OUT_DIR/timings.tsv"
echo -e "genome\ttool\tsubcommand\taligner\twall_s\trss_kb\texit" > "$RESULTS"

# A single timed run. stdout → .log, stderr → .err. Appends one row to
# $RESULTS. Resilient to the run's own failure — we still record wall/rss.
time_run() {
  local label="$1"  # genome_tool_subcommand_aligner
  local genome="$2" tool="$3" subcmd="$4" aligner="$5"
  shift 5
  local log_dir="$OUT_DIR/logs"
  mkdir -p "$log_dir"
  local time_file="$log_dir/${label}.time"
  local log="$log_dir/${label}.log"
  echo ">>> $label" >&2
  /usr/bin/time -v -o "$time_file" "$@" > "$log" 2>&1
  local exit_code=$?
  local wall_s rss_kb
  # /usr/bin/time emits "Elapsed (wall clock) time (h:mm:ss or m:ss): H:MM:SS.ss" — grab the last token and convert.
  local raw
  raw=$(grep "Elapsed" "$time_file" | awk '{print $NF}')
  wall_s=$(python3 -c "
import sys
t = '$raw'.strip()
if not t: print('0'); sys.exit()
parts = t.split(':')
s = 0.0
for p in parts:
    s = s * 60 + float(p)
print(f'{s:.2f}')
")
  rss_kb=$(grep "Maximum resident set size" "$time_file" | awk '{print $NF}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$(basename "$genome")" "$tool" "$subcmd" "$aligner" "$wall_s" "$rss_kb" "$exit_code" >> "$RESULTS"
}

for genome_path in "$GENOME_DIR"/*.faa.gz "$GENOME_DIR"/*.faa; do
  [ -f "$genome_path" ] || continue
  genome=$(basename "$genome_path")
  stem="${genome%.faa.gz}"; stem="${stem%.faa}"

  for aligner in $ALIGNERS; do
    workdir="$OUT_DIR/workdirs/${stem}_${aligner}"
    rm -rf "$workdir"; mkdir -p "$workdir/rs" "$workdir/r"

    # R uses "blast", Rust uses "blastp".
    rust_aligner="$aligner"
    [ "$aligner" = "blast" ] && rust_aligner="blastp"

    # gapsmith find
    time_run "${stem}_${aligner}_find_rs" "$genome_path" gapsmith find "$aligner" \
      "$GAPSMITH" --data-dir "$DATA_DIR" find -A "$rust_aligner" -p all -b 200 -o "$workdir/rs" "$genome_path"

    # R gapsmith find
    time_run "${stem}_${aligner}_find_r" "$genome_path" gapsmith find "$aligner" \
      "$GAPSEQ_R" find -A "$aligner" -p all -b 200 -f "$workdir/r" "$genome_path"

    # gapsmith find-transport
    time_run "${stem}_${aligner}_ft_rs" "$genome_path" gapsmith find-transport "$aligner" \
      "$GAPSMITH" --data-dir "$DATA_DIR" find-transport -A "$rust_aligner" -b 50 -o "$workdir/rs" "$genome_path"

    # R gapsmith find-transport
    time_run "${stem}_${aligner}_ft_r" "$genome_path" gapsmith find-transport "$aligner" \
      "$GAPSEQ_R" find-transport -A "$aligner" -b 50 -f "$workdir/r" "$genome_path"

    # Assemble input paths for draft+fill stages.
    rxns="$workdir/rs/${stem}-all-Reactions.tbl"
    tr="$workdir/rs/${stem}-Transporter.tbl"
    pwys="$workdir/rs/${stem}-all-Pathways.tbl"

    if [ -s "$rxns" ] && [ -s "$tr" ]; then
      # gapsmith draft
      time_run "${stem}_${aligner}_draft_rs" "$genome_path" gapsmith draft "$aligner" \
        "$GAPSMITH" --data-dir "$DATA_DIR" draft -r "$rxns" -t "$tr" -b neg -o "$workdir/rs"
      # R gapsmith draft — only runs if cobrar is installed.
      if command -v Rscript >/dev/null 2>&1 && Rscript -e 'stopifnot("cobrar" %in% installed.packages()[,1])' >/dev/null 2>&1; then
        time_run "${stem}_${aligner}_draft_r" "$genome_path" gapsmith draft "$aligner" \
          "$GAPSEQ_R" draft -r "$rxns" -t "$tr" -b neg -f "$workdir/r"
      fi
    fi

    # medium + fill: rs-only for now (cobrar gating)
    draft_cbor="$workdir/rs/${stem}-draft.gmod.cbor"
    medium_csv="$workdir/rs/${stem}-medium.csv"
    if [ -s "$draft_cbor" ] && [ -s "$pwys" ]; then
      time_run "${stem}_${aligner}_medium_rs" "$genome_path" gapsmith medium "$aligner" \
        "$GAPSMITH" --data-dir "$DATA_DIR" medium -m "$draft_cbor" -p "$pwys" -o "$medium_csv"
    fi
    if [ -s "$draft_cbor" ] && [ -s "$medium_csv" ]; then
      time_run "${stem}_${aligner}_fill_rs" "$genome_path" gapsmith fill "$aligner" \
        "$GAPSMITH" --data-dir "$DATA_DIR" fill "$draft_cbor" -n "$medium_csv" -r "$rxns" -k 0.01 -o "$workdir/rs" --no-sbml
    fi
  done
done

echo ""
echo "Results written to: $RESULTS"
cat "$RESULTS"
