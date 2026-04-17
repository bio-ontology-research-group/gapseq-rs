#!/usr/bin/env python3
"""Aggregate per-run /usr/bin/time output into a markdown comparison table.

Usage:
  aggregate.py <results_dir>

Reads <results_dir>/logs/*.time files whose stems follow
  <genome>_<aligner>_<stage>_<tool>
(e.g. GCF_000005845.2_..._blast_find_rs / _find_r). Produces a markdown
table comparing rs vs r wall-time + peak-RSS per (genome, stage,
aligner), plus a speedup column.
"""

from __future__ import annotations
import glob
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_time_file(path: Path) -> tuple[float, int, int]:
    """Returns (wall_sec, rss_kb, exit_status)."""
    wall = 0.0
    rss = 0
    exit_status = 0
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith("Elapsed (wall clock) time"):
                raw = line.split(":", 1)[1].strip()
                # "H:MM:SS.ss" or "M:SS.ss"
                parts = raw.split(":")
                # Last token is the clock spec — drop the "(h:mm:ss or m:ss)" annotation.
                # Grab the final timestamp part.
                # /usr/bin/time always prints "Elapsed (wall clock) time (h:mm:ss or m:ss): H:MM:SS.ss"
                # so the time is just the last whitespace-separated token.
                ts = line.split()[-1]
                tparts = ts.split(":")
                t = 0.0
                for p in tparts:
                    t = t * 60 + float(p)
                wall = t
            elif line.startswith("Maximum resident set size"):
                rss = int(line.split()[-1])
            elif line.startswith("Exit status"):
                exit_status = int(line.split()[-1])
    return wall, rss, exit_status


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/mnt/data1/leechuck/gapseq-bench/results")
    logs = root / "logs"
    # Key: (genome, stage, aligner). Value: dict tool -> (wall, rss, exit)
    table: dict = defaultdict(dict)
    for tf in logs.glob("*.time"):
        stem = tf.stem  # <genome>_<aligner>_<stage>_<tool>
        # Parse right-to-left: last token is tool (rs/r), second-to-last is stage, etc.
        parts = stem.rsplit("_", 3)
        if len(parts) != 4:
            continue
        genome, aligner, stage, tool = parts
        try:
            wall, rss, exit_status = parse_time_file(tf)
        except Exception as e:
            print(f"skip {tf}: {e}", file=sys.stderr)
            continue
        table[(genome, stage, aligner)][tool] = (wall, rss, exit_status)

    # Order: genomes in insertion order (sort alphabetically); stages in pipeline order.
    stage_order = ["find", "ft", "draft", "medium", "fill"]
    stage_label = {
        "find": "find -p all",
        "ft": "find-transport",
        "draft": "draft",
        "medium": "medium",
        "fill": "fill",
    }

    # Group rows by genome then stage.
    genomes = sorted({k[0] for k in table})

    # Per-genome size: count proteins from the genomes dir if available.
    # Accept either `<root>/../genomes/*.faa.gz` (when results/ is a sibling
    # of genomes/) or `<root>/genomes/*.faa.gz` (when they're nested).
    geno_sizes = {}
    for cand in (root.parent / "genomes", root / "genomes"):
        faas = glob.glob(str(cand / "*.faa.gz"))
        if not faas:
            continue
        import gzip
        for fa in faas:
            n = 0
            with gzip.open(fa, "rt") as f:
                for line in f:
                    if line.startswith(">"):
                        n += 1
            geno_sizes[Path(fa).name] = n
        break

    print()
    print("| Genome | Proteins | Stage | upstream R (s) | gapsmith (s) | Speedup | Rs RSS (MB) | R RSS (MB) |")
    print("|---|---:|---|---:|---:|---:|---:|---:|")
    organism_names = {
        "GCF_000005845.2_ASM584v2_protein": "*E. coli* K-12",
        "GCF_000006945.2_ASM694v2_protein": "*S. Typhimurium* LT2",
        "GCF_000009045.1_ASM904v1_protein": "*B. subtilis* 168",
        "GCF_000007725.1_ASM772v1_protein": "*B. floridanus*",
    }
    # Normalise geno_sizes keys to match the filename stem used in time-log filenames.
    size_by_stem = {
        k.replace(".faa.gz", "").replace(".faa", ""): v for k, v in geno_sizes.items()
    }
    for g in genomes:
        nprot = size_by_stem.get(g, "—")
        short = organism_names.get(g, g[:30])
        if isinstance(nprot, int):
            nprot_str = f"{nprot:,}"
        else:
            nprot_str = str(nprot)
        for stage in stage_order:
            key = None
            for aligner in ("blast", "diamond", "mmseqs2"):
                k = (g, stage, aligner)
                if k in table:
                    key = k
                    break
            if key is None:
                continue
            tools = table[key]
            r_wall = tools.get("r", (None, None, None))[0]
            rs_wall = tools.get("rs", (None, None, None))[0]
            r_rss = tools.get("r", (None, None, None))[1]
            rs_rss = tools.get("rs", (None, None, None))[1]
            speedup = ""
            if r_wall and rs_wall:
                if rs_wall > 0:
                    speedup = f"**{r_wall / rs_wall:.1f}×**"
            rw = f"{r_wall:.1f}" if r_wall else "—"
            rsw = f"{rs_wall:.1f}" if rs_wall else "—"
            rrss = f"{r_rss / 1024:.0f}" if r_rss else "—"
            rsr = f"{rs_rss / 1024:.0f}" if rs_rss else "—"
            if stage in ("draft", "medium", "fill") and not r_wall:
                rw = "(R deps missing)"
                rrss = "—"
            print(f"| {short} | {nprot_str} | {stage_label[stage]} | {rw} | {rsw} | {speedup} | {rsr} | {rrss} |")
    print()


if __name__ == "__main__":
    main()
