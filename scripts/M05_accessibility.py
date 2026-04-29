#!/usr/bin/env python3
"""
M05_accessibility.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M05 · mRNA Accessibility (RNAplfold)

PURPOSE
    Compute position-wise mRNA accessibility profiles for each
    target gene alignment using RNAplfold (ViennaRNA package).

    RNAplfold computes the probability that each position is
    single-stranded (unpaired) — essential for Cas13a crRNA design,
    since LwCas13a requires accessible (single-stranded) target sites.

    For each target, a consensus sequence is derived from the MSA
    and used as input to RNAplfold. The output is a per-position
    accessibility profile (.lunp file).

    High accessibility (>= min_acc threshold) marks candidate
    regions for crRNA design in M06.

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M05_accessibility.py \
        --config ~/sherlock/config.yaml

OUTPUT
    main/data/05_accessibility/{target}_accessibility.tsv
    main/data/05_accessibility/{target}_consensus.fasta
    main/logs/M05_accessibility.log
    main/reports/M05_accessibility/summary.tsv
    main/reports/M05_accessibility/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG = "~/sherlock/config.yaml"
THREADS        = 8     # change to 16 for maximum

# RNAplfold parameters
WINDOW_SIZE    = 80    # nucleotide window (W parameter)
SPAN_SIZE      = 40    # maximum base-pair span (L parameter)
MIN_ACC        = 0.5   # minimum unpaired probability for crRNA design

# Consensus sequence: minimum fraction of non-gap bases to include position
MIN_BASE_FREQ  = 0.5

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from pipeline_utils import (
    load_config,
    get_logger,
    print_logo,
    save_versions,
    write_tsv,
    write_checkpoint,
)

# =============================================================================
# CONSENSUS SEQUENCE
# =============================================================================

def derive_consensus(aln_path: Path, min_base_freq: float = MIN_BASE_FREQ) -> str:
    """
    Derive a consensus sequence from a MAFFT alignment.

    For each position:
      - Take the most frequent non-gap base
      - If gap frequency > (1 - min_base_freq), skip position
      - Convert T → U (RNA sequence for RNAplfold)

    Returns consensus RNA sequence (U instead of T).
    """
    sequences = []
    current   = []

    with open(aln_path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
            else:
                current.append(line)
        if current:
            sequences.append("".join(current).upper())

    if not sequences:
        return ""

    length    = len(sequences[0])
    consensus = []

    for i in range(length):
        col   = [seq[i] for seq in sequences if i < len(seq)]
        bases = [b for b in col if b not in ("-", "N", "?")]

        if not bases:
            continue
        gap_freq = 1 - len(bases) / len(col)
        if gap_freq > (1 - min_base_freq):
            continue   # skip mostly-gap positions

        most_common = Counter(bases).most_common(1)[0][0]
        consensus.append(most_common)

    # Convert to RNA
    rna = "".join(consensus).replace("T", "U")
    return rna


# =============================================================================
# RNAplfold
# =============================================================================

def run_rnaplfold(sequence: str,
                  target_name: str,
                  out_dir: Path,
                  window: int,
                  span: int,
                  log) -> Path | None:
    """
    Run RNAplfold on a sequence.
    Writes {target_name}_lunp file in out_dir.
    Returns path to .lunp file or None on failure.
    """
    if len(sequence) < window:
        log.warning(f"  {target_name}: sequence too short "
                    f"({len(sequence)} < {window}) — using full length")
        window = len(sequence) // 2
        span   = window // 2

    lunp_file = out_dir / f"{target_name}_lunp"

    # RNAplfold writes output to current directory
    # Use temp dir and move output
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write input fasta
        in_fasta = Path(tmpdir) / "input.fa"
        in_fasta.write_text(f">{target_name}\n{sequence}\n")

        cmd = [
            "RNAplfold",
            "-W", str(window),
            "-L", str(span),
            "-u", "28",   # R3: window = spacer length (Lorenz 2011)           # compute unpaired probabilities
            "--noLP",            # no lonely base pairs
        ]

        try:
            with open(in_fasta) as stdin_f:
                r = subprocess.run(
                    cmd,
                    stdin   = stdin_f,
                    capture_output = True,
                    text    = True,
                    cwd     = tmpdir,
                    timeout = 300,
                )

            if r.returncode != 0:
                log.error(f"  RNAplfold failed: {r.stderr[:200]}")
                return None

            # Find the _lunp file in tmpdir
            lunp_tmp = Path(tmpdir) / f"{target_name}_lunp"
            if not lunp_tmp.exists():
                # Try alternative name
                lunp_files = list(Path(tmpdir).glob("*_lunp"))
                if lunp_files:
                    lunp_tmp = lunp_files[0]
                else:
                    log.error(f"  No _lunp file found for {target_name}")
                    return None

            # Copy to output dir
            import shutil
            shutil.copy(lunp_tmp, lunp_file)
            return lunp_file

        except subprocess.TimeoutExpired:
            log.error(f"  RNAplfold timeout for {target_name}")
            return None
        except Exception as e:
            log.error(f"  RNAplfold error: {e}")
            return None


# =============================================================================
# PARSE LUNP FILE
# =============================================================================

def parse_lunp(lunp_path: Path) -> pd.DataFrame:
    """
    Parse RNAplfold _lunp output file.
    Returns DataFrame with columns: position, unpaired_prob

    _lunp format:
      #i    l=1   l=2   ...
       1    0.95  0.82  ...
    Column l=1 is the probability of position i being unpaired
    in a window of length 1 (single-nucleotide accessibility).
    """
    rows = []
    with open(lunp_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                pos  = int(parts[0])
                prob = float(parts[28])   # R3: read l=28 column   # l=1 column
                rows.append({"position": pos, "unpaired_prob": prob})
            except ValueError:
                continue

    if not rows:
        return pd.DataFrame(columns=["position", "unpaired_prob"])

    return pd.DataFrame(rows)


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config", default=DEFAULT_CONFIG,
                   help="Path to config.yaml")
    return p.parse_args()


def main():
    args   = parse_args()
    cfg    = load_config(str(Path(args.config).expanduser()))
    org    = cfg["organism"]["display"]
    paths  = cfg["paths"]
    acc_cfg= cfg.get("accessibility", {})
    window = acc_cfg.get("window",  WINDOW_SIZE)
    span   = acc_cfg.get("span",    SPAN_SIZE)
    min_acc= acc_cfg.get("min_acc", MIN_ACC)

    # Dirs
    main_dir   = Path(paths["main"])
    msa_dir    = main_dir / "data" / "04_alignment" / "msa"
    out_dir    = main_dir / "data" / "05_accessibility"
    report_dir = main_dir / "reports" / "M05_accessibility"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M05_accessibility", cfg)
    print_logo("M05 · mRNA Accessibility", organism=org)

    save_versions(
        tools       = {"RNAplfold": "RNAplfold"},
        python_pkgs = ["numpy", "pandas"],
        report_dir  = report_dir,
        log         = log,
    )

    # Find all alignments from M04
    aln_files = sorted(msa_dir.glob("*.aln"))
    if not aln_files:
        log.error("No alignment files found. Run M04 first.")
        sys.exit(1)

    log.info(f"Found {len(aln_files)} alignments")
    log.info(f"RNAplfold: W={window}  L={span}  min_acc={min_acc}")

    summary_rows = []

    for aln in tqdm(aln_files, desc="RNAplfold", unit="target"):
        target = aln.stem
        log.info("─" * 56)
        log.info(f"Target: {target}")

        # Derive consensus sequence from alignment
        consensus = derive_consensus(aln)
        if not consensus:
            log.warning(f"  {target}: empty consensus — skipping")
            continue

        log.info(f"  Consensus length: {len(consensus)} nt")

        # Save consensus FASTA
        cons_fasta = out_dir / f"{target}_consensus.fasta"
        cons_fasta.write_text(
            f">{target}_consensus\n{consensus}\n", encoding="utf-8")

        # Run RNAplfold
        lunp_path = run_rnaplfold(
            sequence    = consensus,
            target_name = target,
            out_dir     = out_dir,
            window      = window,
            span        = span,
            log         = log,
        )

        if lunp_path is None:
            summary_rows.append({
                "target":          target,
                "consensus_len":   len(consensus),
                "n_accessible":    0,
                "pct_accessible":  0.0,
                "mean_acc":        0.0,
                "status":          "FAILED",
            })
            continue

        # Parse and save accessibility profile
        df = parse_lunp(lunp_path)
        if df.empty:
            log.warning(f"  {target}: empty lunp file")
            continue

        # Add accessible flag
        df["accessible"] = df["unpaired_prob"] >= min_acc
        df["target"]     = target

        # Save as TSV
        acc_tsv = out_dir / f"{target}_accessibility.tsv"
        df.to_csv(acc_tsv, sep="\t", index=False)

        n_acc   = df["accessible"].sum()
        pct_acc = 100 * n_acc / len(df)
        mean_acc= df["unpaired_prob"].mean()

        log.info(f"  Accessible positions: {n_acc}/{len(df)} "
                 f"({pct_acc:.1f}%)  mean={mean_acc:.3f}")

        summary_rows.append({
            "target":         target,
            "consensus_len":  len(consensus),
            "n_positions":    len(df),
            "n_accessible":   int(n_acc),
            "pct_accessible": round(pct_acc, 1),
            "mean_acc":       round(float(mean_acc), 3),
            "status":         "OK",
        })

    # Write summary
    write_tsv(summary_rows, report_dir / "summary.tsv", log=log)

    log.info("─" * 56)
    log.info("M05 complete.")
    log.info(f"  Profiles → {out_dir}")
    for r in summary_rows:
        if r["status"] == "OK":
            log.info(f"  {r['target']:<20} "
                     f"{r['pct_accessible']}% accessible  "
                     f"mean={r['mean_acc']}")

    write_checkpoint("M05_accessibility", report_dir, log=log)


if __name__ == "__main__":
    main()
