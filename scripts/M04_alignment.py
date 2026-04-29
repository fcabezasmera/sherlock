#!/usr/bin/env python3
"""
M04_alignment.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M04 · Concatenation + MSA + Conservation

PURPOSE
    1. Concatenate per-group FASTA files into per-target files
       (e.g. tcdA_groupA + tcdA_groupB → tcdA_all.fasta)
    2. Align concatenated files with MAFFT
    3. Compute per-position conservation metrics:
         - Shannon entropy  (variability; 0=conserved)
         - Trident index    (physicochemical conservation)
         - Wilcoxon test    (window vs background significance)

    NOTE: TrimAl is NOT used. Removing alignment columns would
    eliminate potential crRNA target sites. Gap filtering is
    handled position-wise in M06 (CaSilico crRNA design).

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M04_alignment.py --config ~/sherlock/config.yaml

    --skip-mafft   reuse existing alignments

OUTPUT
    main/data/04_alignment/msa/{target}.aln
    main/data/04_alignment/conservation/{target}_conservation.tsv
    main/logs/M04_alignment.log
    main/reports/M04_alignment/summary.tsv
    main/reports/M04_alignment/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG = "~/sherlock/config.yaml"
THREADS        = 8     # change to 16 for maximum
MAFFT_MODE     = "auto"
CONSERVATION_WINDOW = 28   # sliding window = crRNA spacer length
MIN_SEQS            = 5    # minimum sequences to attempt alignment

# Concatenation map: {output_name: [input_stems_from_M03]}
# Each output → one MAFFT alignment → one crRNA target
CONCAT_MAP = {
    "tcdA_all":       ["tcdA_groupA",   "tcdA_groupB"],
    "tcdB_clade2":    ["tcdB_groupA"],                    # RT027-like (hypervirulent)
    "tcdB_clade1":    ["tcdB_groupB"],                    # RT012-like (classic)
    "tcdC_wt":        ["tcdC_groupB"],
    "tcdC_junction":  ["tcdC_junction_groupA"],
    "cdtA_groupA":    ["cdtA_groupA"],
    "cdtB_groupA":    ["cdtB_groupA"],
    "tpiA_all":       ["tpiA_groupA", "tpiA_groupB", "tpiA_groupC"],
    "rpoB_all":       ["rpoB_groupA", "rpoB_groupB", "rpoB_groupC"],
}

# Nucleotide physicochemical groups for Trident index
NT_GROUPS = {
    "A": "purine", "G": "purine",
    "C": "pyrimidine", "T": "pyrimidine", "U": "pyrimidine",
}

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from pipeline_utils import (
    load_config,
    get_logger,
    print_logo,
    save_versions,
    write_tsv,
    write_checkpoint,
    count_seqs,
)

# =============================================================================
# STEP 1 — Concatenation
# =============================================================================

def concatenate_fastas(input_stems: list,
                       extract_dir: Path,
                       out_fasta: Path,
                       log) -> int:
    """
    Concatenate multiple FASTA files into one.
    Returns total number of sequences in output.
    """
    records = []
    for stem in input_stems:
        fasta = extract_dir / f"{stem}.fasta"
        if not fasta.exists():
            log.warning(f"  Not found: {fasta.name} — skipping")
            continue
        with open(fasta, encoding="utf-8") as f:
            records.append(f.read())

    if not records:
        log.warning(f"  No input files found for {out_fasta.name}")
        return 0

    out_fasta.write_text("\n".join(records), encoding="utf-8")
    n = count_seqs(out_fasta)
    return n


# =============================================================================
# STEP 2 — MAFFT alignment
# =============================================================================

def run_mafft(input_fasta: Path, output_aln: Path,
              mode: str, threads: int, log) -> bool:
    """Run MAFFT. Returns True if successful."""
    n = count_seqs(input_fasta)
    if n < MIN_SEQS:
        log.warning(f"  {input_fasta.name}: {n} sequences "
                    f"(min={MIN_SEQS}) — skipping")
        return False

    log.info(f"  {input_fasta.name}  ({n} sequences)")

    # Select strategy based on number of sequences
    if mode == "auto":
        if n <= 200:
            strategy = ["--localpair", "--maxiterate", "1000"]
        elif n <= 500:
            strategy = ["--globalpair", "--maxiterate", "100"]
        else:
            strategy = ["--auto"]
    else:
        strategy = [f"--{mode}"]

    cmd = (["mafft"] + strategy +
           ["--thread", str(threads), "--quiet", str(input_fasta)])

    try:
        r = subprocess.run(cmd, capture_output=True, text=True, check=True)
        output_aln.write_text(r.stdout, encoding="utf-8")
        log.info(f"    → {output_aln.name}")
        return True
    except subprocess.CalledProcessError as e:
        log.error(f"  MAFFT failed: {e.stderr[:200]}")
        return False


# =============================================================================
# STEP 3 — Conservation metrics
# =============================================================================

def parse_alignment(aln_path: Path) -> tuple:
    """Parse FASTA alignment. Returns (sequences, length)."""
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
    length = len(sequences[0]) if sequences else 0
    return sequences, length


def shannon_entropy(column: list) -> float:
    """Shannon entropy for one alignment column. 0=conserved."""
    bases = [b for b in column if b not in ("-", "N", "?")]
    if not bases:
        return float("nan")
    counts = Counter(bases)
    total  = sum(counts.values())
    h = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            h -= p * math.log2(p)
    return h


def trident_index(column: list) -> float:
    """Physicochemical conservation (purine/pyrimidine). 1=conserved."""
    bases = [b for b in column if b not in ("-", "N", "?")]
    if not bases:
        return float("nan")
    groups = [NT_GROUPS.get(b, "other") for b in bases]
    majority = max(Counter(groups).values()) / len(groups)
    return majority


def compute_conservation(sequences: list, length: int,
                          window: int) -> pd.DataFrame:
    """
    Per-position conservation metrics + sliding window Wilcoxon test.
    """
    shannons  = []
    tridents  = []
    gap_fracs = []

    for i in range(length):
        col = [seq[i] for seq in sequences if i < len(seq)]
        shannons.append(shannon_entropy(col))
        tridents.append(trident_index(col))
        n_gaps = sum(1 for b in col if b == "-")
        gap_fracs.append(n_gaps / len(col) if col else float("nan"))

    sh_arr     = np.array(shannons, dtype=float)
    background = sh_arr[~np.isnan(sh_arr)]

    win_means = []
    wilcox_p  = []
    wilcox_sig= []

    for i in range(length):
        s = max(0, i - window // 2)
        e = min(length, i + window // 2)
        wv = sh_arr[s:e]
        wv = wv[~np.isnan(wv)]

        if len(wv) < 3:
            win_means.append(float("nan"))
            wilcox_p.append(float("nan"))
            wilcox_sig.append(False)
            continue

        win_means.append(float(np.nanmean(wv)))
        try:
            _, p = stats.mannwhitneyu(wv, background, alternative="less")
            wilcox_p.append(float(p))
            wilcox_sig.append(False)  # placeholder, corrected below
        except Exception:
            wilcox_p.append(float("nan"))
            wilcox_sig.append(False)

    # R2: FDR correction (Benjamini & Hochberg 1995)
    pvals_clean = [p if p == p else 1.0 for p in wilcox_p]  # replace NaN
    if any(p < 1.0 for p in pvals_clean):
        _, pvals_corr, _, _ = multipletests(pvals_clean, method="fdr_bh")
        wilcox_sig = [bool(p < 0.05) for p in pvals_corr]
    return pd.DataFrame({
        "position":             range(1, length + 1),
        "shannon":              shannons,
        "trident":              tridents,
        "gap_fraction":         gap_fracs,
        "window_shannon_mean":  win_means,
        "wilcoxon_pvalue":      wilcox_p,
        "wilcoxon_significant": wilcox_sig,
    })


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
    p.add_argument("--skip-mafft", action="store_true",
                   help="Reuse existing alignments")
    return p.parse_args()


def main():
    args    = parse_args()
    cfg     = load_config(str(Path(args.config).expanduser()))
    org     = cfg["organism"]["display"]
    paths   = cfg["paths"]
    aln_cfg = cfg.get("alignment", {})
    mode    = aln_cfg.get("mode", MAFFT_MODE)
    threads = aln_cfg.get("threads",
                           cfg.get("threads", {}).get("default", THREADS))

    # Dirs
    main_dir    = Path(paths["main"])
    extract_dir = main_dir / "data" / "03_extract"
    concat_dir  = main_dir / "data" / "04_alignment" / "concat"
    msa_dir     = main_dir / "data" / "04_alignment" / "msa"
    cons_dir    = main_dir / "data" / "04_alignment" / "conservation"
    report_dir  = main_dir / "reports" / "M04_alignment"

    for d in (concat_dir, msa_dir, cons_dir, report_dir):
        d.mkdir(parents=True, exist_ok=True)

    log = get_logger("M04_alignment", cfg)
    print_logo("M04 · MSA + Conservation", organism=org)

    save_versions(
        tools       = {"mafft": "mafft"},
        python_pkgs = ["numpy", "pandas", "scipy"],
        report_dir  = report_dir,
        log         = log,
    )

    summary_rows = []

    # ------------------------------------------------------------------
    # STEP 1 — Concatenate per-group FASTAs
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 1 · Concatenating FASTA files")
    log.info("─" * 56)

    for target, sources in CONCAT_MAP.items():
        out_fasta = concat_dir / f"{target}.fasta"
        n = concatenate_fastas(sources, extract_dir, out_fasta, log)
        log.info(f"  {target:<20}: {n} sequences  ({' + '.join(sources)})")
        summary_rows.append({
            "target":  target,
            "sources": " + ".join(sources),
            "n_seqs":  n,
            "aligned": False,
            "n_cons_positions": 0,
        })

    # ------------------------------------------------------------------
    # STEP 2 — MAFFT alignment
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 2 · MAFFT alignment")
    log.info("─" * 56)

    for row in tqdm(summary_rows, desc="MAFFT", unit="target"):
        target    = row["target"]
        in_fasta  = concat_dir / f"{target}.fasta"
        out_aln   = msa_dir / f"{target}.aln"

        if not in_fasta.exists() or row["n_seqs"] == 0:
            continue

        if args.skip_mafft and out_aln.exists():
            log.info(f"  {target}: skipping (--skip-mafft)")
            row["aligned"] = True
            continue

        ok = run_mafft(in_fasta, out_aln, mode, threads, log)
        row["aligned"] = ok

    # ------------------------------------------------------------------
    # STEP 3 — Conservation
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 3 · Conservation metrics")
    log.info("─" * 56)

    for row in tqdm(summary_rows, desc="Conservation", unit="target"):
        target  = row["target"]
        aln     = msa_dir / f"{target}.aln"
        cons_out= cons_dir / f"{target}_conservation.tsv"

        if not aln.exists():
            continue

        sequences, length = parse_alignment(aln)
        if not sequences:
            log.warning(f"  {target}: empty alignment")
            continue

        log.info(f"  {target}: {len(sequences)} seqs × {length} positions")
        df = compute_conservation(sequences, length, CONSERVATION_WINDOW)
        df.to_csv(cons_out, sep="\t", index=False)

        n_sig = df["wilcoxon_significant"].sum()
        mean_sh = df["shannon"].mean()
        log.info(f"    mean Shannon={mean_sh:.3f}  "
                 f"significant windows={n_sig}")
        row["n_cons_positions"] = length

    # ------------------------------------------------------------------
    # Reports
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("Writing reports")
    log.info("─" * 56)

    write_tsv(summary_rows, report_dir / "summary.tsv", log=log)

    log.info("─" * 56)
    log.info("M04 complete.")
    log.info(f"  Targets    : {len(CONCAT_MAP)}")
    log.info(f"  Concat     → {concat_dir}")
    log.info(f"  MSA        → {msa_dir}")
    log.info(f"  Conservation → {cons_dir}")

    write_checkpoint("M04_alignment", report_dir, log=log)


if __name__ == "__main__":
    main()
