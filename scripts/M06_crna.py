#!/usr/bin/env python3
"""
M06_crna.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M06 · crRNA Design (ADAPT)

PURPOSE
    Design crRNA candidates for each of the 9 targets using ADAPT
    (maximize-activity mode) with the pre-trained LwCas13a model.

    ADAPT scores every possible 28-nt guide across the alignment
    using a deep learning model trained on ~19,000 Cas13a guide-target
    pairs (Metsky et al. 2022). This is the gold-standard computational
    approach for SHERLOCK guide design.

    BADGERS was evaluated but excluded: its advantage (artificial
    sequence exploration) is negligible for C. difficile toxin genes
    which are highly conserved (mean Shannon entropy tcdB = 0.024).

    Per target:
      1. ADAPT maximize-activity  → activity-scored guides per window
      2. Biological filters       → GC 30-70%, poly-U ≤ 3, len = 28
      3. Accessibility filter     → mean unpaired_prob ≥ 0.5 (M05)
      4. Composite ranking        → ADAPT activity × accessibility
      5. Top 10 candidates        → output per target + combined

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M06_crna.py --config ~/sherlock/config.yaml

    --skip-adapt    reuse existing ADAPT TSV outputs

OUTPUT
    main/data/06_crna/{target}_adapt.tsv       raw ADAPT output
    main/data/06_crna/{target}_candidates.tsv  top 10 filtered candidates
    main/data/06_crna/all_candidates.tsv        all targets combined
    main/logs/M06_crna.log
    main/reports/M06_crna/summary.tsv
    main/reports/M06_crna/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG  = "~/sherlock/config.yaml"
THREADS         = 8
TOP_N           = 10
SPACER_LEN      = 28
GC_MIN          = 0.30
GC_MAX          = 0.70
POLYU_MAX       = 3       # max consecutive U in spacer
MIN_ACC         = 0.5     # minimum RNAplfold unpaired probability
ADAPT_WINDOW    = 250     # sliding window size

# The 9 MSA targets (from M04)
TARGETS = [
    "tcdA_all",
    "tcdB_clade2",
    "tcdB_clade1",
    "tcdC_wt",
    "tcdC_junction",
    "cdtA_groupA",
    "cdtB_groupA",
    "tpiA_all",
    "rpoB_all",
]

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import subprocess
import sys
from pathlib import Path

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
# BIOLOGICAL FILTERS
# =============================================================================

def gc_content(seq: str) -> float:
    seq = seq.upper().replace("U", "T")
    return sum(1 for b in seq if b in "GC") / len(seq) if seq else 0.0


def max_polyu(seq: str) -> int:
    """Max consecutive U (or T) run in sequence."""
    seq = seq.upper().replace("T", "U")
    max_run = current = 0
    for b in seq:
        if b == "U":
            current += 1
            max_run = max(max_run, current)
        else:
            current = 0
    return max_run


def passes_filters(seq: str, gc_min: float, gc_max: float,
                   polyu_max: int) -> tuple:
    """
    Apply biological filters for LwCas13a.
    Returns (passes: bool, reason: str).
    """
    if len(seq) != SPACER_LEN:
        return False, f"length={len(seq)}≠{SPACER_LEN}"
    gc = gc_content(seq)
    if gc < gc_min or gc > gc_max:
        return False, f"GC={gc:.2f} out of [{gc_min},{gc_max}]"
    pu = max_polyu(seq)
    if pu > polyu_max:
        return False, f"poly-U={pu}>{polyu_max}"
    return True, "OK"


# =============================================================================
# STEP 1 — Run ADAPT
# =============================================================================

def run_adapt(aln_path: Path, out_prefix: Path,
              adapt_env: str, log) -> Path | None:
    """
    Run ADAPT sliding-window maximize-activity on an alignment.
    Returns path to output TSV or None on failure.
    """
    out_tsv = Path(str(out_prefix) + ".tsv")

    cmd = [
        "conda", "run", "--no-capture-output", "-n", adapt_env,
        "design.py", "sliding-window", "fasta",
        str(aln_path),
        "-o", str(out_prefix),
        "--obj", "maximize-activity",
        "-w",   str(ADAPT_WINDOW),
        "-gl",  str(SPACER_LEN),
        "--predict-cas13a-activity-model",
    ]

    log.info("CMD: " + " ".join(cmd))
    try:
        r = subprocess.run(cmd, timeout=3600)
        if r.returncode != 0:
            log.error(f"ADAPT failed (exit {r.returncode})")
            return None
        if out_tsv.exists():
            log.info(f"  → {out_tsv.name}")
            return out_tsv
        # Try alternate naming
        candidates = list(out_prefix.parent.glob(
            f"{out_prefix.name}*.tsv"))
        if candidates:
            return candidates[0]
        log.error(f"ADAPT output TSV not found")
        return None
    except subprocess.TimeoutExpired:
        log.error(f"ADAPT timeout (>1h) for {aln_path.name}")
        return None
    except Exception as e:
        log.error(f"ADAPT error: {e}")
        return None


# =============================================================================
# STEP 2 — Parse ADAPT output
# =============================================================================

def parse_adapt(tsv_path: Path) -> pd.DataFrame:
    """
    Parse ADAPT maximize-activity TSV.

    Key columns used:
      window-start               → guide start position in alignment
      target-sequences           → 28-nt guide sequence
      guide-set-expected-activity → predicted Cas13a activity score
      total-frac-bound           → fraction of sequences detected
    """
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        return pd.DataFrame()

    # Rename for clarity
    rename = {
        "window-start":               "window_start",
        "window-end":                 "window_end",
        "target-sequences":           "guide_seq",
        "guide-set-expected-activity":"adapt_activity",
        "guide-expected-activities":  "adapt_activity_raw",
        "total-frac-bound":           "adapt_frac_bound",
        "target-sequence-positions":  "guide_pos_in_window",
    }
    df = df.rename(columns={k: v for k, v in rename.items()
                             if k in df.columns})

    if "adapt_activity" not in df.columns:
        return pd.DataFrame()

    df["adapt_activity"]   = pd.to_numeric(df["adapt_activity"],   errors="coerce")
    df["adapt_frac_bound"] = pd.to_numeric(df.get("adapt_frac_bound",
                                                    pd.Series()), errors="coerce")

    # Extract real crRNA position: window_start + position_within_window
    # ADAPT reports position as {N} e.g. {188}
    def extract_pos(s):
        try:
            return int(str(s).strip("{}").split(",")[0])
        except Exception:
            return 0

    if "guide_pos_in_window" in df.columns:
        df["guide_pos_in_window"] = df["guide_pos_in_window"].apply(extract_pos)
        df["crna_position"] = (pd.to_numeric(df["window_start"], errors="coerce")
                               + df["guide_pos_in_window"])
    else:
        df["crna_position"] = pd.to_numeric(df["window_start"], errors="coerce")

    return df


# =============================================================================
# STEP 3 — Load accessibility from M05
# =============================================================================

def load_accessibility(acc_dir: Path, target: str) -> pd.DataFrame:
    """Load RNAplfold accessibility profile for a target."""
    f = acc_dir / f"{target}_accessibility.tsv"
    if not f.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(f, sep="\t")
    except Exception:
        return pd.DataFrame()


def window_accessibility(acc_df: pd.DataFrame,
                          window_start: int) -> float:
    """
    Mean unpaired probability for guide window [window_start, window_start+28].
    acc_df positions are 1-based; window_start is 0-based.
    """
    if acc_df.empty:
        return float("nan")
    mask = ((acc_df["position"] >= window_start + 1) &
            (acc_df["position"] <= window_start + SPACER_LEN))
    vals = acc_df.loc[mask, "unpaired_prob"]
    return float(vals.mean()) if not vals.empty else float("nan")


# =============================================================================
# STEP 4 — Filter, rank, top-N
# =============================================================================

def rank_candidates(adapt_df: pd.DataFrame,
                    acc_df: pd.DataFrame,
                    target: str,
                    top_n: int,
                    gc_min: float, gc_max: float,
                    polyu_max: int, min_acc: float,
                    log) -> pd.DataFrame:
    """
    Apply filters, add accessibility, compute composite score, return top-N.

    Composite score = 0.60 × adapt_norm + 0.40 × acc_norm
    (ADAPT activity weighted higher as it directly predicts Cas13a performance)
    """
    if adapt_df.empty or "guide_seq" not in adapt_df.columns:
        log.warning(f"  {target}: no ADAPT candidates")
        return pd.DataFrame()

    rows = []
    for _, row in adapt_df.iterrows():
        seq = str(row.get("guide_seq", "")).strip()
        # R7: total U fraction filter (Milligan 1987)
        u_frac_val = seq.upper().replace("T","U").count("U") / max(len(seq),1)
        if u_frac_val > 0.50:
            rows.append({"target":target,"guide_seq":seq,"window_start":int(row.get("window_start",0)),
                "window_end":int(row.get("window_end",0)),"crna_position":int(row.get("crna_position",row.get("window_start",0))),
                "adapt_activity":float(row.get("adapt_activity",float("nan"))),"adapt_frac_bound":float(row.get("adapt_frac_bound",float("nan"))),
                "accessibility":0.0,"gc_content":round(gc_content(seq),3),"max_polyu":max_polyu(seq),
                "filter_pass":False,"filter_reason":"high_U_frac"})
            continue
        ok, reason = passes_filters(seq, gc_min, gc_max, polyu_max)
        acc = window_accessibility(acc_df, int(row.get("window_start", 0)))

        rows.append({
            "target":          target,
            "guide_seq":       seq,
            "window_start":    int(row.get("window_start", 0)),
            "window_end":      int(row.get("window_end", 0)),
            "crna_position":   int(row.get("crna_position",
                               row.get("window_start", 0))),
            "adapt_activity":  float(row.get("adapt_activity", float("nan"))),
            "adapt_frac_bound":float(row.get("adapt_frac_bound", float("nan"))),
            "accessibility":   acc,
            "gc_content":      round(gc_content(seq), 3),
            "max_polyu":       max_polyu(seq),
            "filter_pass":     ok,
            "filter_reason":   reason,
        })

    df = pd.DataFrame(rows)

    # Biological filters
    df_pass = df[df["filter_pass"]].copy()

    # Accessibility filter (skip if no profile available)
    if not acc_df.empty:
        df_pass = df_pass[
            df_pass["accessibility"].isna() |
            (df_pass["accessibility"] >= min_acc)
        ].copy()

    if df_pass.empty:
        log.warning(f"  {target}: no candidates passed filters")
        return pd.DataFrame()

    # Normalize scores
    def norm01(s):
        v = s.dropna()
        if v.empty or v.max() == v.min():
            return pd.Series(0.5, index=s.index)
        return (s - v.min()) / (v.max() - v.min())

    df_pass["adapt_norm"] = norm01(df_pass["adapt_activity"])
    df_pass["acc_norm"]   = norm01(df_pass["accessibility"])
    df_pass["adapt_norm"] = df_pass["adapt_norm"].fillna(0.5)
    df_pass["acc_norm"]   = df_pass["acc_norm"].fillna(0.5)

    # Composite score
    # R6: G at position 1 penalty (Wessels 2020, Mol Cell)
    df_pass["g_pos1_pen"] = df_pass["guide_seq"].apply(
        lambda s: -0.05 if s.upper().replace("T","U")[:1] == "G" else 0.0)
    # R8: GC score for LwCas13a (25-65% optimal; C. diff genome 29% GC)
    def gc_score(seq):
        gc = sum(1 for b in seq.upper() if b in "GC") / max(len(seq),1)
        if 0.25 <= gc <= 0.65: return 1.0
        if gc < 0.20 or gc > 0.75: return 0.0
        return 0.5
    df_pass["gc_score"] = df_pass["guide_seq"].apply(gc_score)
    # Composite score: ADAPT 55% + accessibility 35% + GC 10% + G-pos1 penalty
    df_pass["score"] = (0.55 * df_pass["adapt_norm"] +
                        0.35 * df_pass["acc_norm"] +
                        0.10 * df_pass["gc_score"] +
                        df_pass["g_pos1_pen"])

    # Deduplicate → top N
    df_pass = (df_pass
               .sort_values("score", ascending=False)
               .drop_duplicates(subset="guide_seq")
               .head(top_n)
               .reset_index(drop=True))
    df_pass["rank"] = df_pass.index + 1

    log.info(f"  {target}: {len(df_pass)} candidates selected (top {top_n})")
    return df_pass


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
    p.add_argument("--skip-adapt", action="store_true",
                   help="Reuse existing ADAPT TSV outputs")
    return p.parse_args()


def main():
    args     = parse_args()
    cfg      = load_config(str(Path(args.config).expanduser()))
    org      = cfg["organism"]["display"]
    paths    = cfg["paths"]
    conda    = cfg.get("conda", {})
    crna_cfg = cfg.get("crna", {})

    adapt_env  = conda.get("adapt", "adapt")
    top_n      = crna_cfg.get("top_n",      TOP_N)
    gc_min     = crna_cfg.get("gc_min",     GC_MIN)
    gc_max     = crna_cfg.get("gc_max",     GC_MAX)
    polyu_max  = crna_cfg.get("polyu_max",  POLYU_MAX)
    min_acc    = crna_cfg.get("min_acc",    MIN_ACC)

    # Dirs
    main_dir   = Path(paths["main"])
    msa_dir    = main_dir / "data" / "04_alignment" / "msa"
    acc_dir    = main_dir / "data" / "05_accessibility"
    out_dir    = main_dir / "data" / "06_crna"
    report_dir = main_dir / "reports" / "M06_crna"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M06_crna", cfg)
    print_logo("M06 · crRNA Design (ADAPT)", organism=org)

    save_versions(
        tools       = {},
        python_pkgs = ["pandas", "numpy"],
        report_dir  = report_dir,
        log         = log,
    )

    all_candidates = []
    summary_rows   = []

    for target in tqdm(TARGETS, desc="Targets", unit="target"):
        log.info("─" * 56)
        log.info(f"TARGET: {target}")

        aln_file = msa_dir / f"{target}.aln"
        if not aln_file.exists():
            log.warning(f"  Alignment not found: {aln_file} — skipping")
            continue

        # ----------------------------------------------------------
        # STEP 1 — ADAPT
        # ----------------------------------------------------------
        adapt_prefix = out_dir / f"{target}_adapt"
        adapt_tsv    = Path(str(adapt_prefix) + ".tsv")

        if args.skip_adapt and adapt_tsv.exists():
            log.info(f"  --skip-adapt: reusing {adapt_tsv.name}")
        else:
            log.info(f"  Running ADAPT maximize-activity...")
            result = run_adapt(aln_file, adapt_prefix, adapt_env, log)
            if result:
                adapt_tsv = result

        adapt_df = pd.DataFrame()
        if adapt_tsv.exists():
            adapt_df = parse_adapt(adapt_tsv)
            log.info(f"  ADAPT: {len(adapt_df)} windows")

        # ----------------------------------------------------------
        # STEP 2 — Accessibility
        # ----------------------------------------------------------
        acc_df = load_accessibility(acc_dir, target)
        if acc_df.empty:
            log.warning(f"  No accessibility profile for {target}")

        # ----------------------------------------------------------
        # STEP 3 — Filter + rank
        # ----------------------------------------------------------
        candidates = rank_candidates(
            adapt_df  = adapt_df,
            acc_df    = acc_df,
            target    = target,
            top_n     = top_n,
            gc_min    = gc_min,
            gc_max    = gc_max,
            polyu_max = polyu_max,
            min_acc   = min_acc,
            log       = log,
        )

        n_cands = 0
        if not candidates.empty:
            target_out = out_dir / f"{target}_candidates.tsv"
            candidates.to_csv(target_out, sep="\t", index=False)
            all_candidates.append(candidates)
            n_cands = len(candidates)

        summary_rows.append({
            "target":          target,
            "adapt_windows":   len(adapt_df),
            "candidates":      n_cands,
        })

    # ----------------------------------------------------------
    # Combined output
    # ----------------------------------------------------------
    log.info("─" * 56)
    log.info("Writing combined output")

    if all_candidates:
        combined = pd.concat(all_candidates, ignore_index=True)
        combined_out = out_dir / "all_candidates.tsv"
        combined.to_csv(combined_out, sep="\t", index=False)
        log.info(f"All candidates → {combined_out}  ({len(combined)} total)")
    else:
        log.warning("No candidates generated")

    write_tsv(summary_rows, report_dir / "summary.tsv", log=log)

    log.info("─" * 56)
    log.info("M06 complete.")
    for r in summary_rows:
        log.info(f"  {r['target']:<20}  candidates={r['candidates']}")

    write_checkpoint("M06_crna", report_dir, log=log)


if __name__ == "__main__":
    main()
