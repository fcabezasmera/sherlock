#!/usr/bin/env python3
"""
M07_primers.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M07 · RT-RPA Primer Co-Design (PrimedRPA)

PURPOSE
    Design RT-RPA primer pairs co-localized with crRNA candidates.

    Strategy (inverse co-design):
      1. PrimedRPA on full MSA → all possible primer sets genome-wide
      2. Load full ADAPT TSV (all windows, not just top-10) → all crRNAs
      3. For each amplicon: find best crRNA falling within boundaries
      4. Rank co-designs by crRNA activity × primer quality
      5. Top sets per target exported

    This inverse approach is robust because:
    - PrimedRPA evaluates conservation over full alignment
    - ADAPT activity covers all positions (5,000+ unique crRNA positions)
    - Co-selection is done in post-processing, not by constraining either tool

    T7 promoter prepended to FP for RT-NASBA transcription.

    Conda env: RPA (PrimedRPA)
    Input: MSA alignments from M04, ADAPT TSVs from M06

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M07_primers.py --config ~/sherlock/config.yaml

    --skip-primedrpa   reuse existing PrimedRPA outputs

OUTPUT
    main/data/07_primers/{target}/{target}_Output_Sets.csv  raw PrimedRPA
    main/data/07_primers/{target}_codesign.tsv              co-designed sets
    main/data/07_primers/all_codesign.tsv                   all targets
    main/logs/M07_primers.log
    main/reports/M07_primers/summary.tsv
    main/reports/M07_primers/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG    = "~/sherlock/config.yaml"
THREADS           = 8
PRIMEDRPA_ENV     = "RPA"

PRIMER_LEN        = 32
IDENTITY_THRESH   = 0.95
AMPLICON_MAX      = 250
AMPLICON_MIN      = 100
NT_REPEAT_LIMIT   = 5
MIN_GC            = 40    # integer (%)
MAX_GC            = 60    # integer (%)
MAX_SETS          = 100
CRNA_MARGIN       = 5     # nt inside amplicon boundary

T7_PROMOTER       = "AATTCTAATACGACTCACTATAGG"
SPACER_LEN        = 28

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
import re
import subprocess
import sys
from pathlib import Path

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
# STEP 1 — Run PrimedRPA on full MSA
# =============================================================================

def run_primedrpa(aln_path: Path, run_id: str, out_dir: Path,
                  rpa_env: str, threads: int, log) -> bool:
    """Run PrimedRPA on full MSA alignment."""
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "conda", "run", "--no-capture-output", "-n", rpa_env,
        "PrimedRPA",
        "--RunID",               run_id,
        "--InputFile",           str(aln_path),
        "--InputFileType",       "MSA",
        "--IdentityThreshold",   str(IDENTITY_THRESH),
        "--PrimerLength",        str(PRIMER_LEN),
        "--ProbeRequired",       "NO",
        "--AmpliconSizeLimit",   str(AMPLICON_MAX),
        "--NucleotideRepeatLimit", str(NT_REPEAT_LIMIT),
        "--MinGC",               str(MIN_GC),
        "--MaxGC",               str(MAX_GC),
        "--MaxSets",             str(MAX_SETS),
        "--Threads",             str(threads),
    ]

    log.info("CMD: " + " ".join(cmd))
    try:
        r = subprocess.run(cmd, cwd=str(out_dir), timeout=3600)
        return r.returncode == 0
    except subprocess.TimeoutExpired:
        log.error("PrimedRPA timeout (>1h)")
        return False
    except Exception as e:
        log.error(f"PrimedRPA error: {e}")
        return False


def find_primedrpa_output(out_dir: Path, run_id: str) -> dict:
    """Find PrimedRPA output CSV files."""
    files = {}
    for f in out_dir.glob(f"{run_id}*.csv"):
        name = f.name.lower()
        if "output_sets" in name:
            files["sets"] = f
        elif "oligo_binding" in name:
            files["oligos"] = f
        elif "alignment_summary" in name:
            files["summary"] = f
    return files


def parse_primer_sets(sets_csv: Path) -> pd.DataFrame:
    """Parse PrimedRPA Output_Sets.csv."""
    try:
        df = pd.read_csv(sets_csv)
    except Exception:
        return pd.DataFrame()

    rename = {
        "Forward Primer (FP)":    "fp_seq",
        "FP GC%":                 "fp_gc",
        "FP Binding Start Site":  "fp_start",
        "Reverse Primer (RP)":    "rp_seq",
        "RP GC%":                 "rp_gc",
        "RP Binding Start Site":  "rp_start",
        "Amplicon Size":          "amplicon_size",
        "Max Dimerisation Percentage Score": "dimerisation_pct",
        "Minimum Primer 3' Identity Anchor": "min_3prime_anchor",
    }
    df = df.rename(columns={k: v for k, v in rename.items()
                             if k in df.columns})
    for col in ["fp_start", "rp_start", "amplicon_size", "fp_gc", "rp_gc"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# =============================================================================
# STEP 2 — Load full ADAPT TSV (all windows)
# =============================================================================

def load_adapt_full(adapt_tsv: Path) -> pd.DataFrame:
    """
    Load complete ADAPT TSV and compute real crRNA alignment position.
    crna_position = window_start + position_within_window
    """
    if not adapt_tsv.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(adapt_tsv, sep="\t")
    except Exception:
        return pd.DataFrame()

    rename = {
        "window-start":               "window_start",
        "target-sequences":           "guide_seq",
        "guide-set-expected-activity":"adapt_activity",
        "total-frac-bound":           "adapt_frac_bound",
        "target-sequence-positions":  "pos_in_window_raw",
    }
    df = df.rename(columns={k: v for k, v in rename.items()
                             if k in df.columns})

    if "guide_seq" not in df.columns:
        return pd.DataFrame()

    def extract_pos(s):
        try:
            return int(str(s).strip("{}").split(",")[0])
        except Exception:
            return 0

    df["pos_in_window"] = df["pos_in_window_raw"].apply(extract_pos) \
        if "pos_in_window_raw" in df.columns else 0

    df["window_start"]  = pd.to_numeric(df.get("window_start", 0),
                                         errors="coerce").fillna(0).astype(int)
    df["crna_position"] = df["window_start"] + df["pos_in_window"]
    df["adapt_activity"]= pd.to_numeric(df.get("adapt_activity", 0),
                                         errors="coerce")
    df["adapt_frac_bound"] = pd.to_numeric(df.get("adapt_frac_bound", 0),
                                            errors="coerce")

    # Deduplicate by guide_seq: keep highest activity
    df = (df.sort_values("adapt_activity", ascending=False)
            .drop_duplicates(subset="guide_seq")
            .reset_index(drop=True))
    return df


# =============================================================================
# STEP 3 — Co-design
# =============================================================================

def codesign(primer_df: pd.DataFrame, adapt_df: pd.DataFrame,
             target: str, t7: str, log) -> pd.DataFrame:
    """
    For each primer set, find the best crRNA within the amplicon.
    Returns ranked co-design table.
    """
    if primer_df.empty:
        log.warning(f"  {target}: no primer sets")
        return pd.DataFrame()
    if adapt_df.empty:
        log.warning(f"  {target}: no ADAPT data")
        return pd.DataFrame()

    # Filter amplicons by size
    primer_df = primer_df[
        (primer_df["amplicon_size"] >= AMPLICON_MIN) &
        (primer_df["amplicon_size"] <= AMPLICON_MAX)
    ].copy()

    if primer_df.empty:
        log.warning(f"  {target}: no primer sets with "
                    f"amplicon {AMPLICON_MIN}-{AMPLICON_MAX} bp")
        return pd.DataFrame()

    rows = []
    for _, primer in primer_df.iterrows():
        fp_start = int(primer.get("fp_start", 0))
        rp_start = int(primer.get("rp_start", 0))
        amp_start = fp_start + CRNA_MARGIN
        amp_end   = rp_start - CRNA_MARGIN - SPACER_LEN

        # Find all crRNAs within this amplicon
        mask = ((adapt_df["crna_position"] >= amp_start) &
                (adapt_df["crna_position"] <= amp_end))
        inside = adapt_df[mask]

        if inside.empty:
            continue

        best = inside.sort_values("adapt_activity", ascending=False).iloc[0]

        fp_seq = str(primer.get("fp_seq", ""))
        rows.append({
            "target":             target,
            "fp_seq":             fp_seq,
            "fp_t7":              t7 + fp_seq,
            "fp_gc":              primer.get("fp_gc", ""),
            "fp_start":           fp_start,
            "rp_seq":             primer.get("rp_seq", ""),
            "rp_gc":              primer.get("rp_gc", ""),
            "rp_start":           rp_start,
            "amplicon_size":      primer.get("amplicon_size", ""),
            "dimerisation_pct":   primer.get("dimerisation_pct", ""),
            "crna_seq":           best["guide_seq"],
            "crna_position":      int(best["crna_position"]),
            "crna_activity":      round(float(best["adapt_activity"]), 4),
            "crna_frac_bound":    round(float(best.get("adapt_frac_bound", 0)), 4),
            "n_crna_in_amplicon": int(mask.sum()),
        })

    if not rows:
        log.warning(f"  {target}: no co-designs found")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = (df.sort_values("crna_activity", ascending=False)
            .reset_index(drop=True))
    df.insert(0, "rank", df.index + 1)
    log.info(f"  {target}: {len(df)} co-designed sets")
    return df


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
    p.add_argument("--skip-primedrpa", action="store_true",
                   help="Reuse existing PrimedRPA outputs")
    return p.parse_args()


def main():
    args    = parse_args()
    cfg     = load_config(str(Path(args.config).expanduser()))
    org     = cfg["organism"]["display"]
    paths   = cfg["paths"]
    pr_cfg  = cfg.get("primers", {})
    crna_cfg= cfg.get("crna", {})

    rpa_env = pr_cfg.get("rpa_env",    PRIMEDRPA_ENV)
    threads = cfg.get("threads", {}).get("default", THREADS)
    t7      = crna_cfg.get("t7",       T7_PROMOTER)

    # Dirs
    main_dir   = Path(paths["main"])
    msa_dir    = main_dir / "data" / "04_alignment" / "msa"
    crna_dir   = main_dir / "data" / "06_crna"
    out_dir    = main_dir / "data" / "07_primers"
    report_dir = main_dir / "reports" / "M07_primers"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M07_primers", cfg)
    print_logo("M07 · RT-RPA Primer Co-Design", organism=org)

    save_versions(
        tools       = {"PrimedRPA": "PrimedRPA"},
        python_pkgs = ["pandas"],
        report_dir  = report_dir,
        log         = log,
    )

    all_codesign = []
    summary_rows = []

    for target in tqdm(TARGETS, desc="Targets", unit="target"):
        log.info("─" * 56)
        log.info(f"TARGET: {target}")

        aln_file  = msa_dir / f"{target}.aln"
        adapt_tsv = crna_dir / f"{target}_adapt.tsv"

        if not aln_file.exists():
            log.warning(f"  Alignment not found — skipping")
            continue

        target_dir = out_dir / target
        target_dir.mkdir(parents=True, exist_ok=True)
        run_id = target

        # ----------------------------------------------------------
        # STEP 1 — PrimedRPA
        # ----------------------------------------------------------
        pfiles = find_primedrpa_output(target_dir, run_id)

        if args.skip_primedrpa and pfiles.get("sets"):
            log.info(f"  --skip-primedrpa: reusing output")
        else:
            log.info(f"  Running PrimedRPA on full MSA...")
            ok = run_primedrpa(aln_file, run_id, target_dir,
                               rpa_env, threads, log)
            if ok:
                pfiles = find_primedrpa_output(target_dir, run_id)

        primer_df = pd.DataFrame()
        if pfiles.get("sets"):
            primer_df = parse_primer_sets(pfiles["sets"])
            log.info(f"  PrimedRPA: {len(primer_df)} primer sets")
        else:
            log.warning(f"  PrimedRPA output not found")

        # ----------------------------------------------------------
        # STEP 2 — Load full ADAPT data
        # ----------------------------------------------------------
        adapt_df = load_adapt_full(adapt_tsv)
        if adapt_df.empty:
            log.warning(f"  ADAPT TSV not found: {adapt_tsv.name}")
        else:
            log.info(f"  ADAPT: {len(adapt_df)} unique crRNA sequences "
                     f"(pos {adapt_df['crna_position'].min()}"
                     f"-{adapt_df['crna_position'].max()})")

        # ----------------------------------------------------------
        # STEP 3 — Co-design
        # ----------------------------------------------------------
        codesign_df = codesign(primer_df, adapt_df, target, t7, log)

        n_codesign = 0
        if not codesign_df.empty:
            out_tsv = out_dir / f"{target}_codesign.tsv"
            codesign_df.to_csv(out_tsv, sep="\t", index=False)
            all_codesign.append(codesign_df)
            n_codesign = len(codesign_df)

        summary_rows.append({
            "target":       target,
            "primer_sets":  len(primer_df),
            "adapt_crnas":  len(adapt_df),
            "codesigns":    n_codesign,
        })

    # ----------------------------------------------------------
    # Combined output
    # ----------------------------------------------------------
    log.info("─" * 56)
    log.info("Writing combined output")

    if all_codesign:
        combined = pd.concat(all_codesign, ignore_index=True)
        combined_out = out_dir / "all_codesign.tsv"
        combined.to_csv(combined_out, sep="\t", index=False)
        log.info(f"All co-designs → {combined_out}  ({len(combined)} total)")

    write_tsv(summary_rows, report_dir / "summary.tsv", log=log)

    log.info("─" * 56)
    log.info("M07 complete.")
    for r in summary_rows:
        log.info(f"  {r['target']:<20}  "
                 f"primers={r['primer_sets']}  "
                 f"codesigns={r['codesigns']}")

    write_checkpoint("M07_primers", report_dir, log=log)


if __name__ == "__main__":
    main()
