#!/usr/bin/env python3
"""
M11_rtqpcr.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M11 · RT-qPCR Reference Primers

PURPOSE
    Design RT-qPCR reference primer pairs for the same target genes
    used in SHERLOCK. These serve as:
      1. Gold-standard validation method for SHERLOCK results
      2. Positive/negative control confirmation in clinical samples
      3. Quantification reference for sensitivity/specificity comparison

    Uses Primer3 on consensus sequences from M05.
    Targets: tcdA, tcdB, tcdC, cdtA, cdtB, tpiA, sodA

    Amplicon: 80-200 bp (optimal for qPCR)
    Tm: 58-62°C (standard qPCR conditions)
    GC: 45-65%

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M11_rtqpcr.py --config ~/sherlock/config.yaml

OUTPUT
    main/data/11_rtqpcr/{gene}_primers.tsv    primer pairs per gene
    main/data/11_rtqpcr/all_rtqpcr_primers.tsv  combined
    main/logs/M11_rtqpcr.log
    main/reports/M11_rtqpcr/summary.tsv
    main/reports/M11_rtqpcr/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG    = "~/sherlock/config.yaml"

PRIMER_OPT_SIZE   = 20
PRIMER_MIN_SIZE   = 18
PRIMER_MAX_SIZE   = 25
PRIMER_OPT_TM     = 60.0
PRIMER_MIN_TM     = 58.0
PRIMER_MAX_TM     = 62.0
PRIMER_MIN_GC     = 45.0
PRIMER_MAX_GC     = 65.0
PRODUCT_SIZE_MIN  = 80
PRODUCT_SIZE_MAX  = 200
NUM_RETURN        = 10

# Gene → MSA target mapping (consensus from M05)
RTQPCR_TARGETS = {
    "tcdA": "tcdA_all",
    "tcdB": "tcdB_all",
    "tcdC": "tcdC_wt",
    "cdtA": "cdtA_groupA",
    "cdtB": "cdtB_groupA",   # added
    "tpiA": "tpiA_all",
    "sodA": "sodA_all",
}

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from pipeline_utils import (
    load_config, get_logger, print_logo,
    save_versions, write_tsv, write_checkpoint,
)

# =============================================================================
# PRIMER3
# =============================================================================

def load_consensus(acc_dir: Path, target: str) -> str:
    f = acc_dir / f"{target}_consensus.fasta"
    if not f.exists():
        return ""
    lines = f.read_text(encoding="utf-8").splitlines()
    return "".join(l for l in lines if not l.startswith(">")).upper().replace("U","T")


def run_primer3(sequence: str, gene_name: str, primer3_exe: str,
                amp_min: int, amp_max: int,
                tm_min: float, tm_max: float, tm_opt: float) -> list:
    if not sequence or len(sequence) < amp_max + 50:
        return []

    p3_input = f"""SEQUENCE_ID={gene_name}
SEQUENCE_TEMPLATE={sequence}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE={PRIMER_OPT_SIZE}
PRIMER_MIN_SIZE={PRIMER_MIN_SIZE}
PRIMER_MAX_SIZE={PRIMER_MAX_SIZE}
PRIMER_OPT_TM={tm_opt}
PRIMER_MIN_TM={tm_min}
PRIMER_MAX_TM={tm_max}
PRIMER_MIN_GC={PRIMER_MIN_GC}
PRIMER_MAX_GC={PRIMER_MAX_GC}
PRIMER_PRODUCT_SIZE_RANGE={amp_min}-{amp_max}
PRIMER_NUM_RETURN={NUM_RETURN}
PRIMER_EXPLAIN_FLAG=1
="""

    try:
        r = subprocess.run([primer3_exe], input=p3_input,
                           capture_output=True, text=True, timeout=60)
        return parse_primer3_output(r.stdout, gene_name)
    except Exception:
        return []


def parse_primer3_output(output: str, gene_name: str) -> list:
    data = {}
    for line in output.splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            data[k.strip()] = v.strip()

    pairs = []
    i = 0
    while True:
        if f"PRIMER_LEFT_{i}_SEQUENCE" not in data:
            break
        fp_seq  = data.get(f"PRIMER_LEFT_{i}_SEQUENCE",  "")
        rp_seq  = data.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
        fp_tm   = float(data.get(f"PRIMER_LEFT_{i}_TM",  0))
        rp_tm   = float(data.get(f"PRIMER_RIGHT_{i}_TM", 0))
        fp_gc   = float(data.get(f"PRIMER_LEFT_{i}_GC_PERCENT",  0))
        rp_gc   = float(data.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0))
        prod    = int(data.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE",   0))
        penalty = float(data.get(f"PRIMER_PAIR_{i}_PENALTY",      0))
        if fp_seq and rp_seq:
            pairs.append({
                "gene":         gene_name,
                "pair_rank":    i + 1,
                "fp_seq":       fp_seq,
                "rp_seq":       rp_seq,
                "fp_tm":        round(fp_tm, 2),
                "rp_tm":        round(rp_tm, 2),
                "fp_gc":        round(fp_gc, 1),
                "rp_gc":        round(rp_gc, 1),
                "amplicon_size":prod,
                "penalty":      round(penalty, 4),
            })
        i += 1
    return pairs


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=DEFAULT_CONFIG)
    return p.parse_args()


def main():
    args     = parse_args()
    cfg      = load_config(str(Path(args.config).expanduser()))
    org      = cfg["organism"]["display"]
    paths    = cfg["paths"]
    tools    = cfg["tools"]
    qpcr_cfg = cfg.get("rtqpcr", {})

    primer3_exe = tools.get("primer3", "primer3_core")
    amp_min  = qpcr_cfg.get("amplicon_min", PRODUCT_SIZE_MIN)
    amp_max  = qpcr_cfg.get("amplicon_max", PRODUCT_SIZE_MAX)
    tm_tgt   = qpcr_cfg.get("tm_target",   PRIMER_OPT_TM)
    tm_tol   = qpcr_cfg.get("tm_tol",      2.0)
    # Override targets from config if set, else use all
    genes    = qpcr_cfg.get("targets", list(RTQPCR_TARGETS.keys()))
    # Ensure cdtB is always included
    if "cdtB" not in genes:
        genes.append("cdtB")

    # Dirs
    main_dir   = Path(paths["main"])
    acc_dir    = main_dir / "data" / "05_accessibility"
    out_dir    = main_dir / "data" / "11_rtqpcr"
    rep_dir    = main_dir / "reports" / "M11_rtqpcr"
    out_dir.mkdir(parents=True, exist_ok=True)
    rep_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M11_rtqpcr", cfg)
    print_logo("M11 · RT-qPCR Reference Primers", organism=org)

    save_versions(tools={"primer3": primer3_exe},
                  python_pkgs=["pandas"],
                  report_dir=rep_dir, log=log)

    all_pairs    = []
    summary_rows = []

    for gene in genes:
        msa_target = RTQPCR_TARGETS.get(gene, f"{gene}_all")
        log.info(f"  {gene} → consensus from {msa_target}")

        consensus = load_consensus(acc_dir, msa_target)
        if not consensus:
            log.warning(f"  No consensus for {msa_target} — skipping")
            summary_rows.append({
                "gene": gene, "msa_target": msa_target,
                "n_pairs": 0, "status": "NO_CONSENSUS"})
            continue

        log.info(f"  Consensus: {len(consensus)} nt")

        pairs = run_primer3(consensus, gene, primer3_exe,
                            amp_min, amp_max,
                            tm_tgt - tm_tol, tm_tgt + tm_tol, tm_tgt)
        log.info(f"  Primer3: {len(pairs)} pairs")

        if pairs:
            df = pd.DataFrame(pairs)
            df.to_csv(out_dir / f"{gene}_primers.tsv", sep="\t", index=False)
            all_pairs.append(df)
            best = pairs[0]
            log.info(f"  Top: FP={best['fp_seq']} Tm={best['fp_tm']}°C | "
                     f"RP={best['rp_seq']} Tm={best['rp_tm']}°C | "
                     f"Amp={best['amplicon_size']}bp")

        summary_rows.append({
            "gene":       gene,
            "msa_target": msa_target,
            "n_pairs":    len(pairs),
            "status":     "OK" if pairs else "NO_PAIRS",
        })

    if all_pairs:
        combined = pd.concat(all_pairs, ignore_index=True)
        combined.to_csv(out_dir / "all_rtqpcr_primers.tsv",
                        sep="\t", index=False)
        log.info(f"All primers → {out_dir/'all_rtqpcr_primers.tsv'} "
                 f"({len(combined)} pairs)")

    write_tsv(summary_rows, rep_dir / "summary.tsv", log=log)

    log.info("─" * 56)
    log.info("M11 complete.")
    for r in summary_rows:
        log.info(f"  {r['gene']:<8}  pairs={r['n_pairs']}  {r['status']}")

    write_checkpoint("M11_rtqpcr", rep_dir, log=log)


if __name__ == "__main__":
    main()
