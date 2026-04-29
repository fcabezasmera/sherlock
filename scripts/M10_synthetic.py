#!/usr/bin/env python3
"""
M10_synthetic.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M10 · Synthetic Control Sequences

PURPOSE
    Generate synthetic positive control sequences for SHERLOCK validation.
    For each final crRNA candidate, produces:

    1. crRNA sequence with LwCas13a direct repeat (DR + spacer)
       Format: T7 + DR + spacer (for in vitro transcription)

    2. Synthetic target RNA template
       Format: T7 + flank + target + flank
       This is the synthetic RNA transcript that activates the crRNA

    3. IDT ordering format
       All sequences formatted for IDT gBlock or Ultramer ordering

    LwCas13a DR: GGGGAUUUAGACUACCCCAAAAACGAAGGGGGGACUAAAAC
    T7 promoter: AATTCTAATACGACTCACTATAGG

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M10_synthetic.py --config ~/sherlock/config.yaml

OUTPUT
    main/data/10_synthetic/crRNA_sequences.tsv    crRNA + DR sequences
    main/data/10_synthetic/target_templates.tsv   synthetic target templates
    main/data/10_synthetic/IDT_order.tsv          IDT ordering sheet
    main/logs/M10_synthetic.log
    main/reports/M10_synthetic/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG = "~/sherlock/config.yaml"
FLANK_NT       = 100   # flanking nt on each side of target in synthetic template
IDT_MAX_LEN    = 200   # max length for IDT Ultramer (longer = gBlock)

# LwCas13a DR (RNA sequence, 5'→3')
DR_RNA  = "GGGGAUUUAGACUACCCCAAAAACGAAGGGGGGACUAAAAC"
# T7 promoter (DNA, added upstream for IVT)
T7_DNA  = "AATTCTAATACGACTCACTATAGG"

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import sys
from pathlib import Path

import pandas as pd

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
# SEQUENCE UTILITIES
# =============================================================================

def dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def rna_to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def reverse_complement(seq: str) -> str:
    comp = {"A":"T","T":"A","G":"C","C":"G","U":"A","N":"N"}
    return "".join(comp.get(b, "N") for b in seq.upper()[::-1])


def make_crna_sequence(spacer_rna: str, dr: str) -> str:
    """
    Build full crRNA: DR + spacer (RNA)
    LwCas13a: direct repeat is 5' of spacer.
    """
    return dr + spacer_rna


def make_ivt_template(crna_rna: str, t7: str) -> str:
    """
    Build IVT template DNA for crRNA:
    T7 promoter + reverse complement of crRNA
    (T7 reads 3'→5' so template is RC of crRNA)
    """
    crna_dna = rna_to_dna(crna_rna)
    return t7 + reverse_complement(crna_dna)


def make_synthetic_target(target_seq: str, flank: int = FLANK_NT) -> str:
    """
    Build synthetic target template for IVT:
    T7 + flank + target + flank
    target_seq should be the crRNA target region (DNA)
    """
    # Use poly-A flanks for simple synthetic template
    flank_seq = "A" * flank
    return T7_DNA + flank_seq + rna_to_dna(target_seq) + flank_seq


def idt_format(name: str, seq: str) -> str:
    """Format sequence for IDT order."""
    seq_type = "Ultramer" if len(seq) <= IDT_MAX_LEN else "gBlock"
    return seq_type


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=DEFAULT_CONFIG)
    return p.parse_args()


def main():
    args   = parse_args()
    cfg    = load_config(str(Path(args.config).expanduser()))
    org    = cfg["organism"]["display"]
    paths  = cfg["paths"]
    crna_cfg = cfg.get("crna", {})

    dr  = crna_cfg.get("dr",  DR_RNA).replace("T", "U")
    t7  = crna_cfg.get("t7",  T7_DNA)
    flank = cfg.get("synthetic", {}).get("flank_nt", FLANK_NT)

    # Dirs
    main_dir   = Path(paths["main"])
    report_dir = main_dir / "data" / "09_report"
    out_dir    = main_dir / "data" / "10_synthetic"
    rep_dir    = main_dir / "reports" / "M10_synthetic"
    out_dir.mkdir(parents=True, exist_ok=True)
    rep_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M10_synthetic", cfg)
    print_logo("M10 · Synthetic Controls", organism=org)

    # Load final ranking from M09
    ranking_tsv = report_dir / "final_ranking.tsv"
    if not ranking_tsv.exists():
        log.error("final_ranking.tsv not found. Run M09 first.")
        sys.exit(1)

    df = pd.read_csv(ranking_tsv, sep="\t")
    log.info(f"Loaded {len(df)} final candidates from M09")

    crna_rows   = []
    target_rows = []
    idt_rows    = []

    for _, row in df.iterrows():
        target    = str(row["target"])
        guide_seq = str(row.get("guide_seq_rna", row.get("guide_seq", "")))  # RNA spacer (U not T)
        rank      = int(row["rank"])

        # Ensure RNA format
        spacer_rna = guide_seq.upper().replace("T", "U")
        spacer_dna = spacer_rna.replace("U", "T")

        # 1. Full crRNA (DR + spacer)
        crna_rna  = make_crna_sequence(spacer_rna, dr)
        crna_name = f"{target}_rank{rank}_crRNA"

        # 2. IVT template for crRNA
        ivt_crna  = make_ivt_template(crna_rna, t7)
        ivt_name  = f"{target}_rank{rank}_IVT_crRNA"

        # 3. Synthetic target template
        syn_target= make_synthetic_target(spacer_dna, flank)
        syn_name  = f"{target}_rank{rank}_target_template"

        crna_rows.append({
            "name":         crna_name,
            "target":       target,
            "rank":         rank,
            "spacer_rna":   spacer_rna,
            "spacer_len":   len(spacer_rna),
            "dr":           dr,
            "full_crna_rna":crna_rna,
            "full_crna_len":len(crna_rna),
            "ivt_template_dna": ivt_crna,
            "ivt_len":      len(ivt_crna),
        })

        target_rows.append({
            "name":         syn_name,
            "target":       target,
            "rank":         rank,
            "spacer_rna":   spacer_rna,
            "template_dna": syn_target,
            "template_len": len(syn_target),
        })

        # IDT ordering rows
        idt_rows.append({
            "Name":         ivt_name,
            "Sequence":     ivt_crna,
            "Scale":        "25nm",
            "Purification": "STD",
            "Type":         idt_format(ivt_name, ivt_crna),
            "Notes":        f"IVT template for {target} crRNA rank {rank}",
        })
        idt_rows.append({
            "Name":         syn_name,
            "Sequence":     syn_target,
            "Scale":        "25nm",
            "Purification": "STD",
            "Type":         idt_format(syn_name, syn_target),
            "Notes":        f"Synthetic target template for {target} rank {rank}",
        })

    # Write outputs
    crna_out   = out_dir / "crRNA_sequences.tsv"
    target_out = out_dir / "target_templates.tsv"
    idt_out    = out_dir / "IDT_order.tsv"

    pd.DataFrame(crna_rows).to_csv(crna_out,   sep="\t", index=False)
    pd.DataFrame(target_rows).to_csv(target_out,sep="\t", index=False)
    pd.DataFrame(idt_rows).to_csv(idt_out,     sep="\t", index=False)

    log.info(f"crRNA sequences  → {crna_out}   ({len(crna_rows)} sequences)")
    log.info(f"Target templates → {target_out}  ({len(target_rows)} templates)")
    log.info(f"IDT order sheet  → {idt_out}     ({len(idt_rows)} oligos)")

    write_tsv([
        {"item": "n_candidates",    "value": len(df)},
        {"item": "n_crna_seqs",     "value": len(crna_rows)},
        {"item": "n_target_templates","value": len(target_rows)},
        {"item": "n_idt_oligos",    "value": len(idt_rows)},
        {"item": "dr_sequence",     "value": dr},
        {"item": "t7_promoter",     "value": t7},
        {"item": "flank_nt",        "value": flank},
    ], rep_dir / "summary.tsv", log=log)

    log.info("M10 complete.")
    write_checkpoint("M10_synthetic", rep_dir, log=log)


if __name__ == "__main__":
    main()
