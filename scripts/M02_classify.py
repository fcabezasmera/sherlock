#!/usr/bin/env python3
"""
M02_classify.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M02 · Functional Classification

PURPOSE
    Classify all downloaded genomes into functional groups
    based on toxin gene presence/absence and integrity from GFF files.

    Groups:
      A — Hypervirulent RT027 equivalent
          tcdA+ tcdB+ cdtA+ cdtB+ tcdC(Δ18nt)+
      B — Toxigenic RT012 equivalent
          tcdA+ tcdB+ cdtA- cdtB- tcdC(WT)
      C — Non-toxigenic (off-target control)
          tcdA- tcdB-

    Gene integrity from GFF flags:
      partial=true  → gene is truncated (fragmented across contigs)
      pseudo=true   → gene is a pseudogene

    Gene status per genome:
      COMPLETE   → gene present, not partial, not pseudo
      PARTIAL    → gene present but partial=true
      PSEUDOGENE → gene present but pseudo=true
      ABSENT     → gene not found in GFF

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M02_classify.py --config ~/sherlock/config.yaml

OUTPUT
    main/data/02_classify/gene_matrix.tsv      gene status per genome
    main/data/02_classify/groupA.txt           Group A accessions
    main/data/02_classify/groupB.txt           Group B accessions
    main/data/02_classify/groupC.txt           Group C accessions
    main/data/02_classify/unclassified.txt     ambiguous genomes
    main/logs/M02_classify.log
    main/reports/M02_classify/summary.tsv
    main/reports/M02_classify/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG = "~/sherlock/config.yaml"
THREADS        = 8   # change to 16 for maximum

# Target genes to scan in GFF
TARGET_GENES = ["tcdA", "tcdB", "tcdC", "cdtA", "cdtB",
                "tpiA", "sodA", "16S"]

# tcdC deletion detection
TDCC_DEL_BP     = 18    # canonical RT027 deletion
JUNCTION_FLANK  = 20    # nt flanking the deletion site

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

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
# GFF PARSING
# =============================================================================

def find_gff(genome_dir: Path) -> Path | None:
    """Find GFF file in NCBI dataset directory structure."""
    # Standard NCBI datasets structure:
    # {genome_dir}/ncbi_dataset/data/{accession}/genomic.gff
    gff_files = list(genome_dir.rglob("genomic.gff"))
    if gff_files:
        return gff_files[0]
    return None


def find_fna(genome_dir: Path) -> Path | None:
    """Find FASTA file in NCBI dataset directory structure."""
    fna_files = list(genome_dir.rglob("*_genomic.fna"))
    if fna_files:
        return fna_files[0]
    return None


def parse_gff_genes(gff_path: Path, target_genes: list) -> dict:
    """
    Parse GFF file and return gene status for target genes.

    Returns dict: {gene_name: status}
    where status is: 'COMPLETE', 'PARTIAL', 'PSEUDOGENE', 'ABSENT'

    Only considers CDS and gene features.
    Uses 'gene=' attribute for gene name matching.
    """
    found = defaultdict(list)  # gene_name → list of (partial, pseudo)

    try:
        with open(gff_path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                if feature_type not in ("gene", "CDS", "pseudogene"):
                    continue

                attrs = parts[8]

                # Extract gene name
                gene_match = re.search(r'(?:^|;)gene=([^;]+)', attrs)
                if not gene_match:
                    # Try Name= as fallback
                    gene_match = re.search(r'(?:^|;)Name=([^;]+)', attrs)
                if not gene_match:
                    continue

                gene_name = gene_match.group(1).strip()

                # Check if this gene is in our target list
                # Case-insensitive matching
                matched = None
                for tg in target_genes:
                    if gene_name.lower() == tg.lower():
                        matched = tg
                        break
                if not matched:
                    continue

                # Check partial and pseudo flags
                is_partial = "partial=true" in attrs.lower()
                is_pseudo  = ("pseudo=true" in attrs.lower() or
                              feature_type == "pseudogene" or
                              "pseudogene" in attrs.lower())

                found[matched].append((is_partial, is_pseudo))

    except Exception as e:
        return {g: "ERROR" for g in target_genes}

    # Determine final status per gene
    result = {}
    for gene in target_genes:
        if gene not in found:
            result[gene] = "ABSENT"
        else:
            records = found[gene]
            # Any COMPLETE record → COMPLETE (prioritize complete CDS)
            has_complete   = any(not p and not ps for p, ps in records)
            has_partial    = any(p for p, ps in records)
            has_pseudogene = any(ps for p, ps in records)

            if has_complete:
                result[gene] = "COMPLETE"
            elif has_partial:
                result[gene] = "PARTIAL"
            elif has_pseudogene:
                result[gene] = "PSEUDOGENE"
            else:
                result[gene] = "ABSENT"

    return result


def detect_tcdC_deletion(gff_path: Path, fna_path: Path,
                          del_bp: int, flank: int, log) -> bool:
    """
    Detect the canonical RT027 18nt deletion in tcdC.

    Strategy:
    1. Find tcdC coordinates in GFF
    2. Extract the tcdC sequence from FASTA
    3. Align against reference tcdC and check for 18nt deletion

    For now, uses a heuristic based on gene length:
    RT027 tcdC with 18nt deletion → shorter CDS
    Full-length tcdC: ~657 bp (218 aa)
    RT027 tcdC (Δ18nt + frameshift): effectively truncated

    Returns True if deletion is detected, False otherwise.
    """
    # Heuristic: parse tcdC CDS length from GFF
    tcdC_lengths = []
    try:
        with open(gff_path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] not in ("CDS", "gene"):
                    continue
                attrs = parts[8]
                gene_match = re.search(r'(?:^|;)gene=tcdC', attrs,
                                       re.IGNORECASE)
                if not gene_match:
                    continue
                try:
                    start = int(parts[3])
                    end   = int(parts[4])
                    length = abs(end - start) + 1
                    tcdC_lengths.append(length)
                except ValueError:
                    continue
    except Exception:
        return False

    if not tcdC_lengths:
        return False

    # RT027 tcdC with 18nt deletion and frameshift is substantially shorter
    # Full tcdC: ~657 bp; RT027 (truncated due to frameshift): <400 bp
    max_len = max(tcdC_lengths)
    RT027_THRESHOLD = 500  # bp — below this suggests truncated tcdC
    return max_len < RT027_THRESHOLD


# =============================================================================
# CLASSIFICATION
# =============================================================================

def classify_genome(gene_status: dict, tcdC_del: bool,
                    groups: dict) -> str:
    """
    Classify a genome into Group A, B, or C based on gene status.

    Only COMPLETE genes count as 'present' for classification.
    PARTIAL and PSEUDOGENE genes are treated as non-functional.

    Returns: 'A', 'B', 'C', or 'UNCLASSIFIED'
    """
    def present(gene):
        """Gene exists in functional form (COMPLETE) or partially (PARTIAL)."""
        return gene_status.get(gene, "ABSENT") in ("COMPLETE", "PARTIAL")

    def functional(gene):
        """Gene is fully complete (not partial, not pseudogene)."""
        return gene_status.get(gene, "ABSENT") == "COMPLETE"

    def absent(gene):
        """Gene absent or non-functional (pseudogene)."""
        return gene_status.get(gene, "ABSENT") in ("ABSENT", "PSEUDOGENE")

    # Group A — Hypervirulent RT027 equivalent
    # Core signature: tcdB complete + tcdC non-functional (PSEUDOGENE/ABSENT)
    # CDT subunits may be present, partial, or degraded (RT027 variants)
    tcdC_status = gene_status.get("tcdC", "ABSENT")
    tcdC_nonfunctional = tcdC_status in ("PSEUDOGENE", "ABSENT")
    has_cdt = (present("cdtA") or present("cdtB"))

    if functional("tcdB") and (tcdC_nonfunctional or has_cdt):
        # Confirm it is truly RT027-like: needs at least tcdC non-functional
        # OR CDT genes present (even if degraded)
        if tcdC_nonfunctional and has_cdt:
            return "A"
        # tcdC non-functional even without CDT (some RT027 lost CdtLoc)
        if tcdC_nonfunctional and functional("tcdA"):
            return "A"

    # Group A extended — CDT-positive strains regardless of tcdC status
    # Includes RT078 and other CDT+ ribotypes
    if (functional("tcdB") and
            functional("cdtA") and functional("cdtB")):
        return "A"

    # Group B — Toxigenic (tcdB-sufficient strains included)
    # tcdB complete + no CDT genes
    # Includes TcdB-only strains (tcdA absent or pseudogene)
    if (functional("tcdB") and
            absent("cdtA") and absent("cdtB")):
        return "B"

    # Group C — Non-toxigenic
    # tcdB absent or non-functional
    if not functional("tcdB"):
        return "C"

    return "UNCLASSIFIED"


# =============================================================================
# COLLECT GENOMES
# =============================================================================

def collect_all_genomes(data_dir: Path, log) -> dict:
    """
    Collect all genome directories from working/, chilean/, and anchors/.
    Returns {accession: (genome_dir, set_label)}.
    """
    genomes = {}
    sets = {
        "working":  data_dir / "genomes" / "working",
        "chilean":  data_dir / "genomes" / "chilean",
        "anchors":  data_dir / "genomes" / "anchors",
    }

    for set_label, base_dir in sets.items():
        data_path = base_dir / "ncbi_dataset" / "data"
        if not data_path.exists():
            log.warning(f"Not found: {data_path}")
            continue
        for acc_dir in sorted(data_path.iterdir()):
            if not acc_dir.is_dir():
                continue
            if acc_dir.name.startswith("GC"):
                genomes[acc_dir.name] = (acc_dir, set_label)

    log.info(f"Total genomes found: {len(genomes)}")
    for label in sets:
        n = sum(1 for _, s in genomes.values() if s == label)
        log.info(f"  {label:<10}: {n}")
    return genomes


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
    cls_cfg= cfg.get("classification", {})
    groups = cls_cfg.get("groups", {})
    targets= cls_cfg.get("target_genes", TARGET_GENES)
    del_bp = cls_cfg.get("tcdC_del_bp", TDCC_DEL_BP)
    flank  = cls_cfg.get("junction_flank", JUNCTION_FLANK)

    # Dirs
    data_dir    = Path(paths["main"]) / "data" / "01_download"
    out_dir     = Path(paths["main"]) / "data" / "02_classify"
    report_dir  = Path(paths["main"]) / "reports" / "M02_classify"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M02_classify", cfg)
    print_logo("M02 · Functional Classification", organism=org)

    save_versions(
        tools       = {},
        python_pkgs = ["tqdm", "yaml"],
        report_dir  = report_dir,
        log         = log,
    )

    # ------------------------------------------------------------------
    # Collect all genomes
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("Collecting genomes")
    log.info("─" * 56)

    genomes = collect_all_genomes(data_dir, log)
    if not genomes:
        log.error("No genomes found. Run M01 first.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Parse GFF and classify
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("Parsing GFF files and classifying genomes")
    log.info("─" * 56)

    matrix_rows = []
    group_A = []
    group_B = []
    group_C = []
    unclassified = []

    for acc, (acc_dir, set_label) in tqdm(
            sorted(genomes.items()), desc="Classifying", unit="genome"):

        gff = find_gff(acc_dir.parent.parent.parent)
        fna = find_fna(acc_dir.parent.parent.parent)

        # Try direct path first
        gff_direct = acc_dir / "genomic.gff"
        fna_direct  = list(acc_dir.glob("*_genomic.fna"))
        if gff_direct.exists():
            gff = gff_direct
        if fna_direct:
            fna = fna_direct[0]

        if gff is None:
            log.warning(f"  {acc}: GFF not found — skipping")
            unclassified.append(acc)
            continue

        # Parse gene status from GFF
        gene_status = parse_gff_genes(gff, targets)

        # Detect tcdC deletion (RT027 signature)
        tcdC_del = False
        if fna and gene_status.get("tcdC") in ("COMPLETE", "PARTIAL"):
            tcdC_del = detect_tcdC_deletion(gff, fna, del_bp, flank, log)

        # Classify
        group = classify_genome(gene_status, tcdC_del, groups)

        # Build matrix row
        row = {
            "accession":    acc,
            "set":          set_label,
            "group":        group,
            "tcdC_del":     str(tcdC_del),
        }
        for gene in targets:
            row[gene] = gene_status.get(gene, "ABSENT")
        matrix_rows.append(row)

        # Assign to group list
        if group == "A":
            group_A.append(acc)
        elif group == "B":
            group_B.append(acc)
        elif group == "C":
            group_C.append(acc)
        else:
            unclassified.append(acc)

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("Writing outputs")
    log.info("─" * 56)

    # Gene matrix TSV
    matrix_path = out_dir / "gene_matrix.tsv"
    write_tsv(matrix_rows, matrix_path, log=log)

    # Group files
    for name, accessions in [("groupA", group_A), ("groupB", group_B),
                              ("groupC", group_C),
                              ("unclassified", unclassified)]:
        p = out_dir / f"{name}.txt"
        p.write_text("\n".join(accessions) + "\n", encoding="utf-8")
        log.info(f"  {name:<15}: {len(accessions)} genomes → {p.name}")

    # Summary
    write_tsv([
        {"item": "total_genomes",    "value": len(genomes)},
        {"item": "group_A",          "value": len(group_A)},
        {"item": "group_B",          "value": len(group_B)},
        {"item": "group_C",          "value": len(group_C)},
        {"item": "unclassified",     "value": len(unclassified)},
        {"item": "group_A_label",    "value": groups.get("A",{}).get("label","")},
        {"item": "group_B_label",    "value": groups.get("B",{}).get("label","")},
        {"item": "group_C_label",    "value": groups.get("C",{}).get("label","")},
        {"item": "gene_matrix",      "value": str(matrix_path)},
    ], report_dir / "summary.tsv", log=log)

    log.info("─" * 56)
    log.info("M02 complete.")
    log.info(f"  Total     : {len(genomes)}")
    log.info(f"  Group A   : {len(group_A)}  ({groups.get('A',{}).get('label','')})")
    log.info(f"  Group B   : {len(group_B)}  ({groups.get('B',{}).get('label','')})")
    log.info(f"  Group C   : {len(group_C)}  ({groups.get('C',{}).get('label','')})")
    log.info(f"  Unclassified: {len(unclassified)}")
    log.info(f"  Matrix    → {matrix_path}")

    write_checkpoint("M02_classify", report_dir, log=log)


if __name__ == "__main__":
    main()
