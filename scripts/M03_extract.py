#!/usr/bin/env python3
"""
M03_extract.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M03 · CDS Sequence Extraction

PURPOSE
    Extract CDS sequences for target genes from GFF + FASTA files.
    Only extracts genes marked as COMPLETE in M02 gene matrix.
    PARTIAL and PSEUDOGENE sequences are excluded.

    For each target gene, produces per-group FASTA files:
      tcdA_groupA.fasta, tcdA_groupB.fasta, etc.

    Special handling for tcdC:
      - groupA → tcdC_groupA.fasta  (RT027: PSEUDOGENE, junction region)
      - groupB → tcdC_groupB.fasta  (RT012: COMPLETE, WT sequence)

    The extracted sequences are the CDS (mRNA-equivalent in bacteria).
    UTR flanks (configurable) are added to approximate full transcript.

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M03_extract.py --config ~/sherlock/config.yaml

OUTPUT
    main/data/03_extract/{gene}_group{X}.fasta
    main/data/03_extract/tcdC_junction_groupA.fasta   (RT027 tcdC region)
    main/logs/M03_extract.log
    main/reports/M03_extract/summary.tsv
    main/reports/M03_extract/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG  = "~/sherlock/config.yaml"
THREADS         = 8    # change to 16 for maximum
UTR_FLANK_NT    = 50   # nt added upstream and downstream of CDS

# Target genes and which groups to extract for each
# Format: {gene: [groups]}
GENE_GROUP_MAP = {
    "tcdA":  ["groupA", "groupB"],
    "tcdB":  ["groupA", "groupB"],
    "tcdC":  ["groupB"],              # WT tcdC from Group B only
    "cdtA":  ["groupA"],
    "cdtB":  ["groupA"],
    "tpiA":  ["groupA", "groupB", "groupC"],
    "rpoB":  ["groupA", "groupB", "groupC"],
}

# Alternative search terms for genes not annotated with gene= attribute
# Format: {gene_name: [product_keywords, ...]}
GENE_PRODUCT_ALIASES = {
    "rpoB": ["DNA-directed RNA polymerase subunit beta",
             "RNA polymerase subunit beta"],
}

# Feature types to accept per gene (default: CDS)
GENE_FEATURE_TYPES = {

}

# Alternative search terms for genes not annotated with gene= attribute
# Format: {gene_name: [product_keywords, ...]}
GENE_PRODUCT_ALIASES = {
    "rpoB": ["DNA-directed RNA polymerase subunit beta",
             "RNA polymerase subunit beta"],
}

# Feature types to accept per gene (default: CDS)
GENE_FEATURE_TYPES = {

}

# tcdC junction extraction for RT027 (Group A)
# Extract region around where tcdC exists (even if pseudogene)
EXTRACT_TDCC_JUNCTION = True
JUNCTION_FLANK_NT     = 200  # nt on each side of tcdC locus

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import re
import sys
from pathlib import Path

from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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
# GFF PARSING
# =============================================================================

def parse_gene_coords(gff_path: Path, gene_name: str,
                      complete_only: bool = True) -> list:
    """
    Extract coordinates for a specific gene from GFF.
    Searches by gene= attribute first, then product= aliases.

    Returns list of:
        (contig, start, end, strand, is_partial, is_pseudo)
    Coordinates are 1-based (GFF format).
    """
    # Get accepted feature types for this gene
    accepted_types = GENE_FEATURE_TYPES.get(gene_name, ["CDS", "gene", "pseudogene"])
    # Get product aliases for fallback matching
    aliases = GENE_PRODUCT_ALIASES.get(gene_name, [])

    results = []
    try:
        with open(gff_path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                if feature_type not in accepted_types:
                    continue

                attrs = parts[8]

                # Try gene= attribute first
                matched = False
                gene_match = re.search(r'(?:^|;)gene=([^;]+)', attrs)
                if gene_match:
                    if gene_match.group(1).strip().lower() == gene_name.lower():
                        matched = True

                # Try product= aliases as fallback
                if not matched and aliases:
                    product_match = re.search(r'product=([^;]+)', attrs)
                    if product_match:
                        product = product_match.group(1).strip().lower()
                        if any(alias.lower() in product for alias in aliases):
                            matched = True

                if not matched:
                    continue

                is_partial = "partial=true" in attrs.lower()
                is_pseudo  = ("pseudo=true" in attrs.lower() or
                              feature_type == "pseudogene" or
                              "pseudogene" in attrs.lower())

                if complete_only and (is_partial or is_pseudo):
                    continue

                try:
                    start  = int(parts[3])
                    end    = int(parts[4])
                    strand = parts[6]
                    contig = parts[0]
                    results.append((contig, start, end, strand,
                                    is_partial, is_pseudo))
                except ValueError:
                    continue

    except Exception:
        pass

    return results


def extract_sequence(fna_path: Path, contig: str,
                     start: int, end: int, strand: str,
                     flank: int = 0) -> Seq | None:
    """
    Extract sequence from FASTA by contig name and coordinates.
    Applies strand complement if needed.
    Adds UTR flanks (clamped to contig boundaries).

    GFF coords are 1-based inclusive → convert to 0-based for Python.
    """
    try:
        genome = SeqIO.to_dict(SeqIO.parse(fna_path, "fasta"))
    except Exception:
        return None

    if contig not in genome:
        return None

    seq    = genome[contig].seq
    length = len(seq)

    # Apply flanks (0-based)
    start0 = max(0, start - 1 - flank)
    end0   = min(length, end + flank)

    subseq = seq[start0:end0]

    if strand == "-":
        subseq = subseq.reverse_complement()

    return subseq


# =============================================================================
# LOAD CLASSIFICATION DATA
# =============================================================================

def load_group_accessions(classify_dir: Path) -> dict:
    """
    Load group membership from M02 output files.
    Returns {group_name: set of accessions}.
    """
    groups = {}
    for group in ("groupA", "groupB", "groupC"):
        f = classify_dir / f"{group}.txt"
        if f.exists():
            accs = {l.strip() for l in f.read_text().splitlines()
                    if l.strip()}
            groups[group] = accs
        else:
            groups[group] = set()
    return groups


def load_gene_matrix(classify_dir: Path) -> dict:
    """
    Load gene status matrix from M02.
    Returns {accession: {gene: status}}.
    """
    matrix = {}
    f = classify_dir / "gene_matrix.tsv"
    if not f.exists():
        return matrix
    with open(f, encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            row   = dict(zip(header, parts))
            acc   = row.get("accession", "")
            if acc:
                matrix[acc] = row
    return matrix


# =============================================================================
# FIND GENOME FILES
# =============================================================================

def find_genome_files(data_dir: Path, accession: str) -> tuple:
    """
    Find GFF and FNA files for a given accession.
    Searches working/, chilean/, and anchors/ directories.
    Returns (gff_path, fna_path) or (None, None).
    """
    for subdir in ("working", "chilean", "anchors"):
        base = (data_dir / "01_download" / "genomes" /
                subdir / "ncbi_dataset" / "data" / accession)
        if base.exists():
            gff = base / "genomic.gff"
            fna_files = list(base.glob("*_genomic.fna"))
            if gff.exists() and fna_files:
                return gff, fna_files[0]
    return None, None


# =============================================================================
# MAIN EXTRACTION
# =============================================================================

def extract_gene_for_group(gene: str, group: str,
                           accessions: set,
                           gene_matrix: dict,
                           data_dir: Path,
                           out_dir: Path,
                           flank: int,
                           log) -> int:
    """
    Extract CDS sequences for a gene from all accessions in a group.
    Writes to {gene}_{group}.fasta.
    Returns number of sequences written.
    """
    out_file = out_dir / f"{gene}_{group}.fasta"
    records  = []

    for acc in tqdm(sorted(accessions),
                    desc=f"{gene}/{group}", unit="genome", leave=False):

        # Check gene status in matrix
        row = gene_matrix.get(acc, {})
        status = row.get(gene, "ABSENT")
        if status != "COMPLETE":
            continue  # skip PARTIAL, PSEUDOGENE, ABSENT

        gff, fna = find_genome_files(data_dir, acc)
        if gff is None:
            log.warning(f"  {acc}: files not found")
            continue

        coords = parse_gene_coords(gff, gene, complete_only=True)
        if not coords:
            continue

        # Take the first complete CDS entry
        contig, start, end, strand, _, _ = coords[0]
        seq = extract_sequence(fna, contig, start, end, strand, flank)
        if seq is None or len(seq) < 100:
            continue

        records.append(SeqRecord(
            seq,
            id=f"{acc}_{gene}",
            description=f"group={group} contig={contig} "
                        f"coords={start}-{end} strand={strand}",
        ))

    if records:
        with open(out_file, "w") as f:
            SeqIO.write(records, f, "fasta")
        log.info(f"  {gene}/{group}: {len(records)} sequences → {out_file.name}")
    else:
        log.warning(f"  {gene}/{group}: no sequences extracted")

    return len(records)


def extract_tcdC_junction_groupA(accessions: set,
                                  gene_matrix: dict,
                                  data_dir: Path,
                                  out_dir: Path,
                                  junction_flank: int,
                                  log) -> int:
    """
    Extract tcdC locus region from Group A genomes.
    Even if tcdC is PSEUDOGENE, the surrounding region contains
    the RT027-specific junction sequence used for crRNA design.
    """
    out_file = out_dir / "tcdC_junction_groupA.fasta"
    records  = []

    for acc in tqdm(sorted(accessions),
                    desc="tcdC_junction/groupA", unit="genome", leave=False):

        row    = gene_matrix.get(acc, {})
        status = row.get("tcdC", "ABSENT")

        # Include even PSEUDOGENE — we want the junction sequence
        if status == "ABSENT":
            continue

        gff, fna = find_genome_files(data_dir, acc)
        if gff is None:
            continue

        # Get tcdC coords (including pseudogenes)
        coords = parse_gene_coords(gff, "tcdC", complete_only=False)
        if not coords:
            continue

        contig, start, end, strand, is_partial, is_pseudo = coords[0]
        # Extract with large flanks to capture junction region
        seq = extract_sequence(fna, contig, start, end, strand,
                               junction_flank)
        if seq is None or len(seq) < 50:
            continue

        records.append(SeqRecord(
            seq,
            id=f"{acc}_tcdC_junction",
            description=f"group=A pseudo={is_pseudo} partial={is_partial} "
                        f"contig={contig} coords={start}-{end}",
        ))

    if records:
        with open(out_file, "w") as f:
            SeqIO.write(records, f, "fasta")
        log.info(f"  tcdC_junction/groupA: {len(records)} sequences → {out_file.name}")

    return len(records)


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
    args    = parse_args()
    cfg     = load_config(str(Path(args.config).expanduser()))
    org     = cfg["organism"]["display"]
    paths   = cfg["paths"]
    ext_cfg = cfg.get("extraction", {})
    flank   = ext_cfg.get("utr_flank_nt", UTR_FLANK_NT)
    targets = ext_cfg.get("targets", list(GENE_GROUP_MAP.keys()))

    # Dirs
    main_dir     = Path(paths["main"])
    data_dir     = main_dir / "data"
    classify_dir = data_dir / "02_classify"
    out_dir      = data_dir / "03_extract"
    report_dir   = main_dir / "reports" / "M03_extract"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M03_extract", cfg)
    print_logo("M03 · CDS Extraction", organism=org)

    save_versions(
        tools       = {},
        python_pkgs = ["Bio", "tqdm"],
        report_dir  = report_dir,
        log         = log,
    )

    # ------------------------------------------------------------------
    # Load classification data
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("Loading classification data from M02")
    log.info("─" * 56)

    groups     = load_group_accessions(classify_dir)
    gene_matrix= load_gene_matrix(classify_dir)

    for grp, accs in groups.items():
        log.info(f"  {grp}: {len(accs)} genomes")

    if not gene_matrix:
        log.error("Gene matrix not found. Run M02 first.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Extract CDS per gene per group
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info(f"Extracting CDS sequences  (flank={flank}nt)")
    log.info("─" * 56)

    summary_rows = []

    for gene in targets:
        gene_groups = GENE_GROUP_MAP.get(gene, [])
        for group in gene_groups:
            if group not in groups:
                continue
            n = extract_gene_for_group(
                gene       = gene,
                group      = group,
                accessions = groups[group],
                gene_matrix= gene_matrix,
                data_dir   = data_dir,
                out_dir    = out_dir,
                flank      = flank,
                log        = log,
            )
            summary_rows.append({
                "gene":  gene,
                "group": group,
                "n_seqs": n,
                "file":  f"{gene}_{group}.fasta",
            })

    # ------------------------------------------------------------------
    # tcdC junction for Group A (RT027 signature)
    # ------------------------------------------------------------------
    if EXTRACT_TDCC_JUNCTION:
        log.info("─" * 56)
        log.info("Extracting tcdC junction region (Group A / RT027)")
        log.info("─" * 56)

        n = extract_tcdC_junction_groupA(
            accessions   = groups.get("groupA", set()),
            gene_matrix  = gene_matrix,
            data_dir     = data_dir,
            out_dir      = out_dir,
            junction_flank= JUNCTION_FLANK_NT,
            log          = log,
        )
        summary_rows.append({
            "gene":  "tcdC_junction",
            "group": "groupA",
            "n_seqs": n,
            "file":  "tcdC_junction_groupA.fasta",
        })

    # ------------------------------------------------------------------
    # Reports
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("Writing reports")
    log.info("─" * 56)

    write_tsv(summary_rows, report_dir / "summary.tsv", log=log)

    # List all output files
    fasta_files = sorted(out_dir.glob("*.fasta"))
    log.info(f"Output files: {len(fasta_files)}")
    for f in fasta_files:
        n = count_seqs(f)
        log.info(f"  {f.name:<40} {n} sequences")

    log.info("─" * 56)
    log.info("M03 complete.")
    log.info(f"  Sequences → {out_dir}")

    write_checkpoint("M03_extract", report_dir, log=log)


if __name__ == "__main__":
    main()
