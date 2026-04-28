#!/usr/bin/env python3
"""
M01_download.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M01 · Genome Download

PURPOSE
    1. Download full C. difficile genome metadata from NCBI
    2. Filter by quality criteria (Complete Genome, GCF, scaffolds <= 2)
    3. Download filtered genomes (FASTA + GFF3)
    4. Download anchor strains (Chile-relevant references)
    5. Download reference gene sequences via efetch

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M01_download.py --config ~/sherlock/config.yaml

    # Optional flags to skip completed steps:
    --skip-metadata    reuse existing ncbi_dataset.tsv
    --skip-genomes     metadata + filter only, no downloads
    --skip-ref-genes   skip individual gene downloads

OUTPUT
    main/data/01_download/ncbi_dataset.tsv
    main/data/01_download/filtered_genomes.tsv
    main/data/01_download/genomes/working/
    main/data/01_download/genomes/anchors/
    main/data/01_download/ref_genes/
    main/logs/M01_download.log
    main/reports/M01_download/summary.tsv
    main/reports/M01_download/genome_status.tsv
    main/reports/M01_download/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG  = "~/sherlock/config.yaml"
TAXON           = "Clostridioides difficile"  # change for different organism
ASSEMBLY_LEVEL  = "Complete Genome"
MAX_SCAFFOLDS   = 2
REQUIRE_GCF     = True
BATCH_SIZE      = 50    # accessions per datasets download call
NCBI_SLEEP      = 1.5   # seconds between batches (NCBI rate limit)
EFETCH_SLEEP    = 0.5   # seconds between efetch calls

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import subprocess
import sys
import time
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
    run_cmd,
)

# =============================================================================
# STEP 1 — Download metadata
# =============================================================================

def download_metadata(taxon: str,
                      fields: str,
                      out_tsv: Path,
                      datasets_exe: str,
                      log) -> int:
    """
    Download genome metadata for a taxon from NCBI and save as TSV.
    Returns total number of assemblies found.
    """
    log.info(f"Querying NCBI for: {taxon}")

    # Get JSON lines
    cmd_summary = [
        datasets_exe, "summary", "genome",
        "taxon", taxon,
        "--as-json-lines",
    ]
    log.info("CMD: " + " ".join(cmd_summary))
    r = subprocess.run(cmd_summary, capture_output=True, text=True)

    if r.returncode != 0:
        log.error(f"datasets summary failed:\n{r.stderr}")
        raise RuntimeError("datasets summary failed")

    lines = [l for l in r.stdout.splitlines() if l.strip()]
    log.info(f"JSON lines received: {len(lines)}")

    # Save raw JSON lines (keep for debugging)
    jsonl = out_tsv.parent / "ncbi_raw.jsonl"
    jsonl.write_text(r.stdout, encoding="utf-8")

    # Convert to TSV with selected fields
    cmd_tsv = [
        "dataformat", "tsv", "genome",
        "--inputfile", str(jsonl),
        "--fields", fields,
    ]
    log.info("CMD: " + " ".join(cmd_tsv))
    r2 = subprocess.run(cmd_tsv, capture_output=True, text=True)

    if r2.returncode != 0:
        log.error(f"dataformat failed:\n{r2.stderr}")
        raise RuntimeError("dataformat conversion failed")

    out_tsv.write_text(r2.stdout, encoding="utf-8")
    n = len(r2.stdout.splitlines()) - 1  # subtract header
    log.info(f"Metadata saved → {out_tsv}  ({n} assemblies)")
    return n


# =============================================================================
# STEP 2 — Filter metadata
# =============================================================================

def filter_metadata(tsv_path: Path, cols: dict, log) -> tuple:
    """
    Parse and filter metadata TSV.
    Returns (all_rows, filtered_rows) as lists of dicts.

    Filters applied:
      - Assembly Level == Complete Genome
      - Accession starts with GCF (RefSeq-paired)
      - Number of scaffolds <= MAX_SCAFFOLDS
    """
    col_acc      = cols.get("accession", "Assembly Accession")
    col_level    = cols.get("level",     "Assembly Level")
    col_scaffolds= cols.get("scaffolds", "Assembly Stats Number of Scaffolds")

    all_rows = []
    with open(tsv_path, encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            if not line.strip():
                continue
            values = line.strip().split("\t")
            while len(values) < len(header):
                values.append("")
            all_rows.append(dict(zip(header, values)))

    log.info(f"Total assemblies in NCBI: {len(all_rows)}")

    filtered = []
    n_level = n_gcf = n_scaff = 0

    for row in all_rows:
        # Assembly level filter
        if row.get(col_level, "").strip() != ASSEMBLY_LEVEL:
            n_level += 1
            continue
        # GCF filter
        acc = row.get(col_acc, "").strip()
        if REQUIRE_GCF and not acc.startswith("GCF"):
            n_gcf += 1
            continue
        # Scaffold filter
        try:
            n_sc = int(row.get(col_scaffolds, "999").strip() or "999")
        except ValueError:
            n_sc = 999
        if n_sc > MAX_SCAFFOLDS:
            n_scaff += 1
            continue
        filtered.append(row)

    log.info("Filter results:")
    log.info(f"  Removed (not Complete Genome)  : {n_level}")
    log.info(f"  Removed (not GCF)              : {n_gcf}")
    log.info(f"  Removed (scaffolds > {MAX_SCAFFOLDS})         : {n_scaff}")
    log.info(f"  Passing filter                 : {len(filtered)}")

    return all_rows, filtered


# =============================================================================
# STEPS 3 & 4 — Download genomes
# =============================================================================

def download_genomes(accessions: list,
                     out_dir: Path,
                     datasets_exe: str,
                     retries: int,
                     label: str,
                     log) -> dict:
    """
    Download genome FASTA + GFF3 for a list of accessions in batches.
    Returns {accession: "OK" | "FAILED"} dict.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    status  = {}
    batches = [accessions[i:i+BATCH_SIZE]
               for i in range(0, len(accessions), BATCH_SIZE)]

    log.info(f"Downloading {len(accessions)} genomes [{label}] "
             f"in {len(batches)} batch(es) of {BATCH_SIZE}")

    with tqdm(total=len(accessions),
              desc=label, unit="genome") as pbar:

        for batch in batches:
            zip_file = out_dir / "_batch_tmp.zip"

            for attempt in range(1, retries + 1):
                try:
                    cmd = [
                        datasets_exe, "download", "genome",
                        "accession", *batch,
                        "--include", "genome,gff3",
                        "--filename", str(zip_file),
                    ]
                    run_cmd(cmd, log, capture=False, check=True)

                    # Unzip
                    run_cmd(
                        ["unzip", "-o", str(zip_file), "-d", str(out_dir)],
                        log, capture=True, check=False,
                    )
                    zip_file.unlink(missing_ok=True)

                    for acc in batch:
                        status[acc] = "OK"
                    break

                except Exception as e:
                    zip_file.unlink(missing_ok=True)
                    log.warning(f"  Attempt {attempt}/{retries} failed: {e}")
                    if attempt == retries:
                        for acc in batch:
                            status[acc] = "FAILED"
                    else:
                        time.sleep(5 * attempt)

            pbar.update(len(batch))
            time.sleep(NCBI_SLEEP)

    n_ok   = sum(1 for s in status.values() if s == "OK")
    n_fail = sum(1 for s in status.values() if s == "FAILED")
    log.info(f"  [{label}] OK={n_ok}  FAILED={n_fail}")
    return status


# =============================================================================
# STEP 5 — Download reference genes
# =============================================================================

def download_ref_genes(ref_genes: dict,
                       out_dir: Path,
                       efetch_exe: str,
                       log) -> None:
    """
    Download individual gene reference sequences via efetch.
    ref_genes: {gene_name: ncbi_accession}
    Output: {gene}_reference.fasta per gene
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Downloading {len(ref_genes)} reference gene sequences")

    for gene, accession in tqdm(ref_genes.items(),
                                desc="Ref genes", unit="gene"):
        out_file = out_dir / f"{gene}_reference.fasta"

        if out_file.exists():
            log.info(f"  {gene}: exists — skipping")
            continue

        try:
            cmd = [
                efetch_exe,
                "-db", "nuccore",
                "-id", accession,
                "-format", "fasta",
            ]
            r = subprocess.run(
                cmd, capture_output=True, text=True, timeout=60)

            if r.returncode == 0 and r.stdout.strip():
                out_file.write_text(r.stdout, encoding="utf-8")
                log.info(f"  {gene} ({accession}) → {out_file.name}")
            else:
                log.warning(f"  {gene} ({accession}): FAILED")

        except Exception as e:
            log.warning(f"  {gene} ({accession}): {e}")

        time.sleep(EFETCH_SLEEP)


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
    p.add_argument("--skip-metadata", action="store_true",
                   help="Reuse existing ncbi_dataset.tsv")
    p.add_argument("--skip-genomes", action="store_true",
                   help="Skip genome downloads (metadata + filter only)")
    p.add_argument("--skip-ref-genes", action="store_true",
                   help="Skip reference gene downloads")
    return p.parse_args()


def main():
    args   = parse_args()
    cfg    = load_config(str(Path(args.config).expanduser()))
    org    = cfg["organism"]["display"]
    paths  = cfg["paths"]
    tools  = cfg["tools"]
    dl_cfg = cfg["download"]
    fields = cfg.get("metadata_fields",
                     "accession,assminfo-name,assminfo-level,"
                     "assmstats-number-of-scaffolds")
    cols   = cfg.get("metadata_cols", {})

    # Dirs
    data_dir    = Path(paths["main"]) / "data" / "01_download"
    working_dir = data_dir / "genomes" / "working"
    anchors_dir = data_dir / "genomes" / "anchors"
    ref_dir     = data_dir / "ref_genes"
    meta_tsv    = data_dir / "ncbi_dataset.tsv"
    filt_tsv    = data_dir / "filtered_genomes.tsv"
    report_dir  = Path(paths["main"]) / "reports" / "M01_download"

    for d in (data_dir, working_dir, anchors_dir, ref_dir, report_dir):
        d.mkdir(parents=True, exist_ok=True)

    log = get_logger("M01_download", cfg)
    print_logo("M01 · Genome Download", organism=org)

    save_versions(
        tools       = {"datasets": tools["datasets"],
                       "efetch":   tools["efetch"]},
        python_pkgs = ["tqdm", "yaml"],
        report_dir  = report_dir,
        log         = log,
    )

    # ------------------------------------------------------------------
    # STEP 1 — Metadata
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 1 · Download metadata from NCBI")
    log.info("─" * 56)

    if args.skip_metadata and meta_tsv.exists():
        log.info(f"--skip-metadata: reusing {meta_tsv}")
        n_total = sum(1 for _ in open(meta_tsv)) - 1
        log.info(f"Existing entries: {n_total}")
    else:
        n_total = download_metadata(
            taxon       = TAXON,
            fields      = fields,
            out_tsv     = meta_tsv,
            datasets_exe= tools["datasets"],
            log         = log,
        )

    # ------------------------------------------------------------------
    # STEP 2 — Filter
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 2 · Filter metadata")
    log.info("─" * 56)

    all_rows, filtered_rows = filter_metadata(meta_tsv, cols, log)
    write_tsv(filtered_rows, filt_tsv, log=log)

    col_acc    = cols.get("accession", "Assembly Accession")
    accessions = [r[col_acc] for r in filtered_rows
                  if r.get(col_acc, "").startswith("GCF")]
    log.info(f"Accessions to download: {len(accessions)}")

    # ------------------------------------------------------------------
    # STEP 3 — Filtered genomes
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 3 · Download filtered genomes")
    log.info("─" * 56)

    if args.skip_genomes:
        log.info("--skip-genomes: skipping")
        genome_status = {acc: "SKIPPED" for acc in accessions}
        n_ok = n_fail = 0
    else:
        genome_status = download_genomes(
            accessions  = accessions,
            out_dir     = working_dir,
            datasets_exe= tools["datasets"],
            retries     = dl_cfg.get("retries", 3),
            label       = "working",
            log         = log,
        )
        n_ok   = sum(1 for s in genome_status.values() if s == "OK")
        n_fail = sum(1 for s in genome_status.values() if s == "FAILED")

    # ------------------------------------------------------------------
    # STEP 4 — Anchors
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 4 · Download anchor strains  (Chile-relevant)")
    log.info("─" * 56)

    anchors = dl_cfg.get("anchors", {})
    if args.skip_genomes:
        log.info("--skip-genomes: skipping anchors")
        anchor_status = {v: "SKIPPED" for v in anchors.values()}
    else:
        anchor_status = download_genomes(
            accessions  = list(anchors.values()),
            out_dir     = anchors_dir,
            datasets_exe= tools["datasets"],
            retries     = dl_cfg.get("retries", 3),
            label       = "anchors",
            log         = log,
        )
        for name, acc in anchors.items():
            log.info(f"  {name} ({acc}): {anchor_status.get(acc, '?')}")

    # ------------------------------------------------------------------
    # STEP 5 — Reference genes
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 5 · Download reference genes")
    log.info("─" * 56)

    ref_genes = dl_cfg.get("ref_genes", {})
    if args.skip_ref_genes:
        log.info("--skip-ref-genes: skipping")
    else:
        download_ref_genes(
            ref_genes  = ref_genes,
            out_dir    = ref_dir,
            efetch_exe = tools["efetch"],
            log        = log,
        )

    # ------------------------------------------------------------------
    # Reports
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 6 · Writing reports")
    log.info("─" * 56)

    write_tsv([
        {"item": "taxon",           "value": TAXON},
        {"item": "assembly_level",  "value": ASSEMBLY_LEVEL},
        {"item": "max_scaffolds",   "value": MAX_SCAFFOLDS},
        {"item": "total_ncbi",      "value": n_total},
        {"item": "passing_filter",  "value": len(filtered_rows)},
        {"item": "genomes_ok",      "value": n_ok},
        {"item": "genomes_failed",  "value": n_fail},
        {"item": "anchors",         "value": len(anchors)},
        {"item": "ref_genes",       "value": len(ref_genes)},
        {"item": "metadata_tsv",    "value": str(meta_tsv)},
        {"item": "filtered_tsv",    "value": str(filt_tsv)},
        {"item": "working_dir",     "value": str(working_dir)},
        {"item": "anchors_dir",     "value": str(anchors_dir)},
        {"item": "ref_genes_dir",   "value": str(ref_dir)},
    ], report_dir / "summary.tsv", log=log)

    write_tsv(
        [{"accession": acc, "status": st}
         for acc, st in {**genome_status, **anchor_status}.items()],
        report_dir / "genome_status.tsv", log=log,
    )

    log.info("─" * 56)
    log.info("M01 complete.")
    log.info(f"  Metadata        → {meta_tsv}")
    log.info(f"  Filtered        → {filt_tsv}  ({len(filtered_rows)} genomes)")
    log.info(f"  Working genomes → {working_dir}")
    log.info(f"  Anchors         → {anchors_dir}")
    log.info(f"  Ref genes       → {ref_dir}")

    write_checkpoint("M01_download", report_dir, log=log)


if __name__ == "__main__":
    main()
