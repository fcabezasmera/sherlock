#!/usr/bin/env python3
"""
M01_download.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M01 · Genome Download

PURPOSE
    1. Download genome metadata from NCBI (all C. difficile assemblies)
    2. Filter by assembly level (Complete Genome + Chromosome) + GCF
    3. Download filtered genomes (FASTA + GFF3)
    4. Download anchor strains (Chile-relevant references)
    5. Download Chilean strains (all 12 available in NCBI)

    Gene integrity is verified from GFF flags (partial=true, pseudo=true)
    in M02. No individual gene downloads needed.

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M01_download.py --config ~/sherlock/config.yaml

    # Re-run flags:
    --skip-metadata   reuse existing ncbi_dataset.tsv
    --skip-genomes    filter only, no downloads

OUTPUT
    main/data/01_download/ncbi_dataset.tsv
    main/data/01_download/filtered_genomes.tsv
    main/data/01_download/genomes/working/    global set
    main/data/01_download/genomes/anchors/    Chile references
    main/data/01_download/genomes/chilean/    Chilean clinical strains
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
BATCH_SIZE      = 50    # accessions per datasets call
NCBI_SLEEP      = 1.5   # seconds between batches
THREADS         = 8     # default (change to 16 for maximum)

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

def download_metadata(taxon, fields, out_tsv, datasets_exe, log):
    """Download genome metadata from NCBI and save as TSV."""
    log.info(f"Querying NCBI: {taxon}")

    # Get JSON lines
    cmd = [datasets_exe, "summary", "genome",
           "taxon", taxon, "--as-json-lines"]
    log.info("CMD: " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)

    if r.returncode != 0:
        log.error(f"datasets summary failed:\n{r.stderr}")
        raise RuntimeError("datasets summary failed")

    lines = [l for l in r.stdout.splitlines() if l.strip()]
    log.info(f"JSON lines received: {len(lines)}")

    # Save raw JSON (keep for debugging)
    jsonl = out_tsv.parent / "ncbi_raw.jsonl"
    jsonl.write_text(r.stdout, encoding="utf-8")

    # Convert to TSV
    cmd2 = ["dataformat", "tsv", "genome",
            "--inputfile", str(jsonl),
            "--fields", fields]
    log.info("CMD: " + " ".join(cmd2))
    r2 = subprocess.run(cmd2, capture_output=True, text=True)

    if r2.returncode != 0:
        log.error(f"dataformat failed:\n{r2.stderr}")
        raise RuntimeError("dataformat failed")

    out_tsv.write_text(r2.stdout, encoding="utf-8")
    n = len(r2.stdout.splitlines()) - 1
    log.info(f"Metadata → {out_tsv}  ({n} assemblies)")
    return n


# =============================================================================
# STEP 2 — Filter metadata
# =============================================================================

def filter_metadata(tsv_path, assembly_levels, cols, log):
    """
    Filter metadata by assembly level and GCF requirement.
    Returns (all_rows, filtered_rows).
    """
    col_acc   = cols.get("accession", "Assembly Accession")
    col_level = cols.get("level",     "Assembly Level")

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

    # Count by level for report
    from collections import Counter
    level_counts = Counter(r.get(col_level, "").strip() for r in all_rows)
    for level, count in sorted(level_counts.items(), key=lambda x: -x[1]):
        log.info(f"  {level:<25}: {count}")

    filtered = []
    n_level = n_gcf = 0

    for row in all_rows:
        if row.get(col_level, "").strip() not in assembly_levels:
            n_level += 1
            continue
        acc = row.get(col_acc, "").strip()
        if not acc.startswith("GCF"):
            n_gcf += 1
            continue
        filtered.append(row)

    log.info(f"Filter results:")
    log.info(f"  Removed (level not in {assembly_levels}): {n_level}")
    log.info(f"  Removed (not GCF)                      : {n_gcf}")
    log.info(f"  Passing filter                          : {len(filtered)}")
    return all_rows, filtered


# =============================================================================
# DOWNLOAD HELPER
# =============================================================================

def download_genomes(accessions, out_dir, datasets_exe,
                     retries, label, log):
    """
    Download FASTA + GFF3 for a list of accessions in batches.
    Returns {accession: 'OK' | 'FAILED'}.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    status  = {}
    batches = [accessions[i:i+BATCH_SIZE]
               for i in range(0, len(accessions), BATCH_SIZE)]

    log.info(f"Downloading {len(accessions)} genomes [{label}] "
             f"in {len(batches)} batch(es)")

    with tqdm(total=len(accessions), desc=label, unit="genome") as pbar:
        for batch in batches:
            zip_file = out_dir / "_tmp.zip"
            for attempt in range(1, retries + 1):
                try:
                    cmd = [datasets_exe, "download", "genome",
                           "accession", *batch,
                           "--include", "genome,gff3",
                           "--filename", str(zip_file)]
                    run_cmd(cmd, log, capture=False, check=True)
                    run_cmd(["unzip", "-o", str(zip_file),
                             "-d", str(out_dir)],
                            log, capture=True, check=False)
                    zip_file.unlink(missing_ok=True)
                    for acc in batch:
                        status[acc] = "OK"
                    break
                except Exception as e:
                    zip_file.unlink(missing_ok=True)
                    log.warning(f"  Attempt {attempt}/{retries}: {e}")
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
                   help="Filter only, skip downloads")
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

    # Assembly levels to include
    assembly_levels = dl_cfg.get("assembly_levels",
                                  ["Complete Genome", "Chromosome"])

    # Chilean strains and anchors
    chilean = dl_cfg.get("chilean_strains", {})
    anchors = dl_cfg.get("anchors", {})
    retries = dl_cfg.get("retries", 3)

    # Dirs
    data_dir    = Path(paths["main"]) / "data" / "01_download"
    working_dir = data_dir / "genomes" / "working"
    anchors_dir = data_dir / "genomes" / "anchors"
    chilean_dir = data_dir / "genomes" / "chilean"
    meta_tsv    = data_dir / "ncbi_dataset.tsv"
    filt_tsv    = data_dir / "filtered_genomes.tsv"
    report_dir  = Path(paths["main"]) / "reports" / "M01_download"

    for d in (data_dir, working_dir, anchors_dir, chilean_dir, report_dir):
        d.mkdir(parents=True, exist_ok=True)

    log = get_logger("M01_download", cfg)
    print_logo("M01 · Genome Download", organism=org)

    save_versions(
        tools       = {"datasets": tools["datasets"]},
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
            taxon        = TAXON,
            fields       = fields,
            out_tsv      = meta_tsv,
            datasets_exe = tools["datasets"],
            log          = log,
        )

    # ------------------------------------------------------------------
    # STEP 2 — Filter
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 2 · Filter metadata")
    log.info(f"         Levels: {assembly_levels}")
    log.info("─" * 56)

    all_rows, filtered_rows = filter_metadata(
        meta_tsv, assembly_levels, cols, log)
    write_tsv(filtered_rows, filt_tsv, log=log)

    col_acc    = cols.get("accession", "Assembly Accession")
    accessions = [r[col_acc] for r in filtered_rows
                  if r.get(col_acc, "").startswith("GCF")]
    log.info(f"Accessions to download: {len(accessions)}")

    # ------------------------------------------------------------------
    # STEP 3 — Global genomes (working/)
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 3 · Download global genomes → working/")
    log.info("─" * 56)

    if args.skip_genomes:
        log.info("--skip-genomes: skipping")
        genome_status = {acc: "SKIPPED" for acc in accessions}
        n_ok = n_fail = 0
    else:
        genome_status = download_genomes(
            accessions, working_dir, tools["datasets"],
            retries, "working", log)
        n_ok   = sum(1 for s in genome_status.values() if s == "OK")
        n_fail = sum(1 for s in genome_status.values() if s == "FAILED")

    # ------------------------------------------------------------------
    # STEP 4 — Chilean strains (chilean/)
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 4 · Download Chilean strains → chilean/")
    log.info(f"         {len(chilean)} strains (all available in NCBI)")
    log.info("─" * 56)

    if args.skip_genomes:
        log.info("--skip-genomes: skipping")
        chilean_status = {v: "SKIPPED" for v in chilean.values()}
    else:
        chilean_status = download_genomes(
            list(chilean.values()), chilean_dir, tools["datasets"],
            retries, "chilean", log)
        for name, acc in chilean.items():
            log.info(f"  {name:<10} ({acc}): {chilean_status.get(acc,'?')}")

    # ------------------------------------------------------------------
    # STEP 5 — Anchor strains (anchors/)
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 5 · Download anchor strains → anchors/")
    log.info("─" * 56)

    if args.skip_genomes:
        log.info("--skip-genomes: skipping")
        anchor_status = {v: "SKIPPED" for v in anchors.values()}
    else:
        anchor_status = download_genomes(
            list(anchors.values()), anchors_dir, tools["datasets"],
            retries, "anchors", log)
        for name, acc in anchors.items():
            log.info(f"  {name:<20} ({acc}): {anchor_status.get(acc,'?')}")

    # ------------------------------------------------------------------
    # Reports
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 6 · Writing reports")
    log.info("─" * 56)

    ch_ok   = sum(1 for s in chilean_status.values() if s == "OK")
    ch_fail = sum(1 for s in chilean_status.values() if s == "FAILED")

    write_tsv([
        {"item": "taxon",            "value": TAXON},
        {"item": "assembly_levels",  "value": str(assembly_levels)},
        {"item": "total_ncbi",       "value": n_total},
        {"item": "passing_filter",   "value": len(filtered_rows)},
        {"item": "global_ok",        "value": n_ok},
        {"item": "global_failed",    "value": n_fail},
        {"item": "chilean_total",    "value": len(chilean)},
        {"item": "chilean_ok",       "value": ch_ok},
        {"item": "chilean_failed",   "value": ch_fail},
        {"item": "anchors_total",    "value": len(anchors)},
        {"item": "total_working",    "value": n_ok + ch_ok},
        {"item": "metadata_tsv",     "value": str(meta_tsv)},
        {"item": "filtered_tsv",     "value": str(filt_tsv)},
        {"item": "working_dir",      "value": str(working_dir)},
        {"item": "chilean_dir",      "value": str(chilean_dir)},
        {"item": "anchors_dir",      "value": str(anchors_dir)},
    ], report_dir / "summary.tsv", log=log)

    write_tsv(
        [{"accession": acc, "status": st, "set": "global"}
         for acc, st in genome_status.items()] +
        [{"accession": acc, "status": st, "set": "chilean"}
         for acc, st in chilean_status.items()] +
        [{"accession": acc, "status": st, "set": "anchor"}
         for acc, st in anchor_status.items()],
        report_dir / "genome_status.tsv", log=log,
    )

    log.info("─" * 56)
    log.info("M01 complete.")
    log.info(f"  Metadata   → {meta_tsv}")
    log.info(f"  Filtered   → {filt_tsv}  ({len(filtered_rows)} genomes)")
    log.info(f"  Working    → {working_dir}  ({n_ok} global)")
    log.info(f"  Chilean    → {chilean_dir}  ({ch_ok}/{len(chilean)} strains)")
    log.info(f"  Anchors    → {anchors_dir}  ({len(anchors)})")

    write_checkpoint("M01_download", report_dir, log=log)


if __name__ == "__main__":
    main()
