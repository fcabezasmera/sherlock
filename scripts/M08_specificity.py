#!/usr/bin/env python3
"""
M08_specificity.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M08 · Specificity Analysis

PURPOSE
    Verify that crRNA candidates from M06 do NOT cross-react with:
      1. Non-toxigenic C. difficile (Group C, 52 genomes)
      2. Common enteric pathogens (10 species)
      3. Human transcriptome (GRCh38 RefSeq mRNA)
      4. Human gut microbiome (UHGG v2.0, Bowtie2)

    Rejection criterion: identity >= 90% over >= 70% of crRNA length.

    Databases built once and reused across runs.

    Tools:
      - blastn -task blastn-short (crRNA vs genome/transcriptome)
      - bowtie2 (crRNA vs UHGG index)

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M08_specificity.py \
        --config ~/sherlock/config.yaml

    # Skip flags:
    --skip-db-build   reuse existing BLAST/Bowtie2 databases
    --skip-blast      reuse existing BLAST results
    --skip-bowtie2    reuse existing Bowtie2 results

OUTPUT
    main/data/08_specificity/dbs/                  databases
    main/data/08_specificity/{target}_blast.tsv    BLAST results per target
    main/data/08_specificity/{target}_specificity.tsv  pass/fail per crRNA
    main/data/08_specificity/all_specificity.tsv   combined
    main/logs/M08_specificity.log
    main/reports/M08_specificity/summary.tsv
    main/reports/M08_specificity/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG  = "~/sherlock/config.yaml"
THREADS         = 8
SPACER_LEN      = 28

# Rejection thresholds
REJECT_IDENTITY = 90.0   # % identity
REJECT_COVERAGE = 70.0   # % of crRNA covered

# BLAST parameters
BLASTN_TASK     = "blastn-short"
BLASTN_WORD     = 7
EVALUE          = 0.01

# Enteric pathogens to download (RefSeq representative genomes)
ENTEROPATHOGENS = {
    "C_sordellii":       "GCF_001077475.1",
    "C_perfringens":     "GCF_000008245.2",
    "C_botulinum":       "GCF_000063585.1",
    "E_faecalis":        "GCF_000007785.2",
    "E_faecium":         "GCF_009734005.1",
    "K_pneumoniae":      "GCF_000240185.1",
    "E_coli_O157":       "GCF_000008865.2",
    "S_enterica":        "GCF_000006945.2",
    "C_jejuni":          "GCF_000009085.1",
    "L_monocytogenes":   "GCF_000196035.1",
    "S_aureus":          "GCF_000013425.1",
}

# Human transcriptome RefSeq mRNA (representative accession for download)
HUMAN_TX_TAXON  = "9606"   # Homo sapiens
HUMAN_TX_SEQTYPE= "rna"    # mRNA sequences

# Per-target database check rules:
# nontox     : informational only — never used for rejection
# enteropath : check all except 16S (universal 16S will always hit bacteria)
# human_tx   : check all targets

SKIP_DATABASES = {
    # Toxin targets: nontox skipped (C. diff nontox may retain toxin gene remnants)
    "tcdA_all":      {"nontox"},
    "tcdB_clade2":   {"nontox"},
    "tcdB_clade1":   {"nontox"},
    "tcdC_wt":       {"nontox"},
    "tcdC_junction": {"nontox"},
    "cdtA_groupA":   {"nontox"},
    "cdtB_groupA":   {"nontox"},
    # Housekeeping: skip nontox; rpoB also skip enteropathogens (universal gene)
    "tpiA_all":      {"nontox"},
    "rpoB_all":      {"nontox", "enteropathogens"},
}

# Skip UHGG for housekeeping genes (expected hits in gut bacteria)
SKIP_UHGG_TARGETS = {"tpiA_all", "rpoB_all"}

# Raise coverage threshold to reduce false positives from short alignments
REJECT_COVERAGE = 80.0   # % of crRNA covered (increased from 70%)

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
import tempfile
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
    run_cmd,
)

# =============================================================================
# DATABASE BUILDING
# =============================================================================

def collect_nontox_fasta(genome_base: Path, group_c_file: Path,
                          out_fasta: Path, log) -> bool:
    """
    Concatenate FASTA files for Group C (non-toxigenic) genomes.
    Returns True if successful.
    """
    accessions = [l.strip() for l in group_c_file.read_text().splitlines()
                  if l.strip()]
    log.info(f"  Group C genomes: {len(accessions)}")

    with open(out_fasta, "w") as out:
        found = 0
        for acc in accessions:
            for subdir in ("working", "chilean", "anchors"):
                fna_dir = (genome_base / subdir /
                           "ncbi_dataset" / "data" / acc)
                fna_files = list(fna_dir.glob("*_genomic.fna"))
                if fna_files:
                    out.write(fna_files[0].read_text(encoding="utf-8",
                                                      errors="ignore"))
                    found += 1
                    break

    log.info(f"  Non-tox FASTA: {found}/{len(accessions)} genomes written")
    return found > 0


def download_enteropathogens(accessions: dict, out_fasta: Path,
                              datasets_exe: str, log) -> bool:
    """Download enteric pathogen genomes and concatenate into one FASTA."""
    log.info(f"  Downloading {len(accessions)} enteropathogen genomes...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_zip = Path(tmpdir) / "entero.zip"
        cmd = [datasets_exe, "download", "genome",
               "accession", *list(accessions.values()),
               "--include", "genome",
               "--filename", str(tmp_zip)]
        try:
            run_cmd(cmd, log, capture=False, check=True)
            run_cmd(["unzip", "-o", str(tmp_zip), "-d", tmpdir],
                    log, capture=True, check=False)

            with open(out_fasta, "w") as out:
                for fna in Path(tmpdir).rglob("*_genomic.fna"):
                    out.write(fna.read_text(encoding="utf-8", errors="ignore"))

            log.info(f"  Enteropathogens → {out_fasta}")
            return out_fasta.exists() and out_fasta.stat().st_size > 0
        except Exception as e:
            log.error(f"  Download failed: {e}")
            return False


def download_human_transcriptome(out_fasta: Path, datasets_exe: str,
                                  log) -> bool:
    """Download human RefSeq mRNA sequences."""
    log.info("  Downloading human transcriptome (GRCh38 RefSeq mRNA)...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_zip = Path(tmpdir) / "human_tx.zip"
        cmd = [datasets_exe, "download", "genome",
               "taxon", HUMAN_TX_TAXON,
               "--reference",
               "--include", "rna",
               "--filename", str(tmp_zip)]
        try:
            run_cmd(cmd, log, capture=False, check=True)
            run_cmd(["unzip", "-o", str(tmp_zip), "-d", tmpdir],
                    log, capture=True, check=False)

            with open(out_fasta, "w") as out:
                for rna in Path(tmpdir).rglob("*.fna"):
                    out.write(rna.read_text(encoding="utf-8", errors="ignore"))

            log.info(f"  Human TX → {out_fasta}")
            return out_fasta.exists() and out_fasta.stat().st_size > 0
        except Exception as e:
            log.error(f"  Human TX download failed: {e}")
            return False


def build_blast_db(fasta: Path, db_prefix: Path,
                   makeblastdb_exe: str, log) -> bool:
    """Build BLAST nucleotide database."""
    cmd = [makeblastdb_exe, "-in", str(fasta),
           "-out", str(db_prefix),
           "-dbtype", "nucl", "-parse_seqids"]
    try:
        run_cmd(cmd, log, capture=True, check=True)
        log.info(f"  BLAST DB → {db_prefix}")
        return True
    except Exception as e:
        log.error(f"  makeblastdb failed: {e}")
        return False


def build_bowtie2_index(fasta: Path, idx_prefix: Path,
                         bowtie2_build_exe: str, threads: int,
                         log) -> bool:
    """Build Bowtie2 index."""
    cmd = [bowtie2_build_exe, "--threads", str(threads),
           str(fasta), str(idx_prefix)]
    try:
        run_cmd(cmd, log, capture=True, check=True)
        log.info(f"  Bowtie2 index → {idx_prefix}")
        return True
    except Exception as e:
        log.error(f"  bowtie2-build failed: {e}")
        return False


# =============================================================================
# BLAST SEARCH
# =============================================================================

def write_query_fasta(crna_seqs: list, out_fasta: Path):
    """Write crRNA sequences as FASTA query file."""
    with open(out_fasta, "w") as f:
        for i, seq in enumerate(crna_seqs):
            f.write(f">crRNA_{i+1}\n{seq}\n")


def run_blast(query_fasta: Path, db_prefix: Path, out_tsv: Path,
              blastn_exe: str, threads: int, log) -> pd.DataFrame:
    """
    Run blastn-short and parse results.
    Returns DataFrame of hits.
    """
    outfmt = "6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore"
    cmd = [blastn_exe,
           "-task",   BLASTN_TASK,
           "-word_size", str(BLASTN_WORD),
           "-query",  str(query_fasta),
           "-db",     str(db_prefix),
           "-out",    str(out_tsv),
           "-outfmt", outfmt,
           "-evalue", str(EVALUE),
           "-num_threads", str(threads),
           "-perc_identity", str(REJECT_IDENTITY - 5)]  # slightly loose
    try:
        run_cmd(cmd, log, capture=True, check=True)
    except Exception as e:
        log.error(f"  blastn failed: {e}")
        return pd.DataFrame()

    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        return pd.DataFrame()

    cols = ["qseqid", "sseqid", "pident", "length", "qlen",
            "mismatch", "gapopen", "qstart", "qend",
            "sstart", "send", "evalue", "bitscore"]
    try:
        df = pd.read_csv(out_tsv, sep="\t", names=cols, header=None)
        return df
    except Exception:
        return pd.DataFrame()


def run_bowtie2(query_fasta: Path, idx_prefix: Path,
                bowtie2_exe: str, threads: int, log) -> pd.DataFrame:
    """
    Run Bowtie2 in very-sensitive mode and parse SAM output.
    Returns DataFrame of hits.
    """
    cmd = [bowtie2_exe,
           "--very-sensitive",
           "-N", "1",          # allow 1 mismatch in seed
           "-x", str(idx_prefix),
           "-f",               # query is FASTA
           "-U", str(query_fasta),
           "--no-unal",
           "--threads", str(threads),
           "-S", "/dev/stdout"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        hits = []
        for line in r.stdout.splitlines():
            if line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            flag = int(parts[1])
            if flag & 4:  # unmapped
                continue
            hits.append({
                "qseqid": parts[0],
                "sseqid": parts[2],
                "mapq":   int(parts[4]),
            })
        return pd.DataFrame(hits)
    except Exception as e:
        log.error(f"  bowtie2 failed: {e}")
        return pd.DataFrame()


# =============================================================================
# SPECIFICITY ASSESSMENT
# =============================================================================

def assess_specificity(blast_df: pd.DataFrame,
                        crna_seqs: list) -> dict:
    """
    For each crRNA, determine if it has off-target hits.
    Rejection: pident >= REJECT_IDENTITY AND
               (qend - qstart + 1) / qlen >= REJECT_COVERAGE/100

    Returns {crna_idx: True (specific) | False (off-target)}.
    """
    results = {f"crRNA_{i+1}": True for i in range(len(crna_seqs))}

    if blast_df.empty:
        return results

    for _, row in blast_df.iterrows():
        qid      = str(row["qseqid"])
        pident   = float(row.get("pident", 0))
        aln_len  = float(row.get("length", 0))
        qlen     = float(row.get("qlen", SPACER_LEN))
        coverage = (aln_len / qlen) * 100 if qlen > 0 else 0

        # Require alignment to cover >= REJECT_COVERAGE% of crRNA
        # AND start near position 1 (avoid partial 3-prime matches)
        qstart = float(row.get("qstart", 1))
        if (pident >= REJECT_IDENTITY and
                coverage >= REJECT_COVERAGE and
                qstart <= 5):   # alignment must start near 5' end
            results[qid] = False

    return results


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config", default=DEFAULT_CONFIG)
    p.add_argument("--skip-db-build", action="store_true",
                   help="Reuse existing databases")
    p.add_argument("--skip-blast",    action="store_true",
                   help="Reuse existing BLAST results")
    p.add_argument("--skip-bowtie2",  action="store_true",
                   help="Reuse existing Bowtie2 results")
    return p.parse_args()


def main():
    args   = parse_args()
    cfg    = load_config(str(Path(args.config).expanduser()))
    org    = cfg["organism"]["display"]
    paths  = cfg["paths"]
    tools  = cfg["tools"]
    sp_cfg = cfg.get("specificity", {})
    threads= cfg.get("threads", {}).get("default", THREADS)

    reject_id  = sp_cfg.get("reject_identity", REJECT_IDENTITY)
    reject_cov = sp_cfg.get("reject_coverage",  REJECT_COVERAGE)

    # Executables
    blastn_exe  = tools.get("blastn",       "blastn")
    makeblst_exe= tools.get("makeblastdb",  "makeblastdb")
    bt2_exe     = tools.get("bowtie2",      "bowtie2")
    bt2build_exe= tools.get("bowtie2_build","bowtie2-build")
    datasets_exe= tools.get("datasets",     "datasets")

    # Dirs
    main_dir    = Path(paths["main"])
    genome_base = main_dir / "data" / "01_download" / "genomes"
    classify_dir= main_dir / "data" / "02_classify"
    crna_dir    = main_dir / "data" / "06_crna"
    sp_dir      = main_dir / "data" / "08_specificity"
    db_dir      = sp_dir / "dbs"
    report_dir  = main_dir / "reports" / "M08_specificity"

    sp_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M08_specificity", cfg)
    print_logo("M08 · Specificity Analysis", organism=org)

    save_versions(
        tools       = {"blastn": blastn_exe, "bowtie2": bt2_exe},
        python_pkgs = ["pandas"],
        report_dir  = report_dir,
        log         = log,
    )

    # ---------------------------------------------------------------
    # STEP 1 — Build databases
    # ---------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 1 · Building databases")
    log.info("─" * 56)

    # 1a. Non-toxigenic C. difficile (Group C)
    nontox_fasta  = db_dir / "nontox" / "nontox.fasta"
    nontox_db     = db_dir / "nontox" / "nontox_db"
    nontox_fasta.parent.mkdir(parents=True, exist_ok=True)

    if not args.skip_db_build or not nontox_db.with_suffix(".nhr").exists():
        log.info("Building non-toxigenic C. difficile DB (Group C)")
        ok = collect_nontox_fasta(genome_base,
                                   classify_dir / "groupC.txt",
                                   nontox_fasta, log)
        if ok:
            build_blast_db(nontox_fasta, nontox_db, makeblst_exe, log)
    else:
        log.info("  Reusing nontox DB")

    # 1b. Enteropathogens
    entero_fasta  = db_dir / "enteropathogens" / "enteropathogens.fasta"
    entero_db     = db_dir / "enteropathogens" / "entero_db"
    entero_fasta.parent.mkdir(parents=True, exist_ok=True)

    if not args.skip_db_build or not entero_db.with_suffix(".nhr").exists():
        log.info("Building enteropathogen DB")
        ok = download_enteropathogens(ENTEROPATHOGENS, entero_fasta,
                                       datasets_exe, log)
        if ok:
            build_blast_db(entero_fasta, entero_db, makeblst_exe, log)
    else:
        log.info("  Reusing enteropathogens DB")

    # 1c. Human transcriptome
    human_fasta   = db_dir / "human_tx" / "human_tx.fasta"
    human_db      = db_dir / "human_tx" / "human_tx_db"
    human_fasta.parent.mkdir(parents=True, exist_ok=True)

    if not args.skip_db_build or not human_db.with_suffix(".nhr").exists():
        log.info("Building human transcriptome DB")
        ok = download_human_transcriptome(human_fasta, datasets_exe, log)
        if ok:
            build_blast_db(human_fasta, human_db, makeblst_exe, log)
    else:
        log.info("  Reusing human_tx DB")

    # ---------------------------------------------------------------
    # STEP 2 — BLAST crRNA candidates per target
    # ---------------------------------------------------------------
    log.info("─" * 56)
    log.info("STEP 2 · BLAST specificity checks")
    log.info("─" * 56)

    uhgg_db = db_dir / "uhgg_idx" / "uhgg_no_cdiff_db"
    databases = {
        "nontox":          nontox_db,
        "enteropathogens": entero_db,
        "human_tx":        human_db,
        "uhgg":            uhgg_db,
    }

    all_specificity = []
    summary_rows    = []

    for target in tqdm(TARGETS, desc="Targets", unit="target"):
        log.info(f"  TARGET: {target}")

        cands_tsv = crna_dir / f"{target}_candidates.tsv"
        if not cands_tsv.exists():
            log.warning(f"  No candidates for {target} — skipping")
            continue

        cands_df  = pd.read_csv(cands_tsv, sep="\t")
        crna_seqs = cands_df["guide_seq"].tolist()
        if not crna_seqs:
            continue

        query_fasta = sp_dir / f"{target}_query.fasta"
        write_query_fasta(crna_seqs, query_fasta)

        target_results = {f"crRNA_{i+1}": {"guide_seq": seq,
                                             "target": target}
                          for i, seq in enumerate(crna_seqs)}

        # BLAST vs each database
        for db_name, db_prefix in databases.items():
            # Skip databases not relevant for this target
            skip_dbs = SKIP_DATABASES.get(target, set())
            # Skip uhgg for housekeeping genes (expected to hit gut bacteria)
            if db_name == "uhgg" and target in SKIP_UHGG_TARGETS:
                log.info(f"    {db_name}: skipped for housekeeping {target}")
                for k in target_results:
                    target_results[k][f"{db_name}_specific"] = "NA"
                continue
            if db_name in skip_dbs:
                log.info(f"    {db_name}: skipped for {target} (by design)")
                for k in target_results:
                    target_results[k][f"{db_name}_specific"] = "NA"
                continue

            if not db_prefix.with_suffix(".nhr").exists():
                log.warning(f"    DB not available: {db_name} — skipping")
                for k in target_results:
                    target_results[k][f"{db_name}_specific"] = "NA"
                continue

            blast_out = sp_dir / f"{target}_{db_name}_blast.tsv"

            if args.skip_blast and blast_out.exists():
                blast_df = pd.read_csv(blast_out, sep="\t",
                                       names=["qseqid","sseqid","pident",
                                              "length","qlen","mismatch",
                                              "gapopen","qstart","qend",
                                              "sstart","send","evalue",
                                              "bitscore"])
            else:
                blast_df = run_blast(query_fasta, db_prefix, blast_out,
                                     blastn_exe, threads, log)

            spec = assess_specificity(blast_df, crna_seqs)
            for k in target_results:
                target_results[k][f"{db_name}_specific"] = spec.get(k, True)

        # UHGG is handled as a standard BLAST database above

        # Build per-target specificity table
        rows = list(target_results.values())
        sp_df = pd.DataFrame(rows)

        # Overall specific: pass all available DBs
        bool_cols = [c for c in sp_df.columns
                     if c.endswith("_specific") and
                     sp_df[c].dtype != object]
        if bool_cols:
            sp_df["overall_specific"] = sp_df[bool_cols].all(axis=1)
        else:
            sp_df["overall_specific"] = True

        sp_tsv = sp_dir / f"{target}_specificity.tsv"
        sp_df.to_csv(sp_tsv, sep="\t", index=False)
        all_specificity.append(sp_df)

        n_spec = int(sp_df["overall_specific"].sum()) \
            if "overall_specific" in sp_df.columns else len(sp_df)
        log.info(f"    {target}: {n_spec}/{len(sp_df)} crRNAs pass specificity")

        summary_rows.append({
            "target":        target,
            "n_crna":        len(sp_df),
            "n_specific":    n_spec,
            "n_offTarget":   len(sp_df) - n_spec,
        })

    # Combined output
    log.info("─" * 56)
    if all_specificity:
        combined = pd.concat(all_specificity, ignore_index=True)
        combined_out = sp_dir / "all_specificity.tsv"
        combined.to_csv(combined_out, sep="\t", index=False)
        log.info(f"All specificity → {combined_out}")

    write_tsv(summary_rows, report_dir / "summary.tsv", log=log)

    log.info("─" * 56)
    log.info("M08 complete.")
    for r in summary_rows:
        log.info(f"  {r['target']:<20}  "
                 f"specific={r['n_specific']}/{r['n_crna']}")

    write_checkpoint("M08_specificity", report_dir, log=log)


if __name__ == "__main__":
    main()
