#!/usr/bin/env python3
"""
M00_verify.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M00 · Installation Verification

PURPOSE
    Verifies all tools, Python packages, paths, and databases
    required by the pipeline.
    Generates a TSV + TXT report and a DONE.txt checkpoint.

USAGE
    conda activate sherlock_main
    python ~/scripts/M00_verify.py --config ~/scripts/config.yaml

OUTPUT
    ~/main/00_setup/reports/versions.tsv
    ~/main/00_setup/reports/verification.tsv
    ~/main/00_setup/reports/verification.txt
    ~/main/00_setup/reports/DONE.txt
    ~/main/00_setup/logs/M00_verify_YYYYMMDD_HHMMSS.log
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# Edit to add / remove tools for a different organism.
# =============================================================================

DEFAULT_CONFIG = "~/scripts/config.yaml"

# CLI tools to verify: {display_name: executable}
# Tools that live in sherlock env
TOOLS = {
    # NCBI retrieval
    "datasets":       "datasets",
    "efetch":         "efetch",
    "esearch":        "esearch",
    # QC and typing
    "mlst":           "mlst",
    "BUSCO":          "busco",
    # Annotation
    # Sequence analysis
    "MAFFT":          "mafft",
    "blastn":         "blastn",
    "makeblastdb":    "makeblastdb",
    "Bowtie2":        "bowtie2",
    "bowtie2-build":  "bowtie2-build",
    "seqkit":         "seqkit",
    "Primer3":        "primer3_core",
    # RNA structure
    "RNAplfold":      "RNAplfold",
    # Phylogenetics
    "IQ-TREE":        "iqtree",
    "ModelTest-NG":   "modeltest-ng",
    # Workflow
    "Snakemake":      "snakemake",
}

# Python packages to verify: {display_name: import_name}
PACKAGES = {
    "Biopython":  "Bio",
    "pandas":     "pandas",
    "numpy":      "numpy",
    "matplotlib": "matplotlib",
    "seaborn":    "seaborn",
    "PyYAML":     "yaml",
    "tqdm":       "tqdm",
    "click":      "click",
    "scipy":      "scipy",
}

# Bakta DB expected files (subset check)
BAKTA_DB_FILES = ["bakta.db", "pfam.h3m"]

# BUSCO lineage expected for this organism
BUSCO_LINEAGES = ["clostridia_odb10", "clostridiales_odb10"]

# =============================================================================
# imports
# =============================================================================

import argparse
import importlib
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# Add ~/scripts to path so pipeline_utils is importable
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
# helpers
# =============================================================================

def check_tool(exe: str):
    """Return (found: bool, info: str)."""
    path = shutil.which(exe)
    if path is None:
        return False, "NOT FOUND"
    # Try common version flags
    for flag in ["--version", "-version", "-v", "version"]:
        try:
            r = subprocess.run(
                [exe, flag], capture_output=True, text=True, timeout=8)
            out = (r.stdout + r.stderr).strip()
            if out:
                lines = [l for l in out.splitlines() if l.strip()]
                if lines:
                    return True, f"{path}  |  {lines[0][:70]}"
        except Exception:
            continue
    return True, f"{path}  |  version unknown"


def check_package(import_name: str):
    """Return (found: bool, version: str)."""
    try:
        mod = importlib.import_module(import_name)
        ver = getattr(mod, "__version__", "?")
        return True, f"v{ver}"
    except ImportError:
        return False, "NOT FOUND"


def check_dir(path: str):
    """Return (exists: bool | None, info: str)."""
    if not path:
        return None, "not configured"
    p = Path(path).expanduser()
    if p.exists():
        n = sum(1 for _ in p.iterdir())
        return True, f"{p}  ({n} items)"
    return False, f"NOT FOUND — {p}"


def check_bakta_db(db_path: str):
    """Return (ok: bool | None, info: str)."""
    if not db_path:
        return None, "paths.bakta_db not set — update config.yaml"
    p = Path(db_path).expanduser()
    if not p.exists():
        return False, f"directory not found: {p}"
    missing = [f for f in BAKTA_DB_FILES if not (p / f).exists()]
    if missing:
        return False, f"incomplete DB — missing: {', '.join(missing)}"
    return True, str(p)


def check_busco_lineages(lineages_dir: str):
    """Return list of (lineage, found: bool, info: str)."""
    results = []
    if not lineages_dir:
        for lin in BUSCO_LINEAGES:
            results.append((lin, None, "paths.busco_lineages not set"))
        return results
    base = Path(lineages_dir).expanduser()
    for lin in BUSCO_LINEAGES:
        p = base / lin
        if p.exists():
            results.append((lin, True, str(p)))
        else:
            results.append((lin, False, f"NOT FOUND: {p}"))
    return results


def check_badgers(badgers_dir: str):
    """Return (found: bool, info: str)."""
    if not badgers_dir:
        return False, "tools.badgers_dir not set in config.yaml"
    p = Path(badgers_dir).expanduser()
    if p.exists():
        return True, str(p)
    return False, f"not found: {p}"


def check_conda_env(env_name: str):
    """Return (found: bool, info: str)."""
    try:
        r = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True, text=True, timeout=30)
        if env_name in r.stdout:
            return True, f"'{env_name}' exists"
        return False, f"'{env_name}' NOT FOUND"
    except Exception:
        return False, "conda unavailable"


def check_metadata(meta_path: str):
    """Return (found: bool | None, info: str)."""
    if not meta_path:
        return None, "paths.metadata not set"
    p = Path(meta_path).expanduser()
    if p.exists():
        n = sum(1 for _ in open(p)) - 1
        return True, f"{p}  ({n} genomes)"
    return None, f"NOT FOUND — {p}\nCopy with:  cp ncbi_dataset.tsv {p}"


# =============================================================================
# main
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--config",
        default=DEFAULT_CONFIG,
        help=f"Path to config.yaml (default: {DEFAULT_CONFIG})",
    )
    return p.parse_args()


def main():
    args = parse_args()
    cfg  = load_config(str(Path(args.config).expanduser()))

    org   = cfg.get("organism", {}).get("display", "Unknown")
    paths = cfg.get("paths", {})
    tools = cfg.get("tools", {})
    conda = cfg.get("conda", {})

    # Dirs
    log        = get_logger("M00_verify", cfg)
    report_dir = Path(cfg["paths"]["main"]) / "reports" / "M00_verify"
    report_dir.mkdir(parents=True, exist_ok=True)

    # Logger

    # Logo
    print_logo("M00 · Installation Verification", organism=org)

    # Collect results
    rows = []   # {section, item, status, info}

    def record(section, item, ok, info):
        if ok is True:
            status = "PASS"
            log.info(f"[{status}]  {item:<32} {info}")
        elif ok is False:
            status = "FAIL"
            log.error(f"[{status}]  {item:<32} {info}")
        else:
            status = "WARN"
            log.warning(f"[{status}]  {item:<32} {info}")
        rows.append({"section": section, "item": item,
                     "status": status, "info": info})

    # ------------------------------------------------------------------
    # 1. Conda environments
    # ------------------------------------------------------------------
    log.info("─" * 55)
    log.info("1 · CONDA ENVIRONMENTS")
    log.info("─" * 55)

    for key in ("main_env", "badgers_env"):
        name = conda.get(key, "")
        ok, info = check_conda_env(name)
        record("conda", name, ok, info)

    # ------------------------------------------------------------------
    # 2. CLI tools
    # ------------------------------------------------------------------
    log.info("─" * 55)
    log.info("2 · COMMAND-LINE TOOLS")
    log.info("─" * 55)

    for display, exe in TOOLS.items():
        ok, info = check_tool(exe)
        record("cli", display, ok, info)

    # ------------------------------------------------------------------
    # 3. Python packages
    # ------------------------------------------------------------------
    log.info("─" * 55)
    log.info("3 · PYTHON PACKAGES")
    log.info("─" * 55)

    for display, imp in PACKAGES.items():
        ok, info = check_package(imp)
        record("python", display, ok, info)

    # ------------------------------------------------------------------
    # 4. Paths
    # ------------------------------------------------------------------
    log.info("─" * 55)
    log.info("4 · PATHS")
    log.info("─" * 55)

    for key in ("scripts", "main"):
        ok, info = check_dir(paths.get(key, ""))
        record("paths", key, ok, info)

    ok, info = check_bakta_db(paths.get("bakta_db", ""))
    record("paths", "bakta_db", ok, info)

    # ------------------------------------------------------------------
    # 5. BUSCO lineages
    # ------------------------------------------------------------------
    log.info("─" * 55)
    log.info("5 · BUSCO LINEAGES")
    log.info("─" * 55)

    for lin, ok, info in check_busco_lineages(paths.get("busco_lineages", "")):
        record("busco", lin, None if not ok else ok, info)

    # ------------------------------------------------------------------
    # 6. BADGERS
    # ------------------------------------------------------------------
    log.info("─" * 55)
    log.info("6 · BADGERS")
    log.info("─" * 55)

    ok, info = check_badgers(tools.get("badgers_dir", cfg.get("conda", {}).get("badgers_dir", "")))
    record("badgers", "design_guides.py", ok, info)

    # ------------------------------------------------------------------
    # 7. Input metadata file
    # ------------------------------------------------------------------
    log.info("─" * 55)
    log.info("7 · INPUT FILES")
    log.info("─" * 55)

    meta = str(Path(cfg["paths"]["main"]) / "01_download" / "ncbi_dataset.tsv")
    ok, info = check_metadata(meta)
    record("inputs", "ncbi_dataset.tsv", None, info)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    log.info("─" * 55)

    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    n_warn = sum(1 for r in rows if r["status"] == "WARN")
    n_fail = sum(1 for r in rows if r["status"] == "FAIL")

    log.info(f"RESULT   PASS={n_pass}  WARN={n_warn}  FAIL={n_fail}")

    if n_fail == 0:
        log.info("All critical checks passed. Ready for M01.")
    else:
        log.error(f"{n_fail} check(s) failed. Fix before running M01.")

    # ------------------------------------------------------------------
    # Write reports
    # ------------------------------------------------------------------

    # versions.tsv (reproducibility snapshot)
    save_versions(
        tools       = {k: v for k, v in tools.items()
                       if isinstance(v, str) and "/" not in v},
        python_pkgs = list(PACKAGES.values()),
        report_dir  = report_dir,
        log         = log,
    )

    # verification.tsv
    write_tsv(rows, report_dir / "verification.tsv", log=log)

    # verification.txt  (human-readable)
    txt = report_dir / "verification.txt"
    with open(txt, "w") as f:
        f.write(f"SHERLOCK Pipeline — M00 Verification\n")
        f.write(f"Generated : {datetime.now()}\n")
        f.write(f"Organism  : {org}\n\n")
        f.write(f"{'ITEM':<32} {'STATUS':<6}  INFO\n")
        f.write("─" * 80 + "\n")
        for r in rows:
            f.write(f"{r['item']:<32} {r['status']:<6}  {r['info']}\n")
        f.write("\n")
        f.write(f"PASS={n_pass}  WARN={n_warn}  FAIL={n_fail}\n")
    log.info(f"Report → {txt}")

    # DONE.txt checkpoint (only if no failures)
    if n_fail == 0:
        write_checkpoint("M00_verify", Path(report_dir), log=log)

    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
