#!/usr/bin/env python3
"""
M00_verify.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M00 · Installation Verification

PURPOSE
    Verifies all tools, Python packages, conda environments,
    and paths required by the pipeline.
    Reports are overwritten on each run.

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M00_verify.py --config ~/sherlock/config.yaml

OUTPUT
    ~/sherlock/main/logs/M00_verify.log
    ~/sherlock/main/reports/M00_verify/verification.tsv
    ~/sherlock/main/reports/M00_verify/verification.txt
    ~/sherlock/main/reports/M00_verify/versions.tsv
    ~/sherlock/main/reports/M00_verify/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG = "~/sherlock/config.yaml"

# CLI tools in sherlock env: {display_name: executable}
TOOLS = {
    "datasets":      "datasets",
    "efetch":        "efetch",
    "MAFFT":         "mafft",
    "blastn":        "blastn",
    "makeblastdb":   "makeblastdb",
    "Bowtie2":       "bowtie2",
    "bowtie2-build": "bowtie2-build",
    "seqkit":        "seqkit",
    "mlst":          "mlst",
    "IQ-TREE":       "iqtree",
    "RNAplfold":     "RNAplfold",
    "Primer3":       "primer3_core",
    "PrimedRPA":     "PrimedRPA",
    "Snakemake":     "snakemake",
}

# Tools in external conda envs: {display_name: (env, executable)}
EXTERNAL_TOOLS = {
    "bakta":    ("bakta",    "bakta"),
    "checkm2":  ("checkm2", "checkm2"),
}

# Conda environments to verify
CONDA_ENVS = [
    "sherlock",
    "bakta",
    "checkm2",
    "adapt",
    "badgers-cas13",
]

# Python packages: {display_name: import_name}
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

# Bakta DB expected files
BAKTA_DB_FILES = ["bakta.db", "pfam.h3m"]

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import importlib
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

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
# CHECK HELPERS
# =============================================================================

def check_tool(exe: str):
    path = shutil.which(exe)
    if not path:
        return False, "NOT FOUND"
    for flag in ["--version", "-version", "-v", "version"]:
        try:
            r = subprocess.run(
                [exe, flag], capture_output=True, text=True, timeout=8)
            out = (r.stdout + r.stderr).strip()
            if out:
                lines = [l for l in out.splitlines() if l.strip()]
                if lines:
                    return True, lines[0][:60]
        except Exception:
            continue
    return True, "found"


def check_external_tool(env: str, exe: str):
    try:
        r = subprocess.run(
            ["conda", "run", "-n", env, exe, "--version"],
            capture_output=True, text=True, timeout=30)
        out = (r.stdout + r.stderr).strip()
        if out:
            lines = [l for l in out.splitlines() if l.strip()]
            if lines:
                return True, f"[{env}] {lines[0][:55]}"
        return True, f"[{env}] found"
    except Exception as e:
        return False, str(e)


def check_conda_env(env: str):
    try:
        r = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True, text=True, timeout=30)
        if env in r.stdout:
            return True, "exists"
        return False, "NOT FOUND"
    except Exception:
        return False, "conda unavailable"


def check_package(import_name: str):
    try:
        mod = importlib.import_module(import_name)
        ver = getattr(mod, "__version__", "?")
        return True, f"v{ver}"
    except ImportError:
        return False, "NOT FOUND"


def check_dir(path: str):
    if not path:
        return None, "not configured"
    p = Path(path).expanduser()
    if p.exists():
        n = sum(1 for _ in p.iterdir())
        return True, f"{p}  ({n} items)"
    return False, f"NOT FOUND: {p}"


def check_bakta_db(db_path: str):
    if not db_path:
        return None, "paths.bakta_db not set"
    p = Path(db_path).expanduser()
    if not p.exists():
        return False, f"NOT FOUND: {p}"
    missing = [f for f in BAKTA_DB_FILES if not (p / f).exists()]
    if missing:
        return False, f"incomplete — missing: {', '.join(missing)}"
    return True, str(p)


def check_checkm2_db(db_path: str):
    if not db_path:
        return None, "paths.checkm2_db not set"
    p = Path(db_path).expanduser()
    if p.exists():
        return True, str(p)
    return False, f"NOT FOUND: {p}"


def check_casilico():
    try:
        r = subprocess.run(
            ["R", "--vanilla", "-e", "library(CaSilico); cat('OK')"],
            capture_output=True, text=True, timeout=30)
        if "OK" in r.stdout + r.stderr:
            return True, "CaSilico loaded in R"
        return False, "CaSilico not found in R"
    except Exception as e:
        return False, str(e)


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
    args  = parse_args()
    cfg   = load_config(str(Path(args.config).expanduser()))
    org   = cfg["organism"]["display"]
    paths = cfg["paths"]
    tools = cfg["tools"]

    # Dirs
    logs_dir   = Path(paths["main"]) / "logs"
    report_dir = Path(paths["main"]) / "reports" / "M00_verify"
    logs_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M00_verify", cfg)
    print_logo("M00 · Installation Verification", organism=org)

    rows = []

    def record(section, item, ok, info):
        if ok is True:
            status = "PASS"
            log.info(f"  [PASS]  {item:<30}  {info}")
        elif ok is False:
            status = "FAIL"
            log.error(f"  [FAIL]  {item:<30}  {info}")
        else:
            status = "WARN"
            log.warning(f"  [WARN]  {item:<30}  {info}")
        rows.append({"section": section, "item": item,
                     "status": status, "info": info})

    # ------------------------------------------------------------------
    # 1. Conda environments
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("1 · CONDA ENVIRONMENTS")
    log.info("─" * 56)
    for env in CONDA_ENVS:
        ok, info = check_conda_env(env)
        record("conda", env, ok, info)

    # ------------------------------------------------------------------
    # 2. Tools in sherlock env
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("2 · TOOLS  (sherlock env)")
    log.info("─" * 56)
    for display, exe in TOOLS.items():
        ok, info = check_tool(exe)
        record("tools", display, ok, info)

    # ------------------------------------------------------------------
    # 3. Tools in external envs
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("3 · TOOLS  (external envs)")
    log.info("─" * 56)
    for display, (env, exe) in EXTERNAL_TOOLS.items():
        ok, info = check_external_tool(env, exe)
        record("external", display, ok, info)

    # ------------------------------------------------------------------
    # 4. Python packages
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("4 · PYTHON PACKAGES")
    log.info("─" * 56)
    for display, imp in PACKAGES.items():
        ok, info = check_package(imp)
        record("python", display, ok, info)

    # ------------------------------------------------------------------
    # 5. Paths and databases
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("5 · PATHS AND DATABASES")
    log.info("─" * 56)
    for key in ("project", "scripts", "main"):
        ok, info = check_dir(paths.get(key, ""))
        record("paths", key, ok, info)
    ok, info = check_bakta_db(paths.get("bakta_db", ""))
    record("paths", "bakta_db", ok, info)
    ok, info = check_checkm2_db(paths.get("checkm2_db", ""))
    record("paths", "checkm2_db", ok, info)

    # ------------------------------------------------------------------
    # 6. CaSilico (R)
    # ------------------------------------------------------------------
    log.info("─" * 56)
    log.info("6 · CaSilico (R)")
    log.info("─" * 56)
    ok, info = check_casilico()
    record("casilico", "CaSilico", ok, info)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    log.info("─" * 56)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    n_warn = sum(1 for r in rows if r["status"] == "WARN")
    n_fail = sum(1 for r in rows if r["status"] == "FAIL")
    log.info(f"RESULT   PASS={n_pass}  WARN={n_warn}  FAIL={n_fail}")
    if n_fail == 0:
        log.info("All critical checks passed. Ready for M01.")
    else:
        log.error(f"{n_fail} check(s) failed. Fix before running M01.")

    # ------------------------------------------------------------------
    # Reports
    # ------------------------------------------------------------------
    save_versions(
        tools       = {k: v for k, v in tools.items()
                       if isinstance(v, str) and "/" not in v},
        python_pkgs = list(PACKAGES.values()),
        report_dir  = report_dir,
        log         = log,
    )
    write_tsv(rows, report_dir / "verification.tsv", log=log)

    txt = report_dir / "verification.txt"
    with open(txt, "w") as f:
        f.write(f"SHERLOCK Pipeline — M00 Verification\n")
        f.write(f"Generated : {datetime.now()}\n")
        f.write(f"Organism  : {org}\n\n")
        f.write(f"{'ITEM':<30}  {'STATUS':<6}  INFO\n")
        f.write("─" * 70 + "\n")
        for r in rows:
            f.write(f"{r['item']:<30}  {r['status']:<6}  {r['info']}\n")
        f.write(f"\nPASS={n_pass}  WARN={n_warn}  FAIL={n_fail}\n")
    log.info(f"Report → {txt}")

    if n_fail == 0:
        write_checkpoint("M00_verify", report_dir, log=log)

    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
