"""
pipeline_utils.py
SHERLOCK crRNA Design Pipeline · v1.0

Shared utilities imported by every module.
"""

import csv
import importlib
import logging
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

PIPELINE_VERSION = "1.0"

# =============================================================================
# LOGO
# =============================================================================

_LOGO = r"""
 _____ _  _ ___ ___ _    ___   ___ _  _
/ ____| || | __| _ \ |  / _ \ / __| |/ /
\__ \| __ | _||   / |_| (_) | (__| ' 
|___/|_||_|___|_|_\____\___/ \___|_|\_\

  CRISPR-Cas13a  ·  crRNA Design  ·  v{version}
  Organism : {organism}
  Module   : {module}
  Run date : {ts}
{sep}"""


def print_logo(module_name: str,
               organism: str = "Clostridioides difficile") -> None:
    print(_LOGO.format(
        version  = PIPELINE_VERSION,
        organism = organism,
        module   = module_name,
        ts       = datetime.now().strftime("%Y-%m-%d  %H:%M:%S"),
        sep      = "─" * 56,
    ))


# =============================================================================
# CONFIG
# =============================================================================

def load_config(config_path: str) -> Dict[str, Any]:
    """Load and validate config.yaml."""
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(
            f"config.yaml not found: {p}\n"
            f"Expected at: ~/sherlock/config.yaml"
        )
    with open(p, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    required = ["organism", "paths", "tools", "conda"]
    missing  = [k for k in required if k not in cfg]
    if missing:
        raise ValueError(f"config.yaml missing keys: {missing}")
    return cfg


# =============================================================================
# LOGGER
# =============================================================================

def get_logger(module_id: str,
               cfg: Dict[str, Any],
               level: int = logging.INFO) -> logging.Logger:
    """
    Logger writing to console AND ~/sherlock/main/logs/{module_id}.log
    Log file is overwritten on each run.
    """
    logs_dir = Path(cfg["paths"]["main"]) / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = logs_dir / f"{module_id}.log"

    fmt = logging.Formatter(
        fmt     = "[%(asctime)s]  %(levelname)-8s  %(message)s",
        datefmt = "%H:%M:%S",
    )
    log = logging.getLogger(module_id)
    log.setLevel(level)
    if log.handlers:
        log.handlers.clear()

    for h in (
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file, mode="w", encoding="utf-8"),
    ):
        h.setLevel(level)
        h.setFormatter(fmt)
        log.addHandler(h)

    log.info(f"Log → {log_file}")
    return log


# =============================================================================
# REPORTS
# =============================================================================

def get_report_dir(module_id: str,
                   cfg: Dict[str, Any]) -> Path:
    """
    Return ~/sherlock/main/reports/{module_id}/
    Created if absent. Overwrites existing reports on re-run.
    """
    report_dir = Path(cfg["paths"]["main"]) / "reports" / module_id
    report_dir.mkdir(parents=True, exist_ok=True)
    return report_dir


# =============================================================================
# REPRODUCIBILITY
# =============================================================================

def save_versions(tools: Dict[str, str],
                  python_pkgs: List[str],
                  report_dir: Path,
                  log: Optional[logging.Logger] = None) -> None:
    """Save tool + package versions to versions.tsv in report_dir."""
    import platform
    rows: List[Dict] = []

    rows += [
        {"category": "system", "name": "OS",
         "version": platform.platform(), "path": ""},
        {"category": "system", "name": "Python",
         "version": sys.version.split()[0], "path": sys.executable},
    ]
    for name, exe in tools.items():
        rows.append({
            "category": "cli",
            "name":     name,
            "version":  _cli_version(exe),
            "path":     shutil.which(exe) or "not found",
        })
    for pkg in python_pkgs:
        try:
            mod = importlib.import_module(pkg)
            ver = getattr(mod, "__version__", "?")
        except ImportError:
            ver = "NOT FOUND"
        rows.append({"category": "python_pkg", "name": pkg,
                     "version": ver, "path": ""})
    rows.append({"category": "git", "name": "commit",
                 "version": _git_hash(), "path": ""})

    write_tsv(rows, report_dir / "versions.tsv", log=log)


def _cli_version(exe: str) -> str:
    for flag in ["--version", "-version", "-v", "version"]:
        try:
            r = subprocess.run(
                [exe, flag], capture_output=True, text=True, timeout=8)
            out = (r.stdout + r.stderr).strip()
            if out:
                lines = [l for l in out.splitlines() if l.strip()]
                if lines:
                    return lines[0][:80]
        except Exception:
            continue
    return "unknown"


def _git_hash() -> str:
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5)
        return r.stdout.strip() or "not-a-git-repo"
    except Exception:
        return "git-unavailable"


# =============================================================================
# SUBPROCESS
# =============================================================================

def run_cmd(cmd: List[str],
            log: logging.Logger,
            cwd: Optional[str] = None,
            check: bool = True,
            capture: bool = False) -> subprocess.CompletedProcess:
    """
    Run external command with logging.
    External tools keep their own output (not suppressed).
    Use capture=True only when output must be parsed.
    """
    log.info("CMD: " + " ".join(str(c) for c in cmd))
    result = subprocess.run(
        [str(c) for c in cmd],
        cwd            = cwd,
        check          = False,
        capture_output = capture,
        text           = capture,
    )
    if check and result.returncode != 0:
        msg = f"Failed (exit {result.returncode}): {' '.join(str(c) for c in cmd)}"
        if capture:
            msg += f"\nSTDERR: {result.stderr}"
        log.error(msg)
        raise RuntimeError(msg)
    return result


# =============================================================================
# CHECKPOINT
# =============================================================================

def write_checkpoint(module_id: str,
                     report_dir: Path,
                     log: Optional[logging.Logger] = None) -> None:
    """Write DONE.txt — used by Snakemake (M00b) to skip completed modules."""
    cp = report_dir / "DONE.txt"
    cp.write_text(
        f"module    : {module_id}\n"
        f"completed : {datetime.now().isoformat()}\n"
        f"pipeline  : SHERLOCK v{PIPELINE_VERSION}\n"
    )
    if log:
        log.info(f"Checkpoint → {cp}")


# =============================================================================
# I/O
# =============================================================================

def write_tsv(rows: List[Dict],
              out: Path,
              log: Optional[logging.Logger] = None) -> None:
    """Write list of dicts to TSV. Overwrites if exists."""
    if not rows:
        if log:
            log.warning(f"No data — skipping {out}")
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    if log:
        log.info(f"Wrote {len(rows)} rows → {out}")


# =============================================================================
# FASTA
# =============================================================================

def count_seqs(fasta: Path) -> int:
    """Count sequences in FASTA without loading into memory."""
    n = 0
    with open(fasta, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n
