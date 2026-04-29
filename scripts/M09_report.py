#!/usr/bin/env python3
"""
M09_report.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M09 · Ranking + HTML Report

PURPOSE
    Integrate results from M06 (crRNA design), M07 (RT-RPA primers),
    and M08 (specificity) into a final ranked candidate table.

    Ranking composite score:
      - ADAPT activity score     (0.35) — predicted Cas13a activity
      - Accessibility score      (0.20) — mRNA unpaired probability
      - Conservation score       (0.20) — Shannon entropy (inverted)
      - Specificity bonus        (0.15) — pass all specificity checks
      - Primer co-design bonus   (0.10) — has matched RT-RPA primer set

    Outputs:
      - TSV ranking per target + combined
      - HTML report with interactive table

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M09_report.py --config ~/sherlock/config.yaml

OUTPUT
    main/data/09_report/final_ranking.tsv
    main/data/09_report/final_ranking.html
    main/logs/M09_report.log
    main/reports/M09_report/summary.tsv
    main/reports/M09_report/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG = "~/sherlock/config.yaml"
TOP_N          = 5   # top candidates per target in final ranking

# Weights for composite score
WEIGHTS = {
    "adapt":        0.35,
    "accessibility":0.20,
    "conservation": 0.20,
    "specificity":  0.15,
    "primer":       0.10,
}

TARGETS = [
    "tcdA_all",
    "tcdB_all",
    "tcdC_wt",
    "tcdC_junction",
    "cdtA_groupA",
    "cdtB_groupA",
    "tpiA_all",
    "sodA_all",
    "16S_all",
]

# Reaction assignment per target
REACTION_MAP = {
    "tcdB_all":       "A",
    "tcdA_all":       "A",
    "16S_all":        "A",
    "cdtA_groupA":    "B",
    "cdtB_groupA":    "B",
    "tcdC_wt":        "B",
    "tcdC_junction":  "B",
    "tpiA_all":       "B",
    "sodA_all":       "B",
}

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
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
# LOAD DATA
# =============================================================================

def load_candidates(crna_dir: Path, target: str) -> pd.DataFrame:
    """Load M06 crRNA candidates for a target."""
    f = crna_dir / f"{target}_candidates.tsv"
    if not f.exists():
        return pd.DataFrame()
    return pd.read_csv(f, sep="\t")


def load_specificity(sp_dir: Path, target: str) -> pd.DataFrame:
    """Load M08 specificity results for a target."""
    f = sp_dir / f"{target}_specificity.tsv"
    if not f.exists():
        return pd.DataFrame()
    return pd.read_csv(f, sep="\t")


def load_conservation(cons_dir: Path, target: str) -> pd.DataFrame:
    """Load M04 conservation metrics for a target."""
    f = cons_dir / f"{target}_conservation.tsv"
    if not f.exists():
        return pd.DataFrame()
    return pd.read_csv(f, sep="\t")


def load_codesign(primers_dir: Path, target: str) -> pd.DataFrame:
    """Load M07 co-design results for a target."""
    f = primers_dir / f"{target}_codesign.tsv"
    if not f.exists():
        return pd.DataFrame()
    return pd.read_csv(f, sep="\t")


# =============================================================================
# RANKING
# =============================================================================

def get_conservation_score(cons_df: pd.DataFrame,
                            crna_position: int) -> float:
    """
    Get conservation score at crRNA position.
    Returns 1 - mean_shannon (higher = more conserved).
    """
    if cons_df.empty or "position" not in cons_df.columns:
        return float("nan")
    mask = ((cons_df["position"] >= crna_position + 1) &
            (cons_df["position"] <= crna_position + 28))
    subset = cons_df.loc[mask, "shannon"]
    if subset.empty:
        return float("nan")
    mean_sh = float(subset.mean())
    # Invert: 0 Shannon = perfectly conserved = score 1.0
    return max(0.0, 1.0 - mean_sh / 2.0)


def norm01(series: pd.Series) -> pd.Series:
    """Normalize series to [0, 1]."""
    v = series.dropna()
    if v.empty or v.max() == v.min():
        return pd.Series(0.5, index=series.index)
    return (series - v.min()) / (v.max() - v.min())


def rank_target(target: str,
                cands_df: pd.DataFrame,
                spec_df: pd.DataFrame,
                cons_df: pd.DataFrame,
                codesign_df: pd.DataFrame,
                top_n: int) -> pd.DataFrame:
    """
    Build final ranking for one target integrating all data.
    """
    if cands_df.empty:
        return pd.DataFrame()

    rows = []
    # Build specificity lookup
    spec_lookup = {}
    if not spec_df.empty and "guide_seq" in spec_df.columns:
        for _, r in spec_df.iterrows():
            spec_lookup[r["guide_seq"]] = bool(r.get("overall_specific", True))

    # Build codesign lookup (guide_seq → primer info)
    codesign_lookup = {}
    if not codesign_df.empty and "crna_seq" in codesign_df.columns:
        for _, r in codesign_df.iterrows():
            codesign_lookup[r["crna_seq"]] = {
                "fp_seq":       r.get("fp_seq", ""),
                "fp_t7":        r.get("fp_t7", ""),
                "rp_seq":       r.get("rp_seq", ""),
                "amplicon_size":r.get("amplicon_size", ""),
                "crna_activity":r.get("crna_activity", ""),
            }

    for _, cand in cands_df.iterrows():
        seq      = str(cand.get("guide_seq", ""))
        crna_pos = int(cand.get("crna_position",
                       cand.get("window_start", 0)))

        is_specific  = spec_lookup.get(seq, True)
        has_primer   = seq in codesign_lookup
        cons_score   = get_conservation_score(cons_df, crna_pos)
        primer_info  = codesign_lookup.get(seq, {})

        rows.append({
            "target":          target,
            "reaction":        REACTION_MAP.get(target, "?"),
            "guide_seq":       seq,
            "crna_position":   crna_pos,
            "adapt_activity":  float(cand.get("adapt_activity", 0)),
            "accessibility":   float(cand.get("accessibility", 0)),
            "conservation":    cons_score,
            "specific":        is_specific,
            "has_primer":      has_primer,
            "fp_seq":          primer_info.get("fp_seq", ""),
            "fp_t7":           primer_info.get("fp_t7", ""),
            "rp_seq":          primer_info.get("rp_seq", ""),
            "amplicon_size":   primer_info.get("amplicon_size", ""),
            "gc_content":      float(cand.get("gc_content", 0)),
            "max_polyu":       int(cand.get("max_polyu", 0)),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Normalize all score components
    df["adapt_norm"] = norm01(df["adapt_activity"])
    df["acc_norm"]   = norm01(df["accessibility"])
    df["cons_norm"]  = norm01(df["conservation"])
    df["spec_score"] = df["specific"].astype(float)
    df["primer_score"]= df["has_primer"].astype(float)

    # Fill NaN
    for col in ["adapt_norm", "acc_norm", "cons_norm"]:
        df[col] = df[col].fillna(0.5)

    # Composite score
    df["final_score"] = (
        WEIGHTS["adapt"]        * df["adapt_norm"] +
        WEIGHTS["accessibility"]* df["acc_norm"] +
        WEIGHTS["conservation"] * df["cons_norm"] +
        WEIGHTS["specificity"]  * df["spec_score"] +
        WEIGHTS["primer"]       * df["primer_score"]
    )

    df = (df.sort_values("final_score", ascending=False)
            .drop_duplicates(subset="guide_seq")
            .head(top_n)
            .reset_index(drop=True))
    df.insert(0, "rank", df.index + 1)

    return df


# =============================================================================
# HTML REPORT
# =============================================================================

def make_html(df: pd.DataFrame, out_path: Path, org: str):
    """Generate interactive HTML report."""

    # Color coding per reaction
    reaction_colors = {"A": "#2196F3", "B": "#4CAF50", "?": "#9E9E9E"}

    rows_html = []
    for _, r in df.iterrows():
        react = r.get("reaction", "?")
        color = reaction_colors.get(react, "#9E9E9E")
        spec  = "✅" if r.get("specific", True) else "⚠️"
        primer= "✅" if r.get("has_primer", False) else "—"
        score = f"{r.get('final_score', 0):.3f}"
        act   = f"{r.get('adapt_activity', 0):.3f}"
        acc   = f"{r.get('accessibility', 0):.3f}"
        cons  = f"{r.get('conservation', 0):.3f}"
        amp   = str(r.get("amplicon_size", "—"))

        rows_html.append(f"""
        <tr>
          <td>{r['rank']}</td>
          <td><span style="background:{color};color:white;padding:2px 6px;
              border-radius:3px;font-size:11px">{react}</span></td>
          <td><b>{r['target']}</b></td>
          <td><code style="font-size:11px">{r['guide_seq']}</code></td>
          <td>{r['crna_position']}</td>
          <td>{act}</td>
          <td>{acc}</td>
          <td>{cons}</td>
          <td>{spec}</td>
          <td>{primer}</td>
          <td>{amp}</td>
          <td><b>{score}</b></td>
        </tr>""")

    rows_str = "\n".join(rows_html)
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>SHERLOCK crRNA Candidates — {org}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 20px; color: #333; }}
    h1 {{ color: #1565C0; }}
    h2 {{ color: #424242; font-size: 14px; font-weight: normal; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
    th {{ background: #1565C0; color: white; padding: 8px 10px;
          text-align: left; position: sticky; top: 0; }}
    td {{ padding: 6px 10px; border-bottom: 1px solid #eee; }}
    tr:hover {{ background: #f5f5f5; }}
    tr:nth-child(even) {{ background: #fafafa; }}
    .badge {{ display: inline-block; padding: 2px 8px; border-radius: 10px;
              font-size: 11px; font-weight: bold; }}
    .filter {{ margin: 10px 0; }}
    input {{ padding: 6px; border: 1px solid #ccc; border-radius: 4px; width: 200px; }}
    .summary {{ background: #E3F2FD; padding: 10px 15px; border-radius: 6px;
                margin: 15px 0; font-size: 13px; }}
  </style>
</head>
<body>
  <h1>🔬 SHERLOCK crRNA Candidate Report</h1>
  <h2>Organism: {org} | Generated: {now} | Pipeline v1.0</h2>

  <div class="summary">
    <b>Reaction A</b> (LwCas13a): tcdB + tcdA + 16S &nbsp;|&nbsp;
    <b>Reaction B</b> (LwCas13a): cdtA + cdtB + tcdC_wt + tcdC_junction + tpiA + sodA<br>
    Ranking weights: ADAPT activity 35% · Accessibility 20% · Conservation 20% ·
    Specificity 15% · Primer co-design 10%
  </div>

  <div class="filter">
    🔍 Filter: <input type="text" id="filterInput"
      onkeyup="filterTable()" placeholder="Search guide or target...">
  </div>

  <table id="rankTable">
    <thead>
      <tr>
        <th>#</th>
        <th>Rxn</th>
        <th>Target</th>
        <th>Guide Sequence (5'→3')</th>
        <th>Position</th>
        <th>ADAPT Activity</th>
        <th>Accessibility</th>
        <th>Conservation</th>
        <th>Specific</th>
        <th>Primer</th>
        <th>Amplicon (bp)</th>
        <th>Final Score</th>
      </tr>
    </thead>
    <tbody>
{rows_str}
    </tbody>
  </table>

  <script>
    function filterTable() {{
      var input = document.getElementById("filterInput").value.toLowerCase();
      var rows  = document.getElementById("rankTable").getElementsByTagName("tr");
      for (var i = 1; i < rows.length; i++) {{
        var txt = rows[i].textContent.toLowerCase();
        rows[i].style.display = txt.includes(input) ? "" : "none";
      }}
    }}
  </script>
</body>
</html>"""

    out_path.write_text(html, encoding="utf-8")


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
    rank_cfg = cfg.get("ranking", {})
    top_n  = rank_cfg.get("top_n", TOP_N)

    # Dirs
    main_dir    = Path(paths["main"])
    crna_dir    = main_dir / "data" / "06_crna"
    primers_dir = main_dir / "data" / "07_primers"
    sp_dir      = main_dir / "data" / "08_specificity"
    cons_dir    = main_dir / "data" / "04_alignment" / "conservation"
    out_dir     = main_dir / "data" / "09_report"
    report_dir  = main_dir / "reports" / "M09_report"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M09_report", cfg)
    print_logo("M09 · Ranking + Report", organism=org)

    save_versions(
        tools={}, python_pkgs=["pandas", "numpy"],
        report_dir=report_dir, log=log)

    all_ranked   = []
    summary_rows = []

    for target in TARGETS:
        log.info(f"  {target}")

        cands_df    = load_candidates(crna_dir, target)
        spec_df     = load_specificity(sp_dir, target)
        cons_df     = load_conservation(cons_dir, target)
        codesign_df = load_codesign(primers_dir, target)

        ranked = rank_target(target, cands_df, spec_df,
                              cons_df, codesign_df, top_n)

        if not ranked.empty:
            ranked.to_csv(out_dir / f"{target}_ranked.tsv",
                          sep="\t", index=False)
            all_ranked.append(ranked)
            log.info(f"    → {len(ranked)} candidates ranked")

        summary_rows.append({
            "target":    target,
            "reaction":  REACTION_MAP.get(target, "?"),
            "n_ranked":  len(ranked),
            "top1_score":round(float(ranked.iloc[0]["final_score"]), 3)
                          if not ranked.empty else 0,
            "top1_seq":  ranked.iloc[0]["guide_seq"]
                          if not ranked.empty else "",
        })

    # Combined ranking
    if all_ranked:
        combined = pd.concat(all_ranked, ignore_index=True)
        tsv_out  = out_dir / "final_ranking.tsv"
        html_out = out_dir / "final_ranking.html"
        combined.to_csv(tsv_out, sep="\t", index=False)
        make_html(combined, html_out, org)
        log.info(f"TSV  → {tsv_out}")
        log.info(f"HTML → {html_out}")
        log.info(f"Total candidates in report: {len(combined)}")

    write_tsv(summary_rows, report_dir / "summary.tsv", log=log)

    log.info("─" * 56)
    log.info("M09 complete.")
    for r in summary_rows:
        log.info(f"  {r['target']:<20}  "
                 f"Rxn={r['reaction']}  "
                 f"top1={r['top1_score']}  "
                 f"{r['top1_seq']}")

    write_checkpoint("M09_report", report_dir, log=log)


if __name__ == "__main__":
    main()
