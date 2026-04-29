#!/usr/bin/env python3
"""
M12_figures.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M12 · Publication-Ready Figures (Nature style)

PURPOSE
    Generate all figures for publication in a high-impact journal
    (Nature, Nature Methods, Nature Communications, PLOS Pathogens).

    Figure specifications:
      - DPI: 300 (print quality)
      - Format: PNG
      - Font: Arial (Nature standard)
      - Color palette: colorblind-safe
      - Panel labels: bold uppercase A, B, C...
      - Line width: 0.75 pt minimum
      - Figure width: 89mm (1-column) or 183mm (2-column)

    Figures produced:
      Fig 1. Pipeline overview (schematic)
      Fig 2. Genome classification heatmap (Groups A/B/C)
      Fig 3. Conservation analysis (Shannon entropy per target)
      Fig 4. mRNA accessibility profiles (RNAplfold)
      Fig 5. crRNA candidate ranking (dot plot)
      Fig 6. Specificity matrix (databases × crRNAs)
      Fig 7. RT-RPA amplicon map (co-design schematic)
      Fig 8. Summary table figure (top candidates per reaction)

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M12_figures.py --config ~/sherlock/config.yaml

OUTPUT
    main/data/12_figures/fig1_pipeline.png
    main/data/12_figures/fig2_classification.png
    main/data/12_figures/fig3_conservation.png
    main/data/12_figures/fig4_accessibility.png
    main/data/12_figures/fig5_ranking.png
    main/data/12_figures/fig6_specificity.png
    main/data/12_figures/fig7_amplicon_map.png
    main/data/12_figures/fig8_summary_table.png
    main/logs/M12_figures.log
    main/reports/M12_figures/DONE.txt
"""

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

DEFAULT_CONFIG = "~/sherlock/config.yaml"
DPI            = 300
FONT_FAMILY    = "Arial"
FONT_SIZE      = 7       # Nature: 5-7 pt for labels
TITLE_SIZE     = 8
LABEL_SIZE     = 8
LEGEND_SIZE    = 6
LINE_WIDTH     = 0.75
COL1_WIDTH     = 89 / 25.4   # 89mm in inches (1-column)
COL2_WIDTH     = 183 / 25.4  # 183mm in inches (2-column)

# Colorblind-safe palette (Wong 2011, Nature Methods)
COLORS = {
    "groupA":    "#E69F00",   # orange  — Hypervirulent RT027
    "groupB":    "#56B4E9",   # sky blue — Toxigenic RT012
    "groupC":    "#009E73",   # green   — Non-toxigenic
    "specific":  "#0072B2",   # blue
    "offTarget": "#D55E00",   # vermillion
    "rxnA":      "#0072B2",   # blue
    "rxnB":      "#009E73",   # green
    "neutral":   "#999999",   # grey
    "highlight": "#CC79A7",   # pink
}

TARGETS_ORDER = [
    "tcdA_all", "tcdB_all", "tcdC_wt", "tcdC_junction",
    "cdtA_groupA", "cdtB_groupA", "tpiA_all", "sodA_all", "16S_all",
]

TARGET_LABELS = {
    "tcdA_all":       "tcdA",
    "tcdB_all":       "tcdB",
    "tcdC_wt":        "tcdC (WT)",
    "tcdC_junction":  "tcdC (RT027 jct)",
    "cdtA_groupA":    "cdtA",
    "cdtB_groupA":    "cdtB",
    "tpiA_all":       "tpiA",
    "sodA_all":       "sodA",
    "16S_all":        "16S rRNA",
}

REACTION_MAP = {
    "tcdA_all": "A", "tcdB_all": "A", "16S_all": "A",
    "cdtA_groupA": "B", "cdtB_groupA": "B",
    "tcdC_wt": "B", "tcdC_junction": "B",
    "tpiA_all": "B", "sodA_all": "B",
}

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent))
from pipeline_utils import (
    load_config, get_logger, print_logo,
    save_versions, write_tsv, write_checkpoint,
)

# =============================================================================
# STYLE SETUP
# =============================================================================

def setup_style():
    """Apply Nature journal matplotlib style."""
    plt.rcParams.update({
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Arial", "DejaVu Sans"],
        "font.size":         FONT_SIZE,
        "axes.titlesize":    TITLE_SIZE,
        "axes.labelsize":    LABEL_SIZE,
        "xtick.labelsize":   FONT_SIZE,
        "ytick.labelsize":   FONT_SIZE,
        "legend.fontsize":   LEGEND_SIZE,
        "axes.linewidth":    LINE_WIDTH,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "xtick.major.width": LINE_WIDTH,
        "ytick.major.width": LINE_WIDTH,
        "xtick.minor.width": LINE_WIDTH * 0.5,
        "ytick.minor.width": LINE_WIDTH * 0.5,
        "lines.linewidth":   LINE_WIDTH,
        "patch.linewidth":   LINE_WIDTH,
        "figure.dpi":        DPI,
        "savefig.dpi":       DPI,
        "savefig.bbox":      "tight",
        "savefig.pad_inches":0.05,
    })


def add_panel_label(ax, label, x=-0.18, y=1.05):
    """Add bold panel label (A, B, C...) in Nature style."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=LABEL_SIZE + 2, fontweight="bold",
            va="top", ha="left")


def save_fig(fig, path: Path, log):
    fig.savefig(path, dpi=DPI, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    log.info(f"  → {path.name}")


# =============================================================================
# FIG 2 — GENOME CLASSIFICATION
# =============================================================================

def fig_classification(classify_dir: Path, out_dir: Path, log):
    """
    Bar chart showing genome counts per group with gene matrix summary.
    """
    matrix_f = classify_dir / "gene_matrix.tsv"
    if not matrix_f.exists():
        log.warning("gene_matrix.tsv not found — skipping Fig 2")
        return

    df = pd.read_csv(matrix_f, sep="\t")

    group_counts = df["group"].value_counts()
    groups   = ["A", "B", "C"]
    labels   = ["Hypervirulent\n(RT027-like)", "Toxigenic\n(RT012-like)", "Non-toxigenic"]
    counts   = [group_counts.get(g, 0) for g in groups]
    colors   = [COLORS["groupA"], COLORS["groupB"], COLORS["groupC"]]

    fig, axes = plt.subplots(1, 2, figsize=(COL2_WIDTH, 2.5))

    # Panel A — bar chart
    ax = axes[0]
    bars = ax.bar(range(3), counts, color=colors, width=0.6,
                  edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=FONT_SIZE)
    ax.set_ylabel("Number of genomes")
    ax.set_title("Genome classification")
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                str(count), ha="center", va="bottom", fontsize=FONT_SIZE)
    add_panel_label(ax, "A")

    # Panel B — gene presence heatmap per group
    ax2 = axes[1]
    genes    = ["tcdA", "tcdB", "tcdC", "cdtA", "cdtB", "tpiA", "sodA", "16S"]
    statuses = ["COMPLETE", "PARTIAL", "PSEUDOGENE", "ABSENT"]
    status_colors = {"COMPLETE": "#0072B2", "PARTIAL": "#56B4E9",
                     "PSEUDOGENE": "#E69F00", "ABSENT": "#EEEEEE"}

    group_order = ["A", "B", "C"]
    x_ticks = []
    x_labels = []
    x = 0
    for gi, grp in enumerate(group_order):
        gdf = df[df["group"] == grp]
        for gene in genes:
            if gene not in gdf.columns:
                x += 1
                continue
            counts_s = gdf[gene].value_counts()
            total    = len(gdf)
            bottom   = 0
            for st in statuses:
                n = counts_s.get(st, 0)
                if n > 0:
                    ax2.bar(x, n/total, bottom=bottom,
                            color=status_colors[st],
                            width=0.8, edgecolor="none")
                    bottom += n/total
            x += 1
        x_ticks.append(gi * len(genes) + len(genes)/2 - 0.5)
        x_labels.append(f"Group {grp}")

    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_labels)
    ax2.set_ylabel("Fraction of genomes")
    ax2.set_title("Gene status by group")
    ax2.set_ylim(0, 1)

    # Gene labels at bottom
    gene_x = list(range(len(genes)*3))
    ax2.set_xticks(gene_x, minor=True)
    ax2.set_xticklabels(genes*3, minor=True, fontsize=5, rotation=90)
    ax2.tick_params(axis="x", which="minor", length=0)

    legend_patches = [mpatches.Patch(color=v, label=k)
                      for k, v in status_colors.items()]
    ax2.legend(handles=legend_patches, loc="upper right",
               fontsize=FONT_SIZE - 1, frameon=False)
    add_panel_label(ax2, "B")

    fig.tight_layout()
    save_fig(fig, out_dir / "fig2_classification.png", log)


# =============================================================================
# FIG 3 — CONSERVATION
# =============================================================================

def fig_conservation(cons_dir: Path, out_dir: Path, log):
    """Shannon entropy profiles for all 9 targets."""
    fig, axes = plt.subplots(3, 3, figsize=(COL2_WIDTH, 5.5))
    axes = axes.flatten()

    for i, target in enumerate(TARGETS_ORDER):
        f = cons_dir / f"{target}_conservation.tsv"
        ax = axes[i]
        if not f.exists():
            ax.set_visible(False)
            continue

        df   = pd.read_csv(f, sep="\t")
        pos  = df["position"].values
        sh   = df["shannon"].values
        rxn  = REACTION_MAP.get(target, "?")
        col  = COLORS["rxnA"] if rxn == "A" else COLORS["rxnB"]

        # Rolling mean for clarity
        win = max(1, len(pos)//100)
        sh_smooth = pd.Series(sh).rolling(win, center=True,
                                           min_periods=1).mean().values

        ax.fill_between(pos, sh_smooth, alpha=0.3, color=col)
        ax.plot(pos, sh_smooth, color=col, lw=0.8)
        ax.axhline(0, color="#cccccc", lw=0.5, ls="--")

        mean_sh = np.nanmean(sh)
        ax.set_title(f"{TARGET_LABELS[target]}\n"
                     f"(H̄={mean_sh:.3f})",
                     fontsize=FONT_SIZE, pad=2)
        ax.set_xlabel("Position (nt)", fontsize=FONT_SIZE - 1)
        ax.set_ylabel("Shannon H", fontsize=FONT_SIZE - 1)
        ax.set_ylim(-0.05, 2.1)

        # Reaction label
        ax.text(0.98, 0.95, f"Rxn {rxn}",
                transform=ax.transAxes, fontsize=FONT_SIZE - 1,
                ha="right", va="top",
                color=col, fontweight="bold")

        if i == 0:
            add_panel_label(ax, "A", x=-0.22)

    fig.suptitle("Conservation analysis — Shannon entropy",
                 fontsize=TITLE_SIZE, y=1.01)
    fig.tight_layout()
    save_fig(fig, out_dir / "fig3_conservation.png", log)


# =============================================================================
# FIG 4 — ACCESSIBILITY
# =============================================================================

def fig_accessibility(acc_dir: Path, out_dir: Path, log):
    """mRNA accessibility profiles from RNAplfold."""
    fig, axes = plt.subplots(3, 3, figsize=(COL2_WIDTH, 5.5))
    axes = axes.flatten()

    for i, target in enumerate(TARGETS_ORDER):
        f   = acc_dir / f"{target}_accessibility.tsv"
        ax  = axes[i]
        if not f.exists():
            ax.set_visible(False)
            continue

        df  = pd.read_csv(f, sep="\t")
        pos = df["position"].values
        acc = df["unpaired_prob"].values
        rxn = REACTION_MAP.get(target, "?")
        col = COLORS["rxnA"] if rxn == "A" else COLORS["rxnB"]

        ax.fill_between(pos, acc, alpha=0.25, color=col)
        ax.plot(pos, acc, color=col, lw=0.7, alpha=0.8)
        ax.axhline(0.5, color=COLORS["highlight"], lw=0.8,
                   ls="--", label="min acc = 0.5")

        pct = (acc >= 0.5).mean() * 100
        ax.set_title(f"{TARGET_LABELS[target]}\n"
                     f"({pct:.0f}% accessible)",
                     fontsize=FONT_SIZE, pad=2)
        ax.set_xlabel("Position (nt)", fontsize=FONT_SIZE - 1)
        ax.set_ylabel("P(unpaired)", fontsize=FONT_SIZE - 1)
        ax.set_ylim(-0.02, 1.05)

        if i == 0:
            add_panel_label(ax, "A", x=-0.22)
            ax.legend(fontsize=FONT_SIZE - 2, frameon=False,
                      loc="upper right")

    fig.suptitle("mRNA accessibility (RNAplfold, W=80, L=40)",
                 fontsize=TITLE_SIZE, y=1.01)
    fig.tight_layout()
    save_fig(fig, out_dir / "fig4_accessibility.png", log)


# =============================================================================
# FIG 5 — crRNA RANKING
# =============================================================================

def fig_ranking(report_dir: Path, out_dir: Path, log):
    """Dot plot of final crRNA candidates colored by score components."""
    f = report_dir / "final_ranking.tsv"
    if not f.exists():
        log.warning("final_ranking.tsv not found — skipping Fig 5")
        return

    df = pd.read_csv(f, sep="\t")
    df["target_label"] = df["target"].map(TARGET_LABELS).fillna(df["target"])

    fig, ax = plt.subplots(figsize=(COL2_WIDTH, 3.0))

    rxn_colors = {"A": COLORS["rxnA"], "B": COLORS["rxnB"], "?": COLORS["neutral"]}

    for rxn, gdf in df.groupby("reaction"):
        col = rxn_colors.get(rxn, COLORS["neutral"])
        sc  = ax.scatter(
            gdf["adapt_activity"],
            gdf["accessibility"],
            s=gdf["final_score"] * 80,
            c=col,
            alpha=0.75,
            edgecolors="white",
            linewidths=0.4,
            label=f"Reaction {rxn}",
            zorder=3,
        )

    ax.set_xlabel("ADAPT predicted activity (LwCas13a)")
    ax.set_ylabel("mRNA accessibility (P unpaired)")
    ax.set_title("crRNA candidate landscape")

    # Size legend
    for sz, lbl in [(0.6, "0.60"), (0.75, "0.75"), (0.85, "0.85")]:
        ax.scatter([], [], s=sz*80, c=COLORS["neutral"],
                   alpha=0.75, label=f"Score {lbl}")

    ax.legend(fontsize=FONT_SIZE, frameon=False,
              loc="lower right", ncol=2)
    ax.axhline(0.5, color=COLORS["highlight"], lw=0.8,
               ls="--", alpha=0.7, label="min acc")

    add_panel_label(ax, "A")
    fig.tight_layout()
    save_fig(fig, out_dir / "fig5_ranking.png", log)


# =============================================================================
# FIG 6 — SPECIFICITY MATRIX
# =============================================================================

def fig_specificity(sp_dir: Path, out_dir: Path, log):
    """Heatmap of specificity results per target × database."""
    f = sp_dir / "all_specificity.tsv"
    if not f.exists():
        log.warning("all_specificity.tsv not found — skipping Fig 6")
        return

    df = pd.read_csv(f, sep="\t")

    db_cols = [c for c in df.columns if c.endswith("_specific")
               and c != "overall_specific"]

    db_labels = {
        "nontox_specific":          "Non-tox\nC. difficile",
        "enteropathogens_specific":  "Enteric\npathogens",
        "human_tx_specific":         "Human\ntranscriptome",
        "uhgg_specific":             "Gut\nmicrobiome\n(UHGG)",
    }

    # Build matrix: rows = crRNA × target, cols = databases
    rows = []
    row_labels = []
    for target in TARGETS_ORDER:
        tdf = df[df["target"] == target]
        for j, (_, r) in enumerate(tdf.iterrows()):
            vals = []
            for db in db_cols:
                v = r.get(db, "NA")
                if v == "NA" or str(v) == "nan":
                    vals.append(0.5)   # grey = NA
                elif v is True or str(v).lower() == "true":
                    vals.append(1.0)   # blue = specific
                else:
                    vals.append(0.0)   # red = off-target
            rows.append(vals)
            row_labels.append(f"{TARGET_LABELS.get(target, target)} #{j+1}")

    matrix = np.array(rows)

    fig, ax = plt.subplots(figsize=(COL1_WIDTH + 0.5, len(rows) * 0.22 + 1.5))

    cmap = matplotlib.colors.ListedColormap(
        [COLORS["offTarget"], "#CCCCCC", COLORS["specific"]])
    bounds = [0, 0.33, 0.67, 1.01]
    norm   = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(db_cols)))
    ax.set_xticklabels([db_labels.get(c, c) for c in db_cols],
                       fontsize=FONT_SIZE)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=FONT_SIZE - 1)
    ax.set_title("Specificity analysis", fontsize=TITLE_SIZE, pad=6)

    # Legend
    legend_elements = [
        mpatches.Patch(color=COLORS["specific"],  label="Specific"),
        mpatches.Patch(color="#CCCCCC",           label="Not checked"),
        mpatches.Patch(color=COLORS["offTarget"], label="Off-target"),
    ]
    ax.legend(handles=legend_elements, loc="upper right",
              bbox_to_anchor=(1.0, -0.05), ncol=3,
              fontsize=FONT_SIZE, frameon=False)

    # Target group separators
    count = 0
    for target in TARGETS_ORDER:
        tdf = df[df["target"] == target]
        count += len(tdf)
        ax.axhline(count - 0.5, color="white", lw=1.5)

    add_panel_label(ax, "A")
    fig.tight_layout()
    save_fig(fig, out_dir / "fig6_specificity.png", log)


# =============================================================================
# FIG 7 — AMPLICON MAP
# =============================================================================

def fig_amplicon_map(primers_dir: Path, out_dir: Path, log):
    """
    Schematic of RT-RPA amplicon with crRNA position for each target.
    """
    all_f = primers_dir / "all_codesign.tsv"
    if not all_f.exists():
        log.warning("all_codesign.tsv not found — skipping Fig 7")
        return

    df = pd.read_csv(all_f, sep="\t")
    # Top co-design per target
    df = df.sort_values("crna_activity", ascending=False)
    df = df.drop_duplicates(subset="target").set_index("target")

    fig, axes = plt.subplots(3, 3, figsize=(COL2_WIDTH, 5.5))
    axes = axes.flatten()

    for i, target in enumerate(TARGETS_ORDER):
        ax = axes[i]
        if target not in df.index:
            ax.set_visible(False)
            continue

        row = df.loc[target]
        fp  = int(row.get("fp_start",  0))
        rp  = int(row.get("rp_start",  fp + 150))
        cp  = int(row.get("crna_position", (fp + rp)//2))
        amp = int(row.get("amplicon_size", rp - fp))
        rxn = REACTION_MAP.get(target, "?")
        col = COLORS["rxnA"] if rxn == "A" else COLORS["rxnB"]

        # Normalize to 0-100 for display
        length = max(rp - fp, 1)
        fp_n, rp_n = 0, 100
        cp_n = min(100, max(0, (cp - fp) / length * 100))
        cr_n = min(100, cp_n + 28/length*100)

        # Amplicon bar
        ax.barh(0, 100, height=0.3, color="#E0E0E0",
                edgecolor="#BDBDBD", linewidth=0.5)
        # crRNA region
        ax.barh(0, cr_n - cp_n, left=cp_n, height=0.3,
                color=col, alpha=0.85, edgecolor="white", linewidth=0.5)
        # FP arrow
        ax.annotate("", xy=(8, 0.2), xytext=(0, 0.2),
                    arrowprops=dict(arrowstyle="->", color=COLORS["groupA"],
                                   lw=0.8))
        # RP arrow
        ax.annotate("", xy=(92, 0.2), xytext=(100, 0.2),
                    arrowprops=dict(arrowstyle="->", color=COLORS["groupA"],
                                   lw=0.8))

        ax.set_xlim(-5, 105)
        ax.set_ylim(-0.5, 0.8)
        ax.set_yticks([])
        ax.set_xticks([0, 50, 100])
        ax.set_xticklabels([str(fp), str(fp + amp//2), str(rp)],
                           fontsize=FONT_SIZE - 1)
        ax.set_xlabel("Alignment position (nt)", fontsize=FONT_SIZE - 1)
        ax.set_title(f"{TARGET_LABELS[target]}\n"
                     f"Amplicon={amp}bp | crRNA pos={cp}",
                     fontsize=FONT_SIZE, pad=2)

        if i == 0:
            add_panel_label(ax, "A", x=-0.22)

    fig.suptitle("RT-RPA amplicon + crRNA co-design",
                 fontsize=TITLE_SIZE, y=1.01)
    fig.tight_layout()
    save_fig(fig, out_dir / "fig7_amplicon_map.png", log)


# =============================================================================
# FIG 8 — SUMMARY TABLE
# =============================================================================

def fig_summary_table(report_dir: Path, out_dir: Path, log):
    """Publication-ready summary table of top candidates."""
    f = report_dir / "final_ranking.tsv"
    if not f.exists():
        log.warning("final_ranking.tsv not found — skipping Fig 8")
        return

    df = pd.read_csv(f, sep="\t")
    df["target_label"] = df["target"].map(TARGET_LABELS).fillna(df["target"])

    # Select columns for display
    show_cols = ["reaction", "target_label", "guide_seq_rna",
                 "adapt_activity", "accessibility", "conservation",
                 "specific", "amplicon_size", "final_score"]
    show_labels = ["Rxn", "Target", "Spacer (5'→3' RNA)",
                   "ADAPT\nActivity", "Accessibility",
                   "Conservation", "Specific", "Amplicon\n(bp)", "Score"]

    dfs = df[show_cols].copy()
    dfs["adapt_activity"]  = dfs["adapt_activity"].apply(lambda x: f"{x:.3f}")
    dfs["accessibility"]   = dfs["accessibility"].apply(lambda x: f"{x:.3f}")
    dfs["conservation"]    = dfs["conservation"].apply(
        lambda x: f"{x:.3f}" if not pd.isna(x) else "—")
    dfs["specific"]        = dfs["specific"].apply(
        lambda x: "Yes" if x else "No")
    dfs["amplicon_size"]   = dfs["amplicon_size"].apply(
        lambda x: str(int(x)) if pd.notna(x) and x != "" else "—")
    dfs["final_score"]     = dfs["final_score"].apply(lambda x: f"{x:.3f}")

    fig, ax = plt.subplots(figsize=(COL2_WIDTH, len(dfs) * 0.28 + 1.0))
    ax.axis("off")

    table = ax.table(
        cellText   = dfs.values,
        colLabels  = show_labels,
        loc        = "center",
        cellLoc    = "center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(FONT_SIZE)
    table.scale(1, 1.4)

    # Style header
    for j in range(len(show_labels)):
        table[(0, j)].set_facecolor("#0D47A1")
        table[(0, j)].set_text_props(color="white", fontweight="bold")

    # Color rows by reaction
    rxn_col_map = {"A": "#E3F2FD", "B": "#E8F5E9"}
    for i, (_, row) in enumerate(dfs.iterrows()):
        rxn = row["reaction"]
        bg  = rxn_col_map.get(rxn, "white")
        for j in range(len(show_cols)):
            table[(i+1, j)].set_facecolor(bg)
            # Red for non-specific
            if show_cols[j] == "specific" and row["specific"] == "No":
                table[(i+1, j)].set_facecolor("#FFEBEE")
                table[(i+1, j)].set_text_props(color=COLORS["offTarget"],
                                                fontweight="bold")

    ax.set_title("Top crRNA candidates — SHERLOCK C. difficile",
                 fontsize=TITLE_SIZE, pad=10, fontweight="bold")

    # Legend
    patches = [
        mpatches.Patch(color="#E3F2FD", label="Reaction A"),
        mpatches.Patch(color="#E8F5E9", label="Reaction B"),
    ]
    ax.legend(handles=patches, loc="lower right",
              bbox_to_anchor=(1.0, -0.02), fontsize=FONT_SIZE,
              frameon=False)

    add_panel_label(ax, "A", x=0.0, y=1.02)
    fig.tight_layout()
    save_fig(fig, out_dir / "fig8_summary_table.png", log)


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=DEFAULT_CONFIG)
    return p.parse_args()


def main():
    args = parse_args()
    cfg  = load_config(str(Path(args.config).expanduser()))
    org  = cfg["organism"]["display"]
    paths= cfg["paths"]

    main_dir    = Path(paths["main"])
    classify_dir= main_dir / "data" / "02_classify"
    cons_dir    = main_dir / "data" / "04_alignment" / "conservation"
    acc_dir     = main_dir / "data" / "05_accessibility"
    sp_dir      = main_dir / "data" / "08_specificity"
    report_dir  = main_dir / "data" / "09_report"
    primers_dir = main_dir / "data" / "07_primers"
    out_dir     = main_dir / "data" / "12_figures"
    rep_dir     = main_dir / "reports" / "M12_figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    rep_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M12_figures", cfg)
    print_logo("M12 · Publication Figures (Nature style)", organism=org)

    save_versions(tools={},
                  python_pkgs=["matplotlib", "numpy", "pandas"],
                  report_dir=rep_dir, log=log)

    setup_style()

    log.info("─" * 56)
    log.info("Generating figures...")
    log.info("─" * 56)

    log.info("Fig 2 — Genome classification")
    fig_classification(classify_dir, out_dir, log)

    log.info("Fig 3 — Conservation analysis")
    fig_conservation(cons_dir, out_dir, log)

    log.info("Fig 4 — mRNA accessibility")
    fig_accessibility(acc_dir, out_dir, log)

    log.info("Fig 5 — crRNA ranking")
    fig_ranking(report_dir, out_dir, log)

    log.info("Fig 6 — Specificity matrix")
    fig_specificity(sp_dir, out_dir, log)

    log.info("Fig 7 — Amplicon map")
    fig_amplicon_map(primers_dir, out_dir, log)

    log.info("Fig 8 — Summary table")
    fig_summary_table(report_dir, out_dir, log)

    figs = sorted(out_dir.glob("*.png"))
    log.info("─" * 56)
    log.info(f"M12 complete. {len(figs)} figures generated → {out_dir}")
    for f in figs:
        log.info(f"  {f.name}")

    write_tsv([{"figure": f.name, "size_kb": f.stat().st_size // 1024}
               for f in figs],
              rep_dir / "summary.tsv", log=log)

    write_checkpoint("M12_figures", rep_dir, log=log)


if __name__ == "__main__":
    main()
