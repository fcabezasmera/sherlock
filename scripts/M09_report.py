#!/usr/bin/env python3
# M09_report.py - SHERLOCK Pipeline v1.0
# Uses M07 co-designs as primary source

DEFAULT_CONFIG = "~/sherlock/config.yaml"
TOP_N = 5
WEIGHTS = {"adapt":0.35,"accessibility":0.20,"conservation":0.20,"specificity":0.15,"primer":0.10}
TARGETS = [
    "tcdA_all", "tcdB_clade2", "tcdB_clade1",
    "tcdC_wt", "tcdC_junction",
    "cdtA_groupA", "cdtB_groupA",
    "tpiA_all", "rpoB_all",
]
REACTION_MAP = {
    "tcdA_all":"A","tcdB_clade2":"A","tcdB_clade1":"A","rpoB_all":"A",
    "tcdC_wt":"B","tcdC_junction":"B",
    "cdtA_groupA":"B","cdtB_groupA":"B","tpiA_all":"B",
}
TARGET_LABELS = {
    "tcdA_all":"tcdA",
    "tcdB_clade2":"tcdB (clade 2, RT027)",
    "tcdB_clade1":"tcdB (clade 1, RT012)",
    "tcdC_wt":"tcdC (WT)",
    "tcdC_junction":"tcdC (RT027 jct)",
    "cdtA_groupA":"cdtA",
    "cdtB_groupA":"cdtB",
    "tpiA_all":"tpiA",
    "rpoB_all":"rpoB"}
DR_RNA = "GGGGAUUUAGACUACCCCAAAAACGAAGGGGGGACUAAAAC"
T7_DNA = "AATTCTAATACGACTCACTATAGG"

import argparse, itertools, sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
from pipeline_utils import load_config, get_logger, print_logo, save_versions, write_tsv, write_checkpoint

def reverse_complement(seq):
    comp = {"A":"T","T":"A","G":"C","C":"G","U":"A","N":"N"}
    return "".join(comp.get(b.upper(),"N") for b in seq[::-1])

def rna_to_dna(seq):
    return seq.upper().replace("U","T")

def make_ivt_template(spacer_rna, dr_dna, t7):
    return reverse_complement(t7 + dr_dna + rna_to_dna(spacer_rna))

def make_full_crna(spacer_rna, dr_rna):
    return dr_rna + spacer_rna

def max_polyu(seq):
    seq = seq.upper().replace("T","U")
    return max((len(list(g)) for c,g in itertools.groupby(seq) if c=="U"), default=0)

def load_tsv(path):
    p = Path(path)
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()

def get_conservation_score(cons_df, crna_position):
    if cons_df.empty or "position" not in cons_df.columns:
        return float("nan")
    mask = ((cons_df["position"] >= crna_position + 1) &
            (cons_df["position"] <= crna_position + 28))
    subset = cons_df.loc[mask, "shannon"]
    if subset.empty:
        return float("nan")
    return max(0.0, 1.0 - float(subset.mean()) / 2.0)

def norm01(series):
    v = series.dropna()
    if v.empty or v.max() == v.min():
        return pd.Series(0.5, index=series.index)
    return (series - v.min()) / (v.max() - v.min())

def rank_target(target, codesign_df, adapt_df, spec_df, cons_df, t7, dr_dna, dr_rna, top_n):
    if codesign_df.empty:
        return pd.DataFrame()

    # Build ADAPT lookup from full TSV
    adapt_lookup = {}
    if not adapt_df.empty and "target-sequences" in adapt_df.columns:
        for _, r in adapt_df.iterrows():
            seq = str(r.get("target-sequences","")).upper().replace("U","T")
            try:
                pos_in_win = int(str(r.get("target-sequence-positions","0")).strip("{}").split(",")[0])
            except Exception:
                pos_in_win = 0
            if seq not in adapt_lookup:
                adapt_lookup[seq] = {
                    "adapt_activity":   float(r.get("guide-set-expected-activity", 0) or 0),
                    "adapt_frac_bound": float(r.get("total-frac-bound", 0) or 0),
                    "crna_position":    int(r.get("window-start", 0)) + pos_in_win,
                }

    # Build specificity lookup
    spec_lookup, spec_db_lookup = {}, {}
    if not spec_df.empty and "guide_seq" in spec_df.columns:
        for _, r in spec_df.iterrows():
            seq = str(r["guide_seq"]).upper().replace("U","T")
            spec_lookup[seq] = bool(r.get("overall_specific", True))
            spec_db_lookup[seq] = {
                "nontox": r.get("nontox_specific","NA"),
                "enteropathogens": r.get("enteropathogens_specific","NA"),
                "human_tx": r.get("human_tx_specific","NA"),
                "uhgg": r.get("uhgg_specific","NA"),
            }

    rows = []
    for _, row in codesign_df.iterrows():
        crna_seq   = str(row.get("crna_seq","")).upper().replace("U","T")
        spacer_rna = crna_seq.replace("T","U")
        crna_pos   = int(row.get("crna_position", 0))

        adapt_info  = adapt_lookup.get(crna_seq, {})
        adapt_act   = float(row.get("crna_activity", adapt_info.get("adapt_activity", 0) or 0))
        adapt_frac  = float(adapt_info.get("adapt_frac_bound", 0) or 0)
        is_specific = spec_lookup.get(crna_seq, True)
        spec_dbs    = spec_db_lookup.get(crna_seq, {})
        cons_score  = get_conservation_score(cons_df, crna_pos)
        gc          = sum(1 for b in crna_seq if b in "GC") / max(len(crna_seq), 1)

        rows.append({
            "target":target, "reaction":REACTION_MAP.get(target,"?"),
            "guide_seq_rna":spacer_rna,
            "full_crna_rna":make_full_crna(spacer_rna, dr_rna),
            "ivt_template_dna":make_ivt_template(spacer_rna, dr_dna, t7),
            "crna_position":crna_pos,
            "adapt_activity":adapt_act, "adapt_frac_bound":adapt_frac,
            "conservation":cons_score, "gc_content":round(gc,3),
            "max_polyu":max_polyu(crna_seq),
            "specific":is_specific,
            "nontox_specific":spec_dbs.get("nontox","NA"),
            "entero_specific":spec_dbs.get("enteropathogens","NA"),
            "human_specific":spec_dbs.get("human_tx","NA"),
            "uhgg_specific":spec_dbs.get("uhgg","NA"),
            "has_primer":True,
            "fp_seq":row.get("fp_seq",""), "fp_t7":row.get("fp_t7",""),
            "rp_seq":row.get("rp_seq",""),
            "amplicon_size":row.get("amplicon_size",""),
            "dimerisation_pct":row.get("dimerisation_pct",""),
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df["adapt_norm"] = norm01(df["adapt_activity"]).fillna(0.5)
    df["cons_norm"]  = norm01(df["conservation"]).fillna(0.5)
    df["spec_score"] = df["specific"].astype(float)
    df["final_score"] = (
        WEIGHTS["adapt"]        * df["adapt_norm"] +
        WEIGHTS["conservation"] * df["cons_norm"] +
        WEIGHTS["specificity"]  * df["spec_score"] +
        WEIGHTS["primer"]       * 1.0
    )
    df = (df.sort_values("final_score", ascending=False)
            .drop_duplicates(subset="guide_seq_rna")
            .head(top_n).reset_index(drop=True))
    df.insert(0, "rank", df.index + 1)
    return df

def make_idt_crna_sheet(df):
    rows = []
    for _, r in df.iterrows():
        seq = str(r["ivt_template_dna"])
        rows.append({
            "Name": f"{r['target']}_Rxn{r['reaction']}_rank{r['rank']}_crRNA_IVT",
            "Sequence": seq, "Length": len(seq),
            "Type": "Ultramer" if len(seq) <= 200 else "gBlock",
            "Scale":"25nm","Purification":"STD",
            "Notes": f"IVT template. T7 transcribes: {r['full_crna_rna']} | Spacer: {r['guide_seq_rna']} | pos={r['crna_position']}",
        })
    return pd.DataFrame(rows)

def make_idt_primer_sheet(df, rtqpcr_data):
    rows, seen = [], set()
    for _, r in df.iterrows():
        for ptype, col, note in [
            ("FP_T7","fp_t7",f"RPA FP+T7 | {r['target']} rank{r['rank']} | amp={r['amplicon_size']}bp"),
            ("RP","rp_seq",f"RPA RP | {r['target']} rank{r['rank']} | amp={r['amplicon_size']}bp"),
        ]:
            seq = str(r.get(col,"")).strip()
            if not seq or seq=="nan" or seq in seen:
                continue
            seen.add(seq)
            rows.append({"Name":f"{r['target']}_Rxn{r['reaction']}_rank{r['rank']}_{ptype}",
                "Sequence":seq,"Length":len(seq),
                "Type":"Ultramer" if len(seq)<=200 else "Standard",
                "Scale":"25nm","Purification":"STD","Notes":note})
    for gene, qdf in rtqpcr_data.items():
        if qdf.empty: continue
        best = qdf.iloc[0]
        for ptype, col in [("qPCR_FP","fp_seq"),("qPCR_RP","rp_seq")]:
            seq = str(best.get(col,"")).strip()
            if not seq or seq in seen: continue
            seen.add(seq)
            tm_k = "fp_tm" if "FP" in ptype else "rp_tm"
            rows.append({"Name":f"{gene}_{ptype}_rank1","Sequence":seq,"Length":len(seq),
                "Type":"Standard","Scale":"25nm","Purification":"STD",
                "Notes":f"RT-qPCR | {gene} | Tm={best.get(tm_k,'')}°C | amp={best.get('amplicon_size','')}bp"})
    return pd.DataFrame(rows)

def make_html(df, rtqpcr_data, cfg, out_path, org):
    dr_rna = cfg.get("crna",{}).get("dr",DR_RNA)
    t7     = cfg.get("crna",{}).get("t7",T7_DNA)
    now    = datetime.now().strftime("%Y-%m-%d %H:%M")
    rxn_colors = {"A":"#1565C0","B":"#2E7D32","?":"#757575"}

    rows_html = []
    for _, r in df.iterrows():
        react = r.get("reaction","?")
        color = rxn_colors.get(react,"#757575")
        spec  = "OK" if r.get("specific",True) else "FAIL"
        spec_c= "#1B5E20" if r.get("specific",True) else "#B71C1C"
        score = f"{r.get('final_score',0):.3f}"
        act   = f"{r.get('adapt_activity',0):.3f}"
        frac  = f"{r.get('adapt_frac_bound',0):.4f}"
        cons  = f"{r.get('conservation',0):.3f}" if not pd.isna(r.get('conservation',float('nan'))) else "n/a"
        amp   = str(r.get("amplicon_size","—"))
        dimer = str(r.get("dimerisation_pct","—"))
        def fmt_spec(v):
            if v is None or str(v) in ("nan","NaN","NA",""): return "skipped"
            return "OK" if str(v).lower() in ("true","1") else "FAIL"
        spec_detail = (f"nontox={fmt_spec(r.get('nontox_specific'))} | "
                       f"entero={fmt_spec(r.get('entero_specific'))} | "
                       f"human={fmt_spec(r.get('human_specific'))} | "
                       f"UHGG={fmt_spec(r.get('uhgg_specific'))}")
        rows_html.append(f"""
<tr class="mr" onclick="toggle(this)">
  <td>{r['rank']}</td>
  <td><span style="background:{color};color:white;padding:2px 8px;border-radius:3px;font-size:11px;font-weight:bold">{react}</span></td>
  <td><b>{TARGET_LABELS.get(r['target'],r['target'])}</b></td>
  <td><code style="font-family:monospace;font-size:10px;background:#ECEFF1;padding:1px 4px">{r['guide_seq_rna']}</code></td>
  <td>{r['crna_position']}</td><td>{act}</td><td>{frac}</td><td>{cons}</td>
  <td><b style="color:{spec_c}">{spec}</b></td>
  <td>{amp}</td><td>{dimer}%</td><td><b>{score}</b></td>
</tr>
<tr class="dr" style="display:none">
  <td colspan="12" style="background:#F5F5F5;padding:10px 16px;border-left:4px solid #0D47A1;font-size:11px">
    <b>Full crRNA (RNA 5'→3'):</b> <code style="font-family:monospace;font-size:10px;background:#ECEFF1;padding:1px 4px;word-break:break-all">{r['full_crna_rna']}</code><br>
    <b>IVT template DNA (order to IDT as ssDNA Ultramer):</b> <code style="font-family:monospace;font-size:10px;background:#ECEFF1;padding:1px 4px;word-break:break-all">{r['ivt_template_dna']}</code><br>
    <b>RPA FP (5'→3', T7 prepended):</b> <code style="font-family:monospace;font-size:10px;background:#ECEFF1;padding:1px 4px;word-break:break-all">{r.get('fp_t7','—')}</code><br>
    <b>RPA RP (5'→3'):</b> <code style="font-family:monospace;font-size:10px;background:#ECEFF1;padding:1px 4px;word-break:break-all">{r.get('rp_seq','—')}</code><br>
    <small style="color:#757575">Specificity: {spec_detail}</small>
  </td>
</tr>""")

    qpcr_rows = []
    for gene, qdf in rtqpcr_data.items():
        if qdf.empty: continue
        best = qdf.iloc[0]
        qpcr_rows.append(f"<tr><td><b>{gene}</b></td><td><code style='font-family:monospace;font-size:10px'>{best.get('fp_seq','')}</code></td><td>{best.get('fp_tm','')}°C</td><td><code style='font-family:monospace;font-size:10px'>{best.get('rp_seq','')}</code></td><td>{best.get('rp_tm','')}°C</td><td>{best.get('amplicon_size','')} bp</td></tr>")

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<title>SHERLOCK Report — {org}</title>
<style>
body{{font-family:'Segoe UI',Arial,sans-serif;margin:20px;color:#212121;background:#fafafa}}
h1{{color:#0D47A1;border-bottom:2px solid #0D47A1;padding-bottom:8px}}
h2{{color:#37474F;margin-top:24px}}
table{{border-collapse:collapse;width:100%;font-size:12px;background:white;box-shadow:0 1px 3px rgba(0,0,0,.1);margin-bottom:20px}}
th{{background:#0D47A1;color:white;padding:8px 10px;text-align:left;position:sticky;top:0;z-index:1}}
td{{padding:5px 10px;border-bottom:1px solid #e0e0e0;vertical-align:top}}
.mr{{cursor:pointer}}.mr:hover td{{background:#E3F2FD}}
.ib{{background:#E8F5E9;border-left:4px solid #388E3C;padding:12px 16px;margin:15px 0;border-radius:4px;font-size:13px}}
.wb{{background:#FFF8E1;border-left:4px solid #F9A825;padding:12px 16px;margin:15px 0;border-radius:4px;font-size:13px}}
.pm{{font-family:monospace;font-size:12px;background:#263238;color:#ECEFF1;padding:12px;border-radius:4px;overflow-x:auto}}
input{{padding:6px 10px;border:1px solid #bbb;border-radius:4px;width:280px;font-size:13px}}
code{{font-family:'Courier New',monospace;font-size:10px;background:#ECEFF1;padding:1px 4px;border-radius:2px}}
ol,ul{{font-size:13px;color:#546E7A}}li{{margin:4px 0}}
</style></head><body>
<h1>SHERLOCK crRNA Design Report — {org}</h1>
<p style="color:#546E7A"><b>Generated:</b> {now} &nbsp;|&nbsp; <b>Pipeline v1.0</b></p>
<div class="ib">
<b>Reaction A</b> (LwCas13a): tcdA + tcdB clade2 + tcdB clade1 + rpoB (species control)<br>
    <b>Reaction B</b> (LwCas13a): tcdC_WT + tcdC_RT027_junction + cdtA + cdtB + tpiA<br><br>
Each row = co-designed set: RT-RPA primer pair + crRNA validated to co-localize within the same amplicon.<br>
Ranking weights: ADAPT activity 35% · Conservation 20% · Specificity 15% · Primer co-design 10%
</div>
<div class="wb">
<b>crRNA synthesis (Kellner et al. 2019, Nat Protoc):</b><br>
Order IVT template to IDT as ssDNA Ultramer: <code>RC( T7_promoter + DR_DNA + spacer_DNA )</code><br>
T7 transcription produces functional crRNA: <code>DR_RNA + spacer_RNA</code><br>
<b>LwCas13a DR (RNA):</b> <code>{dr_rna}</code><br>
<b>T7 promoter (DNA):</b> <code>{t7}</code><br>
<b>RPA FP:</b> T7 promoter prepended to 5' end of forward primer (sense strand).<br>
<b>crRNA spacer:</b> complementary to transcribed RNA from RPA amplicon.
</div>
<h2>Final crRNA Candidates — click row to expand sequences</h2>
<div style="margin:10px 0"><input type="text" id="fi" onkeyup="filt()" placeholder="Search sequence or target..."></div>
<table id="rt">
<thead><tr><th>#</th><th>Rxn</th><th>Target</th><th>Spacer (RNA 5'→3')</th><th>Pos</th><th>ADAPT Activity</th><th>Frac Bound</th><th>Conservation</th><th>Specific</th><th>Amplicon (bp)</th><th>Dimer %</th><th>Score</th></tr></thead>
<tbody>{"".join(rows_html)}</tbody></table>
<h2>RT-qPCR Reference Primers (Primer3, Tm=60°C)</h2>
<table><thead><tr><th>Gene</th><th>Forward Primer (5'→3')</th><th>FP Tm</th><th>Reverse Primer (5'→3')</th><th>RP Tm</th><th>Amplicon</th></tr></thead>
<tbody>{"".join(qpcr_rows)}</tbody></table>
<h2>Specificity Databases</h2>
<ul>
<li><b>Non-toxigenic C. difficile</b> (Group C, n=52 genomes, NCBI)</li>
<li><b>Enteric pathogens</b> (n=11 species: C. sordellii, C. perfringens, C. botulinum, E. faecalis, E. faecium, K. pneumoniae, E. coli O157, S. enterica, C. jejuni, L. monocytogenes, S. aureus)</li>
<li><b>Human transcriptome</b> (GRCh38 RefSeq mRNA, NCBI)</li>
<li><b>Human Gut Microbiome</b> (UHGG v2.0.2, n=4,742 species representatives, excl. C. difficile)</li>
</ul>
<p style="color:#757575">Rejection: ≥90% identity · ≥80% crRNA coverage · alignment start ≤pos 5</p>
<h2>Pipeline Parameters</h2>
<div class="pm">ADAPT     : sliding-window, maximize-activity, W=250, gl=28, LwCas13a model (Metsky 2022)
MAFFT     : --localpair --maxiterate 1000 (n≤200) / --auto (n&gt;500)
RNAplfold : W=80, L=40, u=1, --noLP
PrimedRPA : IdentityThreshold=0.95, PrimerLength=32, AmpliconSizeLimit=250, MSA mode
BLAST     : blastn-short, word_size=7, evalue=0.01
Filters   : GC 30-70%, poly-U ≤3, accessibility ≥0.5
Genomes   : Complete+Chromosome (n=288 global) + 12 Chilean strains (NCBI Apr 2026)</div>
<h2>Key References</h2>
<ol>
<li>Kellner MJ et al. SHERLOCK: nucleic acid detection with CRISPR nucleases. <i>Nat Protoc</i> 2019.</li>
<li>Metsky HC et al. Designing sensitive viral diagnostics with efficient machine learning. <i>Nat Biotechnol</i> 2022.</li>
<li>Higgins M et al. PrimedRPA. <i>Bioinformatics</i> 2019.</li>
<li>Almeida A et al. A unified catalog of 204,938 reference genomes from the human gut microbiome. <i>Nat Biotechnol</i> 2021.</li>
</ol>
<script>
function toggle(r){{var n=r.nextElementSibling;if(n&&n.classList.contains("dr"))n.style.display=n.style.display==="none"?"":"none";}}
function filt(){{var q=document.getElementById("fi").value.toLowerCase();var rows=document.querySelectorAll("#rt tbody tr");var i=0;while(i<rows.length){{var mr=rows[i],dr=rows[i+1];if(mr.classList.contains("dr")){{i++;continue;}}var show=mr.textContent.toLowerCase().includes(q);mr.style.display=show?"":"none";if(dr)dr.style.display=show?"":"none";i+=2;}}}}
</script>
</body></html>"""
    Path(out_path).write_text(html, encoding="utf-8")

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=DEFAULT_CONFIG)
    return p.parse_args()

def main():
    args = parse_args()
    cfg  = load_config(str(Path(args.config).expanduser()))
    org  = cfg["organism"]["display"]
    paths= cfg["paths"]
    crna_cfg = cfg.get("crna", {})
    rank_cfg = cfg.get("ranking", {})
    top_n   = rank_cfg.get("top_n", TOP_N)
    dr_rna  = crna_cfg.get("dr", DR_RNA).replace("T","U")
    dr_dna  = rna_to_dna(dr_rna)
    t7      = crna_cfg.get("t7", T7_DNA)

    main_dir    = Path(paths["main"])
    crna_dir    = main_dir / "data" / "06_crna"
    primers_dir = main_dir / "data" / "07_primers"
    sp_dir      = main_dir / "data" / "08_specificity"
    cons_dir    = main_dir / "data" / "04_alignment" / "conservation"
    rtqpcr_dir  = main_dir / "data" / "11_rtqpcr"
    out_dir     = main_dir / "data" / "09_report"
    report_dir  = main_dir / "reports" / "M09_report"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M09_report", cfg)
    print_logo("M09 · Ranking + Report", organism=org)
    save_versions(tools={}, python_pkgs=["pandas","numpy"], report_dir=report_dir, log=log)

    qpcr_genes  = cfg.get("rtqpcr",{}).get("targets",["tcdA","tcdB_clade1","tcdB_clade2","tcdC","cdtA","cdtB","tpiA","rpoB"])
    rtqpcr_data = {g: load_tsv(rtqpcr_dir / f"{g}_primers.tsv") for g in qpcr_genes}

    all_ranked, summary_rows = [], []
    for target in TARGETS:
        log.info(f"  {target}")
        codesign_df = load_tsv(primers_dir / f"{target}_codesign.tsv")
        adapt_df    = load_tsv(crna_dir / f"{target}_adapt.tsv")
        spec_df     = load_tsv(sp_dir / f"{target}_specificity.tsv")
        cons_df     = load_tsv(cons_dir / f"{target}_conservation.tsv")

        ranked = rank_target(target, codesign_df, adapt_df, spec_df, cons_df, t7, dr_dna, dr_rna, top_n)

        if not ranked.empty:
            ranked.to_csv(out_dir / f"{target}_ranked.tsv", sep="\t", index=False)
            all_ranked.append(ranked)
            r0 = ranked.iloc[0]
            log.info(f"    {len(ranked)} candidates | score={r0['final_score']:.3f} | FP={str(r0['fp_seq'])[:20]}...")

        summary_rows.append({"target":target,"reaction":REACTION_MAP.get(target,"?"),
            "n_ranked":len(ranked),
            "top1_score":round(float(ranked.iloc[0]["final_score"]),3) if not ranked.empty else 0,
            "top1_seq":ranked.iloc[0]["guide_seq_rna"] if not ranked.empty else ""})

    if not all_ranked:
        log.error("No candidates."); sys.exit(1)

    combined = pd.concat(all_ranked, ignore_index=True)
    combined.to_csv(out_dir / "final_ranking.tsv", sep="\t", index=False)
    make_html(combined, rtqpcr_data, cfg, out_dir / "final_ranking.html", org)
    make_idt_crna_sheet(combined).to_csv(out_dir / "IDT_order_crRNA.tsv", sep="\t", index=False)
    make_idt_primer_sheet(combined, rtqpcr_data).to_csv(out_dir / "IDT_order_primers.tsv", sep="\t", index=False)

    log.info(f"TSV → {out_dir/'final_ranking.tsv'}  ({len(combined)} candidates)")
    log.info(f"HTML → {out_dir/'final_ranking.html'}")
    write_tsv(summary_rows, report_dir / "summary.tsv", log=log)
    log.info("─"*56)
    log.info("M09 complete.")
    for r in summary_rows:
        log.info(f"  {r['target']:<20} Rxn={r['reaction']} n={r['n_ranked']} score={r['top1_score']} {r['top1_seq']}")
    write_checkpoint("M09_report", report_dir, log=log)

if __name__ == "__main__":
    main()
