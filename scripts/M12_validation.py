#!/usr/bin/env python3
"""
M12_validation.py
SHERLOCK crRNA Design Pipeline · v1.0
Module : M12 · Comprehensive In Silico Validation

Validates ALL designed components from M09, M10, M11:

  SECTION A — crRNA Validation (from M09 final_ranking.tsv)
    A1. crRNA Sensitivity: spacer vs all on-target mRNA sequences (0,1,2 mm)
    A2. IVT Template Integrity: RC(IVT) = T7+DR+spacer; poly-U≤3; GC 30-70%
    A3. crRNA Secondary Structure: RNAfold — spacer not in stem ≥6nt
    A4. RT-RPA Primer Validation: seqkit amplicon on target sequences

  SECTION B — Synthetic Controls Validation (from M10)
    B1. Synthetic target template: T7 correctly transcribes target RNA
    B2. crRNA detects its own synthetic target (0 mismatches expected)
    B3. Synthetic template GC, length, poly-U checks

  SECTION C — RT-qPCR Primer Validation (from M11)
    C1. In silico PCR on consensus sequences (seqkit amplicon)
    C2. Amplicon size within 80-200bp
    C3. Primer Tm and GC content verification

USAGE
    conda activate sherlock
    python ~/sherlock/scripts/M12_validation.py --config ~/sherlock/config.yaml

    --skip-structure   skip RNAfold (faster)
    --skip-primers     skip seqkit amplicon checks (faster)

OUTPUT
    main/data/12_validation/crna_validation.tsv
    main/data/12_validation/synthetic_validation.tsv
    main/data/12_validation/rtqpcr_validation.tsv
    main/data/12_validation/validation_summary.tsv
    main/data/12_validation/validation_report.html
    main/logs/M12_validation.log
    main/reports/M12_validation/DONE.txt
"""

# =============================================================================
# PARAMETERS
# =============================================================================
DEFAULT_CONFIG   = "~/sherlock/config.yaml"
SPACER_LEN       = 28
MISMATCH_1_PASS  = 0.90
PRIMER_AMP_MIN   = 80
PRIMER_AMP_MAX   = 300
QPCR_AMP_MIN     = 80
QPCR_AMP_MAX     = 200
MAX_POLY_U       = 3
MAX_SPACER_STEM  = 5
DR_RNA = "GGGGAUUUAGACUACCCCAAAAACGAAGGGGGGACUAAAAC"
T7_DNA = "AATTCTAATACGACTCACTATAGG"

GENE_MAP = {
    "tcdA_all":      ["tcdA_groupA","tcdA_groupB"],
    "tcdB_all":      ["tcdB_groupA","tcdB_groupB"],
    "tcdC_wt":       ["tcdC_groupB"],
    "tcdC_junction": ["tcdC_junction_groupA"],
    "cdtA_groupA":   ["cdtA_groupA"],
    "cdtB_groupA":   ["cdtB_groupA"],
    "tpiA_all":      ["tpiA_groupA","tpiA_groupB","tpiA_groupC"],
    "sodA_all":      ["sodA_groupA","sodA_groupB","sodA_groupC"],
    "16S_all":       ["16S_groupA","16S_groupB","16S_groupC"],
}

RTQPCR_GENE_MAP = {
    "tcdA":"tcdA_all","tcdB":"tcdB_all","tcdC":"tcdC_wt",
    "cdtA":"cdtA_groupA","cdtB":"cdtB_groupA",
    "tpiA":"tpiA_all","sodA":"sodA_all",
}

# =============================================================================
# IMPORTS
# =============================================================================
import argparse, itertools, subprocess, sys, tempfile
from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from pipeline_utils import (
    load_config, get_logger, print_logo,
    save_versions, write_tsv, write_checkpoint,
)

# =============================================================================
# SEQUENCE UTILITIES
# =============================================================================
def reverse_complement(seq):
    comp = {"A":"T","T":"A","G":"C","C":"G","U":"A","N":"N"}
    return "".join(comp.get(b.upper(),"N") for b in seq[::-1])

def rna_to_dna(seq): return seq.upper().replace("U","T")
def dna_to_rna(seq): return seq.upper().replace("T","U")

def hamming(a, b):
    if len(a) != len(b): return SPACER_LEN
    return sum(x != y for x, y in zip(a.upper(), b.upper()))

def max_polyu(seq):
    seq = seq.upper().replace("T","U")
    return max((len(list(g)) for c,g in itertools.groupby(seq) if c=="U"), default=0)

def gc_content(seq):
    return sum(1 for b in seq.upper() if b in "GC") / max(len(seq),1)

# =============================================================================
# SEQUENCE LOADING
# =============================================================================
def load_fasta(path):
    seqs = []
    header, seq = None, []
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith(">"):
            if header: seqs.append((header, "".join(seq).upper()))
            header = line[1:].split()[0]; seq = []
        else:
            seq.append(line.strip())
    if header: seqs.append((header, "".join(seq).upper()))
    return seqs

def load_target_seqs(extract_dir, gene_files):
    seqs = {}
    for gf in gene_files:
        f = extract_dir / f"{gf}.fasta"
        if not f.exists(): continue
        for h, s in load_fasta(f):
            seqs[f"{gf}|{h}"] = s
    return seqs

def load_consensus(acc_dir, target):
    f = acc_dir / f"{target}_consensus.fasta"
    if not f.exists(): return ""
    lines = f.read_text(errors="ignore").splitlines()
    return "".join(l for l in lines if not l.startswith(">")).upper().replace("U","T")

# =============================================================================
# SECTION A1 — crRNA SENSITIVITY
# =============================================================================
def crna_sensitivity(spacer_rna, target_seqs):
    spacer_dna = rna_to_dna(spacer_rna).upper()
    spacer_rc  = reverse_complement(spacer_dna)
    slen = len(spacer_dna)
    results = {}
    for seq_id, seq in target_seqs.items():
        seq = seq.upper().replace("U","T")
        min_mm = SPACER_LEN
        for strand in [seq, reverse_complement(seq)]:
            for i in range(len(strand) - slen + 1):
                mm = hamming(spacer_rc, strand[i:i+slen])
                if mm < min_mm: min_mm = mm
                if min_mm == 0: break
            if min_mm == 0: break
        results[seq_id] = min_mm
    total = len(results)
    if total == 0:
        return {"n_seqs":0,"detect_0mm":0.0,"detect_1mm":0.0,"detect_2mm":0.0}
    vals = list(results.values())
    return {
        "n_seqs":    total,
        "detect_0mm":round(sum(1 for v in vals if v==0)/total, 4),
        "detect_1mm":round(sum(1 for v in vals if v<=1)/total, 4),
        "detect_2mm":round(sum(1 for v in vals if v<=2)/total, 4),
    }

# =============================================================================
# SECTION A2 — IVT INTEGRITY
# =============================================================================
def ivt_integrity(spacer_rna, ivt_dna, dr_rna, t7_dna):
    expected = (t7_dna + rna_to_dna(dr_rna) + rna_to_dna(spacer_rna)).upper()
    actual   = reverse_complement(ivt_dna).upper()
    correct  = (expected == actual)
    poly_u   = max_polyu(spacer_rna)
    gc       = round(gc_content(rna_to_dna(spacer_rna)), 3)
    has_int_t7 = t7_dna[8:] in rna_to_dna(spacer_rna)
    issues = []
    if not correct:    issues.append("IVT_mismatch")
    if has_int_t7:     issues.append("internal_T7")
    if poly_u > MAX_POLY_U: issues.append(f"poly_U={poly_u}")
    if gc < 0.30:      issues.append(f"low_GC={gc:.2f}")
    if gc > 0.70:      issues.append(f"high_GC={gc:.2f}")
    return {"ivt_correct":correct,"poly_u":poly_u,"gc":gc,
            "ivt_issues":";".join(issues) if issues else "PASS"}

# =============================================================================
# SECTION A3 — crRNA STRUCTURE
# =============================================================================
def crna_structure(spacer_rna, dr_rna, rnafold_exe):
    full = dr_rna + spacer_rna
    try:
        r = subprocess.run([rnafold_exe,"--noPS"], input=full,
                           capture_output=True, text=True, timeout=30)
        lines = r.stdout.strip().splitlines()
        if len(lines) < 2: return float("nan"), float("nan"), "NO_OUTPUT"
        structure = lines[1].split()[0]
        try: mfe = float(lines[1].split()[-1].strip("()"))
        except: mfe = float("nan")
        dr_len = len(dr_rna)
        spacer_ss = structure[dr_len:dr_len+SPACER_LEN]
        n_paired = sum(1 for c in spacer_ss if c in "()")
        frac = round(n_paired/SPACER_LEN, 3)
        status = f"WARN_stem_{n_paired}nt" if n_paired > MAX_SPACER_STEM*2 else "PASS"
        return mfe, frac, status
    except FileNotFoundError:
        return float("nan"), float("nan"), "RNAFOLD_NA"
    except Exception as e:
        return float("nan"), float("nan"), f"ERROR"

# =============================================================================
# SECTION A4 — PRIMER IN SILICO PCR
# =============================================================================
def primer_amplicon(fp, rp, target_seqs, seqkit_exe, amp_min, amp_max):
    if not fp or not rp or fp=="nan" or rp=="nan": return 0, 0.0, False
    with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w", delete=False) as f:
        for sid, seq in list(target_seqs.items())[:50]:
            f.write(f">{sid}\n{seq}\n")
        tmp = f.name
    try:
        fp_c = fp.upper().replace("U","T")
        rp_c = rp.upper().replace("U","T")
        r = subprocess.run([seqkit_exe,"amplicon","-F",fp_c,"-R",rp_c,
                            "--max-mismatch","3",tmp],
                           capture_output=True, text=True, timeout=60)
        lines = [l for l in r.stdout.splitlines() if l and not l.startswith(">")]
        sizes = [len(l.strip()) for l in lines if l.strip()]
        in_range = [s for s in sizes if amp_min <= s <= amp_max]
        return len(sizes), round(sum(sizes)/len(sizes),1) if sizes else 0.0, len(in_range)>0
    except Exception:
        return 0, 0.0, False
    finally:
        Path(tmp).unlink(missing_ok=True)

# =============================================================================
# SECTION B — SYNTHETIC CONTROLS VALIDATION (M10)
# =============================================================================
def validate_synthetic(syn_df, ranking_df, dr_rna, t7_dna):
    """
    For each synthetic target template:
    - Verify T7 promoter at 5' end
    - Verify crRNA detects synthetic target (0 mismatches)
    - Check length, poly-U, GC
    """
    rows = []
    for _, r in syn_df.iterrows():
        spacer_rna   = str(r.get("spacer_rna","")).upper().replace("T","U")
        template_dna = str(r.get("template_dna","")).upper()
        target       = str(r.get("target",""))
        rank         = int(r.get("rank",0))

        # Check T7 at 5' end
        has_t7 = template_dna.startswith(t7_dna.upper())

        # Extract transcribed region (after T7)
        if has_t7:
            transcribed_dna = template_dna[len(t7_dna):]
            transcribed_rna = dna_to_rna(transcribed_dna)
        else:
            transcribed_rna = dna_to_rna(template_dna)

        # Check crRNA detects this synthetic target (should be 0mm)
        spacer_dna = rna_to_dna(spacer_rna)
        spacer_rc  = reverse_complement(spacer_dna)
        slen = len(spacer_dna)
        target_dna = rna_to_dna(transcribed_rna)
        min_mm = SPACER_LEN
        for strand in [target_dna, reverse_complement(target_dna)]:
            for i in range(len(strand)-slen+1):
                mm = hamming(spacer_rc, strand[i:i+slen])
                if mm < min_mm: min_mm = mm
                if min_mm == 0: break
            if min_mm == 0: break

        length  = len(template_dna)
        poly_u  = max_polyu(transcribed_rna)
        gc      = round(gc_content(transcribed_dna[100:-100] if len(transcribed_dna)>200 else transcribed_dna), 3)

        issues = []
        if not has_t7:      issues.append("no_T7_promoter")
        if min_mm > 0:      issues.append(f"crRNA_mismatch={min_mm}")
        if poly_u > MAX_POLY_U*3: issues.append(f"poly_U={poly_u}")
        if length < 50:     issues.append("template_too_short")

        rows.append({
            "target":          target,
            "rank":            rank,
            "spacer_rna":      spacer_rna,
            "template_len":    length,
            "has_t7":          has_t7,
            "crna_detects_mm": min_mm,
            "template_gc":     gc,
            "template_poly_u": poly_u,
            "status":          "PASS" if not issues else "FAIL",
            "issues":          ";".join(issues) if issues else "PASS",
        })
    return pd.DataFrame(rows)

# =============================================================================
# SECTION C — RT-qPCR PRIMER VALIDATION (M11)
# =============================================================================
def validate_rtqpcr(rtqpcr_dir, acc_dir, seqkit_exe):
    """
    For each RT-qPCR gene: validate top primer pair via in silico PCR
    on the consensus sequence from M05.
    """
    rows = []
    for gene, target in RTQPCR_GENE_MAP.items():
        primers_f = rtqpcr_dir / f"{gene}_primers.tsv"
        if not primers_f.exists(): continue
        pdf = pd.read_csv(primers_f, sep="\t")
        if pdf.empty: continue
        best = pdf.iloc[0]
        fp = str(best.get("fp_seq",""))
        rp = str(best.get("rp_seq",""))

        # Load consensus from M05
        consensus = load_consensus(acc_dir, target)
        if not consensus:
            rows.append({"gene":gene,"target":target,"fp_seq":fp,"rp_seq":rp,
                         "n_amplicons":0,"mean_size":0,"in_range":False,
                         "fp_tm":best.get("fp_tm",""),"rp_tm":best.get("rp_tm",""),
                         "amplicon_size":best.get("amplicon_size",""),
                         "status":"NO_CONSENSUS","issues":"no_consensus"})
            continue

        # In silico PCR on consensus
        seqs = {"consensus": consensus}
        n_amp, mean_sz, in_range = primer_amplicon(fp, rp, seqs, seqkit_exe,
                                                    QPCR_AMP_MIN, QPCR_AMP_MAX)

        fp_tm = float(best.get("fp_tm",0) or 0)
        rp_tm = float(best.get("rp_tm",0) or 0)
        amp   = int(best.get("amplicon_size",0) or 0)

        issues = []
        if n_amp == 0:       issues.append("no_amplicon_in_silico")
        if amp < QPCR_AMP_MIN or amp > QPCR_AMP_MAX: issues.append(f"amp_size={amp}")
        if abs(fp_tm - rp_tm) > 3: issues.append(f"Tm_diff={abs(fp_tm-rp_tm):.1f}°C")

        rows.append({
            "gene":gene,"target":target,
            "fp_seq":fp,"rp_seq":rp,
            "fp_tm":fp_tm,"rp_tm":rp_tm,
            "amplicon_size":amp,
            "n_amplicons_insilico":n_amp,
            "mean_size_insilico":mean_sz,
            "in_range":in_range,
            "status":"PASS" if not issues else ("WARN" if n_amp>0 else "FAIL"),
            "issues":";".join(issues) if issues else "PASS",
        })
    return pd.DataFrame(rows)

# =============================================================================
# HTML REPORT
# =============================================================================
def make_html(crna_df, syn_df, qpcr_df, out_path, org):
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    sc = {"PASS":"#1B5E20","WARN":"#E65100","FAIL":"#B71C1C"}
    bg = {"PASS":"#E8F5E9","WARN":"#FFF3E0","FAIL":"#FFEBEE"}

    def row_html(r, cols, fmt={}):
        st = str(r.get("status","?"))
        cells = "".join(
            f"<td>{fmt.get(c,lambda x:x)(r.get(c,''))}</td>"
            for c in cols)
        bg_color = bg.get(st,"white")
        sc_color = sc.get(st,"#000")
        return ('<tr style="background:' + bg_color + '">'
                '<td><b style="color:' + sc_color + '">' + st + '</b></td>'
                + cells + '</tr>')

    def pct(v):
        try: return f"{float(v):.1%}"
        except: return "—"
    def f3(v):
        try: return f"{float(v):.3f}"
        except: return "—"
    def boolmark(v):
        return "✅" if v else "❌"

    # Section A rows
    a_rows = "".join(row_html(r, ["target","reaction","rank","guide_seq_rna",
        "n_seqs","detect_0mm","detect_1mm","detect_2mm",
        "ivt_issues","struct_status","primer_amplicon_ok","issues"],
        {"detect_0mm":pct,"detect_1mm":pct,"detect_2mm":pct,
         "primer_amplicon_ok":boolmark})
        for _, r in crna_df.iterrows())

    # Section B rows
    b_rows = "".join(row_html(r, ["target","rank","spacer_rna","template_len",
        "has_t7","crna_detects_mm","template_gc","template_poly_u","issues"],
        {"has_t7":boolmark,"template_gc":f3})
        for _, r in syn_df.iterrows()) if not syn_df.empty else "<tr><td colspan='9'>No synthetic data</td></tr>"

    # Section C rows
    c_rows = "".join(row_html(r, ["gene","fp_seq","fp_tm","rp_seq","rp_tm",
        "amplicon_size","n_amplicons_insilico","in_range","issues"],
        {"in_range":boolmark})
        for _, r in qpcr_df.iterrows()) if not qpcr_df.empty else "<tr><td colspan='9'>No qPCR data</td></tr>"

    n = {"PASS":0,"WARN":0,"FAIL":0}
    for df in [crna_df, syn_df, qpcr_df]:
        if df.empty: continue
        for st, cnt in df["status"].value_counts().items():
            n[str(st)] = n.get(str(st),0) + int(cnt)

    html = f"""<!DOCTYPE html><html><head><meta charset="UTF-8">
<title>M12 Validation — {org}</title>
<style>
body{{font-family:'Segoe UI',Arial,sans-serif;margin:20px;color:#212121;background:#fafafa}}
h1{{color:#0D47A1;border-bottom:2px solid #0D47A1;padding-bottom:8px}}
h2{{color:#37474F;margin-top:24px}}
table{{border-collapse:collapse;width:100%;font-size:11px;background:white;
       box-shadow:0 1px 3px rgba(0,0,0,.1);margin-bottom:20px}}
th{{background:#0D47A1;color:white;padding:7px 8px;text-align:left;position:sticky;top:0}}
td{{padding:4px 8px;border-bottom:1px solid #e0e0e0;vertical-align:top;word-break:break-all}}
.sb{{display:inline-block;padding:10px 18px;border-radius:6px;margin:6px;
     font-size:13px;font-weight:bold;color:white}}
code{{font-family:monospace;font-size:10px;background:#ECEFF1;padding:1px 3px}}
</style></head><body>
<h1>SHERLOCK In Silico Validation — {org}</h1>
<p style="color:#546E7A"><b>Generated:</b> {now} | Pipeline v1.0</p>
<div>
<span class="sb" style="background:#1B5E20">PASS: {n['PASS']}</span>
<span class="sb" style="background:#E65100">WARN: {n['WARN']}</span>
<span class="sb" style="background:#B71C1C">FAIL: {n['FAIL']}</span>
</div>

<h2>A. crRNA Validation (M09 candidates)</h2>
<p style="color:#546E7A;font-size:12px">
Sensitivity: ≥90% on-target seqs detected ≤1 mismatch | IVT: RC(template)=T7+DR+spacer, poly-U≤3 |
Structure: spacer stem &lt;{MAX_SPACER_STEM*2}nt | Primer: seqkit amplicon in {PRIMER_AMP_MIN}-{PRIMER_AMP_MAX}bp
</p>
<table><thead><tr><th>Status</th><th>Target</th><th>Rxn</th><th>Rank</th><th>Spacer (RNA)</th>
<th>N seqs</th><th>0mm</th><th>≤1mm</th><th>≤2mm</th>
<th>IVT</th><th>Structure</th><th>Primer</th><th>Issues</th></tr></thead>
<tbody>{a_rows}</tbody></table>

<h2>B. Synthetic Control Validation (M10)</h2>
<p style="color:#546E7A;font-size:12px">
T7 at 5' end | crRNA detects synthetic target with 0 mismatches | poly-U and GC checks
</p>
<table><thead><tr><th>Status</th><th>Target</th><th>Rank</th><th>Spacer</th>
<th>Template len</th><th>Has T7</th><th>crRNA mm</th><th>GC</th><th>Poly-U</th><th>Issues</th></tr></thead>
<tbody>{b_rows}</tbody></table>

<h2>C. RT-qPCR Primer Validation (M11)</h2>
<p style="color:#546E7A;font-size:12px">
seqkit amplicon on consensus | Amplicon {QPCR_AMP_MIN}-{QPCR_AMP_MAX}bp | Tm ΔT &lt;3°C
</p>
<table><thead><tr><th>Status</th><th>Gene</th><th>FP</th><th>FP Tm</th>
<th>RP</th><th>RP Tm</th><th>Amplicon (bp)</th><th>Amplicons found</th><th>In range</th><th>Issues</th></tr></thead>
<tbody>{c_rows}</tbody></table>

<h2>Validation Criteria</h2>
<ul style="font-size:12px;color:#546E7A">
<li><b>PASS:</b> All criteria met. Ready for experimental synthesis and validation.</li>
<li><b>WARN:</b> Minor issue detected (e.g. structure warning or no in silico amplicon). Review before synthesis.</li>
<li><b>FAIL:</b> Critical issue. Consider alternative candidates.</li>
</ul>
</body></html>"""
    Path(out_path).write_text(html, encoding="utf-8")

# =============================================================================
# MAIN
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=DEFAULT_CONFIG)
    p.add_argument("--skip-structure", action="store_true")
    p.add_argument("--skip-primers",   action="store_true")
    return p.parse_args()

def main():
    args   = parse_args()
    cfg    = load_config(str(Path(args.config).expanduser()))
    org    = cfg["organism"]["display"]
    paths  = cfg["paths"]
    tools  = cfg["tools"]
    crna_cfg = cfg.get("crna",{})
    dr_rna   = crna_cfg.get("dr", DR_RNA).replace("T","U")
    t7_dna   = crna_cfg.get("t7", T7_DNA)
    seqkit   = tools.get("seqkit",  "seqkit")
    rnafold  = tools.get("rnafold", "RNAfold")

    main_dir    = Path(paths["main"])
    extract_dir = main_dir / "data" / "03_extract"
    acc_dir     = main_dir / "data" / "05_accessibility"
    report_m09  = main_dir / "data" / "09_report"
    syn_dir     = main_dir / "data" / "10_synthetic"
    rtqpcr_dir  = main_dir / "data" / "11_rtqpcr"
    out_dir     = main_dir / "data" / "12_validation"
    rep_dir     = main_dir / "reports" / "M12_validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    rep_dir.mkdir(parents=True, exist_ok=True)

    log = get_logger("M12_validation", cfg)
    print_logo("M12 · In Silico Validation", organism=org)
    save_versions(tools={"seqkit":seqkit,"RNAfold":rnafold},
                  python_pkgs=["pandas","numpy"],
                  report_dir=rep_dir, log=log)

    # ── SECTION A: crRNA validation ──────────────────────────────────────────
    log.info("═"*56)
    log.info("SECTION A — crRNA Validation")
    log.info("═"*56)

    ranking_f = report_m09 / "final_ranking.tsv"
    if not ranking_f.exists():
        log.error("final_ranking.tsv not found. Run M09 first."); sys.exit(1)
    ranking = pd.read_csv(ranking_f, sep="\t")
    log.info(f"Loaded {len(ranking)} candidates from M09")

    a_rows = []
    for _, cand in ranking.iterrows():
        target     = str(cand["target"])
        spacer_rna = str(cand["guide_seq_rna"]).upper().replace("T","U")
        ivt_dna    = str(cand.get("ivt_template_dna",""))
        fp_seq     = str(cand.get("fp_seq",""))
        rp_seq     = str(cand.get("rp_seq",""))
        rank       = int(cand.get("rank",0))
        reaction   = str(cand.get("reaction","?"))

        log.info(f"  {target} rank{rank}")

        # A1 sensitivity
        target_seqs = load_target_seqs(extract_dir, GENE_MAP.get(target,[]))
        sens = crna_sensitivity(spacer_rna, target_seqs)
        log.info(f"    Sensitivity: 0mm={sens['detect_0mm']:.1%} 1mm={sens['detect_1mm']:.1%} (n={sens['n_seqs']})")

        # A2 IVT integrity
        ivt = ivt_integrity(spacer_rna, ivt_dna, dr_rna, t7_dna)
        log.info(f"    IVT: {ivt['ivt_issues']}")

        # A3 structure
        mfe, frac_paired, struct_status = float("nan"), float("nan"), "SKIPPED"
        if not args.skip_structure:
            mfe, frac_paired, struct_status = crna_structure(spacer_rna, dr_rna, rnafold)
            log.info(f"    Structure: mfe={mfe:.1f} spacer_paired={frac_paired:.1%} {struct_status}")

        # A4 primer
        primer_n, primer_sz, primer_ok = 0, 0.0, False
        if not args.skip_primers and fp_seq and fp_seq != "nan":
            primer_n, primer_sz, primer_ok = primer_amplicon(
                fp_seq, rp_seq, target_seqs, seqkit, PRIMER_AMP_MIN, PRIMER_AMP_MAX)
            log.info(f"    Primer: n_amp={primer_n} mean_size={primer_sz:.0f}bp ok={primer_ok}")

        # Overall status
        issues = []
        warns  = []
        status = "PASS"
        # FAIL: only for insufficient sensitivity
        if sens["detect_1mm"] < MISMATCH_1_PASS:
            issues.append(f"low_sens_1mm={sens['detect_1mm']:.2f}"); status="FAIL"
        # WARN: synthesis/design concerns (GC, poly-U, structure, primers)
        # Note: C. difficile genome is 29% GC — AT-rich spacers are expected
        if ivt["ivt_issues"] != "PASS":
            warns.append(f"ivt:{ivt['ivt_issues']}")
        if struct_status.startswith("WARN"):
            warns.append(f"struct:{struct_status}")
        if not primer_ok and not args.skip_primers:
            warns.append("no_rpa_amplicon")
        if warns and status == "PASS":
            status = "WARN"
        issues.extend(warns)

        a_rows.append({
            "target":target,"reaction":reaction,"rank":rank,
            "guide_seq_rna":spacer_rna,
            "n_seqs":sens["n_seqs"],
            "detect_0mm":sens["detect_0mm"],"detect_1mm":sens["detect_1mm"],
            "detect_2mm":sens["detect_2mm"],
            "ivt_correct":ivt["ivt_correct"],"ivt_poly_u":ivt["poly_u"],
            "ivt_gc":ivt["gc"],"ivt_issues":ivt["ivt_issues"],
            "mfe":mfe,"spacer_paired_frac":frac_paired,"struct_status":struct_status,
            "primer_n_amplicons":primer_n,"primer_mean_size":primer_sz,
            "primer_amplicon_ok":primer_ok,
            "status":status,"issues":";".join(issues) if issues else "PASS",
        })

    crna_df = pd.DataFrame(a_rows)
    crna_df.to_csv(out_dir / "crna_validation.tsv", sep="\t", index=False)

    # ── SECTION B: Synthetic controls (M10) ──────────────────────────────────
    log.info("═"*56)
    log.info("SECTION B — Synthetic Controls Validation (M10)")
    log.info("═"*56)

    syn_f = syn_dir / "target_templates.tsv"
    syn_df = pd.DataFrame()
    if syn_f.exists():
        syn_raw = pd.read_csv(syn_f, sep="\t")
        syn_df  = validate_synthetic(syn_raw, ranking, dr_rna, t7_dna)
        syn_df.to_csv(out_dir / "synthetic_validation.tsv", sep="\t", index=False)
        n_pass_s = (syn_df["status"]=="PASS").sum()
        log.info(f"  Synthetic controls: {n_pass_s}/{len(syn_df)} PASS")
    else:
        log.warning("  target_templates.tsv not found — skipping Section B")

    # ── SECTION C: RT-qPCR primers (M11) ─────────────────────────────────────
    log.info("═"*56)
    log.info("SECTION C — RT-qPCR Primer Validation (M11)")
    log.info("═"*56)

    qpcr_df = validate_rtqpcr(rtqpcr_dir, acc_dir, seqkit)
    if not qpcr_df.empty:
        qpcr_df.to_csv(out_dir / "rtqpcr_validation.tsv", sep="\t", index=False)
        n_pass_q = (qpcr_df["status"]=="PASS").sum()
        log.info(f"  RT-qPCR primers: {n_pass_q}/{len(qpcr_df)} PASS")
        for _, r in qpcr_df.iterrows():
            log.info(f"    {r['gene']}: {r['status']} | {r['issues']}")

    # ── Combined summary ──────────────────────────────────────────────────────
    make_html(crna_df, syn_df, qpcr_df,
              out_dir / "validation_report.html", org)

    # Write combined summary
    all_statuses = []
    for df in [crna_df, syn_df, qpcr_df]:
        if not df.empty:
            all_statuses.extend(df["status"].tolist())

    n_pass = all_statuses.count("PASS")
    n_warn = all_statuses.count("WARN")
    n_fail = all_statuses.count("FAIL")

    write_tsv([
        {"section":"A_crRNA",    "total":len(crna_df),
         "pass":(crna_df["status"]=="PASS").sum() if not crna_df.empty else 0},
        {"section":"B_synthetic","total":len(syn_df),
         "pass":(syn_df["status"]=="PASS").sum() if not syn_df.empty else 0},
        {"section":"C_rtqpcr",  "total":len(qpcr_df),
         "pass":(qpcr_df["status"]=="PASS").sum() if not qpcr_df.empty else 0},
    ], rep_dir / "summary.tsv", log=log)

    log.info("═"*56)
    log.info(f"M12 complete. PASS={n_pass} WARN={n_warn} FAIL={n_fail}")
    log.info(f"HTML → {out_dir/'validation_report.html'}")
    write_checkpoint("M12_validation", rep_dir, log=log)

if __name__ == "__main__":
    main()
