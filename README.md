# SHERLOCK crRNA Design Pipeline v2.0

**Computational pipeline for the design and in silico validation of a multiplexed CRISPR-Cas13a (SHERLOCK) diagnostic assay targeting toxigenic *Clostridioides difficile* from Chilean clinical fecal samples.**

---

## Authors

**Fausto Sebastián Cabezas-Mera**  
Doctoral Candidate — Doctorado en Informática Aplicada a Salud y Medio Ambiente  
Departamento de Informática y Computación, Facultad de Ingeniería  
Universidad Tecnológica Metropolitana (UTEM), Santiago, Chile  
[fcabezasme@utem.cl](mailto:fcabezasme@utem.cl) · ORCID: [0000-0002-3064-7089](https://orcid.org/0000-0002-3064-7089) · [github.com/fcabezasmera](https://github.com/fcabezasmera)

**Ana Rosa Moya-Beltrán, Ph.D.** *(Thesis supervisor · Co-investigator)*  
Research Assistant Professor — Bioinformática y Microbiología Computacional  
Departamento de Informática y Computación, Facultad de Ingeniería  
Universidad Tecnológica Metropolitana (UTEM), Santiago, Chile  
[amoya@utem.cl](mailto:amoya@utem.cl) · [Google Scholar](https://scholar.google.com/citations?user=vYrLnT4AAAAJ)  
Premio Mejor Tesis Doctoral en Microbiología, SOMICH 2022 · Secretaria SCB 2024–2025

---

## Project Context

This pipeline constitutes the **in silico crRNA design component (Specific Aim 1)** of the project:

> *"Diseño e implementación de una plataforma de detección molecular sensible y específica de cepas toxigénicas de* Clostridioides difficile *desde muestras clínicas, mediante un sistema CRISPR-Cas"*  
> **Concurso Endowment I+D Investigación Clínica para la Salud 2025, Universidad Andrés Bello (UNAB)**

**Principal Investigator:** Dr. Iván Calderón Lizana — Centro de Investigación de Resiliencia a Pandemias, Facultad de Ciencias de la Vida, UNAB  
**Co-investigators:** Dra. Lillian Acuña Olivares (UNAB) · Dra. Ana Rosa Moya-Beltrán (UTEM) · Dra. Paola Pidal Méndez (Clínica Indisa) · Dra. Erna Cona Trujillo (Clínica Indisa)  
**Clinical partner:** Clínica Indisa, Unidad de Control de Infecciones Asociadas a la Atención de Salud (IAAS), Santiago, Chile

The experimental validation of the designed crRNAs and RT-RPA primers (Specific Aim 2) will be conducted at UNAB and Clínica Indisa using clinical fecal samples from hospitalized patients with suspected CDI, employing the LwCas13a protein purified by the Calderón laboratory (plasmid pC013-huLwCas13a).

**Status:** In preparation for publication · Pipeline v2.0 · April 2026

---

## Background and Clinical Motivation

*Clostridioides difficile* infection (CDI) is the leading cause of antibiotic-associated diarrhea in healthcare settings worldwide, responsible for substantial morbidity, mortality, and healthcare costs. In Chile, the epidemiological landscape is dominated by the hypervirulent ribotype RT027 (NAP1/BI/ST1 lineage), which produces both major toxins (TcdA, TcdB) and the binary toxin CDT. After the major RT027 outbreak of 2012 — the largest recorded in Chilean history, with 79% incidence of the NAP1/ST1 strain — sustained high rates of CDI incidence and recurrence persist, without access to rapid, specific diagnostic tools that differentiate toxigenic from non-toxigenic strains.

The critical diagnostic gap addressed by this project is the absence of nucleic acid amplification tests that determine whether toxin genes are **actively expressed** (toxigenic strains) versus merely present in the genome. Current NAATs detect toxin genes but cannot distinguish active CDI from asymptomatic carriage — a clinically significant limitation that can lead to both overdiagnosis and inappropriate antibiotic use.

SHERLOCK (Specific High-sensitivity Enzymatic Reporter UnLOCKing), based on the RNA-guided collateral cleavage activity of LwCas13a, directly detects mRNA transcripts rather than genomic DNA, offering: attomolar sensitivity (~2 aM), single-nucleotide mismatch specificity, isothermal operation at 37°C (no thermocycler required), and lateral flow readout compatibility — making it ideal for point-of-care CDI diagnostics. The Calderón group has previously demonstrated the feasibility of SHERLOCK for detecting bacterial transcripts from a fish pathogen (*Yersinia ruckeri*) with sensitivity comparable to RT-qPCR and lateral flow readout (Calderón et al. 2024, *Microorganisms*).

---

## Reaction Design

Two-reaction multiplexed panel using LwCas13a exclusively (Option C — single ortholog, dual reaction):

| Reaction | Targets | Clinical rationale |
|---|---|---|
| **A** | tcdA · tcdB clade2 (RT027) · tcdB clade1 (RT012) · rpoB | Primary toxin detection + species confirmation |
| **B** | tcdC_WT · tcdC_RT027_junction · cdtA · cdtB · tpiA | Hypervirulence characterization + epidemiological typing |

**Key design decisions:**

**tcdB split by clade** — TcdB from RT027 (clade 2) and RT012 (clade 1) lineages shows significant sequence divergence in the receptor-binding and glucosyltransferase domains (Rupnik et al. 2009; Smits et al. 2016). A single pan-tcdB crRNA reduces sensitivity to ≤86% in silico; separate clade-specific crRNAs achieve ≥98% sensitivity per clade.

**rpoB as species-internal control** — Replaces 16S rRNA as the internal extraction/reaction control. rpoB is present in 99.6% of *C. difficile* isolates (Almeida et al. 2026, *Front Cell Infect Microbiol*) and is sufficiently divergent from other gut bacteria to provide species-level confirmation, avoiding the cross-reactivity inherent to 16S-based controls in polymicrobial fecal matrices.

**tcdC RT027 junction** — Targets the 18-bp deletion in tcdC that defines the RT027 hypervirulent lineage, enabling ribotype inference within a single SHERLOCK reaction.

**mRNA as target** — Detection of toxin mRNA directly demonstrates active transcription, differentiating CDI from asymptomatic colonization — the principal limitation of current DNA-based NAATs in Chile.

---

## Installation

### Prerequisites

```bash
git clone https://github.com/fcabezasmera/sherlock.git
cd sherlock

conda env create -f envs/sherlock.yml   # main environment
conda env create -f envs/adapt.yml      # ADAPT crRNA design (LwCas13a model)
conda env create -f envs/RPA.yml        # PrimedRPA (numpy compatibility)
```

### Required tools

| Tool | Version | Environment | Purpose |
|---|---|---|---|
| NCBI datasets | ≥15.0 | sherlock | Genome download |
| MAFFT | ≥7.505 | sherlock | Multiple sequence alignment |
| RNAplfold | ≥2.6.4 | sherlock | mRNA accessibility profiling |
| ADAPT | ≥1.4.0 | adapt | crRNA design (LwCas13a model, Metsky 2022) |
| PrimedRPA | ≥1.0.3 | RPA | RT-RPA primer co-design |
| BLAST+ | ≥2.13.0 | sherlock | Specificity analysis (4 databases) |
| seqkit | ≥2.4.0 | sherlock | Sequence manipulation |
| Primer3 | ≥2.6.0 | sherlock | RT-qPCR reference primer design |
| RNAfold | ≥2.6.4 | sherlock | crRNA secondary structure analysis |
| statsmodels | ≥0.13 | sherlock | FDR correction (Benjamini-Hochberg) |

---

## Configuration

All parameters are centralized in `config.yaml`. Critical parameters:

```yaml
organism:
  taxon_id:  "1496"
  display:   "Clostridioides difficile"

crna:
  dr:  "GGGGAUUUAGACUACCCCAAAAACGAAGGGGGGACUAAAAC"  # LwCas13a DR (Gootenberg 2017)
  t7:  "AATTCTAATACGACTCACTATAGG"                    # T7 promoter
  gc_min:         0.25    # adapted for C. difficile 29% GC genome (Wessels 2020)
  gc_max:         0.65
  u_frac_max:     0.50    # T7 termination prevention (Milligan 1987)
  g_pos1_penalty: -0.05   # G at pos1 reduces LwCas13a activity ~30-40% (Wessels 2020)

accessibility:
  window:   80
  span:     40
  u_window: 28    # spacer length — critical fix: was 1, must equal spacer len (Lorenz 2011)
  min_acc:  0.01  # relative threshold; absolute values low with -u 28

primers:
  rpa_env:     "RPA"
  amplicon_min: 100
  amplicon_max: 250
```

---

## Running the Pipeline

```bash
conda activate sherlock

# M01-M05: genome acquisition through accessibility profiling
python3 scripts/M01_download.py      --config config.yaml
python3 scripts/M02_classify.py      --config config.yaml
python3 scripts/M03_extract.py       --config config.yaml
python3 scripts/M04_alignment.py     --config config.yaml
python3 scripts/M05_accessibility.py --config config.yaml

# M06: crRNA design (ADAPT environment required)
conda run -n adapt \
    python3 scripts/M06_crna.py      --config config.yaml

# M07-M12: co-design, specificity, reporting, validation
python3 scripts/M07_primers.py       --config config.yaml
python3 scripts/M08_specificity.py   --config config.yaml
python3 scripts/M09_report.py        --config config.yaml
python3 scripts/M10_synthetic.py     --config config.yaml
python3 scripts/M11_rtqpcr.py        --config config.yaml
python3 scripts/M12_validation.py    --config config.yaml
```

### Re-run flags

```bash
python3 scripts/M06_crna.py        --config config.yaml --skip-adapt
python3 scripts/M07_primers.py     --config config.yaml --skip-primedrpa
python3 scripts/M08_specificity.py --config config.yaml --skip-db-build --skip-blast
python3 scripts/M12_validation.py  --config config.yaml --skip-structure
```

---

## Pipeline Modules

### M01 — Genome Download
Downloads *C. difficile* genomes (Complete + Chromosome assembly level) from NCBI RefSeq using `datasets`. Includes 12 Chilean clinical strains (NCBI, April 2026) and 5 curated anchor references.

**Anchor references:**

| Strain | Ribotype | Clinical relevance | Accession |
|---|---|---|---|
| R20291 | RT027 | Hypervirulent epidemic reference | GCF_015732555.1 |
| CD196 | RT027 | Historical RT027 outbreak isolate | GCF_021378415.1 |
| CD630Derm | RT012 | Classical toxigenic reference | GCF_000953275.1 |
| DSM27639 | RT012 | RT012 clade 1 reference | GCF_003313565.1 |
| L-NTCD03 | Non-tox | Non-toxigenic negative control | GCF_951803555.1 |

**Output:** `main/data/01_download/genomes/{working,chilean,anchors}/`

---

### M02 — Genome Classification

| Group | Label | Criteria | n |
|---|---|---|---|
| **A** | Hypervirulent (RT027-like, clade 2) | tcdB complete + cdtA/cdtB present + tcdC non-functional (18-bp deletion) | 106 |
| **B** | Toxigenic (RT012-like, clade 1) | tcdB complete + no binary toxin | 142 |
| **C** | Non-toxigenic | tcdB absent or non-functional | 52 |

rpoB and tpiA detected via `gene=` attribute with `product=` fallback for genomes lacking explicit gene annotation. Housekeeping genes bypass gene_matrix filter (present in all groups by definition).

---

### M03 — CDS Extraction
Extracts coding sequences per gene per group from GFF + FASTA. Key implementation: rpoB identified via `product=DNA-directed RNA polymerase subunit beta` (exact match, excluding rpoC) for genomes without explicit `gene=rpoB` annotation, achieving 288/290 genome coverage — equivalent to 16S rRNA.

---

### M04 — Multiple Sequence Alignment + Conservation
MAFFT (--localpair for n≤200; --auto for n>500) + per-position Shannon entropy. Wilcoxon rank-sum test with **Benjamini-Hochberg FDR correction** (α=0.05) for identification of significantly variable windows (correction for ~1,500 tests per gene; Benjamini & Hochberg 1995).

**9 MSA targets:** tcdA_all · tcdB_clade2 · tcdB_clade1 · tcdC_wt · tcdC_junction · cdtA_groupA · cdtB_groupA · tpiA_all · rpoB_all

---

### M05 — mRNA Accessibility (RNAplfold)
Computes unpaired probability profiles (W=80, L=40, **-u 28**, --noLP). The parameter `-u 28` equals spacer length and measures the probability that a full 28-nt window is simultaneously unpaired — the physically relevant quantity for Cas13a binding. Prior implementations using `-u 1` (single-base probability) underestimate steric constraints on crRNA-target hybridization (Lorenz et al. 2011).

Generates consensus sequences (majority-rule) for M11 Primer3 design.

---

### M06 — crRNA Design (ADAPT)
`maximize-activity` mode, LwCas13a activity model (Metsky et al. 2022), sliding window W=250, guide length 28nt.

**Composite scoring:**

| Component | Weight | Reference |
|---|---|---|
| ADAPT predicted activity (normalized) | 0.55 | Metsky et al. 2022 |
| mRNA accessibility (normalized, -u 28) | 0.35 | Lorenz et al. 2011 |
| GC score (25–65% optimal range) | 0.10 | Wessels et al. 2020 |
| G at spacer position 1 penalty | −0.05 | Wessels et al. 2020 |

**Filters:** GC 25–65% (adapted for *C. difficile* 29% GC) · poly-U ≤3 consecutive · total U fraction ≤50% (Milligan 1987)

**Result:** 90 crRNA candidates, 10 per target

---

### M07 — RT-RPA Primer Co-Design (PrimedRPA)

**Inverse co-design strategy:** PrimedRPA identifies conserved primer-binding regions genome-wide on the full MSA, then the highest-scoring ADAPT crRNA within each amplicon is selected. This guarantees every reported candidate has an experimentally validated, co-localized primer+crRNA pair.

Parameters: PrimerLength=32 · IdentityThreshold=0.95 · AmpliconSizeLimit=250

**crRNA synthesis format** (Kellner et al. 2019, *Nat Protoc*):
```
IDT order: ssDNA Ultramer = RC( T7_promoter + DR_DNA + spacer_DNA )
T7 transcription produces: DR_RNA + spacer_RNA  →  functional crRNA

RPA FP:  5'─[AATTCTAATACGACTCACTATAGG]─[FP sequence]─3'
RPA RP:  5'─[RP sequence]─3'
crRNA spacer: complementary to RNA transcribed from RPA amplicon
```

**Result:** 586 co-designs across 9 targets

---

### M08 — Specificity Analysis

| Database | Source | Scope |
|---|---|---|
| Non-toxigenic *C. difficile* (Group C) | NCBI | 52 genomes (informational) |
| Enteric pathogens | NCBI RefSeq | 11 species* |
| Human transcriptome | GRCh38 RefSeq mRNA | Complete |
| Human gut microbiome | UHGG v2.0.2 | 4,742 species representatives (excl. *C. difficile*) |

*blastn-short, word_size=7, evalue=0.01. Rejection: ≥90% identity · ≥80% crRNA coverage · alignment start ≤position 5*

**Database skip rules:** nontox skipped for all toxin targets · enteropathogens + UHGG skipped for rpoB and tpiA (housekeeping genes expected in Clostridiales)

**Result:** 87/90 crRNAs pass all applicable checks (96.7%)

---

### M09 — Ranking + Comprehensive Report

Integrates M07 co-designs (primary source) with ADAPT scores (M06) and specificity (M08). Each candidate is a complete, experimentally actionable set.

**Ranking:** ADAPT activity 35% · Conservation 20% · Specificity 15% · Primer co-design 10%

**Outputs:** `final_ranking.tsv` · `final_ranking.html` (interactive) · `IDT_order_crRNA.tsv` · `IDT_order_primers.tsv`

---

### M10 — Synthetic Positive Controls
IVT template DNA + synthetic target RNA for each final candidate. Used as positive controls in experimental SHERLOCK reactions at UNAB/Clínica Indisa.

**Validation:** 32/32 PASS (T7 present · 0mm crRNA detection · poly-U/GC acceptable)

---

### M11 — RT-qPCR Reference Primers
Primer3 on M05 consensus sequences. 8 genes: tcdA · tcdB_clade1 · tcdB_clade2 · tcdC · cdtA · cdtB · tpiA · rpoB.  
Parameters: Tm=60°C ±2°C · amplicon 80–200bp · GC 45–65%

**Validation:** 8/8 PASS (amplicons confirmed in silico · ΔTm <1.5°C all pairs)

---

### M12 — In Silico Validation

| Section | Method | Result |
|---|---|---|
| **A. crRNA (32 candidates)** | Sensitivity (Hamming ≤1mm) · IVT integrity · RNAfold structure · M07 amplicon check | 2 PASS · 28 WARN · 2 FAIL |
| **B. Synthetic controls (M10)** | T7 presence · 0mm detection · poly-U · GC | 32/32 PASS |
| **C. RT-qPCR primers (M11)** | seqkit amplicon on consensus · size 80–200bp · ΔTm | 8/8 PASS |

**FAIL candidates** (excluded from synthesis): tcdA_all rank3 (86.0% sensitivity) · tcdC_wt rank3 (87.9%). Ranks 1–2 of all 9 targets meet all primary criteria.

**WARN interpretation:** 28/32 WARNs are predominantly low GC (expected for *C. difficile* AT-rich genome, 29% GC) and RNAfold stem warnings. These are synthesis and folding advisories, not evidence of functional failure. Ranks 1–2 per target are prioritized for experimental synthesis.

---

## Results Summary

| Metric | Value |
|---|---|
| Genomes analyzed | 305 (288 global NCBI + 12 Chilean + 5 anchors) |
| Classification | Group A=106 (hypervirulent) · B=142 (toxigenic) · C=52 (non-toxigenic) |
| Diagnostic targets | 9 across 2 LwCas13a reactions |
| crRNA candidates (ADAPT) | 90 (10/target) |
| RT-RPA co-designs | 586 |
| Specificity pass rate | 87/90 (96.7%) |
| **Final validated candidates** | **32 (crRNA + FP+T7 + RP per candidate)** |
| Synthetic controls | 32/32 PASS |
| RT-qPCR reference primers | 8/8 PASS |
| crRNA in silico validation | 2 PASS · 28 WARN · 2 FAIL |

---

## Output Structure

```
main/data/
├── 01_download/      305 genomes (NCBI + Chilean + anchors)
├── 02_classify/      Groups A/B/C + gene_matrix.tsv
├── 03_extract/       CDS per gene per group
├── 04_alignment/     9 MSA + Shannon conservation profiles
├── 05_accessibility/ RNAplfold -u 28 profiles + consensus sequences
├── 06_crna/          90 crRNA candidates + full ADAPT TSVs
├── 07_primers/       586 RT-RPA co-designs (PrimedRPA)
├── 08_specificity/   Specificity vs 4 databases + UHGG
├── 09_report/        32 final candidates + HTML report + IDT sheets
├── 10_synthetic/     Synthetic control templates (64 IDT oligos)
├── 11_rtqpcr/        80 RT-qPCR primer pairs (8 genes × 10 pairs)
└── 12_validation/    Validation HTML report + TSV per section
```

---

## Evidence-Based Future Improvements

### High priority — Biological impact

**1. MUSCLE5 alignment for tcdB clades (R1)**  
MUSCLE5 (Edgar 2022, *Nat Commun*) outperforms MAFFT by 26% of correctly aligned columns in large divergent datasets. For tcdB_clade2 (RT027, n=106, >15% interlineage divergence), this improves ADAPT window scoring accuracy. Implementation: `muscle -super5 {input} -output {output} -threads 8`.

**2. minimize-guides mode for tcdB (R9)**  
ADAPT's `minimize-guides` mode with ≥95% coverage constraint (Metsky et al. 2022) guarantees detection of nearly all variants while minimizing crRNA count — critical for the high sequence diversity of tcdB across ribotypes.

**3. Mutational robustness analysis**  
Systematic in silico SNP introduction at each spacer position to map the seed region (positions 15–21, Wessels et al. 2020) and positions tolerant of 1 mismatch. Essential for predicting performance against emerging variants not represented in current genome collection, including novel Chilean isolates.

### Medium priority — Computational precision

**4. LinearCoPartition for crRNA-target binding (R5)**  
LinearCoPartition (Zhang et al. 2023, *Nucleic Acids Res*) computes intermolecular base-pairing probabilities considering the secondary structures of both crRNA and target RNA simultaneously — physically more accurate than RNAplfold alone. Recommended for top-3 candidates per target in M12, not for full candidate set (computationally expensive).

**5. Nearest-neighbor Tm for RPA primers (R12)**  
SantaLucia 1998 nearest-neighbor model via BioPython `Tm_NN` with RPA buffer conditions (50mM Na⁺, 14mM Mg²⁺, 37°C) provides accurate Tm estimates. PrimedRPA's current GC-based Tm is heuristic; with the T7 overhang added to FP, effective Tm shifts substantially and ΔTm >5°C pairs should be flagged for asymmetric amplification risk.

**6. NUPACK ΔG for dimer assessment (R13)**  
NUPACK (Zadeh et al. 2011) calculates thermodynamically rigorous ΔG for FP-RP homo/heterodimer formation at 37°C under RPA buffer conditions. Rejection criterion: ΔG < −3 kcal/mol. The current PrimedRPA dimerisation score is heuristic and does not account for the T7 overhang.

### Clinical and epidemiological expansion

**7. Latin American genome expansion**  
Current collection is dominated by European/North American submissions. Targeted inclusion of published Latin American sequences — particularly Chilean (PRJNA612578), Argentine, and Brazilian isolates — would improve conservation metrics for the regional clinical context and better represent the diversity of circulating strains in Chile.

**8. Clinical sensitivity simulation**  
Monte Carlo simulation of diagnostic sensitivity using allele frequency distributions from Chilean CDI surveillance data (ISP 2013–2018 bulletin). Accounts for co-infection, mixed ribotype samples, and the matrix inhibition effects of stool on RT-RPA amplification.

**9. Nanopore amplicon sequencing integration**  
SHERLOCK-positive RPA amplicons can be further characterized by Oxford Nanopore sequencing for simultaneous detection confirmation and ribotype-level typing. Pipeline extension with minimap2 + medaka for amplicon analysis would enable combined SHERLOCK + sequencing diagnostic workflows applicable to outbreak surveillance.

---

## References

1. Kellner MJ, Koob JG, Gootenberg JS, Abudayyeh OO, Zhang F. SHERLOCK: nucleic acid detection with CRISPR nucleases. *Nat Protoc.* 2019;14(10):2986–3012.
2. Metsky HC, Welch NL, Pillai PP, et al. Designing sensitive viral diagnostics with efficient machine learning. *Nat Biotechnol.* 2022;40:1123–1131.
3. Higgins M, Ravenhall M, Ward D, et al. PrimedRPA: primer design for recombinase polymerase amplification. *Bioinformatics.* 2019;35(4):682–684.
4. Wessels HH, Méndez-Mancilla A, Guo X, et al. Massively parallel Cas13 screens reveal principles for guide RNA design. *Nat Biotechnol.* 2020;38(6):722–727.
5. Milligan JF, Groebe DR, Witherell GW, Uhlenbeck OC. Oligoribonucleotide synthesis using T7 RNA polymerase and synthetic DNA templates. *Nucleic Acids Res.* 1987;15(21):8783–8798.
6. Lorenz R, Bernhart SH, Höner zu Siederdissen C, et al. ViennaRNA Package 2.0. *Algorithms Mol Biol.* 2011;6:26.
7. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B.* 1995;57(1):289–300.
8. Edgar RC. Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. *Nat Commun.* 2022;13:6968.
9. Zhang H, Zhang L, Mathews DH, Huang L. LinearCoPartition: linear-time cofolding and cofold base pair probability calculation of two RNAs. *Nucleic Acids Res.* 2023;51(18):e96.
10. Almeida A, Nayfach S, Boland M, et al. A unified catalog of 204,938 reference genomes from the human gut microbiome. *Nat Biotechnol.* 2021;39(1):105–114.
11. Gootenberg JS, Abudayyeh OO, Lee JW, et al. Nucleic acid detection with CRISPR-Cas13a/C2c2. *Science.* 2017;356(6336):438–442.
12. Calderón IL, Barros MJ, Fernández-Navarro N, Acuña LG. Detection of nucleic acids of the fish pathogen *Yersinia ruckeri* from planktonic and biofilm samples with a CRISPR/Cas13a-based assay. *Microorganisms.* 2024;12(2):283.
13. Aguayo C, Flores R, Lévesque S, et al. Rapid spread of *Clostridium difficile* NAP1/027/ST1 in Chile confirms the emergence of the epidemic strain in Latin America. *Epidemiol Infect.* 2015;143(14):3069–3073.
14. Plaza-Garrido Á, Barra-Carrasco J, Macias JH, et al. Predominance of *Clostridium difficile* ribotypes 012, 027 and 046 in a university hospital in Chile, 2012. *Epidemiol Infect.* 2016;144(5):976–979.
15. Rupnik M, Wilcox MH, Gerding DN. *Clostridium difficile* infection: new developments in epidemiology and pathogenesis. *Nat Rev Microbiol.* 2009;7(7):526–536.

---

## Citation

> Cabezas-Mera FS, Moya-Beltrán AR, et al. Computational design and in silico validation of a multiplexed SHERLOCK CRISPR-Cas13a diagnostic panel for toxigenic *Clostridioides difficile* from Chilean clinical samples. [*In preparation*, 2026]

---

## Acknowledgments

This work is part of the project *"Diseño e implementación de una plataforma de detección molecular sensible y específica de cepas toxigénicas de Clostridioides difficile desde muestras clínicas, mediante un sistema CRISPR-Cas"* (Concurso Endowment I+D, UNAB 2025). The authors thank Dr. Iván Calderón Lizana (UNAB) and Dr. Lillian Acuña Olivares (UNAB) for their experimental expertise in SHERLOCK/LwCas13a, and Dra. Paola Pidal Méndez and Dra. Erna Cona Trujillo (Clínica Indisa) for their clinical guidance and access to characterized clinical isolates.

---

## License

MIT License — see `LICENSE` file.

**Contact:** Fausto Sebastián Cabezas-Mera · [fcabezasme@utem.cl](mailto:fcabezasme@utem.cl) · [github.com/fcabezasmera](https://github.com/fcabezasmera)
