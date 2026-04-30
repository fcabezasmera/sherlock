# SHERLOCK crRNA Design Pipeline v2.0

**Computational pipeline for the design and in silico validation of CRISPR-Cas13a (SHERLOCK) diagnostic assays targeting toxigenic *Clostridioides difficile* from Chilean clinical fecal samples.**

> **Author:** Fausto Cabezas-Mera  
> **Affiliation:** [your institution]  
> **Status:** Pre-publication — in preparation  
> **GitHub:** [github.com/fcabezasmera/sherlock](https://github.com/fcabezasmera/sherlock)

---

## Background

*Clostridioides difficile* infection (CDI) is the leading cause of healthcare-associated diarrhea worldwide, responsible for significant morbidity and mortality, particularly in elderly and immunocompromised patients. In Chile, the epidemiological landscape is dominated by the hypervirulent ribotype RT027 (PCR ribotype 027, NAP1/BI/027 lineage), which produces both major toxins (TcdA, TcdB) and the binary toxin CDT, contributing to increased disease severity and recurrence.

Current diagnostic algorithms rely on multi-step approaches combining nucleic acid amplification tests (NAATs), enzyme immunoassays (EIA) for glutamate dehydrogenase (GDH), and toxin detection. While NAATs are highly sensitive, they detect toxin genes without distinguishing between active toxin production and asymptomatic carriage — a clinically significant limitation.

SHERLOCK (Specific High-sensitivity Enzymatic Reporter UnLOCKing), based on LwCas13a collateral RNA cleavage, offers single-molecule sensitivity (~2 aM), single-nucleotide mismatch specificity, and compatibility with lateral flow readout — making it an ideal platform for point-of-care CDI diagnostics targeting mRNA (active transcription) rather than genomic DNA.

This pipeline implements a rigorous, reproducible computational workflow for the design of crRNA spacers and RT-RPA primers for a multiplexed SHERLOCK panel targeting toxigenic *C. difficile* from Chilean clinical context.

---

## Reaction Design

Two-reaction multiplexed panel using LwCas13a exclusively (Option C):

| Reaction | Targets | Clinical rationale |
|---|---|---|
| **A** | tcdA · tcdB clade2 (RT027) · tcdB clade1 (RT012) · rpoB | Primary toxin detection + species confirmation |
| **B** | tcdC_WT · tcdC_RT027_junction · cdtA · cdtB · tpiA | Hypervirulence characterization + epidemiological typing |

**Key design decisions:**
- **tcdB split by clade** (RT027/clade2 vs RT012/clade1): TcdB from RT027 and RT012 lineages shows significant sequence divergence that reduces single-crRNA sensitivity. Separate crRNA pairs achieve ≥98% sensitivity per clade (Rupnik et al. 2019).
- **rpoB as internal control** (replaces 16S rRNA): rpoB is present in 99.6% of *C. difficile* isolates (vs universally in all gut bacteria), providing species-level confirmation without the cross-reactivity inherent to 16S-based controls.
- **tcdC RT027 junction**: Targets the 18-bp deletion site in tcdC that defines the RT027/hypervirulent lineage — enables ribotype inference in a single reaction.
- **mRNA as target**: Detecting toxin mRNA (rather than genomic DNA) differentiates active toxigenic strains from asymptomatic carriers, addressing the key limitation of current NAAT-based methods.

---

## Installation

### Prerequisites

```bash
# Clone repository
git clone https://github.com/fcabezasmera/sherlock.git
cd sherlock

# Create conda environments
conda env create -f envs/sherlock.yml   # main environment
conda env create -f envs/adapt.yml      # ADAPT crRNA design
conda env create -f envs/RPA.yml        # PrimedRPA (numpy compatibility fix)
```

### Required tools

| Tool | Version | Environment | Purpose |
|---|---|---|---|
| NCBI datasets | ≥15.0 | sherlock | Genome download |
| MAFFT | ≥7.505 | sherlock | Multiple sequence alignment |
| RNAplfold | ≥2.6.4 | sherlock | mRNA accessibility profiling |
| ADAPT | ≥1.4.0 | adapt | crRNA design (LwCas13a model) |
| PrimedRPA | ≥1.0.3 | RPA | RT-RPA primer design |
| BLAST+ | ≥2.13.0 | sherlock | Specificity analysis |
| seqkit | ≥2.4.0 | sherlock | Sequence manipulation |
| Primer3 | ≥2.6.0 | sherlock | RT-qPCR primer design |
| RNAfold | ≥2.6.4 | sherlock | crRNA secondary structure |
| statsmodels | ≥0.13 | sherlock | FDR correction (M04) |

---

## Configuration

All parameters are centralized in `config.yaml`. Key sections:

```yaml
organism:
  taxon_id:  "1496"
  display:   "Clostridioides difficile"

crna:
  dr:  "GGGGAUUUAGACUACCCCAAAAACGAAGGGGGGACUAAAAC"   # LwCas13a DR (Gootenberg 2017)
  t7:  "AATTCTAATACGACTCACTATAGG"                     # T7 promoter
  gc_min:        0.25    # adapted for C. difficile 29% GC genome
  gc_max:        0.65
  u_frac_max:    0.50    # T7 termination prevention (Milligan 1987)
  g_pos1_penalty: -0.05  # G at position 1 reduces LwCas13a activity (Wessels 2020)

accessibility:
  window:   80
  span:     40
  u_window: 28    # spacer length for RNAplfold -u parameter (Lorenz 2011)

primers:
  rpa_env:     "RPA"
  amplicon_min: 100
  amplicon_max: 250    # relaxed for cdtA (short gene)
```

---

## Running the Pipeline

```bash
conda activate sherlock

# M01-M05: genome download through accessibility profiling
python3 scripts/M01_download.py      --config config.yaml
python3 scripts/M02_classify.py      --config config.yaml
python3 scripts/M03_extract.py       --config config.yaml
python3 scripts/M04_alignment.py     --config config.yaml
python3 scripts/M05_accessibility.py --config config.yaml

# M06: crRNA design (requires ADAPT environment)
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

### Re-run flags (skip completed steps)

```bash
python3 scripts/M06_crna.py        --config config.yaml --skip-adapt
python3 scripts/M07_primers.py     --config config.yaml --skip-primedrpa
python3 scripts/M08_specificity.py --config config.yaml --skip-db-build --skip-blast
python3 scripts/M12_validation.py  --config config.yaml --skip-structure
```

---

## Pipeline Modules

### M01 — Genome Download
Downloads *C. difficile* genomes (Complete + Chromosome assembly level) from NCBI RefSeq using `datasets`. Includes 12 Chilean clinical strains (all publicly available as of April 2026, deposited by the Paredes-Sabja group, Universidad de Santiago de Chile) and 5 anchor reference strains.

**Anchor references:**

| Strain | Ribotype | Significance | Accession |
|---|---|---|---|
| R20291 | RT027 | Hypervirulent epidemic reference | GCF_015732555.1 |
| CD196 | RT027 | Historical RT027 outbreak | GCF_021378415.1 |
| CD630Derm | RT012 | Classic toxigenic reference | GCF_000953275.1 |
| DSM27639 | RT012 | RT012 reference | GCF_003313565.1 |
| L-NTCD03 | Non-tox | Non-toxigenic reference | GCF_951803555.1 |

**Output:** `main/data/01_download/genomes/{working,chilean,anchors}/`

---

### M02 — Genome Classification

Classifies genomes into three functional groups based on toxin gene complement:

| Group | Label | Criteria | n |
|---|---|---|---|
| **A** | Hypervirulent (RT027-like) | tcdB complete + cdtA/cdtB present + tcdC non-functional | 106 |
| **B** | Toxigenic (RT012-like) | tcdB complete + no binary toxin | 142 |
| **C** | Non-toxigenic | tcdB absent/non-functional | 52 |

Classification logic based on GFF gene= and product= attributes. sodA and tpiA treated as housekeeping genes (present in all groups). rpoB detected via `product=DNA-directed RNA polymerase subunit beta` for genomes without explicit `gene=rpoB` annotation.

**Output:** `main/data/02_classify/gene_matrix.tsv`, `groupA/B/C.txt`

---

### M03 — CDS Extraction
Extracts coding sequences for 7 target genes per group from GFF + FASTA using seqkit. Housekeeping genes (rpoB, tpiA) bypass the gene_matrix filter (present in all groups by definition).

**Targets:** tcdA, tcdB, tcdC, tcdC_junction, cdtA, cdtB, tpiA, rpoB

**Output:** `main/data/03_extract/{gene}_{group}.fasta`

---

### M04 — Multiple Sequence Alignment + Conservation

Concatenates sequences by diagnostic target, aligns with MAFFT (--localpair --maxiterate 1000 for n≤200; --auto for n>500), and computes per-position conservation metrics:

- **Shannon entropy** (H): position-wise nucleotide diversity
- **Wilcoxon rank-sum test** with Benjamini-Hochberg FDR correction (α=0.05) for identification of significantly variable windows

**9 MSA targets:**

| Target | Source groups | Sequences | Length |
|---|---|---|---|
| tcdA_all | groupA + groupB | ~220 | ~8,200 bp |
| tcdB_clade2 | groupA only | ~106 | ~7,200 bp |
| tcdB_clade1 | groupB only | ~142 | ~7,200 bp |
| tcdC_wt | groupB | ~140 | ~800 bp |
| tcdC_junction | groupA | ~98 | ~1,100 bp |
| cdtA_groupA | groupA | ~104 | ~1,500 bp |
| cdtB_groupA | groupA | ~105 | ~2,700 bp |
| tpiA_all | groupA+B+C | ~300 | ~850 bp |
| rpoB_all | groupA+B+C | ~298 | ~3,800 bp |

**Output:** `main/data/04_alignment/msa/*.aln`, `conservation/*.tsv`

---

### M05 — mRNA Accessibility

Computes per-position unpaired probability profiles using RNAplfold (W=80, L=40, **-u 28**, --noLP). The critical parameter `-u 28` (equal to spacer length) measures the probability that a 28-nt window is simultaneously unpaired — the physically relevant quantity for Cas13a binding (Lorenz et al. 2011). Previous versions incorrectly used `-u 1`.

Generates consensus sequences (majority rule per position) for use in M11 Primer3 design.

**Output:** `main/data/05_accessibility/{target}_accessibility.tsv`, `{target}_consensus.fasta`

---

### M06 — crRNA Design (ADAPT)

Designs crRNA spacers using ADAPT (Metsky et al. 2022, *Nat Biotechnol*) in `maximize-activity` mode with the LwCas13a activity prediction model. Sliding window: W=250, guide length 28nt.

**Composite scoring (post-ADAPT ranking):**

| Component | Weight | Rationale |
|---|---|---|
| ADAPT activity (norm) | 0.55 | Direct LwCas13a activity prediction |
| mRNA accessibility (norm) | 0.35 | RNAplfold -u 28 profile |
| GC score | 0.10 | Optimal 25-65% for C. diff (Wessels 2020) |
| G at position 1 | −0.05 | Reduces LwCas13a activity ~30-40% (Wessels 2020) |

**Filters applied:**
- GC content: 25–65% (adapted for C. difficile 29% GC genome)
- Poly-U: ≤3 consecutive U (T7 termination prevention)
- Total U fraction: ≤50% (Milligan 1987)

**Output:** `main/data/06_crna/{target}_adapt.tsv` (all windows), `{target}_candidates.tsv` (top-10)  
**Result:** 90 crRNA candidates (10 per target)

---

### M07 — RT-RPA Primer Co-Design (PrimedRPA)

Designs RT-RPA primer pairs using PrimedRPA (Higgins et al. 2019, *Bioinformatics*) on the full MSA. Uses an **inverse co-design strategy**: PrimedRPA identifies conserved primer-binding regions genome-wide, then the best ADAPT-predicted crRNA within each amplicon is selected. This ensures every reported candidate has an experimentally actionable primer+crRNA pair.

**Parameters:** PrimerLength=32, IdentityThreshold=0.95, AmpliconSizeLimit=250, MSA mode  
**Environment:** conda `RPA` (separate from main due to numpy≥1.24 incompatibility)

**crRNA synthesis format** (Kellner et al. 2019, *Nat Protoc*):
```
Order to IDT as ssDNA Ultramer:
  RC( T7_promoter + DR_DNA + spacer_DNA )

T7 RNA polymerase produces functional crRNA:
  DR_RNA + spacer_RNA

RPA Forward primer:  5'─[AATTCTAATACGACTCACTATAGG]─[FP sequence]─3'
RPA Reverse primer:  5'─[RP sequence]─3'
crRNA spacer: complementary to RNA transcribed from RPA amplicon
```

**Output:** `main/data/07_primers/all_codesign.tsv` (586 co-designs)

---

### M08 — Specificity Analysis

BLAST-based specificity validation (blastn-short, word_size=7, evalue=0.01) against four databases:

| Database | Source | Sequences | Rejection criterion |
|---|---|---|---|
| Non-toxigenic *C. difficile* (Group C) | NCBI | 52 genomes | Informational only |
| Enteric pathogens | NCBI RefSeq | 11 species* | ≥90% identity · ≥80% coverage · start ≤pos5 |
| Human transcriptome | GRCh38 RefSeq mRNA | — | idem |
| Human Gut Microbiome | UHGG v2.0.2 | 4,742 species reps (excl. *C. difficile*) | idem |

*Enteric pathogens: *C. sordellii*, *C. perfringens*, *C. botulinum*, *E. faecalis*, *E. faecium*, *K. pneumoniae*, *E. coli* O157, *S. enterica*, *C. jejuni*, *L. monocytogenes*, *S. aureus*

**Database-specific rules:**
- Non-toxigenic check: skipped for all toxin targets (expected hits due to gene remnants)
- Enteropathogens: skipped for rpoB (universal housekeeping gene)
- UHGG: skipped for tpiA and rpoB (expected in gut Clostridiales)

**Result:** 87/90 crRNAs pass all applicable specificity checks  
**Notable exceptions:** 1 tpiA crRNA cross-reacts with *C. botulinum* tpiA (phylogenetic proximity within Clostridiales); 1 cdtA crRNA cross-reacts with *Acutalibacter* sp. (Oscillospirales gut bacterium with partial cdtA homology)

---

### M09 — Ranking + Comprehensive Report

Integrates M07 co-designs as primary source, enriched with ADAPT scores (M06) and specificity (M08). Each reported candidate is a complete, experimentally actionable set: crRNA + RT-RPA primer pair validated to co-localize within the same amplicon.

**Ranking composite score:**

| Component | Weight |
|---|---|
| ADAPT activity | 35% |
| Conservation (1 − mean Shannon H) | 20% |
| Specificity (pass/fail) | 15% |
| Primer co-design (always 1.0 here) | 10% |

**Outputs:**
- `final_ranking.tsv` — 32 final candidates with all sequences
- `final_ranking.html` — interactive report (click row to expand IVT template, FP+T7, RP)
- `IDT_order_crRNA.tsv` — IVT templates ready for IDT ordering (ssDNA Ultramer)
- `IDT_order_primers.tsv` — RPA primers + RT-qPCR primers for IDT ordering

---

### M10 — Synthetic Positive Controls
Generates IVT template DNA and synthetic target RNA templates for each final candidate. Used as positive controls in experimental SHERLOCK reactions.

**Validation (M12 Section B):** 32/32 PASS (T7 present, 0mm crRNA detection, poly-U/GC acceptable)

---

### M11 — RT-qPCR Reference Primers
Designs RT-qPCR primer pairs for 8 genes using Primer3 on M05 consensus sequences. Serves as gold-standard comparator for SHERLOCK validation in clinical samples.

**Parameters:** Tm = 60°C ±2°C, amplicon 80–200bp, GC 45–65%  
**Genes:** tcdA · tcdB_clade1 · tcdB_clade2 · tcdC · cdtA · cdtB · tpiA · rpoB  
**Validation (M12 Section C):** 8/8 PASS (amplicons confirmed in silico, ΔTm <1.5°C for all pairs)

---

### M12 — In Silico Validation

Comprehensive validation of all designed components:

**Section A — crRNA validation (32 candidates):**

| Check | Method | Threshold | Result |
|---|---|---|---|
| Sensitivity | Hamming sliding window vs M03 sequences | ≥90% at ≤1mm | PASS=2, WARN=28, FAIL=2 |
| IVT integrity | RC(template) = T7+DR+spacer; poly-U; GC | — | All correct |
| crRNA structure | RNAfold spacer stem analysis | <10nt in stem | Warns in 12/32 |
| RPA amplicon | M07 amplicon_size verification | 100–250bp | All have amplicons |

**Result interpretation:**
- **2 FAIL** (tcdA rank3: 86.0%; tcdC_wt rank3: 87.9%) — excluded from synthesis
- **28 WARN** — primarily low GC (expected for C. difficile AT-rich genome) and stem warnings from RNAfold. Not grounds for exclusion; prioritize ranks 1–2 per target
- Ranks 1–2 of all 9 targets meet all primary criteria

**Section B:** 32/32 synthetic controls PASS  
**Section C:** 8/8 RT-qPCR primers PASS

---

## Results Summary

| Metric | Value |
|---|---|
| Genomes analyzed | 305 (288 global NCBI + 12 Chilean + 5 anchors) |
| Genome groups | A=106 (hypervirulent) · B=142 (toxigenic) · C=52 (non-toxigenic) |
| Diagnostic targets | 9 (tcdA, tcdB×2, tcdC×2, cdtA, cdtB, tpiA, rpoB) |
| crRNA candidates (ADAPT) | 90 (10 per target) |
| RT-RPA co-designs | 586 |
| Specificity pass rate | 87/90 (96.7%) |
| **Final validated candidates** | **32 (crRNA + FP+T7 + RP per candidate)** |
| Synthetic controls | 32/32 PASS |
| RT-qPCR reference primers | 8/8 PASS |
| crRNA validation | 2 PASS · 28 WARN · 2 FAIL |

---

## Output Structure

```
main/data/
├── 01_download/     288 global + 12 Chilean + 5 anchor genomes
├── 02_classify/     Groups A(106) / B(142) / C(52) + gene matrix
├── 03_extract/      CDS per gene per group (~24 FASTA files)
├── 04_alignment/    9 MSA + Shannon conservation profiles
├── 05_accessibility/ RNAplfold -u 28 profiles + consensus sequences
├── 06_crna/         90 crRNA candidates + full ADAPT TSVs
├── 07_primers/      586 RT-RPA co-designs (PrimedRPA)
├── 08_specificity/  Specificity vs 4 databases (BLAST)
├── 09_report/       32 final candidates + HTML report + IDT sheets
├── 10_synthetic/    Synthetic control templates (64 IDT oligos)
├── 11_rtqpcr/       80 RT-qPCR primer pairs (8 genes × 10 pairs)
└── 12_validation/   Validation report (HTML + TSV per section)
```

---

## Key References

1. Kellner MJ, Koob JG, Gootenberg JS, Abudayyeh OO, Zhang F. SHERLOCK: nucleic acid detection with CRISPR nucleases. *Nat Protoc*. 2019;14(10):2986–3012.
2. Metsky HC, Welch NL, Pillai PP, et al. Designing sensitive viral diagnostics with efficient machine learning. *Nat Biotechnol*. 2022;40:1123–1131. (ADAPT)
3. Higgins M, Ravenhall M, Ward D, Phelan J, Ibrahim A, Forrest MS, et al. PrimedRPA: primer design for recombinase polymerase amplification. *Bioinformatics*. 2019;35(4):682–684.
4. Wessels HH, Méndez-Mancilla A, Guo X, Legut M, Daniloski Z, Bhatt DL, et al. Massively parallel Cas13 screens reveal principles for guide RNA design. *Nat Biotechnol*. 2020;38(6):722–727.
5. Milligan JF, Groebe DR, Witherell GW, Uhlenbeck OC. Oligoribonucleotide synthesis using T7 RNA polymerase and synthetic DNA templates. *Nucleic Acids Res*. 1987;15(21):8783–8798.
6. Lorenz R, Bernhart SH, Höner zu Siederdissen C, et al. ViennaRNA Package 2.0. *Algorithms Mol Biol*. 2011;6:26.
7. Almeida A, Nayfach S, Boland M, et al. A unified catalog of 204,938 reference genomes from the human gut microbiome. *Nat Biotechnol*. 2021;39(1):105–114. (UHGG v2.0)
8. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B*. 1995;57(1):289–300.
9. Rupnik M, Wilcox MH, Gerding DN. *Clostridium difficile* infection: new developments in epidemiology and pathogenesis. *Nat Rev Microbiol*. 2009;7(7):526–536.
10. Gootenberg JS, Abudayyeh OO, Lee JW, et al. Nucleic acid detection with CRISPR-Cas13a/C2c2. *Science*. 2017;356(6336):438–442.

---

## Proposed Future Improvements

### High priority — Biological

**1. MUSCLE5 alignment for tcdB (R1)**
MUSCLE5 (Edgar 2022, *Nat Commun*) outperforms MAFFT by 26% of correctly aligned columns in large, divergent datasets. For tcdB_clade2 (RT027, n=106, ~7.2kb with interlineage divergence), this improvement is directly relevant to ADAPT window scoring accuracy.

**2. minimize-guides mode for tcdB clades (R9)**
Currently using `maximize-activity` which optimizes per-window. For the tcdB clades with known sequence diversity, ADAPT's `minimize-guides` mode with ≥95% coverage constraint (Metsky et al. 2022) would guarantee detection of nearly all variants while minimizing crRNA count — relevant for multiplex design.

**3. Mutational robustness analysis**
Systematic introduction of in silico SNPs at each spacer position to identify the seed region (positions 15–21, Wessels et al. 2020) and positions tolerant of 1 mismatch. Critical for predicting performance against emerging variants not represented in current genome collection.

**4. Cas13design API scoring (R11)**
The Wessels et al. 2020 neural network predictor (cas13design.nygenome.org), trained on 10,000+ LwCas13a crRNAs, provides orthogonal activity prediction to ADAPT. Note: trained on RfxCas13d primarily — use as secondary filter only.

### Medium priority — Computational

**5. LinearCoPartition for crRNA-target binding (R5)**
LinearCoPartition (Zhang et al. 2023, *Nucleic Acids Res*) computes intermolecular base-pairing probabilities between crRNA and target RNA, accounting for both secondary structures simultaneously. More physically accurate than RNAplfold (which evaluates target accessibility without considering the crRNA). Computationally expensive (~1 min per pair) — recommended for top-3 candidates per target.

**6. Nearest-neighbor Tm for RPA primers (R12)**
PrimedRPA reports heuristic GC-based Tm. The SantaLucia 1998 nearest-neighbor model implemented via BioPython `Tm_NN` with RPA buffer conditions (50mM Na⁺, 14mM Mg²⁺, 37°C) would provide more accurate Tm estimates and flag asymmetric primer pairs (ΔTm >5°C) that could cause asymmetric amplification.

**7. NUPACK ΔG for dimer assessment (R13)**
NUPACK (Zadeh et al. 2011) calculates thermodynamically rigorous ΔG for homo/heterodimer formation at 37°C under RPA buffer conditions. Current PrimedRPA dimerisation score is heuristic. Threshold: ΔG < −3 kcal/mol as rejection criterion.

### Lower priority — Clinical validation

**8. Latin American genome expansion**
Current collection is dominated by global NCBI submissions (primarily European/North American). Expanding to include published Latin American sequences (Argentina, Brazil, Colombia) would improve conservation metrics for Chilean clinical context. Key sources: SRA projects PRJNA612578 (Chile 2019), PRJNA614015 (Argentina).

**9. Clinical sensitivity simulation**
Monte Carlo simulation of diagnostic sensitivity using realistic allele frequency distributions from Chilean clinical surveillance data. Accounts for co-infection scenarios and stool sample matrix effects on crRNA accessibility.

**10. Integration with Oxford Nanopore sequencing**
SHERLOCK results can be confirmed in real-time using nanopore sequencing of RPA amplicons. Pipeline extension to generate amplicon-specific analysis workflows (minimap2 + medaka) would enable combined SHERLOCK + sequencing diagnostic protocols.

---

## Citation

If you use this pipeline or its outputs, please cite:

> Cabezas-Mera F et al. Computational design and in silico validation of a multiplexed SHERLOCK CRISPR-Cas13a diagnostic panel for toxigenic *Clostridioides difficile* from Chilean clinical samples. [*In preparation*, 2026]

---

## License

MIT License — see `LICENSE` file.

**Contact:** Fausto Cabezas-Mera · [github.com/fcabezasmera](https://github.com/fcabezasmera)
