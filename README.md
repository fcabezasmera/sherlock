# SHERLOCK crRNA Design Pipeline

**Computational pipeline for CRISPR-Cas13a (SHERLOCK) diagnostic assay design targeting toxigenic *Clostridioides difficile* from Chilean clinical fecal samples.**

Designed for the LwCas13a platform (Option C: two multiplexed reactions). All candidate crRNAs are co-designed with RT-RPA primer pairs for experimental implementation.

---

## Overview

```
Genome download → Classification → CDS extraction → Alignment → Accessibility
     M01              M02              M03            M04           M05
      ↓                                                               ↓
 crRNA design (M06, ADAPT) → RT-RPA co-design (M07) → Specificity (M08) → Report (M09)
                                                                              ↓
                              Synthetic controls (M10) + RT-qPCR (M11) + Validation (M12)
```

---

## Reaction Design

| Reaction | Enzyme | Targets |
|---|---|---|
| **A** | LwCas13a | tcdB + tcdA + 16S rRNA (internal control) |
| **B** | LwCas13a | cdtA + cdtB + tcdC_WT + tcdC_RT027_junction + tpiA + sodA |

---

## Installation

```bash
git clone https://github.com/fcabezasmera/sherlock.git
cd sherlock

conda env create -f envs/sherlock.yml && conda activate sherlock
conda env create -f envs/adapt.yml
conda env create -f envs/RPA.yml
```

### Required tools

| Tool | Env | Purpose |
|---|---|---|
| NCBI datasets ≥15 | sherlock | Genome download |
| MAFFT ≥7.5 | sherlock | MSA |
| RNAplfold ≥2.6 | sherlock | mRNA accessibility |
| ADAPT ≥1.4 | adapt | crRNA design |
| PrimedRPA ≥1.0.3 | RPA | RT-RPA primer design |
| BLAST+ ≥2.13 | sherlock | Specificity |
| seqkit ≥2.4 | sherlock | Sequence utilities |
| Primer3 ≥2.6 | sherlock | RT-qPCR primers |
| RNAfold ≥2.6 | sherlock | crRNA structure |

---

## Configuration

All parameters centralized in `config.yaml`:

```yaml
organism:
  taxon_id: "1496"
  display:  "Clostridioides difficile"
paths:
  main: ~/sherlock/main
crna:
  dr:  "GGGGAUUUAGACUACCCCAAAAACGAAGGGGGGACUAAAAC"
  t7:  "AATTCTAATACGACTCACTATAGG"
primers:
  rpa_env:     "RPA"
  amplicon_min: 100
  amplicon_max: 250
threads:
  default: 8
```

---

## Running the Pipeline

```bash
conda activate sherlock

python3 scripts/M01_download.py      --config config.yaml
python3 scripts/M02_classify.py      --config config.yaml
python3 scripts/M03_extract.py       --config config.yaml
python3 scripts/M04_alignment.py     --config config.yaml
python3 scripts/M05_accessibility.py --config config.yaml
conda run -n adapt \
  python3 scripts/M06_crna.py        --config config.yaml
python3 scripts/M07_primers.py       --config config.yaml
python3 scripts/M08_specificity.py   --config config.yaml
python3 scripts/M09_report.py        --config config.yaml
python3 scripts/M10_synthetic.py     --config config.yaml
python3 scripts/M11_rtqpcr.py        --config config.yaml
python3 scripts/M12_validation.py    --config config.yaml
```

### Skip flags (re-runs)

```bash
python3 scripts/M06_crna.py         --config config.yaml --skip-adapt
python3 scripts/M07_primers.py      --config config.yaml --skip-primedrpa
python3 scripts/M08_specificity.py  --config config.yaml --skip-db-build --skip-blast
python3 scripts/M12_validation.py   --config config.yaml --skip-structure --skip-primers
```

---

## Pipeline Modules

### M01 — Genome Download
Downloads *C. difficile* genomes (Complete + Chromosome) from NCBI using `datasets`. Includes 12 Chilean clinical strains (NCBI Apr 2026) and 5 anchor references.

**Output:** `main/data/01_download/genomes/{working,chilean,anchors}/`

---

### M02 — Classification

| Group | Description | n |
|---|---|---|
| **A** | Hypervirulent (RT027-like): tcdB + cdtA/cdtB + non-functional tcdC | 106 |
| **B** | Toxigenic (RT012-like): tcdB, no binary toxin | 142 |
| **C** | Non-toxigenic: no functional tcdB | 52 |

**Output:** `main/data/02_classify/gene_matrix.tsv`, `groupA/B/C.txt`

---

### M03 — CDS Extraction
Extracts coding sequences for 9 target genes from GFF + FASTA using seqkit.

**Targets:** tcdA, tcdB, tcdC, tcdC_junction (RT027), cdtA, cdtB, tpiA, sodA, 16S rRNA

---

### M04 — Alignment + Conservation
Concatenates sequences by target, aligns with MAFFT, computes Shannon entropy per position.

**9 MSA targets:** tcdA_all, tcdB_all, tcdC_wt, tcdC_junction, cdtA_groupA, cdtB_groupA, tpiA_all, sodA_all, 16S_all

---

### M05 — mRNA Accessibility
Unpaired probability profiles using RNAplfold (W=80, L=40, --noLP). Generates consensus sequences for Primer3 (M11).

---

### M06 — crRNA Design (ADAPT)
Maximize-activity mode, LwCas13a model, sliding window W=250, guide length 28nt.

**Result:** 81 crRNA candidates, all frac_bound ≈ 1.0

---

### M07 — RT-RPA Co-design (PrimedRPA)
**Inverse strategy:** PrimedRPA on full MSA → find best crRNA from ADAPT within each amplicon (not the reverse).

- PrimerLength=32, IdentityThreshold=0.95, AmpliconSizeLimit=250
- Conda env: RPA (numpy compatibility)

**Result:** 475 co-designed primer+crRNA sets across 9 targets

---

### M08 — Specificity (4 databases)

| Database | Source | n |
|---|---|---|
| Non-toxigenic *C. difficile* (Group C) | NCBI | 52 genomes |
| Enteric pathogens | NCBI RefSeq | 11 species |
| Human transcriptome | GRCh38 RefSeq mRNA | — |
| UHGG gut microbiome | EBI UHGG v2.0.2 | 4,742 species reps (excl. *C. difficile*) |

**Rejection:** ≥90% identity · ≥80% crRNA coverage · alignment start ≤pos 5

**Result:** 79/81 M06 crRNAs pass specificity

---

### M09 — Ranking + Report

**Primary source:** M07 co-designs (crRNA + primer pair co-localized in amplicon)

**Ranking weights:** ADAPT activity 35% · Conservation 20% · Specificity 15% · Primer co-design 10%

**crRNA synthesis format** (Kellner et al. 2019, Nat Protoc):
```
Order to IDT as ssDNA Ultramer:
  RC( T7_promoter + DR_DNA + spacer_DNA )

T7 transcribes:  DR_RNA + spacer_RNA  →  functional crRNA

RPA FP:  5'─[AATTCTAATACGACTCACTATAGG]─[FP sequence]─3'
RPA RP:  5'─[RP sequence]─3'
crRNA spacer: complementary to RNA transcribed from amplicon
```

**Outputs:**
- `final_ranking.tsv` — 29 final candidates with all sequences
- `final_ranking.html` — interactive report (click row to expand)
- `IDT_order_crRNA.tsv` — IVT templates ready for IDT
- `IDT_order_primers.tsv` — RPA + RT-qPCR primers for IDT

---

### M10 — Synthetic Controls
IVT template DNA per crRNA + synthetic target RNA template for positive controls.

**Result:** 45 synthetic control sets, 90 IDT oligos

---

### M11 — RT-qPCR Reference Primers
Primer3 on M05 consensus sequences. 7 genes: tcdA, tcdB, tcdC, cdtA, cdtB, tpiA, sodA.

**Parameters:** Tm=60°C ±2°C, amplicon 80-200bp, GC 45-65%

**Result:** 70 primer pairs (10 per gene)

---

### M12 — In Silico Validation

| Section | Checks | Result |
|---|---|---|
| **A. crRNA** | Sensitivity (0/1/2mm), IVT integrity, RNAfold structure, seqkit amplicon | 26 WARN, 3 FAIL (tcdB variability) |
| **B. Synthetic (M10)** | T7 presence, crRNA detects target 0mm, poly-U, GC | 45/45 PASS |
| **C. RT-qPCR (M11)** | seqkit amplicon, size 80-200bp, Tm delta | 7/7 PASS |

Note: 3 tcdB FAIL due to low sensitivity (73-86%) — tcdB has high variability between RT027 and RT012 ribotypes. This is biologically expected and documented.

Note: WARNs across all targets are predominantly low GC (spacers are AT-rich) — expected since *C. difficile* genome is 29% GC.

---

## Final Results Summary

| Metric | Value |
|---|---|
| Genomes analyzed | 305 (288 global + 12 Chilean + 5 anchors) |
| crRNA candidates (ADAPT) | 81 across 9 targets |
| RT-RPA co-designs | 475 |
| Specificity pass rate | 79/81 (97.5%) |
| **Final validated sets** | **32 (crRNA + FP+T7 + RP)** |
| Synthetic controls | 32/32 PASS |
| RT-qPCR reference primers | 8/8 PASS |
| crRNA validation | 2 PASS · 28 WARN (GC/struct) · 2 FAIL (sensitivity) |

---

## References

1. Kellner MJ et al. SHERLOCK: nucleic acid detection with CRISPR nucleases. *Nat Protoc* 2019.
2. Metsky HC et al. Designing sensitive viral diagnostics with efficient machine learning. *Nat Biotechnol* 2022.
3. Higgins M et al. PrimedRPA. *Bioinformatics* 2019.
4. Almeida A et al. A unified catalog of 204,938 reference genomes from the human gut microbiome. *Nat Biotechnol* 2021.
5. Lorenz R et al. ViennaRNA Package 2.0. *Algorithms Mol Biol* 2011.

---

## Citation

> Cabezas-Mera F et al. Computational design of SHERLOCK CRISPR-Cas13a diagnostic assay for toxigenic *Clostridioides difficile* from Chilean clinical samples. [In preparation]

---

## License

MIT License — see `LICENSE` file.

**Contact:** Francisco Cabezas-Mera · [github.com/fcabezasmera](https://github.com/fcabezasmera)
