# Genetic Architecture and Clinical Consequences of Adaptive Immune Cell Abundance

## Overview

Adaptive immune cell abundance varies substantially between individuals and plays a critical role in immune function and disease susceptibility. In this study, we applied ImmuneLens to whole-genome sequencing data from 849,071 participants in the All of Us Research Program and UK Biobank to estimate:

- T-cell fraction (TCF)
- B-cell fraction (BCF)
- T-cell count (TCC)
- B-cell count (BCC)

We performed genome-wide association studies, fine-mapping, functional annotation, rare variant analyses, longitudinal modeling, phenome-wide association studies, Mendelian randomization, and experimental validation to characterize the genetic architecture and clinical consequences of adaptive immune cell abundance.

---

## Study Components

### 1. Common Variant Association Analyses

- TCF GWAS
- BCF GWAS
- TCC GWAS
- BCC GWAS
- Fixed-effects meta-analysis of UK Biobank and All of Us

### 2. Fine-Mapping and Functional Annotation

- Statistical fine-mapping
- SCAVENGE analysis
- Variant-to-gene mapping
- Single-cell chromatin accessibility integration
- eQTL integration
- PromoterAI annotation

Code: https://github.com/bicklab/tbcf/tree/main/scdata_integ

### 3. Rare Variant Analyses

- SKAT-O burden testing
- COAST allelic series analyses

### 4. Population-Specific Genetic Architecture

- Ancestry-stratified GWAS
- Age-stratified GWAS
- Sex-stratified GWAS
- Longitudinal lymphocyte decay analyses

### 5. Clinical Consequence Analyses

- Time-to-event PheWAS: [KZ_PheWAS_TCellF_UKB.r](https://github.com/bicklab/tbcf/tree/main/KZ_PheWAS_TCellF_UKB.r)
- Polygenic risk score analyses: [PheWAS_TCF_BCF_PRS.r](https://github.com/bicklab/tbcf/tree/main/PheWAS_TCF_BCF_PRS.r)
- Mendelian randomization
- Gene prioritization and therapeutic target evaluation
- BMI analysis: [TCF_obesity_BMI.r](https://github.com/bicklab/tbcf/tree/main/TCF_obesity_BMI.r) 
  
### 6. Functional Validation

- EBF1 enhancer reporter assays
- CRISPR interference
- Human CD34+ HSPC editing and B-cell differentiation experiments

---

## Data Availability

Individual-level participant data are not included in this repository.

Access to the underlying datasets requires approval through:

- [UK Biobank DNA Nexus Research Analysis Platform](https://ukbiobank.dnanexus.com)
- [All of Us Research Workbench](https://workbench.researchallofus.org/)
- Vanderbilt BioVU (for validation analyses)

Summary statistics and derived results associated with the manuscript will be made available upon publication.

---

## Software

Major software used in this project includes:

- R (≥4.3)
- PLINK 2.0
- Regenie
- METAL
- SCAVENGE
- TwoSampleMR
- COAST
- Bedtools
- samtools
- ImmuneLens

Additional package requirements are listed within individual analysis directories.

---

## Citation

If you use code from this repository, please cite:

Manuscript citation to be added upon publication.

---

## Contact
Kun Zhao: kun.zhao@vumc.org
Alexander G. Bick Lab: alexander.bick@vumc.org
