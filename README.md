# Genetic Architecture and Clinical Consequences of Adaptive Immune Cell Abundance

## Overview

Adaptive immune cell abundance varies substantially between individuals and plays a critical role in immune function and disease susceptibility. In this study, we applied ImmuneLens to whole-genome sequencing data from 849,071 participants in the All of Us Research Program and UK Biobank to estimate:

- T-cell fraction (TCF)
- B-cell fraction (BCF)

We performed genome-wide association studies, fine-mapping, functional annotation, rare variant analyses, longitudinal modeling, phenome-wide association studies, Mendelian randomization, and experimental validation to characterize the genetic architecture and clinical consequences of adaptive immune cell abundance.

---

## Study Components

### 1. Common Variant Association Analyses

- Genome-wide association study
- Fixed-effects meta-analysis of UK Biobank and All of Us

Regenie v.3.3 was ran on the All of Us researcher workbench, UK Biobank research analysis platform, and BioVU Terra.bio platform.

Step 1 options for Regenie v.3.3:

    --step 1 \
    --bed ukb_imp_step1 \
    --phenoFile tcf_bcf.txt \
    --covarFile tcf_bcf.txt \
    --phenoColList tcf \
    --covarColList age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex,smoking \
    --extract /home/dnanexus/PACER_UKB_GWAS_step1QC_plink_mac5000_thinned.snplist \
    --bsize 1000 \
    --use-relative-path \
    --out tcf_ukb

After Step 1 was ran with 500,000 variants, Step 2 was ran separately by chromosome in parallel. 

Step 2 options for Regenie v.3.3 (example of chr1):

    --step 2 \
    --bgen ukb22828_c14_b0_v3.bgen \
    --sample ukb22828_c14_b0_v3.sample \
    --phenoFile tcf_bcf.txt \
    --covarFile tcf_bcf.txt \
    --phenoColList tcf \
    --covarColList age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --catCovarList sex,smoking \
    --pred tcf_ukb_pred.list \
    --extract /home/dnanexus/imputed_UKB_GWAS_step2QC_plink_maf0.001_geno0.1_chr14.snplist \
    --bsize 200 \
    --test additive \
    --ref-first \
    --gz \
    --out tcf_ukb_chr14

Meta-analysis was performed with the METAL version from 2011-03-25 using the standard-error analysis scheme.

### 2. Fine-Mapping and Functional Annotation

- Statistical fine-mapping
- SCAVENGE analysis
- Variant-to-gene mapping
- IEI enrichment

Code: https://github.com/bicklab/tbcf/tree/main/scdata_integ

### 3. Rare Variant Analyses

- SKAT-O burden testing
- COAST allelic series analyses

Code: https://github.com/bicklab/tbcf/tree/main/TCF_BCF_rare_variant_analysis_ukb.py

### 4. Functional Validation

- EBF1 enhancer reporter assays
- CRISPR interference
- Human CD34+ HSPC editing and B-cell differentiation experiments

### 5. Clinical Consequence Analyses

- Time-to-event PheWAS: [KZ_PheWAS_TCellF_UKB.r](https://github.com/bicklab/tbcf/tree/main/KZ_PheWAS_TCellF_UKB.r)
- Effect comparison: [compare_phewas_effects_tcf_bcf_lywbc.py](https://github.com/bicklab/tbcf/tree/main/compare_phewas_effects_tcf_bcf_lywbc.py)
- Polygenic risk score analyses: [PheWAS_TCF_BCF_PRS.r](https://github.com/bicklab/tbcf/tree/main/PheWAS_TCF_BCF_PRS.r)
- Mendelian randomization
- Gene prioritization and therapeutic target evaluation
- BMI analysis: [KZ_TCF_obesity_BMI.r](https://github.com/bicklab/tbcf/tree/main/KZ_TCF_obesity_BMI.r) 
  
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
