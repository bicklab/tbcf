# library
library(tidyverse)
library(ggplot2)
library(data.table)
library(TwoSampleMR)
library(coloc)

#load exposure
# CD40 chr20:46,118,314-46,129,858
exposure_dat <- read_exposure_data(
  filename = "~/sumstats/BCF/cd40_region_BCF_meta_chr20_annotated.txt",
  sep = "\t", 
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  #eaf_col = "eaf",
  pval_col = "P",
  #phenotype_col = "phenotype_name", 
  samplesize_col = "N"    
)
exposure_dat$exposure <- "CD40_BCF"

# Perform LD clumping with customized parameters
exposure_dat_clumped <- clump_data(
  exposure_dat,
  clump_p1 = 5e-6,
  clump_kb = 1000,
  clump_r2 = 0.2
)

#load outcome data (CD40, RA)
outcome <- read_outcome_data(
  snps = exposure_dat_clumped$SNP, 
  filename = "~/sumstats/ra_gwas_NatGen2021/RA_chr20_sumstats.tsv",
  sep = "\t", 
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  #eaf_col = "af_alt",
  pval_col = "p_value" 
)
outcome$outcome <- "RA"

dat <- harmonise_data(
  exposure_dat = exposure_dat_clumped, 
  outcome_dat = outcome)

results <- mr(dat)

#### Coloc CD40 & RA ####
exposure_dat <- read_exposure_data(
  filename = "~/sumstats/BCF/cd40_region_BCF_meta_chr20_annotated.txt",
  sep = "\t", 
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  #eaf_col = "eaf",
  pval_col = "P",
  #phenotype_col = "phenotype_name", 
  samplesize_col = "N"    
)
exposure_dat$exposure <- "CD40_BCF"

outcome <- read_outcome_data(
  snps = exposure_dat$SNP, 
  filename = "~/sumstats/ra_gwas_NatGen2021/RA_chr20_sumstats.tsv",
  sep = "\t", 
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  #eaf_col = "af_alt",
  pval_col = "p_value" 
)
outcome$outcome <- "RA"

# Merge full regional un-clumped datasets
coloc_dat <- merge(exposure_dat, outcome, by="SNP")

n_exposure <- 849071  # Example N for BCF
n_outcome <- 276020  # Example N for RA GWAS
s_proportion <- 0.13 # Example proportion of cases in RA GWAS

coloc_result <- coloc.abf(
  dataset1 = list(
    beta = coloc_dat$beta.exposure,
    varbeta = coloc_dat$se.exposure^2,
    N = n_exposure, 
    sdY = 1,  # Add this line for quantitative traits if MAF is missing
    type = "quant",
    snp = coloc_dat$SNP
  ),
  dataset2 = list(
    beta = coloc_dat$beta.outcome,
    varbeta = coloc_dat$se.outcome^2,
    N = n_outcome,
    type = "cc",
    s = s_proportion,
    snp = coloc_dat$SNP
  )
)
