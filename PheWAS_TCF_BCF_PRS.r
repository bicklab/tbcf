library(broom)
library(scales)
#library(stringr)
library(arrow)
library(survival)
library(tidyr)
library(dplyr)
library(haven)
#library(reshape2)
#library(tidyverse)
library(data.table)
library(lubridate)
library(ggplot2)

system("gsutil cp -r gs://bicklab-main-storage/Users/Hannah_Poisner/biovu_prs_scores .")

test <- fread("biovu_prs_scores/EUR_TCF_sum_scores_0.5.txt")

test

score_files <- list.files("biovu_prs_scores", pattern = "_sum_scores_.*\\.txt$", full.names = TRUE)
score_files

test <- fread("biovu_prs_scores//EUR_TCF_sum_scores_0.5.txt")

head(test)

prs_data <- map_dfr(score_files, function(file) {
  info <- str_match(basename(file), "(AFR|EUR)_(TCF|BCF)_sum_scores_(0\\.\\d+).txt")
  tibble::tibble(
    data = list(read.table(file, header = TRUE)),
    population = info[2],
    phenotype = info[3],
    p_threshold = info[4]
  )
}) %>%
  unnest(data) %>%
  mutate(p_threshold = as.numeric(p_threshold))

p <- ggplot(prs_data, aes(x = SUM_SCORE)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  facet_grid(population + phenotype ~ p_threshold, scales = "free_y") +
  labs(
    title = "PRS Score Distributions by Ancestry, Phenotype and P-value Threshold",
    x = "PRS Sum Score",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)
p

ggsave("PRS_TCF_BCF_byP.png", plot = p, width = 30, height = 30, dpi = 500)

install.packages("pheatmap")
library(pheatmap)

files <- list.files("biovu_prs_scores", pattern = ".*_TCF_sum_scores_.*\\.txt$", full.names = TRUE)
files

get_label <- function(file) {
  ancestry <- str_extract(basename(file), "^(AFR|EUR)")
  pval <- str_extract(file, "(?<=scores_)[\\d\\.]+")
  paste0(ancestry, "_", pval)
}

prs_list <- lapply(files, function(file) {
  df <- read_delim(file, delim = "\t")
  colname <- get_label(file)
  rename(df, !!colname := SUM_SCORE)
})

prs_merged <- reduce(prs_list, full_join, by = "FID") %>%
  column_to_rownames("FID")

cor_mat <- cor(prs_merged, use = "pairwise.complete.obs")^2

p <- pheatmap(cor_mat,
         display_numbers = TRUE,
         main = "Correlation (R²) between TCF PRS across ancestries and thresholds",
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize_number = 8,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
p

ggsave("PRS_TCF_heatmap.jpg", plot = p, width = 20, height = 20, dpi = 500)

files <- list.files("biovu_prs_scores", pattern = ".*_BCF_sum_scores_.*\\.txt$", full.names = TRUE)
files

get_label <- function(file) {
  ancestry <- str_extract(basename(file), "^(AFR|EUR)")
  pval <- str_extract(file, "(?<=scores_)[\\d\\.]+")
  paste0(ancestry, "_", pval)
}

prs_list <- lapply(files, function(file) {
  df <- read_delim(file, delim = "\t")
  colname <- get_label(file)
  rename(df, !!colname := SUM_SCORE)
})

prs_merged <- reduce(prs_list, full_join, by = "FID") %>%
  column_to_rownames("FID")

cor_mat <- cor(prs_merged, use = "pairwise.complete.obs")^2

p <- pheatmap(cor_mat,
         display_numbers = TRUE,
         main = "Correlation (R²) between BCF PRS across ancestries and thresholds",
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize_number = 8,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
p

ggsave("PRS_BCF_heatmap.jpg", plot = p, width = 20, height = 20, dpi = 500)

x <- load("BioVU_PhenoData_250k.Rdata")
x

saveRDS(data_final, file = "BioVU_PhenoData_250k.rds")

system("gsutil cp BioVU_PhenoData_250k.rds gs://bicklab-main-storage/Users/Kun_Zhao")

head(data_final)
dim(data_final)

agd_250k_demos <- fread("agd_250k_demos.tsv")
dim(agd_250k_demos)
head(agd_250k_demos)

table(agd_250k_demos$gender)

data_all_pheno_biovu_250k <- merge(agd_250k_demos, data_final, by = "person_id", all = F)

dim(data_all_pheno_biovu_250k)

agd_250k_obs_v2 <- fread("agd_250k_obs_v2.tsv")
head(agd_250k_obs_v2)
dim(agd_250k_obs_v2)

data_all_pheno_biovu_250k <- merge(data_all_pheno_biovu_250k, agd_250k_obs_v2, by = "person_id", all = F)
dim(data_all_pheno_biovu_250k)

dnadate_biovu <- fread("AGD250kcohortshiftedsampledates.csv")
dim(dnadate_biovu)
head(dnadate_biovu)

data_pheno3 <- merge(data_all_pheno_biovu_250k,dnadate_biovu,by = "GRID", all = F)
dim(data_pheno3)

data_pheno3$SHIFTED_SAMPLE_DATE <- as.Date(data_pheno3$SHIFTED_SAMPLE_DATE, format = "%m/%d/%Y")
summary(data_pheno3$SHIFTED_SAMPLE_DATE)
data_pheno3$birth_datetime <- as.Date(data_pheno3$birth_datetime, format = "%y/%m/%d")
summary(data_pheno3$birth_datetime)

data_pheno3$age <- NA
data_pheno3$age <- as.numeric(difftime(data_pheno3$SHIFTED_SAMPLE_DATE, data_pheno3$birth_datetime, units ="days"))/365.25
summary(data_pheno3$age)

data_all_pheno_biovu_250k_age0 <- filter(data_pheno3, data_pheno3$age >0)

data_all_pheno_biovu_250k_age0$age2 <- data_all_pheno_biovu_250k_age0$age * data_all_pheno_biovu_250k_age0$age

dim(data_all_pheno_biovu_250k_age0)

head(data_all_pheno_biovu_250k_age0$GRID)

saveRDS(data_all_pheno_biovu_250k_age0, file = "BioVU_PhenoData_250k_0811.rds")

system("gsutil cp BioVU_PhenoData_250k_0811.rds gs://bicklab-main-storage/Users/Kun_Zhao")

afr_tcf_0.001 <- fread("biovu_prs_scores/AFR_TCF_sum_scores_0.001.txt")
afr_bcf_0.001 <- fread("biovu_prs_scores/AFR_BCF_sum_scores_0.001.txt")
afr_tcf_0.1 <- fread("biovu_prs_scores/AFR_TCF_sum_scores_0.1.txt")
afr_bcf_0.1 <- fread("biovu_prs_scores/AFR_BCF_sum_scores_0.1.txt")

eur_tcf_0.001 <- fread("EUR_TCF_sum_scores_0.001_clean.txt")
eur_bcf_0.001 <- fread("EUR_BCF_sum_scores_0.001_clean.txt")
eur_tcf_0.1 <- fread("EUR_TCF_sum_scores_0.1_clean.txt")
eur_bcf_0.1 <- fread("EUR_BCF_sum_scores_0.1_clean.txt")
eur_tcf_0.5 <- fread("EUR_TCF_sum_scores_0.05_clean.txt")

eur_nohtn_0.1 <- fread("Users_Hannah_Poisner_EUR_TCF_no_htn_sum_scores_0.1.txt")

system("gsutil cp gs://bicklab-main-storage/Users/Hannah_Poisner/eur_prs_fixed/no_htn_TCF_eur_normalized_p0.1.txt .")

eur_nohtn_nor <- fread("no_htn_TCF_eur_normalized_p0.1.txt")

eur_nor <- fread("TCF_eur_normalized_p0.1.txt")

dim(eur_tcf_0.1)

names(afr_tcf_0.001)[2] <- "afr_tcf_0.001"
names(afr_bcf_0.001)[2] <- "afr_bcf_0.001"
names(afr_tcf_0.1)[2] <- "afr_tcf_0.1"
names(afr_bcf_0.1)[2] <- "afr_bcf_0.1"

afr_tcf <- merge(afr_tcf_0.001, afr_tcf_0.1)
afr_bcf <- merge(afr_bcf_0.001, afr_bcf_0.1)
afr_PRS <- merge(afr_bcf, afr_tcf)

head(afr_PRS)

eur_tcf_0.001 <- eur_tcf_0.001[,c(1,24)]
eur_bcf_0.001 <- eur_bcf_0.001[,c(1,24)]
eur_tcf_0.1 <- eur_tcf_0.1[,c(1,24)]
eur_bcf_0.1 <- eur_bcf_0.1[,c(1,24)]
eur_tcf_0.5 <- eur_tcf_0.5[,c(1,24)]

names(eur_tcf_0.001)[2] <- "eur_tcf_0.001"
names(eur_bcf_0.001)[2] <- "eur_bcf_0.001"
names(eur_tcf_0.1)[2] <- "eur_tcf_0.1"
names(eur_bcf_0.1)[2] <- "eur_bcf_0.1"
names(eur_tcf_0.5)[2] <- "eur_tcf_0.5"

eur_nohtn_0.1 <- eur_nohtn_0.1[,c(1,24)]
names(eur_nohtn_0.1)[2] <- "eur_nohtn_0.1"

eur_nohtn_nor <- eur_nohtn_nor[,c(1,24)]
names(eur_nohtn_nor)[2] <- "eur_nohtn_nor"

eur_nor <- eur_nor[,c(1,24)]
names(eur_nor)[2] <- "eur_nor"

head(eur_nor)

eur_tcf <- merge(eur_tcf_0.001, eur_tcf_0.1)
eur_bcf <- merge(eur_bcf_0.001, eur_bcf_0.1)
eur_PRS <- merge(eur_bcf, eur_tcf)
eur_PRS <- merge(eur_PRS, eur_tcf_0.5)

head(eur_PRS)
dim(eur_PRS)

data_all_pheno_biovu_AFR <- merge(data_all_pheno_biovu_250k_age0,afr_PRS, by.x = "GRID", by.y = "FID", all = F)

data_all_pheno_biovu_EUR <- merge(data_all_pheno_biovu_250k_age0,eur_PRS, by.x = "GRID", by.y = "FID", all = F)

dim(data_all_pheno_biovu_EUR)
dim(data_all_pheno_biovu_AFR)

saveRDS(data_all_pheno_biovu_EUR, "data_all_pheno_biovu_EUR.rds")

saveRDS(data_all_pheno_biovu_AFR, "data_all_pheno_biovu_AFR.rds")

data_all_pheno_biovu_EUR <- readRDS("data_all_pheno_biovu_EUR.rds")

dim(data_all_pheno_biovu_EUR)

data_all_pheno_biovu_EUR <- merge(data_all_pheno_biovu_EUR,eur_nor, by.x = "GRID", by.y = "FID", all = F)

data_all_pheno_biovu_EUR <- merge(data_all_pheno_biovu_EUR,eur_nohtn_nor, by.x = "GRID", by.y = "FID", all = F)

data_all_pheno_biovu_EUR <- merge(data_all_pheno_biovu_EUR,eur_nohtn_0.1, by.x = "GRID", by.y = "FID", all = F)

dim(data_all_pheno_biovu_EUR)

read_auto <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "rds") {
    return(readRDS(file))
  } else if (ext %in% c("txt", "csv", "tsv")) {
    return(fread(file, data.table = FALSE)) 
  } else {
    stop(paste("Unsupported file type:", ext))
  }
}

id_input_file <- "FID"
id_phecode_file <- "GRID"

merged_df <- inner_join(
  test,
  data_all_pheno_biovu_250k_age0,
  by = setNames(id_phecode_file, id_input_file)
)

dim()

disease_cols <- colnames(data_all_pheno_biovu_EUR)[8:3584] #8:3584
length(disease_cols)
disease_cols

run_logistic_phewas <- function(data, disease_col_range, exposure, covariates) {
  disease_cols <- names(data)[disease_col_range]
  cox_results <- list()
  
  for (i in seq_along(disease_cols)) {
    disease_col <- disease_cols[i]
    
    cat("Analyzing", disease_col, "(No", i, ")\n")
    
    unique_values <- unique(data[[disease_col]])
    if (!all(unique_values %in% c(0, 1))) {
      warning(paste("SKIP", disease_col, "is not 0/1"))
      next
    }
    
    formula_str <- paste0("`", disease_col, "` ~ ", exposure, " + ", covariates)
    
    logit_model <- glm(
      formula = as.formula(formula_str),
      data = data,
      family = binomial()
    )
    
    logit_result <- broom::tidy(logit_model) %>%
      dplyr::filter(term == exposure) %>%
      dplyr::mutate(
        total_samples = sum(data[[disease_col]] == 1),
        log_p_value = -log10(p.value)
      )
    
    cox_results[[disease_col]] <- logit_result
  }
  
  results_df <- dplyr::bind_rows(cox_results, .id = "disease") %>%
    dplyr::arrange(p.value)
  
  return(results_df)
}

disease_range <- 8:3584
exposure_var <- "eur_nohtn_nor"
covars <- "age + age2 + gender"

results <- run_logistic_phewas(
  data = data_all_pheno_biovu_EUR,
  disease_col_range = disease_range,
  exposure = exposure_var,
  covariates = covars
)

#Load Phecode "phecode_icd10.csv"
phecodex <- read.csv("phecodeX_ICD_CM_map_flat.csv")
phecodex <- phecodex[,c("phecode","phecode_string","category")]
phecodex <- distinct(phecodex)

results_df <- merge(results,phecodex, by.x ="disease", by.y = "phecode", all.x = T)
results_df

results_df <- results_df %>% arrange(p.value)
results_df

write.csv(results_df, file = "results_prs_eur_tcf_0.1_nonHTN_normalized_0714.csv")

system("gsutil cp results_prs_eur_tcf_0.1_nonHTN_normalized_0714.csv gs://bicklab-main-storage/Users/Hannah_Poisner/tcf_bcf_phewas/")



results_df <- bind_rows(cox_results, .id = "disease")
results_df <- merge(results_df,phecodex, by.x ="disease", by.y = "phecode", all.x = T)
results_df <- results_df %>% arrange(p.value)
results_df

write.csv(results_df, file = "results_prs_EUR_tcf_0.1_0616.csv")

disease_cols <- colnames(data_all_pheno_biovu_EUR)[8:3584] #8:3584
length(disease_cols)

cox_results <- list()

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  
  cat("Analyzing", disease_col, "（No", i, "）\n")

  unique_values <- unique(data_all_pheno_biovu_AFR[[disease_col]])
  if (!all(unique_values %in% c(0, 1))) {
    warning(paste("SKIP", disease_col, "is not 0/1"))
    next
  }
  
  data_surv <- data_all_pheno_biovu_AFR
  
  formula_str <- paste0("`", disease_col, "` ~ afr_tcf_0.1+ age + age2 + gender ")
  
  logit_model <- glm(
    formula = as.formula(formula_str),
    data = data_surv,
    family = binomial()
  )
  
  logit_result <- tidy(logit_model) %>%
    filter(term == "afr_tcf_0.1") %>%
    mutate(
      total_samples = sum(data_surv[[disease_col]] == 1),
      log_p_value = -log10(p.value)  
    )
  
  cox_results[[disease_col]] <- logit_result
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- merge(results_df,phecodex, by.x ="disease", by.y = "phecode", all.x = T)
results_df <- results_df %>% arrange(p.value)
results_df

write.csv(results_df, file = "results_prs_afr_tcf_0.1.csv")

cox_results <- list()

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  
  cat("Analyzing", disease_col, "（No", i, "）\n")

  unique_values <- unique(data_all_pheno_biovu_EUR[[disease_col]])
  if (!all(unique_values %in% c(0, 1))) {
    warning(paste("SKIP", disease_col, "is not 0/1"))
    next
  }
  
  data_surv <- data_all_pheno_biovu_EUR
  
  formula_str <- paste0("`", disease_col, "` ~ eur_tcf_0.1+ age + age2 + gender ")
  
  logit_model <- glm(
    formula = as.formula(formula_str),
    data = data_surv,
    family = binomial()
  )
  
  logit_result <- tidy(logit_model) %>%
    filter(term == "eur_tcf_0.1") %>%
    mutate(
      total_samples = sum(data_surv[[disease_col]] == 1),
      log_p_value = -log10(p.value)  
    )
  
  cox_results[[disease_col]] <- logit_result
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- merge(results_df,phecodex, by.x ="disease", by.y = "phecode", all.x = T)
results_df <- results_df %>% arrange(p.value)
results_df

write.csv(results_df, file = "results_prs_eur_tcf_0.1.csv")
