library(broom)
library(scales)
library(stringr)
library(arrow)
library(survival)
library(tidyr)
library(dplyr)
library(haven)
library(reshape2)
library(tidyverse)
library(data.table)
library(lubridate)
library(ggplot2)

x <- load("data_all_pheno_cox.Rdata")
x

data_cfs <- fread('cfs_v2_with_covs.txt')

dim(data_all_pheno_cox)

data_cfs_all <- merge(data_cfs, data_all_pheno_cox, by.x = "IID" , by.y = "ID_VUMC", all = F)

dim(data_cfs_all)

names(data_cfs_all) <- iconv(names(data_cfs_all), from = "", to = "UTF-8", sub = "byte")

data_cfs_all <- data_cfs_all %>%
  mutate_at(vars(66:2461), ~ replace_na(., 0))

disease_cols <- colnames(data_cfs_all)[66:2461] #66:2461
length(disease_cols)

time_col<- colnames(data_cfs_all)[2462:4857]
length(time_col)

data_cfs_all2 <- data_cfs_all

#cox_results <- list()

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  
  cat("Analyzing", disease_col, "（No", i, "）\n")

  unique_values <- unique(data_cfs_all[[disease_col]])
  if (!all(unique_values %in% c(0, 1))) {
    warning(paste("SKIP", disease_col, "is not 0/1"))
    next
  }
  
  data_surv <- data_cfs_all
  
  formula_str <- paste0("`", disease_col, "` ~ TCF+ Age_At_Blooddraw + Sex + PC1 + PC2 + PC3 + PC4 + PC5+PC6+PC7+PC8+PC9+PC10 ")
  
  logit_model <- glm(
    formula = as.formula(formula_str),
    data = data_surv,
    family = binomial()
  )
  
  logit_result <- tidy(logit_model) %>%
    filter(term == "TCF") %>%
    mutate(
      total_samples = sum(data_surv[[disease_col]] == 1),
      log_p_value = -log10(p.value)  
    )
  
  cox_results[[disease_col]] <- logit_result
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- results_df %>% arrange(p.value)
results_df

write.csv(results_df, file = "PheWAS_TCF_UKB.csv")

phecodex <- read.csv("phecodeX_ICD_CM_map_flat.csv")
phecodex <- phecodex[,c("phecode","phecode_string","category")]
phecodex <- distinct(phecodex)

results_df <- read.csv("PheWAS_TCF_UKB.csv")
results2_df <- filter(results_df,is.na(results_df$estimate) == F)
phecode_definitions3<-phecodex[,c("phecode_string","category")]
phecode_definitions3 <- distinct(phecode_definitions3)

results2_df$disease <- ifelse(results2_df$disease == "BI_170", "Decreased white blood cell count", results2_df$disease)
results2_df$disease <- ifelse(results2_df$disease == "MB_282", "Nicotine dependence (current and history of)", results2_df$disease)

results2_df <- merge(results2_df,phecode_definitions3,by.x = "disease",by.y = "phecode_string", all.x = T)

results2_df <- results2_df %>%
  mutate(outcome_index = row_number(),  
         shape = ifelse(estimate > 0, 24, 25)) 

results2_df$disease <- ifelse(results2_df$disease == "ID_092", "SIRS", results2_df$disease)
results2_df$category <- ifelse(results2_df$disease == "SIRS", "Infections", results2_df$category)

top_10_outcomes <- results2_df %>%
  arrange(desc(log_p_value)) %>%
  slice_head(n = 25) %>%
  pull(disease)

results2_df <- results2_df %>%
  mutate(top_10_label = ifelse(disease %in% top_10_outcomes, disease, NA))

threshold_pvalue_0_05 <- 4.3 

results2_df <- results2_df %>%
  arrange(category, outcome_index) %>%
  mutate(outcome_index = row_number())  

library(viridis)
library(ggrepel)
p <- ggplot(results2_df, aes(x = outcome_index, y = log_p_value)) +
  geom_point(aes(color = category, shape = as.factor(shape), fill = category), size = 3) +
  scale_shape_manual(values = c(24, 25)) +
  geom_hline(yintercept = threshold_pvalue_0_05, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = top_10_label), size = 2.5, max.overlaps = 15, min.segment.length = 0) + 
  theme_minimal() +
  labs(x = "Category", y = "-log10(p.value)", title = "Manhattan Plot of TCF PheWAS in UKB") +
  scale_color_viridis_d(name = "Category", option = "plasma") +  #viridis, magma, plasma, inferno, cividis
  scale_fill_viridis_d(name = "Category", option = "plasma") +  
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = results2_df %>% group_by(category) %>% summarize(mid = mean(outcome_index)) %>% pull(mid),
                     labels = results2_df %>% group_by(category) %>% summarize(cat = first(category)) %>% pull(cat))
p

library(sjPlot)
save_plot("TCF_PheWAS_UKB.svg", fig = p, width=30, height=20)

target_columns <- data_cfs_all[, 66:2461]

freq_counts <- colSums(target_columns != 0, na.rm = TRUE)

low_frequency_columns <- freq_counts[freq_counts <= 10]
cat("Number of columns with frequency < 10:", length(low_frequency_columns), "\n")

low_frequency_columns <- names(low_frequency_columns)
length(low_frequency_columns)

disease_cols <- setdiff(disease_cols, low_frequency_columns)
length(disease_cols)

disease_cols <- disease_cols[1707:1864]
disease_cols

cox_results <- list()

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  time_col <- paste0("time_", disease_col)

  cat("Processing:", time_col, "（Loop", i, "of", length(disease_cols), "）\n")
    
  data_surv <- data_cfs_all %>%
    mutate(
      surv_time = case_when(
        !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), min_date, units = "days")) / 365.25,
        !!sym(disease_col) == 0 ~ as.numeric(difftime(max_date, min_date, units = "days")) / 365.25
      )
    ) %>%
    filter(surv_time > 0)
  
  total_samples <- sum(data_surv[[disease_col]] == 1)

  surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])
  
  cox_model <- coxph(surv_obj ~ TCF+ Age_At_Blooddraw + Sex + smoking_0+ PC1.x + PC2.x + PC3.x + PC4 + PC5.x+PC6.x+PC7.x+PC8.x+PC9.x+PC10.x, data = data_surv)
  
  cox_result <- tidy(cox_model) %>%
    filter(term == "TCF") %>%
    mutate(
      total_samples = total_samples,
      log_p_value = -log10(p.value)
    )
  
  cox_results[[disease_col]] <- cox_result
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- results_df %>% arrange(p.value)
results_df

results_df2 <- fread("TCF_PheWAS_UKB_cox.csv")
results_df2

results_df2 <- results_df2[,c(2:9)]

results_df3 <- rbind(results_df,results_df2)
results_df3 <- results_df3 %>% arrange(p.value)
results_df3

write.csv(results_df3, file = "TCF_PheWAS_UKB_cox_all.csv")

disease_cols <- colnames(data_cfs_all)[66:2461] #66:2461
disease_cols

target_columns <- data_cfs_all[, 66:2461]

freq_counts <- colSums(target_columns != 0, na.rm = TRUE)

low_frequency_columns <- freq_counts[freq_counts <= 10]
cat("Number of columns with frequency < 10:", length(low_frequency_columns), "\n")

low_frequency_columns <- names(low_frequency_columns)
length(low_frequency_columns)

disease_cols <- setdiff(disease_cols, low_frequency_columns)
length(disease_cols)

cox_results <- list()

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  time_col <- paste0("time_", disease_col)

  cat("Processing:", time_col, "（Loop", i, "of", length(disease_cols), "）\n")
    
  data_surv <- data_cfs_all %>%
    mutate(
      surv_time = case_when(
        !!sym(disease_col) == 1 ~ as.numeric(difftime(!!sym(time_col), min_date, units = "days")) / 365.25,
        !!sym(disease_col) == 0 ~ as.numeric(difftime(max_date, min_date, units = "days")) / 365.25
      )
    ) %>%
    filter(surv_time > 0)
  
  total_samples <- sum(data_surv[[disease_col]] == 1)

  surv_obj <- Surv(time = data_surv$surv_time, event = data_surv[[disease_col]])
  
  cox_model <- coxph(surv_obj ~ BCF+ Age_At_Blooddraw + Sex + smoking_0+ PC1.x + PC2.x + PC3.x + PC4 + PC5.x+PC6.x+PC7.x+PC8.x+PC9.x+PC10.x, data = data_surv)
  
  cox_result <- tidy(cox_model) %>%
    filter(term == "BCF") %>%
    mutate(
      total_samples = total_samples,
      log_p_value = -log10(p.value)
    )
  
  cox_results[[disease_col]] <- cox_result
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- results_df %>% arrange(p.value)
results_df
write.csv(results_df, file = "BCF_PheWAS_UKB_cox_all.csv")

for (col in colnames(data_cfs_all)[grep("^time_", colnames(data_cfs_all))]) {
  data_cfs_all[[col]] <- as.Date(data_cfs_all[[col]])
}

data_cfs_all$min_date <- as.Date(as.character(data_cfs_all$min_date))

library(dplyr)
library(lubridate)
library(purrr)

disease_names <- gsub("^time_", "", colnames(data_cfs_all)[2462:4857])
biosample_col <- data_cfs_all$min_date

class(biosample_col)
class(data_cfs_all$time_Sepsis)

for (disease in disease_names) {
  disease_col <- disease
  time_col <- paste0("time_", disease)

  if (disease_col %in% names(data_cfs_all) && time_col %in% names(data_cfs_all)) {
    time_diff <- as.numeric(data_cfs_all[[time_col]] - biosample_col)
    idx_to_zero <- which(is.na(time_diff) | abs(time_diff) > 90)
    data_cfs_all[[disease_col]][idx_to_zero] <- 0
  }
}

table(data_cfs_all$Hypertension)
table(data_cfs_all2$Hypertension)

save(data_cfs_all, file = "data_CFs_UKB_3month.Rdata")

cox_results <- list()

for (i in seq_along(disease_cols)) {
  disease_col <- disease_cols[i]
  
  cat("Analyzing", disease_col, "（No", i, "）\n")

  unique_values <- unique(data_cfs_all[[disease_col]])
  if (!all(unique_values %in% c(0, 1))) {
    warning(paste("SKIP", disease_col, "is not 0/1"))
    next
  }
  
  data_surv <- data_cfs_all
  
  formula_str <- paste0("`", disease_col, "` ~ TCF+ Age_At_Blooddraw + Sex + PC1.x + PC2.x + PC3.x + PC4 + PC5.x+PC6.x+PC7.x+PC8.x+PC9.x+PC10.x ")
  
  logit_model <- glm(
    formula = as.formula(formula_str),
    data = data_surv,
    family = binomial()
  )
  
  logit_result <- tidy(logit_model) %>%
    filter(term == "TCF") %>%
    mutate(
      total_samples = sum(data_surv[[disease_col]] == 1),
      log_p_value = -log10(p.value)  
    )
  
  cox_results[[disease_col]] <- logit_result
}

results_df <- bind_rows(cox_results, .id = "disease")
results_df <- results_df %>% arrange(p.value)
results_df

results_df2 <- filter(results_df, results_df$total_samples >=10)
results_df2

write.csv(results_df, file = "PheWAS_TCF_UKB_90days.csv")
