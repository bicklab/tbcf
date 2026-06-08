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

system("gsutil cp gs://bicklab-main-storage/Users/Hannah_Poisner/eur_prs_fixed/TCF_eur_normalized_p0.1.txt .")

system("gsutil cp gs://bicklab-main-storage/Users/Hannah_Poisner/eur_prs_fixed/EUR_BCF_sum_scores_0.1_clean.txt .")

eur_nor <- fread("TCF_eur_normalized_p0.1.txt")
b_eur_nor <- fread("EUR_BCF_sum_scores_0.1_clean.txt")

head(eur_nor)
head(b_eur_nor)
dim(b_eur_nor)

b_eur_nor <- b_eur_nor %>%
  mutate(
    BCF_PRS = scale(Total_Score)[,1]
  )

eur_nor <- eur_nor[,c("FID","normalized_column")]
names(eur_nor)[2] <- "TCF_PRS"

b_eur_nor <- b_eur_nor[,c("FID","BCF_PRS")]

data_all_pheno_cox7 <- readRDS("BioVU_PhenoData_250k_0922.rds") 
dim(data_all_pheno_cox7)

table(data_all_pheno_cox7$unsupervised_ancestry_cluster_relabel)

data_all_pheno_cox7_eur <- filter(data_all_pheno_cox7, data_all_pheno_cox7$unsupervised_ancestry_cluster_relabel == "k5_(EUR)")

data_all_pheno_cox7_eur <- merge(data_all_pheno_cox7_eur, eur_nor, by.x = "GRID", by.y = "FID", all = F)

dim(data_all_pheno_cox7_eur)

data_all_pheno_cox7_eur <- merge(data_all_pheno_cox7_eur, b_eur_nor, by.x = "GRID", by.y = "FID", all = F)

dim(data_all_pheno_cox7_eur)

bmi_biovu <- readRDS("bmi_all_biovu.rds")

head(bmi_biovu)

library(dplyr)
library(lubridate)

bmi <- bmi_biovu %>%
  mutate(
    measurement_datetime = as.POSIXct(measurement_datetime)
  )

pheno <- data_all_pheno_cox7_eur %>%
  mutate(
    SHIFTED_SAMPLE_DATE = as.POSIXct(SHIFTED_SAMPLE_DATE)
  )

bmi_nearest <- pheno %>%
  select(person_id, SHIFTED_SAMPLE_DATE) %>%
  left_join(
    bmi,
    by = "person_id"
  ) %>%
  mutate(
    time_diff_days = abs(
      as.numeric(
        difftime(
          measurement_datetime,
          SHIFTED_SAMPLE_DATE,
          units = "days"
        )
      )
    )
  ) %>%
  group_by(person_id) %>%
  slice_min(
    order_by = time_diff_days,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

data_all_pheno_cox7_eur <- data_all_pheno_cox7_eur %>%
  left_join(
    bmi_nearest %>%
      select(
        person_id,
        BMI0 = bmi,
        BMI0_date = measurement_datetime,
        BMI0_diff_days = time_diff_days
      ),
    by = "person_id"
  )

summary(data_all_pheno_cox7_eur$BMI0)

summary(data_all_pheno_cox7_eur$age)

hist(data_all_pheno_cox7_eur$BMI0, breaks = 100)

quantile(data_all_pheno_cox7_eur$BMI0,
         c(0.99,0.995,0.999,0.9995),
         na.rm=TRUE)

library(rms)

tmp <- data_all_pheno_cox7_eur %>%
  dplyr::select(
    TCF_PRS,
    BMI0,
    age,
    gender,
    PC1_SUM, PC2_SUM, PC3_SUM, PC4_SUM, PC5_SUM
  ) %>%
  na.omit()

dd <- datadist(tmp)
options(datadist = "dd")

fit <- ols(
  TCF_PRS ~ rcs(BMI0, 4) +
    age+
    gender+
    PC1_SUM+PC2_SUM+PC3_SUM+PC4_SUM+PC5_SUM,
  data = data_all_pheno_cox7_eur
)

anova(fit)

pred <- Predict(fit, BMI0)

p <- ggplot(pred,
       aes(BMI0, yhat)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = lower,
        ymax = upper),
    alpha = 0.2
  ) +
  theme_bw()
p

library(sjPlot)

ggsave(
  "BMI_TCFPRS_spline.svg",
  plot = p,
  width = 6,
  height = 6
)

library(rms)

tmp <- data_all_pheno_cox7_eur %>%
  dplyr::select(
    BCF_PRS,
    BMI0,
    age,
    gender,
    PC1_SUM, PC2_SUM, PC3_SUM, PC4_SUM, PC5_SUM
  ) %>%
  na.omit()

dd <- datadist(tmp)
options(datadist = "dd")

fit <- ols(
  BCF_PRS ~ rcs(BMI0, 4) +
    age+
    gender+
    PC1_SUM+PC2_SUM+PC3_SUM+PC4_SUM+PC5_SUM,
  data = data_all_pheno_cox7_eur
)

anova(fit)

pred <- Predict(fit, BMI0)

p <- ggplot(pred,
       aes(BMI0, yhat)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = lower,
        ymax = upper),
    alpha = 0.2
  ) +
  theme_bw()
p

data_all_pheno_cox7_eur <- data_all_pheno_cox7_eur %>%
  mutate(
    BMI_group = case_when(
      BMI0 < 18.5 ~ "Underweight",
      BMI0 < 25 ~ "Normal",
      BMI0 < 30 ~ "Overweight",
      BMI0 < 40 ~ "Obesity",
      BMI0 >= 40 ~ "Severe_Obesity"
    )
  )

data_all_pheno_cox7_eur$BMI_group = factor(
  data_all_pheno_cox7_eur$BMI_group,
  levels = c(
    "Normal",
    "Underweight",
    "Overweight",
    "Obesity",
    "Severe_Obesity"
  )
)

fit <- lm(
  TCF_PRS ~ BMI_group +
   age+
    gender+
    PC1_SUM+PC2_SUM+PC3_SUM+PC4_SUM+PC5_SUM,
  data = data_all_pheno_cox7_eur
)

summary(fit)

library(broom)

res <- tidy(fit, conf.int = TRUE) %>%
  filter(grepl("^BMI_group", term))

res

fit <- lm(
  BCF_PRS ~ BMI_group +
   age+
    gender+
    PC1_SUM+PC2_SUM+PC3_SUM+PC4_SUM+PC5_SUM,
  data = data_all_pheno_cox7_eur
)

summary(fit)

library(broom)

res <- tidy(fit, conf.int = TRUE) %>%
  filter(grepl("^BMI_group", term))

res


