library(broom)
library(scales)
library(stringr)
library(survival)
library(tidyr)
library(dplyr)
library(haven)
library(reshape2)
library(tidyverse)
library(data.table)
library(lubridate)
library(ggplot2)

v8_cova <- fread("v8_blood_custom_pcs_cfs_covs.txt")
head(v8_cova)

colnames(v8_cova)

data_all_pheno_cox4 <- readRDS("data_stroke_v8_all.rds")

dim(data_all_pheno_cox4)

summary(data_all_pheno_cox4$BMI)

data_all_pheno_cox5 <- merge(v8_cova, data_all_pheno_cox4, by.x = "IID", by.y ="person_id", all = F)
dim(data_all_pheno_cox5)

summary(data_all_pheno_cox5$biosample_collection_date.x)

data_BMI_original_v8 <- readRDS("data_BMI_original_v8.rds")

data_BMI_original_v8 <- filter(data_BMI_original_v8, data_BMI_original_v8$person_id %in% data_all_pheno_cox5$IID)

head(data_BMI_original_v8)

library(dplyr)
library(lubridate)

bmi <- data_BMI_original_v8 %>%
  mutate(
    measurement_datetime = as.POSIXct(measurement_datetime)
  )

pheno <- data_all_pheno_cox5 %>%
  mutate(
    biosample_collection_date.x =
      as.POSIXct(biosample_collection_date.x)
  )

bmi_nearest <- pheno %>%
  select(IID, biosample_collection_date.x) %>%
  left_join(
    bmi,
    by = c("IID" = "person_id")
  ) %>%
  mutate(
    time_diff_days =
      abs(as.numeric(
        difftime(
          measurement_datetime,
          biosample_collection_date.x,
          units = "days"
        )
      ))
  ) %>%
  group_by(IID) %>%
  slice_min(
    order_by = time_diff_days,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

data_all_pheno_cox5 <- data_all_pheno_cox5 %>%
  left_join(
    bmi_nearest %>%
      select(
        IID,
        BMI0 = value_as_number,
        BMI0_date = measurement_datetime,
        BMI0_diff_days = time_diff_days
      ),
    by = "IID"
  )

summary(data_all_pheno_cox5$BMI0)

data_all_pheno_cox5 <- data_all_pheno_cox5 %>%
  mutate(
    BMI0 = ifelse(
      BMI0 < 10 | BMI0 > 60,
      NA,
      BMI0
    )
  )

hist(data_all_pheno_cox5$BMI0, breaks = 100)

quantile(data_all_pheno_cox5$BMI0,
         c(0.99,0.995,0.999,0.9995),
         na.rm=TRUE)

library(ggplot2)

ggplot(data_all_pheno_cox5,
       aes(x = BMI0, y = TCF)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", color = "red") +
  theme_bw()

cor.test(
  data_all_pheno_cox5$BMI,
  data_all_pheno_cox5$TCF,
  method = "spearman",
  use = "complete.obs"
)

data_all_pheno_cox5 %>%
  filter(!is.na(BMI0), !is.na(TCF)) %>%
  mutate(BMI_bin = cut_number(BMI0, 50)) %>%
  group_by(BMI_bin) %>%
  summarise(
    BMI_mean = mean(BMI0),
    TCF_mean = mean(TCF),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(BMI_mean, TCF_mean)) +
  geom_line() +
  geom_point() +
  theme_bw()

tmp <- data_all_pheno_cox5 %>%
  dplyr::select(
    TCF, BCF,
    BMI0,
    Age.x,
    Sex.x,
    PC1.x, PC2.x, PC3.x, PC4.x, PC5.x
  ) %>%
  mutate(
    TCF_z = as.numeric(scale(TCF)),
    BCF_z = as.numeric(scale(BCF))
  ) %>%
  na.omit()

library(rms)
dd <- datadist(tmp)
options(datadist = "dd")

fit <- ols(
  TCF ~ rcs(BMI0, 4) +
    Age.x +
    Sex.x +
    PC1.x + PC2.x + PC3.x + PC4.x + PC5.x,
  data = data_all_pheno_cox5
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

ggsave(
  "BMI_TCF_spline.svg",
  plot = p,
  width = 6,
  height = 6
)

library(rms)

tmp <- data_all_pheno_cox5 %>%
  dplyr::select(
    BCF,
    BMI0,
    Age.x,
    Sex.x,
    PC1.x, PC2.x, PC3.x, PC4.x, PC5.x
  ) %>%
  na.omit()

dd <- datadist(tmp)
options(datadist = "dd")

fit <- ols(
  BCF ~ rcs(BMI0, 4) +
    Age.x +
    Sex.x +
    PC1.x + PC2.x + PC3.x + PC4.x + PC5.x,
  data = data_all_pheno_cox5
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

ggsave(
  "BMI_BCF_spline.svg",
  plot = p,
  width = 6,
  height = 6
)

fit_tcf <- ols(
  TCF_z ~ rcs(BMI0, 4) +
    Age.x +
    Sex.x +
    PC1.x + PC2.x + PC3.x + PC4.x + PC5.x,
  data = tmp
)

fit_bcf <- ols(
  BCF_z ~ rcs(BMI0, 4) +
    Age.x +
    Sex.x +
    PC1.x + PC2.x + PC3.x + PC4.x + PC5.x,
  data = tmp
)

pred_tcf <- as.data.frame(
  Predict(fit_tcf, BMI0)
) %>%
  mutate(Cell = "TCF_z")

pred_bcf <- as.data.frame(
  Predict(fit_bcf, BMI0)
) %>%
  mutate(Cell = "BCF_z")

pred_all <- bind_rows(
  pred_tcf,
  pred_bcf
)

p <- ggplot(
  pred_all,
  aes(
    x = BMI0,
    y = yhat,
    color = Cell,
    fill = Cell
  )
) +
  geom_line(size = 1.2) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15,
    color = NA
  ) +
  scale_color_manual(
    values = c(
      "TCF_z" = "dodgerblue3",
      "BCF_z" = "darkorange2"
    )
  ) +
  scale_fill_manual(
    values = c(
      "TCF_z" = "dodgerblue3",
      "BCF_z" = "darkorange2"
    )
  ) +
  labs(
    x = "BMI",
    y = "Adjusted cell fraction"
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank()
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = 18.5,
    ymin = -Inf, ymax = Inf,
    alpha = 0.1,
    fill = "grey50"
  ) +
  annotate(
    "rect",
    xmin = 25, xmax = 30,
    ymin = -Inf, ymax = Inf,
    alpha = 0.1,
    fill = "gold"
  ) +
  annotate(
    "rect",
    xmin = 30, xmax = 40,
    ymin = -Inf, ymax = Inf,
    alpha = 0.1,
    fill = "orange"
  ) +
  annotate(
    "rect",
    xmin = 40, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    alpha = 0.1,
    fill = "red"
  )

p

ggsave(
  "BMI_TCF_BCF_spline.svg",
  plot = p,
  width = 6,
  height = 3
)

data_all_pheno_cox5 <- data_all_pheno_cox5 %>%
  mutate(
    BMI_group = case_when(
      BMI0 < 18.5 ~ "Underweight",
      BMI0 < 25 ~ "Normal",
      BMI0 < 30 ~ "Overweight",
      BMI0 < 40 ~ "Obesity",
      BMI0 >= 40 ~ "Severe_Obesity"
    )
  )

table(data_all_pheno_cox5$BMI_group, useNA = "always")

data_all_pheno_cox5$BMI_group = factor(
  data_all_pheno_cox5$BMI_group,
  levels = c(
    "Normal",
    "Underweight",
    "Overweight",
    "Obesity",
    "Severe_Obesity"
  )
)

data_all_pheno_cox5 <- data_all_pheno_cox5 %>%
  mutate(
    TCF_z = scale(TCF)[,1]
  )

data_all_pheno_cox5 <- data_all_pheno_cox5 %>%
  mutate(
    BCF_z = scale(BCF)[,1]
  )

fit <- lm(
  TCF_z ~ BMI_group +
    Age.x +
    Sex.x +
    PC1.x + PC2.x + PC3.x + PC4.x + PC5.x,
  data = data_all_pheno_cox5
)

summary(fit)

library(broom)

res <- tidy(fit, conf.int = TRUE) %>%
  filter(grepl("^BMI_group", term))

res

fit <- lm(
  BCF_z ~ BMI_group +
    Age.x +
    Sex.x +
    PC1.x + PC2.x + PC3.x + PC4.x + PC5.x,
  data = data_all_pheno_cox5
)

summary(fit)

library(broom)

res <- tidy(fit, conf.int = TRUE) %>%
  filter(grepl("^BMI_group", term))

res


