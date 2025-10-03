library(tidyverse)
library(cmdstanr)
library(lattice)
library(glmnet)
library(caret)
# library(MatchIt)
# library(JOPS)
library(diversityForest)
library(splines2)
library(scam)
library(CalibrationCurves)
library(betacal)

source("functions.R")
source("new_functions.R")


d <- read_csv("../data_from_ICR/Cancer Model Training Data 20250228.csv")
##d <- read_csv("../data_from_ICR/Cancer Model Training Data 20250602.csv")

d <- d %>%
  mutate_if(is.character, as.factor)

icr <- d %>%
  filter(Vendor == "ICR") %>%
  droplevels()
icr <- janitor::clean_names(icr)

test <- read.csv("cleaned_data/different_test_data.csv") %>%
  mutate_if(is.character, as.factor) %>%
  filter(training_or_blinded != "Z_Do_Not_Analyze_Hemolyzed") %>%
  droplevels()


test <- test %>%
  select(
    vendor, serum_specific_id, collection_date, collection_country,
    age_at_collection, sex, cancer, disease, targeted_disease_for_collection, afp:tgm2, cg_a,
    scc, thrombospondin_2, free_psa_access:ferritin_access, nse_kryptor_corrected, aldh1a1_corrected,
    nse_corrected, percent_free_psa_access,
    billirubin_index_serum
  ) %>%
  select(-c(
    hgf, beta_hcg, he4, s_fasl, total_psa, aldh1a1_not_corrected,
    nse_not_corrected
  )) %>%
  drop_na(total_psa_access, free_psa_access, percent_free_psa_access, ferritin, trail) %>%
  dplyr::rename(bilirubin = billirubin_index_serum) %>%
  mutate(
    nse_kryptor_corrected = as.numeric(as.character(nse_kryptor_corrected)),
    aldh1a1_corrected = as.numeric(as.character(aldh1a1_corrected)),
    nse_corrected = as.numeric(as.character(nse_corrected)),
    percent_free_psa_access = as.numeric(as.character(percent_free_psa_access))
  ) %>%
  droplevels()

## rename first few columns
names(test)[1:5] <- c("vendor", "id", "collection_date", "country", "age")


## tidy disease names
d <- d %>%
  # mutate(Disease=gsub(" cancer", "", Disease, ignore.case = TRUE)) %>%
  filter(Race %in% c("White", "Hispanic")) %>%
  filter(!Stage %in% c("3", "4")) %>%
  filter(Vendor != "ICR") %>%
  droplevels()

d <- janitor::clean_names(d)

d2 <- clean_raw_data(d)


## move this to functions.R
d2$age <- as.numeric(as.character(d2$age))
prostate <- tidy_prostate_data(d2, filter_age = 0)
prostate <- prostate %>%
  mutate(percent_free_psa_access = as.numeric(as.character(percent_free_psa_access))) %>%
  mutate_if(is.character, as.factor) %>%
  filter(!id %in% c("HMN686019S", "HMN456748")) %>%
  as.data.frame()

## ------------------------------------------------------------
## Bayes model
## ------------------------------------------------------------

## Bernouuli likelihood
l <- function(p, y) {
  ## p <- ifelse(p < 1e-5, 1e-5, p)
  ## p <- ifelse(p > 1 - 1e-5, 1 - 1e-5, p)
  p^(y) * (1 - p)^(1 - y)
}

l(0.001, 1)

## prostate$bound_psa <- prostate$total_psa - prostate$free_psa_access
## prostate$bound_psa[prostate$bound_psa < 0] <- 0

## which samples have no PSA data
no_psa_index <- !apply(prostate[, grepl("psa", names(prostate))], 1, function(x) all(is.na(x)))
prostate <- prostate[no_psa_index, ]

## remove samples with too many missing values before imputing
few_missing <- apply(prostate, 2, function(x) sum(is.na(x))) / nrow(prostate)
few_missing <- few_missing[few_missing < 0.2]

prostate <- select(prostate, names(few_missing))

test <- select(test, names(prostate))
## test$vendor <- factor(ifelse(test$vendor == "BioIVT", "BioIVT", "Other"))
## test <- test %>%
##     mutate(country = recode_factor(country,
##       "Russia" = "Ukraine",
##       "Unknown" = "Turkey",
##       "Vietnam" = "Turkey"
##     ))

## test <- test %>%
##     filter(!(country == "Turkey" & vendor == "Other")) %>% dim

train1 <- prostate %>%
  filter(age > 44.5 & age < 85) %>%
  select(names(few_missing)) %>%
  # preproc2(., ratio=2)
  preproc(., test, index = 9)

## train1 %>%
##     mutate(id=paste0("id_", row_number())) %>%
##     select(-c(collection_date, sex, targeted_disease_for_collection, disease)) %>%
##     ## recode vendor
##         mutate(vendor = recode_factor(vendor,
##                 "BioIVT" = "Vendor A",
##                 "Other" = "Vendor B"
##                 )) %>%
##     ## recode country
##         mutate(country = recode_factor(country,
##                 "Ukraine" = "Eastern Europe",
##                 )) %>%
##     write.csv(., "/home/sel/prioris/projects/2025/Imperial_MSc/data/prostate_cancer_data.csv",
##               row.names = FALSE)

#test <- d[!as.character(d$serum_specific_id) %in% as.character(train1$train$id), ]
## ids <- test %>%
##     #filter(sex == "M") %>%
##     #filter(disease %in% c("Healthy", "Prostate Cancer")) %>%
##     #filter(vendor != "ICR") %>%
##     mutate(id=paste0("id_", row_number())) %>%
##     select(-c(collection_date, sex, targeted_disease_for_collection, disease)) %>%
##     ## recode vendor
##         mutate(vendor = recode_factor(vendor,
##                 "BioIVT" = "Vendor A",
##                 "Other" = "Vendor B"
##                 )) %>%
##     ## recode country
##         mutate(country = recode_factor(country,
##                 "Ukraine" = "Eastern Europe",
##                 ))

## write.csv(ids, "cleaned_data/prostate_test_data.csv", row.names=FALSE)

d %>%
  mutate(percent_free_psa_access = as.numeric(as.character(percent_free_psa_access))) %>%
  filter(sex == "M") %>%
  filter(vendor != "ICR") %>%
  filter(disease %in% c("Healthy", "Prostate Cancer")) %>%
  xyplot(percent_free_psa_access ~ factor(cancer) | factor(I(serum_specific_id %in% prostate$id)),
    data = ., # groups=,
    type = c("g", "p", "a"),
    ylab = "", xlab = "", main = "",
    jitter.x = TRUE,
    scales = list(alternating = FALSE), auto.key = list(columns = 2)
  )


mt <- icr %>%
  filter(order_id == "MT")

names(mt)[c(5, 9, 15, 16, 19)] <- c("vendor", "id", "collection_date", "country", "age")

mt <- mt[, names(mt) %in% names(prostate)]
mt2 <- preproc(prostate, mt, index = 9)

plot(log10(train1$total_psa_access + 1) ~ log10(train1$percent_free_psa_access + 1),
  col = train1$cancer + 1,
  pch = ifelse(train1$cancer == 0, 1, 16)
)

points(log10(as.numeric(mt$total_psa_access) + 1) ~ log10(as.numeric(mt$percent_free_psa_access) + 1),
  col = "steelblue", cex = 1, pch = 17
)

plot(log10(train1$ferritin + 1) ~ log10(train1$trail + 1),
  col = train1$cancer + 1,
  pch = ifelse(train1$cancer == 0, 1, 16)
)
points(log10(as.numeric(mt$ferritin) + 1) ~ log10(as.numeric(mt$trail) + 1),
  col = "steelblue", cex = 1, pch = 17
)

## how to account for ICR vendor effect???

## train2 <- train1$train %>%
##   data.frame()

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  lm(trail ~ vend_co + age, data = .) %>%
  car::Anova()

#train1$train %>%

train2 %>%
#  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(total_psa_access ~ train1$train$age,
    data = .,
    type = c("g", "p", "r"), group = train1$train$cancer,
    ylab = "biomarker", xlab = "age", # main = "Total PSA (Access)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right", title = "Cancer"),
  )

#train2 %>%
train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(trail ~  vend_co,
    data = .,
    type = c("g", "p", "a"), group = train1$train$cancer,
    ylab = "biomarker", xlab = "age", # main = "Total PSA (Access)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right", title = "Cancer"),
  )



## cal residuals to remove bath effests
train2 <- correct_batch_effects(train1$train, train1$train)$train_resid %>%
  data.frame()

train2b <- correct_batch_effects(train1$train, train1$train)$train_resid %>%
  data.frame()

plot(train2$ferritin ~ train2b$ferritin,
     col = factor(paste0(train1$train$vendor, "_", train1$train$country)))
abline(0, 1, col = "grey", lwd = 2)

##write.csv(cbind(train1$train[, 1:9], train2), "cleaned_data/prostate.csv", row.names = FALSE)

## change point size and colour for groups
trellis.par.set(
  superpose.symbol = list(pch = c(1, 16, 17), cex = c(0.5, 1, 1.5), col = c("steelblue", "firebrick", "green")),
  plot.symbol = list(pch = c(1, 16, 17), cex = c(0.5, 1, 1.5), col = c("steelblue", "firebrick", "green"))
)

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$percent_free_psa_access ~ train2$total_psa_access,
    data = .,
    type = c("g", "p"), group = cancer,
    ylab = "% Free PSA", xlab = "Total PSA", # main = "Total PSA (Access)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right", title = "Cancer"),
  )


train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$p2psa_access ~ train2$free_psa_access,
    data = .,
    type = c("g", "p"), group = cancer,
    ylab = "P2PSA", xlab = "Free PSA", # main = "Total PSA (Access)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right", title = "Cancer"),
  )


train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$p2psa_access ~ factor(cancer),
    data = .,
    type = c("g", "p"), group = cancer,
    ylab = "P2PSA", xlab = "Free PSA", # main = "Total PSA (Access)",
    jitter.x = TRUE, factor = 1.5, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right", title = "Cancer"),
  )


pca <- prcomp(select(train2, total_psa_access, free_psa_access), center = TRUE, scale. = TRUE)
plot(pca$x[, 1:2], col = train1$train$cancer + 1, pch = ifelse(train1$train$cancer == 0, 1, 16))

## write.csv(data.frame(train1$train[, 1:9], train2), "cleaned_data/cleaned_data.csv",
##           row.names=FALSE)

inds <- which(train1$train$cancer == 1 & train1$train$total_psa_access < -1)
train1$train[inds, 1:9]

par(las = 1)
d2 %>%
  filter(sex == "M") %>%
  filter(vendor != "ICR") %>%
  filter(disease != "Prostate Cancer") %>%
  #    mutate(total_psa_access=log10(total_psa_access + 1)) %>%
  pull(total_psa_access) %>%
  #    round(., 2) %>%
  log10() %>%
  round(., 2) %>%
  table() %>%
  plot(., xlab = "Log10 Total PSA (Access)", ylab = "Frequency of unique values", main = "Total PSA (Access)", col = "steelblue", cex = 0.5)
#    BHH2::dotPlot(.)
mean(tmp < 0.01, na.rm = T)

d2 %>%
  filter(id %in% c("HMN686019S", "HMN456748")) %>%
  # mutate(total_psa_access=log10(total_psa_access + 1)) %>%
  pull(total_psa_access) %>%
  abline(v = ., col = "red", lty = 2, lwd = 2)


## HERE add key variable plots to pres
train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$total_psa_access ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "Total PSA (Access)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )


train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$percent_free_psa_access ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "% Free PSA (Access)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$p2psa_access ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "P2PSA (Access)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$ferritin ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "Ferritin",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$trail ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "TRAIL",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$ca19_9 ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "CA19_9",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$prolactin ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "Prolactin",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$midkine ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "Midkine",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$nse_kryptor_corrected ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "NSE (corrected)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )


train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$percent_free_psa_access ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "NSE (corrected)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )


## HERE: try without chol decomp for psa... maybe tri Mean difference
prostate %>%
  filter(age > 44.5 & age < 85) %>%
  # train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(as.numeric(percent_free_psa_access) ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "NSE (corrected)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )

prostate %>%
  filter(age > 44.5 & age < 85) %>%
  # train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(as.numeric(phi_score_access) ~ factor(cancer),
    data = .,
    type = c("g", "p", "a"), group = vend_co,
    ylab = "Biomarker level", xlab = "Cancer", main = "NSE (corrected)",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 1, space = "right")
  )


prostate %>%
  filter(age > 44.5 & age < 85) %>%
  xyplot(log10(free_psa_access + 1) ~ factor(cancer),
    data = .,
    jitter.x = T
  )

prostate %>%
  filter(age > 44.5 & age < 85) %>%
  xyplot(log10(total_psa_access - free_psa_access + 1) ~ log10(total_psa_access + 1), data = ., groups = cancer)

prostate %>%
  filter(age > 44.5 & age < 85) %>%
  xyplot(log10(free_psa_access / (total_psa_access)) ~ log10(total_psa_access + 1), data = ., groups = cancer)

pairs(log10(select(prostate, total_psa_access, free_psa_access, p2psa_access) + 1), col = prostate$cancer + 1)

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$free_psa_access ~ vend_co,
    data = ., groups = cancer,
    type = c("g", "p", "a"),
    ylab = "", xlab = "Cancer", main = "",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 2, title = "Cancer")
  )

train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(train2$free_psa_access ~ vend_co,
    data = ., groups = cancer,
    type = c("g", "p", "a"),
    ylab = "", xlab = "Cancer", main = "",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
    auto.key = list(columns = 2, title = "Cancer")
  )




class_wt <- 1 / prop.table(table(train1$train$cancer))
obs_wt <- ifelse(train1$train$cancer == 0, class_wt[1], class_wt[2])

set.seed(123)
s1 <- scam(
  train1$train$cancer ~ s(total_psa_access, bs = "mpi") +
    s(percent_free_psa_access, bs = "mpd") +
    s(trail, bs = "mpd") +
    s(ferritin, bs = "mpi"),
  weights = round(obs_wt, 0),
  data = train2,
  family = "binomial",
  gamma = 1.4
)


summary(s1)

par(las = 1)
plot(s1, pages = 1, scale = 0, scheme = 1)

gam_predictions <- rep(NA, nrow(train2))

tmp <- prostate %>%
  filter(age > 44.5 & age < 85)



for (i in 1:nrow(train2)) {
  ## make sure run few missing
  train1 <- preproc(tmp[-i, ], tmp[i, ], index = 9)

  train2 <- correct_batch_effects(train1$train, train1$train)$train_resid %>%
    data.frame()

  class_wt <- 1 / prop.table(table(train1$train$cancer))
  obs_wt <- ifelse(train1$train$cancer == 0, class_wt[1], class_wt[2])

  set.seed(123)
  s1 <- scam(
    train1$train$cancer ~ s(total_psa_access, bs = "mpi") +
      s(percent_free_psa_access, bs = "mpd") +
      s(trail, bs = "mpd") +
      s(ferritin, bs = "mpi"),
    weights = round(obs_wt, 0),
    data = train2,
    family = "binomial",
    gamma = 1.4
  )

  gam_predictions[i] <- predict(s1, newdata = train1$test, type = "response")
}

## train2 %>%
# mutate(vend_co=factor(paste0(vendor, "_", country))) %>%
xyplot(s1$fitted.values ~ factor(train1$train$cancer),
  # groups=train1$train$targeted_disease_for_collection ==  "Benign Prostatic Hyperplasia",
  # groups=train1$train$vendor,
  # groups=train1$train$country,
  # groups=train1$train$age > 65,
  ## data=.,
  type = c("g", "p", "a"),
  xlab = "", ylab = "Log10 Free PSA (access)", main = "",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0), x = list(rot = 45)),
  auto.key = list(columns = 2, title = " Prostate Cancer")
)

xyplot(gam_predictions ~ factor(train1$train$cancer),
  # groups=train1$train$targeted_disease_for_collection ==  "Benign Prostatic Hyperplasia",
  # groups=train1$train$vendor,
  # groups=train1$train$country,
  # groups=train1$train$age > 65,
  groups = factor(paste0(train1$train$vendor, "_", train1$train$country)),
  ## data=.,
  type = c("g", "p", "a"),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "GAM model",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

calc_binomial_metrics(train1$train$cancer, s1$fitted.values, thresh = 0.5, filter = NULL)
calc_binomial_metrics(train1$train$cancer, gam_predictions, thresh = 0.5, filter = NULL)
##      auc brier_mean brier_median logloss_mean logloss_median accuracy sensitivity specificity  PPV
## 1 0.9153     0.1072       0.0025       0.3948         0.0507    0.855        0.84        0.87 0.42
##      NPV    PLR    NLR   F1    MCC     P4 N_healthy N_cancer
## 1 0.9798 6.4593 0.1839 0.56 0.5328 0.3952       223       25

## diverse pred
div_pred <- predict(s1, newdata = train1$test, type = "response")
xyplot(div_pred ~ factor(test$cancer),
  # groups=train1$train$targeted_disease_for_collection ==  "Benign Prostatic Hyperplasia",
  # groups=test$vendor,
  # groups = train1$test$country,
  # groups=test$age > 65,
  # groups = factor(paste0(train1$train$vendor, "_", train1$train$country)),
  ## data=.,
  type = c("g", "p", "a"),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "GAM model",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)
calc_binomial_metrics(test$cancer, div_pred, thresh = 0.5, filter = NULL)
##      auc brier_mean brier_median logloss_mean logloss_median accuracy sensitivity specificity
## 1 0.8152     0.2008        6e-04       1.2261         0.0254   0.7519      0.6667      0.8372
##      PPV    NPV    PLR    NLR     F1    MCC     P4 N_healthy N_cancer
## 1 0.4615 0.9231 4.0952 0.3981 0.5455 0.4402 0.4948        43        9

train1$test$div_pred <- div_pred

train1$test[div_pred < 0.2 & test$cancer == 1, ] %>%
  select(
    id, vendor, country, cancer, targeted_disease_for_collection,
    total_psa_access, percent_free_psa_access, ferritin, trail,
    div_pred
  ) %>%
  View()

train1$test[div_pred > 0.8 & test$cancer == 0, ] %>%
  select(
    id, vendor, country, cancer, targeted_disease_for_collection,
    total_psa_access, percent_free_psa_access, ferritin, trail, div_pred
  ) %>%
  View()


calc_binomial_metrics(train1$train$cancer, gam_predictions,
  thresh = 0.5,
  filter = train1$train$targeted_disease_for_collection != "Benign Prostatic Hyperplasia"
)
## 1 0.916     0.1078       0.0027        0.403         0.0526   0.8565         0.5   -11362347 0.5
##   NPV PLR NLR  F1 MCC P4 N_healthy N_cancer
## 1 0.5   0   0 0.5   0 NA       189       26

calc_binomial_metrics(train1$train$cancer, gam_predictions,
  thresh = 0.5,
  filter = train1$train$vendor == "Other"
)
##      auc brier_mean brier_median logloss_mean logloss_median accuracy sensitivity specificity PPV
## 1 0.9125      0.094       0.0025       0.3599         0.0511   0.8873         0.5   -32537631 0.5
##   NPV PLR NLR  F1 MCC P4 N_healthy N_cancer
## 1 0.5   0   0 0.5   0 NA        66       18

calc_binomial_metrics(train1$train$cancer, gam_predictions,
  thresh = 0.5,
  filter = train1$train$vendor == "BioIVT"
)
##    auc brier_mean brier_median logloss_mean logloss_median accuracy sensitivity specificity
## 0.8607     0.1661       0.0105       0.6184         0.0898   0.7552       0.625      0.8854
##    PPV    NPV    PLR    NLR     F1    MCC     P4 N_healthy N_cancer
## 0.2174 0.9789 5.4514 0.4236 0.3226 0.3165 0.2634       157        8

calc_binomial_metrics(train1$train$cancer, gam_predictions,
  thresh = 0.5,
  filter = train1$train$country == "USA"
)
##      auc brier_mean brier_median logloss_mean logloss_median accuracy sensitivity specificity
## 1 0.9672     0.0546       0.0016        0.188         0.0374   0.9324           1      0.8649
##      PPV NPV PLR NLR     F1    MCC P4 N_healthy N_cancer
## 1 0.7368   1 7.4   0 0.8485 0.7983  1        37       14


calc_binomial_metrics(train1$train$cancer, gam_predictions,
  thresh = 0.5,
  filter = train1$train$country == "Ukraine"
)
## 1 0.9051     0.1255       0.0049       0.4432         0.0689   0.8129         0.5   -27183337 0.5
##   NPV PLR NLR  F1 MCC P4 N_healthy N_cancer
## 1 0.5   0   0 0.5   0 NA        79       10

nd <- expand.grid(
  total_psa_access = seq(min(train2$total_psa_access),
    max(train2$total_psa_access),
    length.out = 20
  ),
  percent_free_psa_access = seq(min(train2$percent_free_psa_access),
    max(train2$percent_free_psa_access),
    length.out = 20
  )
)
nd$trail <- 0.61
nd$ferritin <- 0

nd$p_cancer <- predict(s1, newdata = nd, type = "response")

# x11()
wireframe(p_cancer ~ total_psa_access * percent_free_psa_access,
  data = nd, drape = TRUE,
  colorkey = TRUE, zlim = c(0, 1), zlab = "P(Cancer)",
  ylab = "% Free PSA", xlab = "Total PSA"
) # , screen = list(z = 70, x = -70))

xyplot(trail ~ factor(train1$train$cancer) | cut(total_psa_access, breaks = c(-5, -1, 1, 5)),
  data = train2,
  type = c("g", "p"), jitter.x = TRUE,
  xlab = "Cancer", ylab = "Trail", main = "Trail"
)

xyplot(
  train2$trail ~ train2$total_psa_access | cut(train2$total_psa_access,
    breaks = c(-5, -1, 1, 5)
  ),
  groups = train1$train$cancer
)


## ============================================================
compiled_model <- cmdstan_model("stan_models/prostate_final_model.stan",
  dir = "compiled_models"
)

## train2 <- train1$train
## train2$vend_co <- as.numeric(factor(paste0(train2$vendor, "_", train2$country)))

plot(bound_psa ~ p2psa_access, data = train2, pch = 20, col = train1$train$cancer + 1)

## train2$psa_interaction <- train2$total_psa_access * train2$percent_free_psa_access

B1 <- iSpline(train2$total_psa_access, df = 6, intercept = TRUE)
B2 <- iSpline(train2$percent_free_psa_access, df = 6, intercept = TRUE)
B3 <- iSpline(train2$trail, df = 6, intercept = TRUE)
B4 <- iSpline(train2$ferritin, df = 6, intercept = TRUE)

## X <- cbind(B1, B2, B3, B4)
## H = X %*% MASS::ginv(t(X) %*% X) %*% t(X)
## plot(diag(H))

## diff_order <- 1
## D <- diff(diag(ncol(B1)), diff = diff_order)

## make_weights <- function(y, pred, scale_factor = 1) {
##     abs_diff <- abs(y - pred)
##     abs_diff[abs_diff < 0.2] <- 1L
##     abs_diff[abs_diff >= 0.2 & abs_diff < 0.6] <- 2L
##     abs_diff[abs_diff >= 0.6 & abs_diff < 0.8] <- 3L
##     abs_diff[abs_diff >= 0.8 & abs_diff < 0.999999999] <- 4L
##     return(abs_diff * scale_factor)
## }

etas_base <- etas
sum(make_weights(train1$train$cancer, etas_base, 5))

dat <- list(
  P = ncol(B1),
  N_train = nrow(train2),
  X1_train = B1,
  X2_train = B2,
  X3_train = B3,
  X4_train = B4,
  y_train = train1$train$cancer,
  w = rep(1, nrow(B1))#calc_weights(train1$train$cancer)# * make_weights(train1$train$cancer, etas_base, 5)
)

m1 <- compiled_model$sample(
  data = dat,
  chains = 5,
  parallel_chains = 5,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  seed = 123,
  refresh = 0, 
  fixed_param = FALSE
)

bs <- matrix(NA, nrow=1000, ncol=25)
for(i in 1:1000){
    m1 <- compiled_model$sample(
                             data = dat,
                             chains = 1,
                             parallel_chains = 1,
                             iter_warmup = 0,
                             iter_sampling = 1,
                             refresh = 0, 
                             fixed_param = TRUE
                         )
    dr <- m1$draws(format = "df") %>%
        data.frame()
    
    bs[i, ] <- dr %>% select_at(vars(starts_with("b"))) %>% as.matrix
}


dr <- m1$draws(format = "df") %>%
  data.frame()


etas <- dr %>%
  select_at(vars(starts_with("linpred"))) %>%
  invlogit() %>%
  colMeans()

dr %>% select_at(vars(starts_with("b1"))) %>%
    boxplot(., outline=FALSE)

par(las=1)
boxplot(bs[, 2:7], outline=FALSE, main="Total PSA")

plot(apply(etas, 2, sd) ~ colMeans(etas))

which.max(apply(etas, 2, sd))
plot(density(etas[, 143]), xlim = c(0, 1), main = "Distribution of etas")

xyplot(etas ~ factor(dat$y_train),
  ## groups=train2$targeted_disease_for_collection ==  "Benign Prostatic Hyperplasia",
  ##groups = paste0(train1$train$vendor, "_", train1$train$country),
  ##groups=ifelse(diag(H) > 0.2, "High", "Low"),
  ## groups=make_weights(train1$train$cancer, etas), 
  type = c("g", "p", "a"),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "Bayes GAM model",
  ylim = c(-0.03, 1.03),
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

calc_binomial_metrics(dat$y_train, etas, thresh = 0.5, filter = NULL)
##      auc brier_mean brier_median logloss_mean logloss_median accuracy sensitivity specificity
## 1 0.9372     0.0872       0.0091       0.3018         0.0981   0.8747       0.875      0.8744
##      PPV    NPV    PLR    NLR     F1    MCC     P4 N_healthy N_cancer
## 1 0.4286 0.9848 6.9687 0.1429 0.5753 0.5566 0.7098       223       24


val.prob.ci.2(etas, dat$y_train)
val.prob.ci.2(s1$fitted.values, dat$y_train)

bc <- cal_estimate_isotonic(data.frame(.pred_ = etas, y = dat$y_train), y,
  metrics = cls_met, save_pred = TRUE
)

eta_cal <- beta_predict(etas, bc)

train1$train[etas > 0.8 & train1$train$cancer == 0, ] %>%
  select(
    id, vendor, country, cancer, targeted_disease_for_collection,
    total_psa_access, percent_free_psa_access, ferritin, trail
  ) %>%
  View()




## bayes_gam_predictions_mean <- rep(NA, nrow(train2))
## bayes_gam_predictions_median <- rep(NA, nrow(train2))
bayes_gam_predictions_mean <- matrix(NA, nrow(train2), ncol = 8)

## bayes_gam_predictions_mean <- cbind(bayes_gam_predictions_mean,
##                                     matrix(NA, nrow(train2), ncol = 4))

tmp <- prostate %>%
    filter(age > 44.5 & age < 85)

tempered <- c(1, 2, 5, 10, 50, 100, 500, 1000)

for (w in 1:8) {
  for (i in 1:nrow(train2)) {
    ## make sure run few missing
    t1 <- preproc(tmp[-i, ], tmp[i, ], index = 9)

    t2 <- correct_batch_effects(t1$train, t1$train)$train_resid %>%
      data.frame()

    B1 <- iSpline(t2$total_psa_access, df = 6, intercept = TRUE)
    B2 <- iSpline(t2$percent_free_psa_access, df = 6, intercept = TRUE)
    B3 <- iSpline(t2$trail, df = 6, intercept = TRUE)
    B4 <- iSpline(t2$ferritin, df = 6, intercept = TRUE)

    dat <- list(
      P = ncol(B1),
      N_train = nrow(t2),
      X1_train = B1,
      X2_train = B2,
      X3_train = B3,
      X4_train = B4,
      y_train = t1$train$cancer,
      w = calc_weights(t1$train$cancer) * tempered[w]
    )

    m1 <- compiled_model$sample(
      data = dat,
      chains = 5,
      parallel_chains = 8,
      iter_warmup = 500,
      iter_sampling = 500,
      adapt_delta = 0.95,
      seed = 123,
      refresh = 0
    )


    dr <- m1$draws(format = "df") %>%
      data.frame()


    bayes_gam_predictions_mean[i, w] <- rowMeans(
      predict_prostate(
        dr,
        list(
          predict(B1, newx = t1$test$total_psa_access),
          predict(B2, newx = t1$test$percent_free_psa_access),
          predict(B3, newx = t1$test$trail),
          predict(B4, newx = t1$test$ferritin)
        )
      )
    )
  }
}

write.csv(bayes_gam_predictions_mean, "output/prostate_tempered_outputs.csv", row.names=F)

calc_binomial_metrics(train1$train$cancer, bayes_gam_predictions_mean[, 4], thresh = 0.5, filter = NULL)
##      auc brier_mean brier_median logloss_mean logloss_median accuracy sensitivity specificity
## 1 0.9176     0.0904       0.0061       0.3624          0.081   0.8911      0.9167      0.8655
##      PPV    NPV    PLR    NLR     F1    MCC     P4 N_healthy N_cancer
## 1 0.4231 0.9897 6.8139 0.0963 0.5789 0.5682 0.5683       223       24

res <- matrix(NA, 8, 17)
for (i in 1:8){
    res[i, ] <- calc_binomial_metrics(train1$train$cancer, bayes_gam_predictions_mean[, i]) %>% c() %>% unlist()
}
colnames(res) <- names(calc_binomial_metrics(train1$train$cancer, bayes_gam_predictions_mean[, 1]))
#res <- t(res)
res <- data.frame(res)

par(mfrow=c(3, 3),
    las=1,
    cex=1)
plot(auc ~ log10(tempered), data=res, type="b", main="AUC", ylab="", xlab="Log10(weight)")
plot(brier_mean ~ log10(tempered), data=res, type="b", main="Brier", ylab="", xlab="Log10(weight)")
plot(logloss_mean ~ log10(tempered), data=res, type="b", main="Log-loss", ylab="", xlab="Log10(weight)")
plot(accuracy ~ log10(tempered), data=res, type="b", main="Accuracy", ylab="", xlab="Log10(weight)")
plot(sensitivity ~ log10(tempered), data=res, type="b", main="Sensitivity", ylab="", xlab="Log10(weight)")
plot(specificity ~ log10(tempered), data=res, type="b", main="Specificity", ylab="", xlab="Log10(weight)")
plot(PPV ~ log10(tempered), data=res, type="b", main="PPV", ylab="", xlab="Log10(weight)")
plot(NPV ~ log10(tempered), data=res, type="b", main="NPV", ylab="", xlab="Log10(weight)")
plot(P4 ~ log10(tempered), data=res, type="b", main="P4", ylab="", xlab="Log10(weight)")


pairs(bayes_gam_predictions_mean, col=train1$train$cancer + 1,
      labels=paste0("Weight=", tempered))

xyplot(bayes_gam_predictions_mean[, 4] ~ factor(train1$train$cancer),
  # groups=train1$train$targeted_disease_for_collection ==  "Benign Prostatic Hyperplasia",
  # groups=train1$train$vendor,
  # groups=train1$train$country,
  # groups=train1$train$age > 65,
  groups = factor(paste0(train1$train$vendor, "_", train1$train$country)),
  type = c("g", "p", "a"),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "Weight = 10",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)



## HERE: check which biomarkers need linearity
plot(etas ~ s1$fitted.values, col = train1$train$cancer + 1)
plot(etas ~ s1$fitted.values, col = as.numeric(train1$train$targeted_disease_for_collection))
abline(0, 1)
abline(h = 0.5)

index <- which(etas > 0.4 & s1$fitted.values < 0.05)
index <- train1$train[index, ] %>%
  pull(id) %>%
  as.character()
filter(d, serum_specific_id %in% index) %>%
  View()

dr %>%
  select_at(vars(starts_with("b0"))) %>%
  boxplot(., outline = FALSE)

dr %>%
  select_at(vars(starts_with("b1."))) %>%
  boxplot(., outline = FALSE, ylim = c(-5, 5))
abline(h = 0, col = "red", lty = 2)

dr %>%
  select_at(vars(starts_with("b2."))) %>%
  boxplot(., outline = FALSE, ylim = c(-5, 5))
abline(h = 0, col = "red", lty = 2)

dr %>%
  select_at(vars(starts_with("b3."))) %>%
  boxplot(., outline = FALSE, ylim = c(-5, 5))
abline(h = 0, col = "red", lty = 2)

dr %>%
  select_at(vars(starts_with("b4."))) %>%
  boxplot(., outline = FALSE, ylim = c(-5, 5))
abline(h = 0, col = "red", lty = 2)

## dr %>%
##   select_at(vars(starts_with("b5."))) %>%
##   boxplot(., outline = FALSE, ylim = c(-5, 5))
## abline(h = 0, col = "red", lty = 2)

## dr %>% select_at(vars(starts_with("b6."))) %>%
##     boxplot(., outline=FALSE, ylim=c(-5, 5)); abline(h=0, col="red", lty=2)

dr %>%
  select_at(vars(starts_with("sigma_penalty"))) %>%
  boxplot(., outline = FALSE)

## https://pmc.ncbi.nlm.nih.gov/articles/PMC11835049
## https://pmc.ncbi.nlm.nih.gov/articles/PMC9481733/

par(las = 1,
    mfrow=c(2, 2),
    mar=c(5, 4.5, 1, 1))
## bound_psa
B1_pred <- predict(B1, seq(min(train2$total_psa_access), max(train2$total_psa_access), length.out = 50))
B2_pred <- predict(B2, rep(0, 50))
B3_pred <- predict(B3, rep(0, 50))
B4_pred <- predict(B4, rep(0, 50))

y_pred <- predict_prostate(dr, list(B1_pred, B2_pred, B3_pred, B4_pred))
ci <- apply(y_pred, 1, function(x) quantile(x, c(0.025, 0.975)))

plot(rowMeans(y_pred) ~ seq(min(train2$total_psa_access), max(train2$total_psa_access), length.out = 50),
  type = "l", xlab = "Total PSA", ylab = "Probability of Prostate Cancer", lwd = 2,
  ylim = c(0, 1)
)
lines(ci[1, ] ~ seq(min(train2$total_psa_access), max(train2$total_psa_access), length.out = 50), lty = 2)
lines(ci[2, ] ~ seq(min(train2$total_psa_access), max(train2$total_psa_access), length.out = 50), lty = 2)
rug(train2$total_psa_access)

## percent_free_psa
B1_pred <- predict(B1, rep(0, 50))
B2_pred <- predict(B2, seq(min(train2$percent_free_psa_access), max(train2$percent_free_psa_access), length.out = 50))
B3_pred <- predict(B3, rep(0, 50))
B4_pred <- predict(B4, rep(0, 50))
y_pred <- predict_prostate(dr, list(B1_pred, B2_pred, B3_pred, B4_pred))
ci <- apply(y_pred, 1, function(x) quantile(x, c(0.025, 0.975)))

plot(rowMeans(y_pred) ~ seq(min(train2$percent_free_psa_access), max(train2$percent_free_psa_access), length.out = 50),
  type = "l", xlab = "% Free PSA", ylab = "Probability of Prostate Cancer", lwd = 2,
  ylim = c(0, 1)
)
lines(ci[1, ] ~ seq(min(train2$percent_free_psa_access), max(train2$percent_free_psa_access), length.out = 50), lty = 2)
lines(ci[2, ] ~ seq(min(train2$percent_free_psa_access), max(train2$percent_free_psa_access), length.out = 50), lty = 2)
rug(train2$percent_free_psa_access)


## trail
B1_pred <- predict(B1, rep(0, 50))
B2_pred <- predict(B2, rep(0, 50))
B3_pred <- predict(B3, seq(min(train2$trail), max(train2$trail), length.out = 50))
B4_pred <- predict(B4, rep(0, 50))
y_pred <- predict_prostate(dr, list(B1_pred, B2_pred, B3_pred, B4_pred))
ci <- apply(y_pred, 1, function(x) quantile(x, c(0.025, 0.975)))

plot(rowMeans(y_pred) ~ seq(min(train2$trail), max(train2$trail), length.out = 50),
  type = "l", xlab = "Trail", ylab = "Probability of Prostate Cancer", lwd = 2,
  ylim = c(0, 1)
)
lines(ci[1, ] ~ seq(min(train2$trail), max(train2$trail), length.out = 50), lty = 2)
lines(ci[2, ] ~ seq(min(train2$trail), max(train2$trail), length.out = 50), lty = 2)
rug(train2$trail)



## ferritin
B1_pred <- predict(B1, rep(0, 50))
B2_pred <- predict(B2, rep(0, 50))
B3_pred <- predict(B3, rep(0, 50))
B4_pred <- predict(B4, seq(min(train2$ferritin),
  max(train2$ferritin),
  length.out = 50
))
y_pred <- predict_prostate(dr, list(B1_pred, B2_pred, B3_pred, B4_pred))
ci <- apply(y_pred, 1, function(x) quantile(x, c(0.025, 0.975)))

plot(
  rowMeans(y_pred) ~ seq(min(train2$ferritin),
    max(train2$ferritin),
    length.out = 50
  ),
  type = "l", xlab = "Ferritin", ylab = "Probability of Prostate Cancer",
  ylim = c(0, 1), lwd = 2
)
lines(ci[1, ] ~ seq(min(train2$ferritin), max(train2$ferritin), length.out = 50), lty = 2)
lines(ci[2, ] ~ seq(min(train2$ferritin), max(train2$ferritin), length.out = 50), lty = 2)
rug(train2$ferritin)



## two variables at once
nd <- expand.grid(
  total_psa_access = seq(min(train2$total_psa_access),
    max(train2$total_psa_access),
    length.out = 20
  ),
  percent_free_psa_access = seq(min(train2$percent_free_psa_access),
    max(train2$percent_free_psa_access),
    length.out = 20
  )
)
nd$trail <- 0
nd$ferritin <- 0

B1_pred <- predict(B1, nd$total_psa_access)
B2_pred <- predict(B2, nd$percent_free_psa_access)
B3_pred <- predict(B3, nd$trail)
B4_pred <- predict(B4, nd$ferritin)

nd$p_cancer <- predict_prostate(dr, list(B1_pred, B2_pred, B3_pred, B4_pred)) %>%
  rowMeans()


# x11()
wireframe(p_cancer ~ total_psa_access * percent_free_psa_access,
  data = nd, drape = TRUE, main = "TRAIL = 0\nFerritin = -1",
  colorkey = FALSE, zlim = c(0, 1), zlab = "P(Cancer)",
  ylab = "% Free PSA", xlab = "Total PSA"
) # , screen = list(z = 70, x = -70))




nd <- expand.grid(
  total_psa_access = seq(min(train2$total_psa_access),
    max(train2$total_psa_access),
    length.out = 20
  ),
  percent_free_psa_access = seq(min(train2$percent_free_psa_access),
    max(train2$percent_free_psa_access),
    length.out = 20
  )
)
nd$trail <- -2
nd$ferritin <- 2

B1_pred <- predict(B1, nd$total_psa_access)
B2_pred <- predict(B2, nd$percent_free_psa_access)
B3_pred <- predict(B3, nd$trail)
B4_pred <- predict(B4, nd$ferritin)

nd$p_cancer <- predict_prostate(dr, list(B1_pred, B2_pred, B3_pred, B4_pred)) %>%
  rowMeans()


# x11()
wireframe(p_cancer ~ total_psa_access * percent_free_psa_access,
  data = nd, drape = TRUE, main = "TRAIL = -2\nFerritin = 2",
  colorkey = FALSE, zlim = c(0, 1), zlab = "P(Cancer)",
  ylab = "% Free PSA", xlab = "Total PSA"
) # , screen = list(z = 70, x = -70))





## ============================================================
## Diverse samples
## ============================================================


B1_pred <- predict(B1, train1$test$total_psa_access)
B2_pred <- predict(B2, train1$test$percent_free_psa_access)
B3_pred <- predict(B3, train1$test$trail)
B4_pred <- predict(B4, train1$test$ferritin)

train1$test$p_cancer <- predict_prostate(dr, list(B1_pred, B2_pred, B3_pred, B4_pred)) %>%
  rowMeans()


xyplot(train1$test$p_cancer ~ factor(test$cancer),
  # groups=train1$train$targeted_disease_for_collection ==  "Benign Prostatic Hyperplasia",
  # groups=test$vendor,
  # groups = train1$test$country,
  # groups=test$age > 65,
  # groups = factor(paste0(train1$test$vendor, "_", train1$test$country)),
  ## data=.,
  type = c("g", "p", "a"), ylim = c(-0.03, 1.03),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "Weight = 10",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

calc_binomial_metrics(test$cancer, train1$test$p_cancer, thresh = 0.5, filter = NULL)
##      auc brier_mean brier_median logloss_mean logloss_median accuracy sensitivity specificity
## 1 0.9328      0.112        0.003         0.36         0.0562   0.8398      0.8889      0.7907
##      PPV    NPV    PLR    NLR     F1    MCC     P4 N_healthy N_cancer
## 1 0.4706 0.9714 4.2469 0.1405 0.6154 0.5481 0.7215        43        9

d[match(as.character(train1$test$id), as.character(d$serum_specific_id)), ] %>% View()

## https://stats.stackexchange.com/questions/519465/how-to-correctly-use-i-splines-for-monotone-non-decreasing-increasing-regressio
library(randomForest)
not_normalised <- c("miap", "tnf_alpha")
remove <- c("free_psa_access", "bound_psa", not_normalised)

rf <- randomForest(
  x = select(train2, !all_of(remove)), y = factor(train1$train$cancer),
  classwt = class_wt,
  importance = FALSE, seed = 123
)

plot(rf)
plot(importance(rf))
varImpPlot(rf)



xyplot(predict(rf, type = "prob")[, 2] ~ factor(train1$train$cancer),
  # groups = train2$targeted_disease_for_collection == "Benign Prostatic Hyperplasia",
  type = c("g", "p", "a"),
  xlab = "", ylab = "", main = "",
  ylim = c(-0.03, 1.03),
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

xyplot(train2$nse_kryptor_corrected ~ factor(train1$train$cancer),
  # groups = train2$targeted_disease_for_collection == "Benign Prostatic Hyperplasia",
  type = c("g", "p", "a"),
  xlab = "", ylab = "", main = "",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

xyplot(train2$midkine ~ factor(train1$train$cancer),
  # groups = train2$targeted_disease_for_collection == "Benign Prostatic Hyperplasia",
  type = c("g", "p", "a"),
  xlab = "", ylab = "", main = "",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)


## check the low values in the cancer group... are they imputed?
xyplot(train2$bound_psa ~ factor(dat$y_train),
  type = c("g", "p", "a"),
  xlab = "", ylab = "", main = "",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

ids <- train1$train$id[train1$train$cancer == 1 & train2$bound_psa < 0]

filter(prostate, id %in% ids) %>% View()

filter(d2, id %in% ids) %>%
  select(
    id, vendor, country, age, cancer, total_psa_access, percent_free_psa_access, p2psa_access,
    ferritin, trail
  ) %>%
  View()


plot(log10(total_psa + 1) ~ log10(total_psa_access + 1),
  xlab = "Log10 Total PSA (access)", ylab = "Log10 Total PSA",
  data = d, pch = 20, main = "All male samples",
  col = ifelse(as.character(d$serum_specific_id) %in% as.character(ids), "red", "darkgrey"),
  cex = ifelse(as.character(d$serum_specific_id) %in% as.character(ids), 3, 1)
)


plot(log10(free_psa_millipore + 1) ~ log10(free_psa_access + 1),
  xlab = "Log10 Free PSA (access)", ylab = "Log10 Free PSA",
  data = d, pch = 20, main = "All male samples",
  col = ifelse(as.character(d$serum_specific_id) %in% as.character(ids), "red", "darkgrey"),
  cex = ifelse(as.character(d$serum_specific_id) %in% as.character(ids), 3, 1)
)


pca <- prcomp(train2, scale = TRUE)

plot(pca$x[, 1], pca$x[, 2],
  pch = 20, col = ifelse(as.character(train1$train$id) %in% as.character(ids), "red", "blue"),
  cex = ifelse(as.character(train1$train$id) %in% as.character(ids), 3, 1)
)

plot(pca$x[, 1], pca$x[, 2], pch = 20, cex = 1.5, col = train1$train$cancer + 1)


tmp <- select(train2, !all_of(remove))
rf <- interactionfor(cancer ~ .,
  data = data.frame(tmp,
    cancer = factor(train1$train$cancer)
  ),
  class.weights = class_wt, seed = 123
)

summary(rf)
plot(rf)

plot(train2$ca19_9 ~ train2$ca19_9_access, pch = 20, col = train1$train$cancer + 1)
plot(train1$train$ca19_9 ~ train1$train$ca19_9_access, pch = 20, col = train1$train$cancer + 1)
plot(prostate$ca19_9 ~ prostate$ca19_9_access, pch = 20, col = train1$train$cancer + 1)


train1$train %>%
  mutate(vend_co = factor(paste0(vendor, "_", country))) %>%
  xyplot(prostate$ca19_9_access ~ vend_co,
    groups = cancer,
    data = .,
    type = c("g", "p", "a"),
    xlab = "", ylab = "", main = "",
    jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0), x = list(rot = 45)),
    auto.key = list(columns = 2, title = " Prostate Cancer")
  )

plot(train2$free_psa_access ~ train1$train$free_psa_access, pch = 20, col = train1$train$cancer + 1)
abline(0, 1)

names(train1$train)[-c(1:9)] == names(train2)

nms <- names(train1$train)[-c(1:9)]
## l1cam and miap are not normalised
for (i in nms) {
  print(i)
  print(all.equal(train1$train[, i], train2[, i]))
}

X <- select(
  train2, total_psa_access, percent_free_psa_access,
  trail, ferritin
)


rf <- randomForest(
  x = X, y = factor(train1$train$cancer),
  classwt = class_wt, ntree = 1000,
  importance = FALSE, seed = 123
)

plot(rf)
varImpPlot(rf)

xyplot(predict(rf, type = "prob")[, 2] ~ factor(train1$train$cancer),
  # groups = train2$targeted_disease_for_collection == "Benign Prostatic Hyperplasia",
  type = c("g", "p", "a"),
  xlab = "", ylab = "", main = "",
  ylim = c(-0.03, 1.03),
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

X2 <- X
X2$total_psa_access <- seq(min(train2$total_psa_access),
  max(train2$total_psa_access),
  length.out = nrow(X)
)
X2[, 2:4] <- 0
X2[, 2] <- -1

newpreds <- predict(rf, newdata = X2, type = "prob")

par(las = 1)
plot(newpreds[, 2] ~ X2$total_psa_access,
  type = "l",
  ylim = c(0, 1), lwd = 2, xlab = "Total PSA",
  ylab = "Probability of Prostate Cancer"
)
rug(X$total_psa_access)


X2 <- X
X2$percent_free_psa_access <- seq(min(train2$percent_free_psa_access),
  max(train2$percent_free_psa_access),
  length.out = nrow(X)
)
X2[, c(1, 3, 4)] <- 0


newpreds <- predict(rf, newdata = X2, type = "prob")

plot(newpreds[, 2] ~ X2$percent_free_psa_access,
  type = "l",
  ylim = c(0, 1)
)


X2 <- X
X2$trail <- seq(min(train2$trail),
  max(train2$trail),
  length.out = nrow(X)
)
X2[, c(1, 2, 4)] <- 0


newpreds <- predict(rf, newdata = X2, type = "prob")

plot(newpreds[, 2] ~ X2$trail,
  type = "l",
  ylim = c(0, 1)
)



X2 <- X
X2$ferritin <- seq(min(train2$ferritin),
  max(train2$ferritin),
  length.out = nrow(X)
)
X2[, c(1, 2, 3)] <- 0


newpreds <- predict(rf, newdata = X2, type = "prob")

plot(newpreds[, 2] ~ X2$ferritin,
  type = "l",
  ylim = c(0, 1)
)


library(pre)

ff <- dplyr::select(train2, -c(
  tgm2, leptin, fapa, nse_corrected,
  hepsin, shbg, il8
))

#ff <- train2
ff$y <- factor(train1$train$cancer)

rfit <- pre::pre(y ~ .,
  data = ff, family = "binomial", weights = obs_wt, #maxdepth = 2,
  normalize = FALSE, standardize = FALSE
)

summary(rfit)
importance(rfit)
