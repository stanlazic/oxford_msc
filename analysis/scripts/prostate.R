library(tidyverse)
library(cmdstanr)
library(lattice)
library(caret)
library(splines2)
library(CalibrationCurves)
library(predtools)
library(mice)

source("functions.R")
source("new_functions.R")


# d <- read_csv("../data_from_ICR/Cancer Model Training Data 20250228.csv")
d <- read_csv("../data_from_ICR/Cancer Model Training Data 20250828.csv")

d2 <- select_samples(d) %>%
  clean_raw_data() %>%
  tidy_prostate_data(., remove_many_missing_biomarkers = TRUE) %>%
  filter(age > 52 & age < 83) %>%
  droplevels()


## standardize and normalize
train1 <- preproc(d2, d2, index = 9)

## imput missing values
train1 <- impute_missing(train1$train, train1$test, skip_cols = c(1:4, 6:7, 9), bm_start = 10)

## remove bath effests
train2 <- correct_batch_effects(
  train1$train,
  train1$test,
  index = 9
)$train_resid


## ============================================================
## Bayesian model
## ============================================================


compiled_model <- cmdstan_model("stan_models/prostate_final_model.stan",
  dir = "compiled_models"
)

B1 <- iSpline(train2$bound_psa_access, df = 6)
B2 <- iSpline(train2$free_psa_access, df = 6)
B3 <- iSpline(train2$p2psa_access, df = 6)


dat <- list(
  P = ncol(B1),
  N_train = nrow(train2),
  X1_train = B1,
  X2_train = B2,
  X3_train = B3,
  X4_train = ifelse(train2$hepsin > 0.75, 1, 0),
  y_train = train1$train$cancer
)

m1 <- compiled_model$sample(
  data = dat,
  chains = 5,
  parallel_chains = 5,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  seed = 123,
  refresh = 0
)


dr <- m1$draws(format = "df") %>%
  data.frame()


etas <- dr %>%
  select_at(vars(starts_with("linpred"))) %>%
  invlogit() %>%
  colMeans()


xyplot(etas ~ factor(dat$y_train),
  groups = paste0(train1$train$vendor, "_", train1$train$country),
  type = c("g", "p", "a"),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "Bayes GAM model",
  ylim = c(-0.03, 1.03),
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

xyplot(etas ~ train1$train$targeted_disease_for_collection,
  type = c("g", "p", "a"),
  xlab = "", ylab = "Predicted probability of cancer", main = "Bayes GAM model",
  ylim = c(-0.03, 1.03),
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0), x = list(rot = 45)),
  auto.key = list(columns = 2)
)


calc_binomial_metrics(train1$train$cancer, etas, thresh = 0.5) %>% t()

val.prob.ci.2(etas, dat$y_train)



## ============================================================
## Cross validation
## ============================================================

set.seed(123)
f <- createFolds(paste0(d2$cancer, d2$vendor, d2$country), k = 10, list = TRUE, returnTrain = FALSE)

bayes_gam_post <- matrix(NA, 5000, ncol = nrow(d2))

tmp <- d2
hepsin <- ifelse(tmp$hepsin > 1.65, 1, 0)
hepsin[is.na(hepsin)] <- 0

for (i in 1:length(f)) {
  ## standardize and normalize
  tr1 <- preproc(tmp[-f[[i]], ], tmp[f[[i]], ], index = 9)

  ## imput missing values
  tr1 <- impute_missing(tr1$train, tr1$test, skip_cols = c(1:4, 6:7, 9), bm_start = 10)

  ## remove bath effests
  tr2 <- correct_batch_effects(
    tr1$train,
    tr1$test,
    index = 9
  )

  B1 <- iSpline(tr2$train_resid$bound_psa_access, df = 6)
  B2 <- iSpline(tr2$train_resid$free_psa_access, df = 6)
  B3 <- iSpline(tr2$train_resid$p2psa_access, df = 6)


  dat <- list(
    P = ncol(B1),
    N_train = nrow(tr2$train_resid),
    X1_train = B1,
    X2_train = B2,
    X3_train = B3,
    X4_train = hepsin[-f[[i]]], 
    y_train = tr1$train$cancer
  )

  m1 <- compiled_model$sample(
    data = dat,
    chains = 5,
    parallel_chains = 5,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    seed = 123,
    refresh = 0
  )

  dr <- m1$draws(format = "df") %>%
    data.frame()

  bayes_gam_post[, f[[i]]] <- predict_prostate(
    dr, list(
      predict(B1, newx = tr2$test_resid$bound_psa_access),
      predict(B2, newx = tr2$test_resid$free_psa_access),
      predict(B3, newx = tr2$test_resid$p2psa_access),
      hepsin[f[[i]]]
    )
  ) %>% t()
}


preds <- bayes_gam_post %>%
  data.frame() %>%
  dplyr::reframe(
    mean = apply(., 2, function(x) mean(inv)),
    sd = apply(., 2, sd),
    median = apply(., 2, median),
    p_above = apply(., 2, function(x) mean(x > 0.5))
  )

calc_binomial_metrics(tmp$cancer, preds$mean) %>% t()

calc_binomial_metrics(tmp$cancer, preds$mean, thresh=prop.table(table(tmp$cancer))[2]) %>% t()

val.prob.ci.2(preds$mean, tmp$cancer)

## BSS
1 - (weighted_brier(tmp$cancer, preds$mean) /
  weighted_brier(tmp$cancer, rep(mean(tmp$cancer), nrow(tmp))))


xyplot(preds$mean ~ factor(tmp$cancer),
  #groups = factor(paste0(tmp$vendor, "_", tmp$country)),
  type = c("g", "p", "a"),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)


xyplot(preds$mean ~ factor(tmp$targeted_disease_for_collection),
  type = c("g", "p", "a"),
  xlab = "", ylab = "Predicted probability of cancer", main = "",
  jitter.x = TRUE, scales = list(alternating = FALSE, tck = c(1, 0), x=list(rot=45)),
  auto.key = list(columns = 2)
)

