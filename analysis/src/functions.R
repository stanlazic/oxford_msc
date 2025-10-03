invlogit <- function(x) {
  1 / (1 + exp(-x))
}



select_samples <- function(d) {
  d <- d %>%
    mutate_if(is.character, as.factor) %>%
    filter(Race %in% c("White", "Hispanic")) %>%
    filter(Stage %in% c("1", "2", "N/A")) %>%
    filter(Vendor != "ICR") %>%
    filter(`Collection Country` != "Unknown") %>%
    droplevels()

  return(d)
}


clean_raw_data <- function(d) {
  ## d: data as provided by ICR

  require(tidyverse)
  require(lubridate)
  ## clean names
  d <- janitor::clean_names(d)

  ## fix disease names with lower case "C"
  d$disease <- recode_factor(d$disease,
    "Esophageal cancer" = "Esophageal Cancer",
    "Ovarian cancer" = "Ovarian Cancer"
  )

  ## filter poor quality samples
  d <- d %>%
    filter(training_or_blinded != "Z_Do_Not_Analyze_Hemolyzed") %>%
    filter(training_or_blinded != "Z_Do_Not_Analyze_Conutry_Unknown") %>%
    filter(training_or_blinded != "Z_Do_Not_Analyze_Stage_Unknown")

  ## fix disease cateogry (missing values)
  for (i in 1:nrow(d)) {
    if (is.na(d$disease[i]) & d$targeted_disease_for_collection[i] %in%
      c("Benign Prostatic Hyperplasia", "Chronic Prostatitis", "Microscopic Colitis", "Healthy")) {
      d$disease[i] <- "Healthy"
    }

    if (is.na(d$disease[i]) & !d$targeted_disease_for_collection[i] %in%
      c("Benign Prostatic Hyperplasia", "Chronic Prostatitis", "Microscopic Colitis", "Healthy")) {
      d$disease[i] <- d$targeted_disease_for_collection[i]
    }
  }

  d$cancer <- ifelse(d$disease == "Healthy", 0, 1)


  ## select key variables
  d <- d %>%
    select(
      vendor, serum_specific_id, collection_date, collection_country,
      age_at_collection, sex, cancer, disease, targeted_disease_for_collection, afp:tgm2, cg_a,
      scc, thrombospondin_2, free_psa_access:ferritin_access, nse_kryptor_corrected, aldh1a1_corrected,
      nse_corrected, percent_free_psa_access, thymidine_kinase_1, billirubin_index_serum, lrg1
    ) %>%
    select(-c(
      hgf, beta_hcg, he4, s_fasl, total_psa, aldh1a1_not_corrected,
      nse_not_corrected
    ))
  ## rename first few columns
  names(d)[1:5] <- c("vendor", "id", "collection_date", "country", "age")

  ## rename bilirubin
  d <- d %>%
    dplyr::rename(bilirubin = billirubin_index_serum)

  ## convert to numeric
  d <- d %>%
    mutate(
      nse_kryptor_corrected = as.numeric(as.character(nse_kryptor_corrected)),
      aldh1a1_corrected = as.numeric(as.character(aldh1a1_corrected)),
      nse_corrected = as.numeric(as.character(nse_corrected))
    )


  ## reorder and clean country names
  d <- d %>%
    relocate(cancer, .before = afp) %>%
    relocate(disease, .before = afp) %>%
    relocate(id, .before = vendor) %>%
    mutate(country = recode_factor(country,
      "Russia" = "Ukraine",
      "Georgia" = "Ukraine",
      "Romania" = "Ukraine"
    ))


  ## combine all non BioIVT vendors
  d <- d %>%
    mutate(vendor = recode_factor(vendor,
      "BioOptions" = "Other",
      "Discovery Life Sciences" = "Other",
      "Precision for Medicine" = "Other",
      "PrecisionMed" = "Other",
      "Conversant" = "Other",
      "Reprocell" = "Other"
    )) %>%
    droplevels()


  ## make into date type
  d <- d %>%
    mutate(collection_date = mdy(collection_date)) %>%
    mutate(age = as.numeric(as.character(age)))

  ## remove people with history of other cancers
  d <- d[!d$id %in% c("201652718S", "201652730S", "201707190S"), ]

  ## remove people with prostatectomy
  d <- d[!d$id %in% c("HMN54232", "HMN54233", "HMN54248", "HMN54227"), ]

  return(d)
}


tidy_prostate_data <- function(d, remove_many_missing_biomarkers = TRUE) {
  ## d: cleaned data

  d <- d %>%
    filter(disease %in% c("Healthy", "Prostate Cancer")) %>%
    filter(sex == "M") %>%
    mutate(percent_free_psa_access = as.numeric(as.character(percent_free_psa_access))) %>%
    ## remove bad samples
    filter(!id %in% c("HMN686019S", "HMN456748")) %>%
    mutate_if(is.character, as.factor) %>%
    droplevels() %>%
    as.data.frame()

  ## remove samples with no PSA data
  no_psa_index <- !apply(d[, grepl("psa", names(d))], 1, function(x) all(is.na(x)))
  d <- d[no_psa_index, ]

  if (remove_many_missing_biomarkers) {
    ## remove samples with too many missing values before imputing
    few_missing <- apply(d, 2, function(x) sum(is.na(x))) / nrow(d)
    few_missing <- few_missing[few_missing < 0.2]

    d <- select(d, names(few_missing))
  }

  ## calculate % p2psa
  d <- d %>%
    mutate(percent_p2psa_access = ((p2psa_access / (total_psa_access * 1000)) * 100)) %>%
    mutate(bound_psa_access = (total_psa_access - free_psa_access))

  ## remove some biomarkers with bad distributions or not relevant
  d <- select(d, -c(bilirubin, tgm2, nse_corrected))

  return(d)
}


tidy_breast_data <- function(x, remove_many_missing_biomarkers = TRUE) {
  x <- x %>%
    filter(disease %in% c("Healthy", "Breast Cancer")) %>%
    filter(sex == "F") %>%
    select(-contains("psa")) %>%
    mutate_if(is.character, as.factor) %>%
    droplevels() %>%
    as.data.frame()

  if (remove_many_missing_biomarkers) {
    ## remove samples with too many missing values before imputing
    few_missing <- apply(x, 2, function(y) sum(is.na(y))) / nrow(x)
    few_missing <- few_missing[few_missing < 0.2]

    x <- select(x, names(few_missing))
  }

  ## remove some biomarkers with bad distributions or not relevant
  x <- select(x, -c(tgm2, bilirubin))

  return(x)
}


tidy_generic_test_data <- function(x, filter_age = 0) {
  ## x: cleaned data

  x <- x %>%
    filter(age >= filter_age) %>%
    droplevels()
  return(x)
}



tidy_pancreas_data <- function(x, remove_many_missing_biomarkers = TRUE) {
  x <- x %>%
    filter(disease %in% c("Healthy", "Pancreatic Cancer")) %>%
    select(-contains("psa")) %>%
    mutate_if(is.character, as.factor) %>%
    droplevels() %>%
    as.data.frame()

  if (remove_many_missing_biomarkers) {
    ## remove samples with too many missing values before imputing
    few_missing <- apply(x, 2, function(y) sum(is.na(y))) / nrow(x)
    few_missing <- few_missing[few_missing < 0.2]

    x <- select(x, names(few_missing))
  }

  ## remove some biomarkers with bad distributions or not relevant
  x <- select(x, -c(tgm2, bilirubin))

  return(x)
}


tidy_colorectal_data <- function(x, remove_many_missing_biomarkers = TRUE) {
  x <- x %>%
    filter(disease %in% c("Healthy", "Colorectal Cancer")) %>%
    select(-contains("psa")) %>%
    mutate_if(is.character, as.factor) %>%
    droplevels() %>%
    as.data.frame()

  if (remove_many_missing_biomarkers) {
    ## remove samples with too many missing values before imputing
    few_missing <- apply(x, 2, function(y) sum(is.na(y))) / nrow(x)
    few_missing <- few_missing[few_missing < 0.2]

    x <- select(x, names(few_missing))
  }

  ## remove some biomarkers with bad distributions or not relevant
  x <- select(x, -c(tgm2, bilirubin))

  return(x)
}



tidy_bladder_data <- function(x, remove_many_missing_biomarkers = TRUE) {
  x <- x %>%
    filter(disease %in% c("Healthy", "Bladder Cancer")) %>%
    select(-contains("psa")) %>%
    mutate_if(is.character, as.factor) %>%
    droplevels() %>%
    as.data.frame()

  if (remove_many_missing_biomarkers) {
    ## remove samples with too many missing values before imputing
    few_missing <- apply(x, 2, function(y) sum(is.na(y))) / nrow(x)
    few_missing <- few_missing[few_missing < 0.2]

    x <- select(x, names(few_missing))
  }

  ## remove some biomarkers with bad distributions or not relevant
  x <- select(x, -c(tgm2, bilirubin))


  return(x)
}


tidy_thyroid_data <- function(x, remove_many_missing_biomarkers = TRUE) {
  x <- x %>%
    filter(disease %in% c("Healthy", "Thyroid Cancer")) %>%
    select(-contains("psa")) %>%
    mutate_if(is.character, as.factor) %>%
    droplevels() %>%
    as.data.frame()

  if (remove_many_missing_biomarkers) {
    ## remove samples with too many missing values before imputing
    few_missing <- apply(x, 2, function(y) sum(is.na(y))) / nrow(x)
    few_missing <- few_missing[few_missing < 0.2]

    x <- select(x, names(few_missing))
  }

  ## remove some biomarkers with bad distributions or not relevant
  ## select(-c(vegf, igfbp3, cg_a, tgm2, scc)) %>%
  x <- select(x, -c(tgm2, bilirubin))

  return(x)
}


tidy_liver_data <- function(x, filter_age = 0) {
  ## x: cleaned data

  x <- x %>%
    filter(disease %in% c("Healthy", "Liver Cancer")) %>%
    select(-contains("psa")) %>%
    filter(age > filter_age) %>%
    droplevels()

  return(x)
}


tidy_kidney_data <- function(x, filter_age = 0) {
  ## x: cleaned data

  x <- x %>%
    filter(disease %in% c("Healthy", "Kidney Cancer")) %>%
    select(-contains("psa")) %>%
    filter(age > filter_age) %>%
    droplevels()

  return(x)
}




preproc <- function(train, test, index) {
  require(bestNormalize)
  set.seed(123)

  train2 <- matrix(NA, nrow = nrow(train), ncol = ncol(train) - index)
  colnames(train2) <- names(train)[(index + 1):ncol(train)]

  test2 <- matrix(NA, nrow = nrow(test), ncol = ncol(test) - index)
  colnames(test2) <- names(test)[(index + 1):ncol(test)]

  norm_type <- rep(NA, ncol(train) - index)
  names(norm_type) <- names(train)[(index + 1):ncol(train)]


  for (i in 1:(ncol(train) - index)) {
    ## apply transformation
    tr <- bestNormalize(train[train$cancer == 0, i + index],
      allow_lambert_s = FALSE,
      allow_lambert_h = FALSE,
      allow_exp = FALSE,
      quiet = TRUE,
      warn = TRUE
    )


    train2[train$cancer == 0, i] <- tr$x.t
    train2[train$cancer == 1, i] <- predict(tr, train[train$cancer == 1, i + index])

    test2[, i] <- predict(tr, test[, i + index])

    norm_type[i] <- attributes(tr$chosen_transform)$class
  }

  ## preProcValues <- preProcess(train2, method = "knnImpute")
  ## train2 <- predict(preProcValues, train2)
  ## test2 <- predict(preProcValues, test2)

  train2 <- data.frame(train[, 1:index], train2)
  test2 <- data.frame(test[, 1:index], test2)

  return(list(train = train2, test = test2, norm_type = norm_type))
}



impute_missing <- function(train, test, skip_cols = c(1:4, 6:7, 9), bm_start) {
  d_imps <- mice(
    rbind(train[, -skip_cols], test[, -skip_cols]),
    m = 50,
    maxit = 10,
    method = "pmm",
    ignore = rep(c(FALSE, TRUE), times = c(nrow(train), nrow(test))),
    seed = 123
  )

  ## combine imputed datasets and remove non-biomarker columns
  extra_cols <- max(skip_cols) - length(skip_cols)
  completed_datasets <- lapply(1:50, function(i) complete(d_imps, i)[, -c(1:extra_cols)])


  # Convert list of matrices to a 3D array
  mat_array <- array(unlist(completed_datasets),
    dim = c(
      dim(completed_datasets[[1]]),
      length(completed_datasets)
    )
  )

  # Apply mean across the third dimension (i.e., for each cell)
  cell_means <- apply(mat_array, c(1, 2), mean)
  colnames(cell_means) <- names(train)[bm_start:ncol(train)]

  return(
    list(
      train = data.frame(train[, 1:(bm_start - 1)], cell_means[1:nrow(train), ]),
      test = data.frame(test[, 1:(bm_start - 1)], cell_means[(nrow(train) + 1):nrow(cell_means), ])
    )
  )
}


correct_batch_effects <- function(train, test, index = 9) {
  require(nlme)

  ## remove non-biomarker data
  train2 <- train[, -c(1:index)]
  test2 <- test[, -c(1:index)]

  train$vend_co <- factor(paste(train$vendor, train$country, sep = "_"))
  test$vend_co <- factor(paste(test$vendor, test$country, sep = "_"))

  ## matrix to store results of corrected values
  train_resid <- matrix(NA, nrow = nrow(train2), ncol = ncol(train2))
  test_resid <- matrix(NA, nrow = nrow(test2), ncol = ncol(test2))

  colnames(train_resid) <- names(train2)
  colnames(test_resid) <- names(test2)

  ## loop through each biomarker
  for (i in 1:ncol(train2)) {
    tmp <- data.frame(train, y = train2[, i]) %>%
      filter(cancer == 0)

    m_adj <- lme(y ~ age,
      random = ~ 1 | vend_co,
      data = tmp
    )

    ## calculate residuals
    train_resid[, i] <- winsorize(train2[, i] - predict(m_adj, train, type = "response"))
    test_resid[, i] <- winsorize(test2[, i] - predict(m_adj, test, type = "response"))
  }

  return(list(
    train_resid = data.frame(train_resid),
    test_resid = data.frame(test_resid)
  ))
}

winsorize <- function(x, lower = -4, upper = 4) {
  x[x < lower] <- lower
  x[x > upper] <- upper
  return(x)
}


weighted_brier <- function(actual, predicted) {
  # actual: vector of 0/1 cancer status values (1 = cancer, 0 = healthy)
  # predicted: vector of predicted probabilities

  # Brier score for healthy samples
  brier_0 <- mean((predicted[actual == 0] - actual[actual == 0])^2, na.rm = TRUE)

  # Brier score for cancer samples
  brier_1 <- mean((predicted[actual == 1] - actual[actual == 1])^2, na.rm = TRUE)

  # average Brier score
  (brier_0 + brier_1) / 2
}


weighted_brier_median <- function(actual, predicted) {
  brier_0 <- median((predicted[actual == 0] - actual[actual == 0])^2, na.rm = TRUE)
  brier_1 <- median((predicted[actual == 1] - actual[actual == 1])^2, na.rm = TRUE)

  (brier_0 + brier_1) / 2
}


log_loss <- function(pred, actual) {
  ## calculate logloss
  pred <- pmax(pred, 0.000001)
  pred <- pmin(pred, 0.999999)
  -1 * ((actual) * log(pred) + (1 - actual) * log(1 - pred))
}


weighted_logloss <- function(actual, predicted) {
  ll_0 <- mean(log_loss(predicted[actual == 0], actual[actual == 0]), na.rm = TRUE)
  ll_1 <- mean(log_loss(predicted[actual == 1], actual[actual == 1]), na.rm = TRUE)

  (ll_0 + ll_1) / 2
}

weighted_logloss_median <- function(actual, predicted) {
  ll_0 <- median(log_loss(predicted[actual == 0], actual[actual == 0]), na.rm = TRUE)
  ll_1 <- median(log_loss(predicted[actual == 1], actual[actual == 1]), na.rm = TRUE)

  (ll_0 + ll_1) / 2
}


balanced_acc <- function(actual, predicted) {
  not_missing <- !is.na(predicted)

  acc_0 <- mean(predicted[actual == 0 & not_missing] == actual[actual == 0 & not_missing])
  acc_1 <- mean(predicted[actual == 1 & not_missing] == actual[actual == 1 & not_missing])

  (acc_0 + acc_1) / 2
}


brier_test <- function(pred, guess, actual) {
  ## pred: model prediction
  ## guess: random guess or another prediction
  ## actual: true class lable (0/1)

  ## calculate Brier score for prediction and random guess
  ## for healthy subjects
  brier_pred_0 <- (pred[actual == 0] - actual[actual == 0])^2
  ## for cancer subjects
  brier_pred_1 <- (pred[actual == 1] - actual[actual == 1])^2

  ## for healthy subjects
  brier_guess_0 <- (guess[actual == 0] - actual[actual == 0])^2
  ## for cancer subjects
  brier_guess_1 <- (guess[actual == 1] - actual[actual == 1])^2

  ## compare prediction and guessing
  wt_0 <- wilcox.test(brier_pred_0, brier_guess_0, paired = TRUE, alternative = "less")
  wt_1 <- wilcox.test(brier_pred_1, brier_guess_1, paired = TRUE, alternative = "less")

  ## calculate geometric mean
  mean_p <- sqrt(wt_0$p.value * wt_1$p.value)

  return(p_value = mean_p)
}


ll_test <- function(pred, guess, actual) {
  ## pred: model prediction
  ## guess: random guess or another prediction
  ## actual: true class lable (0/1)

  ## calculate log loss
  ll_pred_0 <- log_loss(actual[actual == 0], pred[actual == 0])
  ll_pred_1 <- log_loss(actual[actual == 1], pred[actual == 1])

  ll_guess_0 <- log_loss(actual[actual == 0], guess[actual == 0])
  ll_guess_1 <- log_loss(actual[actual == 1], guess[actual == 1])

  ## compare prediction and guessing
  wt_0 <- wilcox.test(ll_pred_0, ll_guess_0, paired = TRUE, alternative = "less")
  wt_1 <- wilcox.test(ll_pred_1, ll_guess_1, paired = TRUE, alternative = "less")

  ## calculate geometric mean
  mean_p <- sqrt(wt_0$p.value * wt_1$p.value)

  return(p_value = mean_p)
}

acc_test <- function(pred, guess, actual, thresh = 0.5) {
  ## pred: model prediction
  ## guess: random guess or another prediction
  ## actual: true class lable (0/1)

  ## sample size
  N <- length(actual)

  ## calculate balanced accuracy (proportion correct) for
  ## predictions and guessing
  prop_pred <- balanced_acc(actual, as.numeric(pred > thresh))
  prop_guess <- balanced_acc(actual, as.numeric(guess > thresh))

  ## proportion test
  prop.test(c(prop_pred * N, prop_guess * N),
    c(N, N),
    alternative = "greater"
  )$p.value
}




binomial_gam_model_prostate <- function(train, test, batch_correct = TRUE) {
  ## calculate amount of bound PSA
  train$bound_psa <- train$total_psa_access - train$free_psa_access
  test$bound_psa <- test$total_psa_access - test$free_psa_access

  ## set negative values to zero
  train$bound_psa[train$bound_psa < 0] <- 0
  test$bound_psa[test$bound_psa < 0] <- 0


  ## center, scale, normalise, impute data
  preproc_data <- preproc(train, test)

  ## select matched samples
  m.out <- preproc_data$train %>%
    matchit(cancer ~ age + vendor + country, data = ., method = "nearest", ratio = 2)

  preproc_data$train <- match_data(m.out)

  ## Cholesky decomposition
  correlated_variables <- preproc_data$train[, c(
    "bound_psa",
    "free_psa_access",
    "p2psa_access"
  )] %>%
    as.matrix()

  cor_matrix <- cor(correlated_variables)

  ## Perform Cholesky decomposition
  chol_matrix <- chol(cor_matrix)

  ## Decorrelate variables using Cholesky decomposition
  decorrelated_variables <- correlated_variables %*% solve(chol_matrix)
  decorrelated_variables_test <- as.matrix(preproc_data$test[, c(
    "bound_psa",
    "free_psa_access",
    "p2psa_access"
  )]) %*% solve(chol_matrix)



  train2 <- preproc_data$train
  test2 <- preproc_data$test

  train2[, names(decorrelated_variables)] <- decorrelated_variables
  test2[, names(decorrelated_variables_test)] <- decorrelated_variables_test

  ## calculate class weights
  class_wt <- 1 / prop.table(table(train2$cancer))
  obs_wt <- ifelse(train2$cancer == 0, class_wt[1], class_wt[2])

  set.seed(123)

  s1 <- scam(
    cancer ~ s(bound_psa, bs = "mpi") +
      s(free_psa_access, bs = "mpd") +
      s(prolactin, bs = "mpi") +
      s(ca19_9, bs = "mpi") +
      s(ferritin, bs = "mpi"),
    weights = round(obs_wt, 0),
    data = train2,
    family = "binomial",
    gamma = 1.4
  )


  return(list(
    train_pred = s1$fitted.values,
    y_train = train2$cancer,
    id = test2$id,
    test_pred = predict(s1, newdata = test2, type = "response"),
    y_test = test2$cancer,
    test_data = test2
  ))
}



binomial_gam_model_pancreas <- function(train, test, batch_correct = TRUE) {
  ## center, scale, normalise, impute data
  preproc_data <- preproc(train, test)

  ## remove batch effects
  if (batch_correct) {
    batch_corrected <- correct_batch_effects(
      preproc_data$train,
      preproc_data$test
    )
  } else {
    batch_corrected <- correct_batch_effects(
      preproc_data$train,
      preproc_data$train # stick in training
    )
    ## update with uncorrected-test data
    batch_corrected[[2]] <- preproc_data$test
  }

  train2 <- data.frame(train[, c(1:8)], batch_corrected$train_resid)
  test2 <- data.frame(test[, c(1:8)], batch_corrected$test_resid)

  ## calculate class weights
  class_wt <- 1 / prop.table(table(train2$cancer))
  obs_wt <- ifelse(train2$cancer == 0, class_wt[1], class_wt[2])

  ## raise lower values to close data gap
  train2$cea[train2$cea <= -0.6] <- -0.6
  train2$ca125[train2$ca125 <= 0.5] <- 0.5

  test2$cea[test2$cea <= -0.6] <- -0.6
  test2$ca125[test2$ca125 <= 0.5] <- 0.5

  set.seed(123)

  s1 <- scam(
    cancer ~ s(ca125, bs = "mpi") +
      s(ca19_9, bs = "mpi") +
      s(gdf_15, bs = "mpi") +
      trail,
    weights = round(obs_wt, 0),
    data = train2,
    family = "binomial",
    gamma = 1.4
  )


  return(list(
    train_pred = s1$fitted.values,
    y_train = train2$cancer,
    id = test2$id,
    test_pred = predict(s1, newdata = test2, type = "response"),
    y_test = test2$cancer,
    test_data = test2
  ))
}



binomial_gam_model_colorectal <- function(train, test, batch_correct = TRUE) {
  ## center, scale, normalise, impute data
  preproc_data <- preproc(train, test)

  ## remove batch effects
  if (batch_correct) {
    batch_corrected <- correct_batch_effects(
      preproc_data$train,
      preproc_data$test
    )
  } else {
    batch_corrected <- correct_batch_effects(
      preproc_data$train,
      preproc_data$train # stick in training
    )
    ## update with uncorrected-test data
    batch_corrected[[2]] <- preproc_data$test
  }

  train2 <- data.frame(train[, c(1:8)], batch_corrected$train_resid)
  test2 <- data.frame(test[, c(1:8)], batch_corrected$test_resid)

  ## calculate class weights
  class_wt <- 1 / prop.table(table(train2$cancer))
  obs_wt <- ifelse(train2$cancer == 0, class_wt[1], class_wt[2])

  set.seed(123)

  s1 <- scam(
    cancer ~ s(trail, bs = "mpd", k = 5) +
      s(gdf_15, bs = "mpi", k = 5) +
      s(prolactin, bs = "mpi", k = 5) +
      s(shbg, bs = "mpi", k = 5) +
      s(opn, bs = "mpi", k = 5),
    weights = round(obs_wt, 0),
    data = train2,
    family = "binomial",
    gamma = 1.4
  )

  return(list(
    train_pred = s1$fitted.values,
    y_train = train2$cancer,
    id = test2$id,
    test_pred = predict(s1, newdata = test2, type = "response"),
    y_test = test2$cancer,
    test_data = test2
  ))
}





calc_binomial_metrics <- function(y, y_pred, thresh = 0.5, filter = NULL) {
  if (is.null(filter)) {
    filter <- rep(TRUE, length(y_pred))
  }

  y <- y[filter]
  y_pred <- y_pred[filter]

  N_healthy <- sum(y == 0)
  N_cancer <- sum(y == 1)

  metrics <- data.frame(
    auc = ModelMetrics::auc(y, y_pred),
    brier_mean = ModelMetrics::brier(y, y_pred),
    logloss_mean = ModelMetrics::logLoss(y, y_pred),
    accuracy = 1 - ModelMetrics::ce(y, y_pred > thresh),
    bal_brier_mean = weighted_brier(y, y_pred),
    bal_brier_median = weighted_brier_median(y, y_pred),
    bal_logloss_mean = weighted_logloss(y, y_pred),
    bal_logloss_median = weighted_logloss_median(y, y_pred),
    bal_accuracy = balanced_acc(y, y_pred > thresh),
    sensitivity = ModelMetrics::sensitivity(y, y_pred, cutoff = thresh),
    specificity = ModelMetrics::specificity(y, y_pred, cutoff = thresh),
    PPV = ModelMetrics::ppv(y, y_pred, cutoff = thresh),
    NPV = ModelMetrics::npv(y, y_pred, cutoff = thresh),
    P4 = p4(y, y_pred, thresh = thresh),
    PLR = PLR(y, y_pred, thresh = thresh),
    NLR = NLR(y, y_pred, thresh = thresh),
    mean_0_pred = mean(y_pred[y == 0], na.rm = TRUE),
    mean_1_pred = mean(y_pred[y == 1], na.rm = TRUE),
    mean_diff = mean(y_pred[y == 1], na.rm = TRUE) - mean(y_pred[y == 0], na.rm = TRUE),
    median_0_pred = median(y_pred[y == 0], na.rm = TRUE),
    median_1_pred = median(y_pred[y == 1], na.rm = TRUE),
    median_diff = median(y_pred[y == 1], na.rm = TRUE) - median(y_pred[y == 0], na.rm = TRUE),
    N_healthy = N_healthy,
    N_cancer = N_cancer
  ) %>%
    mutate_if(is.numeric, ~ round(., 4))

  return(metrics)
}


make_empty_metrics <- function(y, y_pred, thresh = 0.5, filter = NULL) {
  if (is.null(filter)) {
    filter <- rep(TRUE, length(y_pred))
  }

  N_healthy <- sum(y[filter] == 0)
  N_cancer <- sum(y[filter] == 1)

  metrics <- data.frame(
    auc = NA,
    brier_mean = NA,
    brier_median = NA,
    logloss_mean = NA,
    logloss_median = NA,
    accuracy = NA,
    sensitivity = NA,
    specificity = NA,
    PPV = NA,
    NPV = NA,
    PLR = NA,
    NLR = NA,
    F1 = NA,
    MCC = NA,
    P4 = NA,
    N_healthy = N_healthy,
    N_cancer = N_cancer
  )

  return(metrics)
}

p4 <- function(actual, predicted, thresh) {
  tp <- sum(actual == 1 & (predicted > thresh))
  tn <- sum(actual == 0 & (predicted < thresh))
  fp <- sum(actual == 0 & (predicted > thresh))
  fn <- sum(actual == 1 & (predicted < thresh))

  (4 * tp * tn) / ((4 * tp * tn) + (tp + tn) * (fp + fn))
}


PLR <- function(y, y_pred, thresh = 0.5) {
  sensitivity <- ModelMetrics::sensitivity(y, y_pred, thresh)
  specificity <- ModelMetrics::specificity(y, y_pred, thresh)

  ## ensure we don't get division by zero
  specificity <- ifelse(specificity == 1, 0.999, specificity)

  return(sensitivity / (1 - specificity))
}

NLR <- function(y, y_pred, thresh = 0.5) {
  sensitivity <- ModelMetrics::sensitivity(y, y_pred, thresh)
  specificity <- ModelMetrics::specificity(y, y_pred, thresh)

  ## ensure we don't get division by zero
  specificity <- ifelse(specificity == 0, 0.001, specificity)

  return((1 - sensitivity) / specificity)
}



test_better_than_random <- function(y, y_pred, thresh) {
  N <- length(y)

  ## always guess P(Cancer) = 0.5
  guess_50_percent <- rep(0.5, N)

  ## always guess P(Cancer) = observed proportion of cancer cases in test set
  guess_observed_percent <- rep(mean(y), N)

  ## random uniform guess between 0.01 and 0.99
  set.seed(123)
  N_sim <- 100
  guess_random <- matrix(runif(N * N_sim, min = 0.01, max = 0.99),
    nrow = N, ncol = N_sim
  )

  ## Brier
  p_brier_50 <- brier_test(y_pred, guess_50_percent, y)
  p_brier_obs <- brier_test(y_pred, guess_observed_percent, y)
  p_brier_rand <- apply(guess_random, 2, function(x) brier_test(y_pred, x, y)) %>%
    mean()


  ## Log loss
  p_ll_50 <- ll_test(y_pred, guess_50_percent, y)
  p_ll_obs <- ll_test(y_pred, guess_observed_percent, y)
  p_ll_rand <- apply(guess_random, 2, function(x) ll_test(y_pred, x, y)) %>%
    mean()


  ## AUC
  p_auc_50 <- roc.test(
    predictor1 = y_pred,
    predictor2 = guess_50_percent,
    response = y,
    alternative = "greater",
    paired = TRUE,
    progress = "none"
  )$p.value

  p_auc_obs <- roc.test(
    predictor1 = y_pred,
    predictor2 = guess_observed_percent,
    response = y,
    alternative = "greater",
    paired = TRUE
  )$p.value

  p_auc_rand <- apply(guess_random, 2, function(x) {
    roc.test(
      predictor1 = y_pred,
      predictor2 = x,
      response = y,
      alternative = "greater",
      paired = TRUE,
      progress = "none"
    )$p.value
  }) %>%
    mean()



  ## Accuracy
  p_acc_50 <- acc_test(y_pred, guess_50_percent, y)
  p_acc_obs <- acc_test(y_pred, guess_observed_percent, y)
  p_acc_rand <- apply(guess_random, 2, function(x) acc_test(y_pred, x, y)) %>%
    mean()

  results <- data.frame(
    comparison = c(
      "Model vs. AUC 50%", "Model vs. AUC observed", "Model vs. AUC random",
      "Model vs. Brier 50%", "Model vs. Brier observed", "Model vs. Brier random",
      "Model vs. Logloss 50%", "Model vs. Logloss observed", "Model vs. Logloss random",
      "Model vs. Accuracy 50%", "Model vs. Accuracy observed", "Model vs. Accuracy random"
    ),
    p_value = c(
      p_auc_50, p_auc_obs, p_auc_rand,
      p_brier_50, p_brier_obs, p_brier_rand,
      p_ll_50, p_ll_obs, p_ll_rand,
      p_acc_50, p_acc_obs, p_acc_rand
    )
  )

  return(results)
}




calc_and_save_all_metrics <- function(test, prediction_folder) {
  if (!any(is.na(test$test_data$disease) | grepl("blind", test$test_data$disease, ignore.case = TRUE))) {
    ## make a folder to save figures if it doesn't
    if (!file.exists(file.path(prediction_folder, "figs"))) {
      dir.create(file.path(prediction_folder, "figs"))
    }

    ## save figures
    png(file.path(prediction_folder, "figs", "binomial_predictions.png"),
      height = 6, width = 5, res = 100, units = "in"
    )
    print(
      xyplot(test$test_pred ~ factor(test$y_test),
        type = c("g", "p", "a"), ylim = c(-0.03, 1.03),
        ylab = "P(Cancer | Data)", xlab = "Cancer", main = "Test",
        jitter.x = TRUE,
        scales = list(alternating = FALSE), auto.key = list(columns = 2)
      )
    )
    dev.off()


    ## calculate metrics
    all_metrics <- calc_binomial_metrics(test$y_test, test$test_pred)

    all_metrics$subset <- "All"

    ## by country
    list_of_countries <- unique(test$test_dat$country)
    for (i in list_of_countries) {
      metrics <- try(
        calc_binomial_metrics(test$y_test, test$test_pred,
          filter = test$test_data$country == i
        )
      )

      if (inherits(metrics, "try-error")) {
        metrics <- make_empty_metrics(test$y_test, test$test_pred,
          filter = test$test_data$country == i
        )
      }

      metrics$subset <- i
      all_metrics <- rbind(all_metrics, metrics)
    }

    ## by vendor
    list_of_vendors <- unique(test$test_data$vendor)
    for (i in list_of_vendors) {
      metrics <- try(
        calc_binomial_metrics(test$y_test, test$test_pred,
          filter = test$test_data$vendor == i
        )
      )

      if (inherits(metrics, "try-error")) {
        metrics <- make_empty_metrics(test$y_test, test$test_pred,
          filter = test$test_data$vendor == i
        )
      }

      metrics$subset <- paste0("Vendor: ", i)
      all_metrics <- rbind(all_metrics, metrics)
    }


    ## by age
    list_of_ages <- list(first = 51:60, 61:100)
    list_of_ages_labels <- c("Age: 51-60", "Age: >60")
    for (i in 1:2) {
      metrics <- try(
        calc_binomial_metrics(test$y_test, test$test_pred,
          filter = test$test_data$age %in% list_of_ages[[i]]
        )
      )

      if (inherits(metrics, "try-error")) {
        metrics <- make_empty_metrics(test$y_test, test$test_pred,
          filter = test$test_data$age %in% list_of_ages[[i]]
        )
      }

      metrics$subset <- list_of_ages_labels[i]
      all_metrics <- rbind(all_metrics, metrics)
    }

    all_metrics <- all_metrics %>%
      relocate(subset, .before = "auc")

    write.csv(all_metrics,
      file.path(prediction_folder, "test_set_all_metrics.csv"),
      row.names = FALSE, quote = FALSE, na = ""
    )



    ## test if the predictions are better than random guessing
    guess <- test_better_than_random(test$y_test, test$test_pred)
    write.csv(guess,
      file.path(prediction_folder, "p_value_for_guessing_all_metrics.csv"),
      row.names = FALSE, quote = FALSE
    )
  }
}





balanced_acc_phi <- function(actual, predicted) {
  not_missing <- !is.na(predicted)

  acc_0 <- mean(predicted[actual == 0 & not_missing] == actual[actual == 0 & not_missing])
  acc_1 <- mean(predicted[actual == 1 & not_missing] == actual[actual == 1 & not_missing])

  (acc_0 + acc_1) / 2
}

tidy_phi_data <- function(d, set) {
  ## tidy names
  d <- janitor::clean_names(d)

  ## remove people with prostatectomy
  d <- d[!d$serum_specific_id %in% c("HMN54232", "HMN54233", "HMN54248", "HMN54227"), ]


  ## combine all non BioIVT vendors
  d <- d %>%
    mutate(vendor = recode_factor(vendor,
      "BioOptions" = "Other",
      "Discovery Life Sciences" = "Other",
      "Precision for Medicine" = "Other",
      "PrecisionMed" = "Other"
    ))

  ## keep only prostate and healthy samples and males
  d <- d %>%
    filter(disease == "Prostate Cancer" | disease == "Healthy") %>%
    filter(sex == "M") %>%
    droplevels()

  ## get location of key columns
  cancer_loc <- grep("^cancer$", names(d), ignore.case = TRUE)
  phi_loc <- grep("^phi", names(d), ignore.case = TRUE)

  ## check uniqueness of column location
  if (!(length(cancer_loc) == 1 & length(phi_loc) == 1)) {
    stop("Cancer and phi columns not found or not unique")
  }

  ## ensure cancer and phi values are numbers
  d[, cancer_loc] <- as.numeric(as.character(d[, cancer_loc]))
  d[, phi_loc] <- as.numeric(as.character(d[, phi_loc]))

  ## remove rows with no phi values
  d <- d %>%
    drop_na(names(d)[phi_loc])

  ## filter blind samples (if any)
  d <- d[!grepl("blind", d[, cancer_loc], ignore.case = TRUE), ] %>%
    droplevels()

  return(list(d = d, cancer_loc = cancer_loc, phi_loc = phi_loc))
}


extract_covariates <- function(d) {
  d <- janitor::clean_names(d)

  ## select key variables
  d <- d %>%
    select(
      serum_specific_id,
      # collection_date, first_assay_date,
      # collection_country, age_at_collection, sex, cancer, disease,
      race, ethnicity, bmi, stage, benign_prostatic_hyperplasia,
      current_co_infections_morbidities, medications,
      tobacco_pack_years, alcohol_history
    ) %>%
    rename(id = serum_specific_id) %>%
    # rename(country = collection_country) %>%
    rename(bph = benign_prostatic_hyperplasia) %>%
    rename(comorbid = current_co_infections_morbidities)

  d <- d %>%
    mutate(stage = factor(gsub("N/A", "0", stage))) %>%
    mutate(bmi = as.numeric(bmi)) %>%
    mutate_if(is.character, as.factor) # %>%
  ## mutate(
  ##   first_assay_date = mdy(first_assay_date),
  ##   collection_date = mdy(collection_date)
  ## )


  ## extract comorbidities
  d$hypertension <- grepl("hypertension|HTN|hyper tension", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$diabetes <- grepl("diabetes", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$copd <- grepl("chronic obstructive pulmonary|COPD|chronic obstruction pulmonary", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$high_lipids <- grepl("hyperlipidemia|HLD|high triglycerides|hypertriglyceridemia|dyslipidemia", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$high_chol <- grepl("hypercholesterolemia|High LDL", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$gerd <- grepl("Acid reflux|GERD|Barrett|gastroesophageal reflux", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$heart_disease <- grepl("coronary artery disease|cardiovascular disease|heart disease|atherosclero|Arteriosclero|myocardial infarction|cardiovascular issues|heart failure|cardiomyopathy", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$hemochromatosis <- grepl("hemochromatosis", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$anemia <- grepl("anemia", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$prostatitis <- grepl("prostatitis", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$prostatectomy <- grepl("prostatectomy", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$thyroid_low <- grepl("thyroid removed|resection of thyroid|hypothyroidism", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$apnea <- grepl("apnea", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$arthritis <- grepl("arthritis", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$pancreatitis <- grepl("pancreatitis", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$osteo <- grepl("osteo", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$bph2 <- grepl("prosatic hyperplasia|prostatic hyerplasia|prostatic hyperplasia|BPH", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$allergy <- grepl("allerg", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  d$gout <- grepl("gout", d$comorbid, ignore.case = TRUE) %>%
    as.numeric()

  ## https://www.mayoclinic.org/diseases-conditions/benign-prostatic-hyperplasia/diagnosis-treatment/drc-20370093

  ## for BPH
  d$tamsulosin <- grepl("tamsulosin|flomax", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## for BPH
  d$cialis <- grepl("cialis|tadalafil", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## for BPH
  d$finasteride <- grepl("finasteride|proscar|propecia", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## for BPH
  d$plavix <- grepl("plavix", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## for BPH
  d$alfuzosin <- grepl("alfuzosin", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## for BPH
  d$doxazosin <- grepl("doxazosin", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## for BPH
  d$rapaflo <- grepl("rapaflo", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## for BPH
  d$avodart <- grepl("avodart", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## check if inflammatory biomarkers differ by nsaid and statin status

  d$nsaid <- grepl("aspirin|ibuprofen|acetaminophen|tylenol|naproxen|lospirin|diclofenac|meloxicam", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  d$statins <- grepl("statin|lipitor|crestor", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## synthetic thyroid hormone
  d$levothyroxine <- grepl("levothyroxine", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  ## for high lipids -> pancreatitis is a side effect
  d$fenofibrate <- grepl("fenofibrate", d$medications, ignore.case = TRUE) %>%
    as.numeric()

  return(d)
}

restrict_range <- function(train, test) {
  test[test < min(train)] <- min(train)
  test[test > max(train)] <- max(train)
  return(test)
}
