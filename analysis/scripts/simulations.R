## https://yahrmason.github.io/bayes/gams-julia/

library(tidyverse)
library(cmdstanr)
library(lattice)
# library(caret)
library(splines2)
# library(CalibrationCurves)
# library(predtools)

##HERE, try if scam can recover the parameters

source("../src/functions.R")


set.seed(8)
N <- 100
x <- sort(runif(N))
B <- iSpline(x, df = 4, intercept=FALSE, degree = 2)


intercept <- -3
coefs <- 1:6 * 0.6## intercept <- -2
## coefs <- rep(1.0, 5)

y_hat <- intercept + B %*% coefs
y <- rbinom(N, 1, plogis(y_hat))


xyplot(factor(y) ~ x, # data=, #groups=,
  type = c("g", "p"),
  xlab = "", ylab = "", main = "",
  jitter.y = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)



compiled_model <- cmdstan_model("../stan_models/single_variable_continuous.stan",
  dir = "../compiled_models"
)


dat <- list(
  P = ncol(B),
  N = N,
  X = B,
  y = y
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
  # invlogit() %>%
  colMeans()

est_coefs <- dr %>%
  select_at(vars(starts_with("b")))

summary(est_coefs)
boxplot(est_coefs, outline = FALSE)
points(c(intercept, coefs), col = "red", pch = 19)


dr %>%
  select_at(vars(starts_with("sigma"))) %>%
  boxplot(., outline = FALSE)


xyplot(y ~ x,
  type = c("g", "p"),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "Bayes GAM model",
  ylim = c(-0.25, 1.25),
  jitter.y = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

xyplot(etas ~ x,
  type = c("g", "p"),
  xlab = "Cancer", ylab = "Predicted probability of cancer", main = "Bayes GAM model",
  # ylim = c(0.75, 2.25),
  jitter.y = TRUE, scales = list(alternating = FALSE, tck = c(1, 0)),
  auto.key = list(columns = 2)
)

plot(jitter(y, factor = 0.1) ~ x)
points(plogis(etas) ~ x, col = "red", pch = 19)
abline(v=knots(B), lty=2, col="steelblue")
abline(v=knots(B, "boundary"), lty=2, col="steelblue")

plot(y_hat ~ etas)
abline(0, 1, lty = 2)

sqrt(mean((y_hat - etas)^2))


## ============================================================
## With diff penalty
## ============================================================


compiled_model <- cmdstan_model("with_D_penalty/prostate_with_diff_penalty_one_variable.stan",
  dir = "../compiled_models"
)


diff_order <- 2
D <- diff(diag(ncol(B)), diff = diff_order)

dat <- list(
  diff_order = diff_order,
  N_D = nrow(D),
  D = D,
  P = ncol(B),
  N = N,
  X = B,
  y = y
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
  # invlogit() %>%
  colMeans()

est_coefs <- dr %>%
  select_at(vars(starts_with("b")))

summary(est_coefs)
boxplot(est_coefs, outline = FALSE)
points(c(intercept, coefs), col = "red", pch = 19)


## ============================================================
## continuous
## ============================================================

set.seed(8)
N <- 100
x <- sort(runif(N))
B <- iSpline(x, df=6, intercept=TRUE, degree = 3)
# df - degree - intercept = number of knots


y_hat <- 5 + 2 * x^2 - 3 * x^3 + 0.5 * sin(5 * x)
y <- rnorm(N, y_hat, 0.3) * -1

plot(y ~ x)

## linear
y_hat <- 5 + 2*x
y <- rnorm(N, y_hat, 0.3)
plot(y ~ x)


compiled_model <- cmdstan_model("../stan_models/single_variable_continuous.stan",
  dir = "../compiled_models"
)


dat <- list(
  P = ncol(B),
  N = N,
  X = B,
  y = y
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
  # invlogit() %>%
  colMeans()

est_coefs <- dr %>%
  select_at(vars(starts_with("b")))

summary(est_coefs)
boxplot(est_coefs, outline = FALSE)
points(c(intercept, coefs), col = "red", pch = 19)


plot(y ~ x)
points(etas ~ x, col = "red", pch = 19)
points(etas ~ x, col = "steelblue", pch = 19)
abline(v=knots(B), lty=2, col="steelblue")
abline(v=knots(B, "boundary"), lty=2, col="steelblue")

plot(I(y_hat * -1)~ etas)
abline(0, 1, lty = 2)
