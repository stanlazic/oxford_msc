library(tidyverse)
library(splines2)
library(cmdstanr)

## nonlinear wiggly function
make_y <- function(x) {
  1 + sin(2 * pi * x) + 0.5 * sin(4 * pi * x)
}

N <- 100
x <- seq(0, 1, length.out = N)

set.seed(123)
y <- make_y(x) + rnorm(length(x), sd = 0.2)

d <- data.frame(x = x, y = y)

knots <- seq(0, 1, length.out = 6)[-c(1, 6)]



cm_direction <- cmdstan_model("../stan_models/single_variable_continuous_direction.stan",
  dir = "../compiled_models"
)

cm_no_direction <- cmdstan_model("../stan_models/single_variable_continuous_no_direction.stan",
  dir = "../compiled_models"
)


B1 <- iSpline(x, df = 6, intercept = TRUE)


dat <- list(
  P = ncol(B1),
  N = N,
  X = B1,
  y = y
)


m1 <- cm_direction$sample(
  data = dat,
  chains = 5,
  parallel_chains = 5,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  seed = 123,
  refresh = 0
)


m2 <- cm_no_direction$sample(
  data = dat,
  chains = 5,
  parallel_chains = 5,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  seed = 123,
  refresh = 0
)


etas1 <- m1$draws(format = "df") %>%
  select_at(vars(starts_with("linpred"))) %>%
  colMeans()


etas2 <- m2$draws(format = "df") %>%
  select_at(vars(starts_with("linpred"))) %>%
  colMeans()


pdf("../../typst/figs/monotonic_example.pdf", height = 6, width = 6, bg = "transparent")
par(
  las = 1,
  mfrow = c(2, 1),
  cex = 1.2,
  mar = c(4.5, 4.5, 2, 1)
)
plot(y ~ x,
  data = d, pch = 19, cex = 0.60, xlab = "X", ylab = "Y",
  main = "No monotonic constraint"
)
abline(v = c(0, knots, 1), lty = 2)
lines(x, etas2, col = "#3C78D8", lwd = 2.5)

plot(y ~ x,
  data = d, pch = 19, cex = 0.60, xlab = "X", ylab = "Y",
  main = "Monotonic constraint enforced"
)
abline(v = c(0, knots, 1), lty = 2)
lines(x, etas1, col = "#3C78D8", lwd = 2.5)

dev.off()


## ============================================================
## Show linear is returned when relationship is linear
## ============================================================

cm_direction <- cmdstan_model("../stan_models/single_variable_continuous_direction.stan",
  dir = "../compiled_models"
)

cm_no_direction <- cmdstan_model("../stan_models/single_variable_continuous_no_direction.stan",
  dir = "../compiled_models"
)

N <- 200
x <- seq(0, 1, length.out = N)

B1 <- iSpline(x, df = 6, intercept = TRUE)


res100 <- data.frame(matrix(0, nrow = 100, ncol = 3))
names(res100) <- c("lm_rmse", "monospline_rmse", "spline")

res200 <- data.frame(matrix(0, nrow = 100, ncol = 2))
names(res200) <- c("lm_rmse", "monospline_rmse", "spline")


## HERE.. do prediction on held out data!!
set.seed(123)
for (i in 1:100) {
  y_hat <- 0 + 1.2 * x
  y <- y_hat + rnorm(length(x), sd = 0.2)
  d <- data.frame(x = x, y = y)

  lm1 <- lm(y ~ x, data = d)

  dat <- list(
    P = ncol(B1),
    N = N,
    X = B1,
    y = y
  )

  m1 <- cm_direction$sample(
    data = dat,
    chains = 5,
    parallel_chains = 5,
    iter_warmup = 1000,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    seed = 123,
    refresh = 0
  )

  m2 <- cm_no_direction$sample(
    data = dat,
    chains = 5,
    parallel_chains = 5,
    iter_warmup = 1000,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    seed = 123,
    refresh = 0
  )


  etas1 <- m1$draws(format = "df") %>%
    select_at(vars(starts_with("linpred"))) %>%
    colMeans()


  etas2 <- m2$draws(format = "df") %>%
    select_at(vars(starts_with("linpred"))) %>%
    colMeans()


  res100[i, 1] <- sqrt(mean((predict(lm1) - y_hat)^2))
  res100[i, 2] <- sqrt(mean((etas1 - y_hat)^2))
  res100[i, 3] <- sqrt(mean((etas2 - y_hat)^2))

  res200[i, 1] <- sqrt(mean((predict(lm1) - y_hat)^2))
  res200[i, 2] <- sqrt(mean((etas1 - y_hat)^2))
  res200[i, 3] <- sqrt(mean((etas2 - y_hat)^2))
}




par(
  las = 1,
  mfrow = c(1, 2),
  cex = 1.2,
  mar = c(4.5, 4.5, 2, 1)
)
boxplot(res100,
  ylim = c(0, 0.8),
  ylab = "RMSE",
  main = "N = 100",
  names = c("Linear\nModel", "Monotonic\nSpline", "Unconstrained\nSpline")
)

boxplot(res200,
  ylim = c(0, 0.8),
  ylab = "RMSE",
  main = "N = 200",
  names = c("Linear model", "Spline model")
)

matplot(t(res100), type = "l")

BHH2::dotPlot(res100[, 1] - res100[, 2],
  main = "N = 100"
)
abline(v=0, lty=2)


par(
  las = 1,
  cex = 1.2,
  mar = c(4.5, 4.5, 2, 1)
)
plot(y ~ x,
  data = d, pch = 19, cex = 0.60, xlab = "X", ylab = "Y",
  main = ""
)
abline(v = c(0, B1$knots, 1), lty = 2)
lines(x, etas2, col = "#3C78D8", lwd = 2.5)
lines(x, etas1, col = "#3C78D8", lwd = 2.5)
lines(x, predict(lm1), col = "red", lwd = 2.5)
lines(x, y_hat, col = "green", lwd = 2.5)
