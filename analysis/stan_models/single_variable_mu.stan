data {
  int P;
  int N;
  matrix[N, P] X;
  array[N] int y;
}

parameters {
  real b0;
  real<lower=0> mu;
  vector[P] b1;
  real<lower=0> sigma_penalty1;
}

transformed parameters {
  vector[N] linpred = b0 + X * b1;
}

model{
  b0 ~ normal(0, 50);  
  b1 ~ normal(mu, sigma_penalty1);
  mu ~ normal(0, 10);
  sigma_penalty1 ~ exponential(0.25);
  
  for (i in 1:N)
    target += bernoulli_logit_lpmf(y[i] | linpred[i]);
}
