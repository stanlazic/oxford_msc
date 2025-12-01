data {
  int P;
  
  int N_train;
  matrix[N_train, P] X1_train;
  matrix[N_train, P] X2_train;
  matrix[N_train, P] X3_train;
  vector[N_train] X4_train;
  array[N_train] int y_train;
}

parameters {
  real b0;
  vector<lower=0>[P] b1;
  vector<upper=0>[P] b2;
  vector<upper=0>[P] b3;
  real b4;
  real<lower=0> sigma_penalty1;
  real<lower=0> sigma_penalty2;
  real<lower=0> sigma_penalty3;
}

transformed parameters {
  vector[N_train] linpred = b0 + X1_train * b1 + X2_train * b2 +
    X3_train * b3 + X4_train * b4;

}

model{
  b0 ~ normal(0, 50);
  
  b1 ~ normal(0, sigma_penalty1);
  b2 ~ normal(0, sigma_penalty2);
  b3 ~ normal(0, sigma_penalty3);
  b4 ~ normal(0, 5);

  // might need different sigma_penalty for each biomarker
  sigma_penalty1 ~ exponential(0.5);
  sigma_penalty2 ~ exponential(0.5);
  sigma_penalty3 ~ exponential(0.5);
  
  for (i in 1:N_train)
    target += bernoulli_logit_lpmf(y_train[i] | linpred[i]);
}

generated quantities {
  vector[N_train] log_lik;

  for (i in 1:N_train)
    log_lik[i] = bernoulli_logit_lpmf(y_train[i] | linpred[i]);
}
