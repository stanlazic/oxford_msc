data {
  int diff_order;
  int P;
  int N_D;
  matrix[N_D, P] D;
  
  int N_train;
  matrix[N_train, P] X1_train;
  matrix[N_train, P] X2_train;
  matrix[N_train, P] X3_train;
  matrix[N_train, P] X4_train;
  matrix[N_train, P] X5_train;
  array[N_train] int y_train;

  /* int N_test; */
  /* matrix[N_test, P] X1_test; */
  /* matrix[N_test, P] X2_test; */
  /* vector[N_test] y_test; */
}

parameters {
  real b0;
  vector[P] b1;
  vector[P] b2;
  vector[P] b3;
  vector[P] b4;
  vector[P] b5;
  // real<lower=0> sigma;
  real<lower=0> sigma_penalty;
}

transformed parameters {
  vector[P - diff_order] b1_diffs = D * b1;
  vector[P - diff_order] b2_diffs = D * b2;
  vector[P - diff_order] b3_diffs = D * b3;
  vector[P - diff_order] b4_diffs = D * b4;
  vector[P - diff_order] b5_diffs = D * b5;
  vector[N_train] eta = b0 + X1_train * b1 + X2_train * b2 +
    X3_train * b3 + X4_train * b4 + X5_train * b5;
}

model{
  b0 ~ normal(0, 10);
  b1 ~ normal(0, 10);
  b2 ~ normal(0, 10);
  b3 ~ normal(0, 10);
  b4 ~ normal(0, 10);
  b5 ~ normal(0, 10);
  /* sigma ~ normal(0, 10); */
  /* b1_diffs ~ normal(0, 0.5); // was 0.05 */
  /* b2_diffs ~ normal(0, 0.5); */
  /* b3_diffs ~ normal(0, 0.5); */
  sigma_penalty ~ normal(0, 1);
  b1_diffs ~ normal(0, sigma_penalty);
  b2_diffs ~ normal(0, sigma_penalty);
  b3_diffs ~ normal(0, sigma_penalty);
  b4_diffs ~ normal(0, sigma_penalty);
  b5_diffs ~ normal(0, sigma_penalty);
  
  y_train ~ bernoulli_logit(eta);
}
