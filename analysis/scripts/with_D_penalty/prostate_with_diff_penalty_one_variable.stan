data {
  int diff_order;
  int P;
  int N_D;
  matrix[N_D, P] D;
  
  int N;
  matrix[N, P] X;
  array[N] int y;
}

parameters {
  real b0;
  vector<lower=0>[P] b1;
  real<lower=0> sigma_penalty;
}

transformed parameters {
  vector[P - diff_order] b1_diffs = D * b1;
  vector[N] linpred = b0 + X * b1;
}

model{
  b0 ~ normal(0, 10);
  b1 ~ normal(0, 10);
  sigma_penalty ~ normal(0, 10);
  b1_diffs ~ normal(0, sigma_penalty);
  
  y ~ bernoulli_logit(linpred);
}
