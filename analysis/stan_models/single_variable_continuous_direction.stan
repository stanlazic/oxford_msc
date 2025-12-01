data {
  int P;
  int N;
  matrix[N, P] X;
  array[N] real y;
}

parameters {
  real b0;
  vector<lower=0>[P] b1;
  real<lower=0> sigma_penalty;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] linpred = b0 + X * b1;
}

model{
  b0 ~ normal(0, 20);  
  b1 ~ normal(0, sigma_penalty);
  sigma ~ normal(0, 10);
  sigma_penalty ~ exponential(1);
  
  for (i in 1:N)
    target += normal_lpdf(y[i] | linpred[i], sigma);
}
