data {
  int P;
  int N;
  matrix[N, P] X;
  array[N] int y;
}

parameters {
  real b0;
  vector<lower=0>[P] b1;
}

transformed parameters {
  vector[N] linpred = b0 + X * b1;
}

model{
  b0 ~ normal(0, 50);  
  b1 ~ normal(0, 50);
  
  for (i in 1:N)
    target += bernoulli_logit_lpmf(y[i] | linpred[i]);
}
