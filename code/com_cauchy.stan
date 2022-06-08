data {
  int<lower=0> N;
  int <lower=0, upper=1> trt[N];
  int <lower=0, upper=1> E[N]; // external data source indicator
  real X1[N]; // age
  real X2[N]; // sex
  real offset[N];
  int<lower=0> nu[N];
  real int_d_1[N];
  real int_d_2[N];
  real a0[N];
}
parameters {
  real beta1;
  real beta2;
  real int1;
  //// Optional: Reduce to regular exponential
  real int2;
  real gamma;
  real delta; // drift
  real<lower=0> delta_sig;
}
transformed parameters  {
  real lp[N];
  real <lower=0> mu[N];
  for (i in 1:N) {
    lp[i] = int1 * int_d_1[i] + int2 * int_d_2[i] + beta1 * X1[i] + beta2 * X2[i] + gamma * trt[i] + delta * E[i] + offset[i]; // linear predictor
    mu[i] = exp(lp[i]); // mean
  }
}
model {
  //// Old Version - Vectorized 
  // gamma ~ normal(0, 100);
  // delta ~ normal(0, delta_sig);
  // delta_sig ~ cauchy(0, 0.3);
  // nu ~ poisson(mu);

  //// New Version - Loop
  // target += normal_lpdf(gamma | 0, 100);
  target += normal_lpdf(delta | 0, delta_sig);
  target += cauchy_lpdf(delta_sig | 0, 0.3);
  for(i in 1:N) {
    target += poisson_lpmf(nu[i] | mu[i]); 
  }
}
