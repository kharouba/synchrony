// Major edits by M Kosmala 
// from file called synchrony1_notype_MK.stan in gelmanhill/synchrony repo

data {
  int N;                                // # data points
  int J;                                // # species
  vector[N] y;                          // DOY of pheno event
  int species[N];                       // species identity, coded as int
  vector[N] year;                       // year of data point
}
parameters {
  vector[J] a;                          // the intercept for each species
  vector[J] b;                          // the slope for each species
  real<lower=0> sigma_y;                // measurement error, noise, etc.

  // hyperparameters
  real mu_a;                            // mean intercept across species
  real<lower=0> sigma_a;                // variation of intercept among species
  real mu_b;                            // mean slope across species
  real<lower=0> sigma_b;                // variation of slope among species

  // we could also model a and b as an array/vector of 2 items
  // vector[J] beta[2];
  // vector[2] mu_beta;
  // vector[2] sigma_beta;

  // technically we should also be modeling the covariance of a and b
  // see stan-reference PDF pp 61-62 for how
}
model {
  real ypred[N];
  for (i in 1:N){
    ypred[i] <- a[species[i]] + b[species[i]] * year[i];

    // this is equivalent to the above one-liner
    // int j;
    // j <- species[i];
    // ypred[i] <- a[j] + b[j] * year[i];
  }
  y ~ normal(ypred, sigma_y);
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
} 
