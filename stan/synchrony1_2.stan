data {
  int N;                                  // # data points
  real y[N];
  int J;                                  // # species
  int species[N];
  real year[N];
  }
//int<lower=1,upper=2> type[J];           // species type
parameters {
  real<lower=0> sigma_y;
  real species_trend[J];         // trend for species, relative to trend for type
  real species_level[J];
  real<lower=0> sigma_species_trend[J];
  real<lower=0> sigma_species_level[J];
}
//real trend[2];                          // avg trend for type
//
//real level[2];                          // avg level for type
//  real<lower=0> sigma_level[2];  

model {
  real ypred[N];
  species_level ~ normal(0,1);
  species_trend ~ normal(0,1);
  for (n in 1:N){
    int s;
    s <- species[n];
    ypred[n] <- sigma_species_level[s]*species_level[s] + (sigma_species_trend[s]*species_trend[s])*year[n];
  }
  y ~ normal(ypred, sigma_y);
} 
//int t;
//t <- type[s];


//level=intercept
//trend=slope
//sigma=variance