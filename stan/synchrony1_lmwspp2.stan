data {
  int<lower=0> N; //#data points. There are N observations, each with predictor x[n] and outcome y[n]
  vector [N] y; // y vector, DOY of pheno event
  vector [N] year; //x vector, year of data point
  int species[N]; //species identity coded as int
    }
parameters {
  real alpha[J];                          // the intercept for each species
  real beta[J];                          // the slope for each species
  real<lower=0> sigma[J];                // measurement error, noise, etc.
//one model per species
}
model {
  for (j in 1:J){
  y[j]~normal(alpha[j]+beta[j]*year[j], sigma[j]);
	}
}
  
 // for (s in 1:J){
//	for (n in 1:N){
//  y[s,n]~normal(alpha[s]+beta[s]*year[s,n], sigma[s]);