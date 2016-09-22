// Major edits by M Kosmala 
// from file called synchrony1_notype_MK.stan in gelmanhill/synchrony repo

data {
  int N;                                // # data points
  int J;                                // # species
  vector[N] y;                          // DOY of pheno event
  int species[N];                       // species identity, coded as int
  vector[N] year;                       // year of data point
  int nVars;					//number of preidictors
  matrix[nVars, nVars] Imat;
}
parameters {
  vector[J] a;                          // the intercept for each species
  vector[J] b;                          // the slope for each species
  real<lower=0> sigma_y;                // measurement error, noise, etc.
  cov_matrix[nVars] Omega3;
 
  // hyperparameters
  real mu_a;                            // mean intercept across species
  real<lower=0> sigma_a;                // variation of intercept among species; implicit uniform prior
  real mu_b;                            // mean slope across species
  real<lower=0> sigma_b;                // variation of slope among species; implicit uniform prior
}

transformed parameters {
   vector[N] log_lik;
for(n in 1:N){
	log_lik[n]<-normal_log(y[n],year[n]*b, sigma_y);	
	}
}

model {
increment_log_prob(-log(sigma_y));    //log prior for p(sigma) proportional to 1/sigma
increment_log_prob(sum(log_lik));   //log likelihood
} 

