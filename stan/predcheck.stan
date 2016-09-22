data {
  int N;                                // # data points
  int J;                                // # species
  int species[N];                       // species identity, coded as int
  vector[N] year;                       // year of data point
  int nVars;					//number of preidictors
  matrix[nVars, nVars] Imat;
  real<lower=0> sigma_y;                // measurement error, noise, etc.
  cov_matrix[nVars] Omega3;
 
  // hyperparameters
  real mu_a;                         // mean intercept across species
  real<lower=0> sigma_a; //variation of intercept among species; implicit uniform prior
  real mu_b;                            // mean slope across species
  real<lower=0> sigma_b;  
}
parameters {
}
model {
	
	//priors for covariance matrix:
	 Omega3 ~ inv_wishart(nVars+1, Imat);  // inverse-Wishart prior for correlations, 
//a ~ normal(mu_a, sigma_a);
  //b ~ normal(mu_b, sigma_b);
}
generated quantities {
real ypred[N];
vector[N] y;
vector[J] a;                          // the intercept for each species
vector[J] b;                          // the slope for each species
 for (i in 1:N){
    ypred[i] <- a[species[i]] + b[species[i]] * year[i];
  }
  y ~ normal(ypred, sigma_y);
}


schools_sim<-extract(fit.notype) 
n_sims<-length(schools_sim$lp__)
y_rep<-array(NA, c(n_sims,J)) 
for(s in 1:n_sims)
 y_rep[s,]<-rnorm(J, schools_sim$theta[s,],sigma)
