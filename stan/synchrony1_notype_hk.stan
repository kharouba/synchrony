// Major edits by M Kosmala 
// from file called synchrony1_notype_MK.stan in gelmanhill/synchrony repo

data {
  int N;                                // # data points
  int J;                                // # species
  vector[N] y;                          // DOY of pheno event
  //int y_prior_mean;
  //int y_prior_sd;
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

model {
	real ypred[N];
	//Priors for covariance matrix:
	 Omega3 ~ inv_wishart(nVars+1, Imat);  // inverse-Wishart prior for correlations,   
	 //Likelihood
	 //y~normal(y_prior_mean, y_prior_sd);	 	
  for (i in 1:N){
    ypred[i] <- a[species[i]] + b[species[i]] * year[i];
  }
  y ~ normal(ypred, sigma_y);
  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
increment_log_prob(-log(sigma_y));    //log prior for p(sigma) proportional to 1/sigma
//increment_log_prob(sum(log_lik));   //log likelihood
} 

generated quantities {
 real tau;
  real residual[N];
 real predicted[N];
  real sq[N];
//Assess model fit using a sums of squares type discrepancy
for (i in 1:N){
residual[i]<-y[i]-mu_b;
predicted[i]<-mu_b;
sq[i]<-pow(residual[i],2); //Squared residuals for observed data
}
tau<-1/(sigma_b*sigma_b);
}

//Extra to get covariance to run:

// for Prior for covariance matrix
//tau (τ)= vector for scales
//omega (Ω)= correlation matrix

//mu_a~normal(0,100);
//mu_b~normal(0,100); 

//transformed parameters{
//  cov_matrix[K] Omega; //variance-covariance matrix
//  corr_matrix[K] phi;  

  //Intercept }
  
  //for covariance matrix
 //corr_matrix[2] Rho_alpha;
    //vector<lower=0>[2] sigma_alpha; 

//model{
//priors for covariance matrix:
 //corr_matrix [K] Omega; //prior correlation; omega (Ω)= correlation matrix for random intercepts and slopes. Include number of predictors in []. 
 //vector<lower=0>[K] tau; //prior scale; tau (τ)= vector for scales
 
 
 //matrix[K,K] Omega_alpha;
	//Sigma_beta<-quad_form_diag(Omega, tau);
	//Omega_alpha<-quad_form_diag(Rho_alpha, sigma_alpha);
 
 //sigma_alpha ~ cauchy(0,2.5); //OR tau
	//Rho_alpha ~ lkj_corr(2);//OR omega

