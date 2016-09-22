data {
  int N;                                // # data points
  vector[N] y;                          // spp diff
  //vector[N] length;                       // length of time series
  vector[N] eco2;
  int nVars;					//number of predictors
  matrix[nVars, nVars] Imat;
  //int<lower=1>K;
  //matrix[N,K] X;
}
parameters {
  real<lower=0> sigma;                // measurement error, noise, etc.
  cov_matrix[nVars] Omega3;
  real alpha; //mean intercept OR real alpha
  real beta; //mean slope   OR real beta                
  }

model {
	//Priors for covariance matrix:
	Omega3 ~ inv_wishart(nVars+1, Imat);  // inverse-Wishart prior for correlations,   
	 //Likelihood	 	
    y~normal(alpha + beta*eco2, sigma);
} 