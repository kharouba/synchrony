// 	three-level random intercept model

## data {
##   // Define variables in data
##   // Number of level-1 observations (an integer)
##   int<lower=0> Ni;
##   // Number of level-2 clusters
##   int<lower=0> Nj;
##   // Number of level-3 clusters
##   int<lower=0> Nk;
##   // Number of fixed effect parameters
##   int<lower=0> p;
## 
##   // Design matrix
##   real desMat[Ni,p];
## 
##   // Cluster IDs
##   int<lower=1> classid[Ni];
##   int<lower=1> schoolid[Ni];
## 
##   // Level 3 look up vector for level 2
##   int<lower=1> schoolLookup[Nj];
## 
##   // Continuous outcome
##   real mathgain[Ni];
##   
##   // Continuous predictor
##   // real X_1ijk[Ni];
## }
## 
## parameters {
##   // Define parameters to estimate
##   // Fixed effects
##   real beta[p];
## 
##   // Level-1 errors
##   real<lower=0> sigma_e0;
## 
##   // Level-2 random effect
##   real u_0jk[Nj];
##   real<lower=0> sigma_u0jk;
## 
##   // Level-3 random effect
##   real u_0k[Nk];
##   real<lower=0> sigma_u0k;
## }
## 
## transformed parameters  {
##   // Varying intercepts
##   real beta_0jk[Nj];
##   real beta_0k[Nk];
## 
##   // Individual mean
##   real mu[Ni];
## 
##   // Varying intercepts definition
##   // Level-3 (10 level-3 random intercepts)
##   for (k in 1:Nk) {
##     beta_0k[k] <- beta[1] + u_0k[k];
##   }
##   // Level-2 (100 level-2 random intercepts)
##   for (j in 1:Nj) {
##     beta_0jk[j] <- beta_0k[schoolLookup[j]] + u_0jk[j];
##   }
##   // Individual mean
##   for (i in 1:Ni) {
##     mu[i] <- beta_0jk[classid[i]] +
##       desMat[i,2]*beta[2] + desMat[i,3]*beta[3] + desMat[i,4]*beta[4] +
##       desMat[i,5]*beta[5] + desMat[i,6]*beta[6];
##   }
## }
## 
## model {
##   // Prior part of Bayesian inference
##   // Flat prior for mu (no need to specify if non-informative)
## 
##   // Random effects distribution
##   u_0k  ~ normal(0, sigma_u0k);
##   u_0jk ~ normal(0, sigma_u0jk);
## 
##   // Likelihood part of Bayesian inference
##   // Outcome model N(mu, sigma^2) (use SD rather than Var)
##   for (i in 1:Ni) {
##     mathgain[i] ~ normal(mu[i], sigma_e0);
##   }
}