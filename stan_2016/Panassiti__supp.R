# Authors: Bernd Panassiti, Florian Hartig, Michael Breuer, Robert Biedermann
# Date created: 10.02.2014; last modified 27.01.2015, prepared as supplement: 28.04.2015

# Supplement to "Bayesian inference of environmental and biotic factors determining the occurrence of the grapevine disease â€˜bois noirâ€™"

# Description
# RStan code to predict prevalences of a phytoplasma disease in grapevines (bois noir) in the Baden region (SW-Germany) using GLMM evaluated in a Bayesian framework.


rm(list=ls(all=TRUE))

# load libraries
library(rstan)
set_cppo("fast")  # for best running speed

# load and define data
read.table("boisnoir.txt")
grapeselection = c(5:13)   # to select grapevines
enviroselection = c(14:26) # to select remaining environmental variables


### Here begins the RStan model
# Data as a list
dataList = list(
  observedPlantsTotal=newdata$Nsum,
  observedDisease=newdata$illsum,
  nObs=length(newdata$ill),
  
  Grapes = newdata[3:(3-1+length(grapeselection))],
  nGrapes = length(grapeselection),
  
  Environment = newdata[(3+length(grapeselection)):(3-1+length(grapeselection)+length(enviroselection))],
  nEnvironment = length(enviroselection),
  
  nRegions = length(unique(newdata$regionData)),
  nOwners = length(unique(newdata$ownerData)),
  ownerData = newdata$ownerData,
  regionData = newdata$regionData
)


# Define model

model.bn.frequentgrapes <- '
data {
// counters
int<lower=0> nObs;
int<lower=0> nGrapes;
int<lower=0> nEnvironment;  
int<lower=0> nRegions; 
int<lower=0> nOwners;

// predictors
matrix[nObs, nGrapes] Grapes;
matrix[nObs, nEnvironment] Environment;

int<lower=1, upper =253> ownerData[nObs];
int<lower=1, upper=6> regionData[nObs];

// response
int<lower=0> observedPlantsTotal[nObs]; // DISEASE - num cases, total number of plants
int<lower=0> observedDisease[nObs];     // num successes
}


parameters {
// parameters that will be used in the model should be defined here
// these will be saved and printed with the output

real b0;                           // intercept
vector[nGrapes] parGrapes;
vector[nEnvironment] parEnviro;

real <lower=0> tauRegions;
real <lower=0> tauOwners;
vector[nRegions] u_raw;
vector[nOwners] v_raw;
}

transformed parameters {
// generation of derived parameters from parameters above should be defined here
// these will be saved and printed with the output

vector[nRegions] u;
vector[nOwners] v;
real realDisease[nObs];
vector[nObs] mu;

for (j in 1:nRegions) {
u[j] <- u_raw[j]*sqrt(tauRegions);
}

for (h in 1:nOwners){
v[h] <- v_raw[h]*sqrt(tauOwners);
}


for( i in 1:nObs ) {
mu[i] <- u[regionData[i]] + v[ownerData[i]];
realDisease[i]  <-(b0  + Environment[i] * parEnviro  + Grapes[i] * parGrapes + mu[i]);
}
}

model{
// This is area for describing the model including the priors and likelihood
// This is the only place where sampling can be done within Stan
// All data and parameters used in this section must be defined above

// Fixed effects priors
b0 ~ normal(-5.5, 100 ); //
parGrapes ~ normal(0 , 10 );
parEnviro ~ normal(0 , 10 );

//Hyperpriors
tauRegions ~ inv_gamma(0.001, 0.001);
tauOwners  ~ inv_gamma(0.001, 0.001);    

//Random effects priors
u_raw ~ normal(0, 1);
v_raw ~ normal(0, 1);

observedDisease ~ binomial_logit(observedPlantsTotal,realDisease);
}

generated quantities{
// predictive inference
real realDiseasePredictions[nObs]; // vector to store predictions

for( i in 1:nObs ) {
realDiseasePredictions[i]  <-binomial_rng(observedPlantsTotal[i],inv_logit(b0  + Environment[i] * parEnviro  + Grapes[i] * parGrapes + mu[i]));
}
}
'




##### Run model - may take some hours

time <- proc.time()
set.seed(2013)
stan.fit <- stan(model_code = model.bn.frequentgrapes, data = dataList, iter = 100000, chains = 4,refresh=1,thin=10)
(proc.time() - time)/60/60

