////////// KEY TO COMMON STAN TERMS
//SIGMA= SD!!!
// sigma^2= variance
//level=intercept
//trend=slope
//ALPHA= SLOPE
//BETA= INTERCEPT
//_0=intercept
//_1=slope

// Three-level (two hierarchical groups) nested random intercept AND slope model
// From Megan & Lizzie in Hawaii, started Jan 6 2017
// 

data{
//counters
	int<lower=0> N; 		//LEVEL 1: number of observations (2000 iterations * number of spp)
	int<lower=0> Nspp; 		//LEVEL 2: number of species (i.e. grouping factor)
	int<lower=0> Nstudy; 	//LEVEL 3: number of studies (i.e. grouping factor)
	
//Group ids
	int<lower=1> species[N];			//spp identity
	int<lower=1> studyid[N]; 	//vector of unique studyids for each species	
	
// predictors
	vector[N] year; 	//year of data point

// response
	real y[N]; 		//mean synch change for each interaction; Continuous
}


parameters{
	vector[Nstudy] a_study;    //estimated intercept for each plot
	vector[Nstudy] b_study;    //estimated slope for each plot
	vector[Nspp] a_spp;    //estimated intercept for each site
	vector[Nspp] b_spp;    //estimated slope for each site

  	real<lower=0> sigma_a_study[Nstudy];  //variance in intercept across plots;
	real<lower=0> sigma_b_study[Nstudy];  //variance in slopes across plots; 
      // the slope for plot j in site s is drawn from a distribution with mean b_study[s] 
      // and standard deviation sig_b_site[s]
  
  	real mu_a;                    //mean intercept across sites; 
	real<lower=0> sigma_a;          //...and standard deviation sig_a
	
  	real mu_b;                    //mean slope across sites; 
      // the site slope are drawn from distribution with mean mu_b...
  
  	real<lower=0> sigma_b;          //...and standard deviation sig_b
  	real<lower=0> sigma_y;          // observation error

}

model{
	real ypred[N]; 
	
	for (i in 1:N) {
    ypred[i] = a_spp[species[i]] + b_spp[species[i]]*year[i];
  	}

  //For estimating a single value for all within-site variances
  for (j in 1:Nspp){    
   	a_spp[j] ~ normal(a_study[studyid[j]], sigma_a_study[studyid[j]]);
  	b_spp[j] ~ normal(b_study[studyid[j]],sigma_b_study);
  }
  
  a_study ~ normal(mu_a,sigma_a);
  b_study ~ normal(mu_b,sigma_b);
  
  sigma_a_study~normal(0, 10);
  sigma_a~normal(0,10);

  sigma_b_study~normal(0, 10);
  sigma_b~normal(0,10);
	
//Likelihood part of Bayesian inference
		y~normal(ypred, sigma_y); //data is distributed normally around predicted (yhat) with s.d. sig_y (this is error of data around predicted values)
}
