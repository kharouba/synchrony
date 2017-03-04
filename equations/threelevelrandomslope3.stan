////////// KEY TO COMMON STAN TERMS
//SIGMA= SD!!!
// sigma^2= variance
//level=intercept
//trend=slope
//ALPHA= SLOPE
//BETA= INTERCEPT
//_0=intercept
//_1=slope

// Three-level (two hierarchical groups) nested random slope model
// From Megan & Lizzie in Hawaii, started Jan 6 2017 (threelevel_plotsinsties_unpooledintercepts.stan)
// 

data{
//counters
	int<lower=0> N; 		//LEVEL 1: number of observations (2000 iterations * number of spp)
	int<lower=0> Nspp; 		//LEVEL 2: number of species (i.e. grouping factor)
	int<lower=0> Nstudy; 	//LEVEL 3: number of studies (i.e. grouping factor)
	
//Group ids
	int<lower=1> species[N];			//spp identity
//	int<lower=1> studyid[N]; 	//vector of unique studyids for each species	//OLD
	int<lower=1> studyid[Nspp]; 	//vector of unique studyids for each species // FOR FAKE DATA
	
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

	//real<lower=0> sigma_a_study[Nstudy];  //variance in intercepts across plots; 
	//real<lower=0> sigma_b_study[Nstudy];  //variance in slopes across plots;
	real<lower=0> sigma_a_study;  //variance in intercepts across plots; 
	real<lower=0> sigma_b_study;  //variance in slopes across plots;  
      // the slope for plot j in site s is drawn from a distribution with mean b_study[s] 
      // and standard deviation sig_b_site[s]
  
  real mu_b;                    //mean slope across sites; 
      // the site slope are drawn from distribution with mean mu_b...
  //real mu_a;                    //mean intercept across sites; 
  
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
  //a_spp[j] ~ normal(a_study[studyid[j]],sigma_a_study[studyid[j]]);  
  //b_spp[j] ~ normal(b_study[studyid[j]],sigma_b_study[studyid[j]]);
  
  a_spp[j] ~ normal(a_study[studyid[j]],sigma_a_study);   
  b_spp[j] ~ normal(b_study[studyid[j]],sigma_b_study); // sigma_b1
  }
  
  b_study ~ normal(mu_b,sigma_b); // sigma_b2
  sigma_a_study~normal(0, 5);
  sigma_b_study~normal(0, 5);
  sigma_b~normal(0,5);
	
//Likelihood part of Bayesian inference
		y~normal(ypred, sigma_y); //data is distributed normally around predicted (yhat) with s.d. sig_y (this is error of data around predicted values)
}
