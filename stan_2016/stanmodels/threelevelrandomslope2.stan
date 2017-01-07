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
// from lday_ind2_indint.stan from Dan Flynn
// Two-level linear model in Stan

data{
//counters
	int<lower=0> N; 		//LEVEL 1: number of observations (2000 iterations * number of spp)
	int<lower=0> Nspp; 		//LEVEL 2: number of species (i.e. grouping factor)
	int<lower=0> Nstudy; 	//LEVEL 3: number of studies (i.e. grouping factor)
//	int<lower=0> p;	 //Number of fixed effect parameters + intercept
	
//Group ids
	int<lower=1> species[N];			//spp identity
	int<lower=1> studyid[N]; 	//vector of unique studyids for each species	
	
// predictors
	vector[N] year; 	//year of data point

// response
	real y[N]; 		//mean synch change for each interaction; Continuous
}


parameters{
	//Fixed effects
	//Level 1:
//	real a_0;  // population/overall intercept (same as beta_0 in classroom ex)
	real b_yr_0;  // overall year effect
	real<lower=0> sigma_y; 		// measurement error, noise etc. (Level 1) population sd

	
	//LEVEL 3: study level (top level)
//	real dev_a_study[Nstudy];	// deviation of intercepts across studies
	real dev_b_yr_study[Nstudy];	//	deviation of slopes across studies
//	real<lower=0> sigma_a_study; /// sd (variation) of deviation of intercept among studies;
	real<lower=0> sigma_b_yr_study; /// sd (variation) of deviation of slope among studies; 
	
	//LEVEL 2: species level
//	real dev_a_spp[Nspp]; // deviation of intercepts across species
	real dev_b_yr_spp[Nspp]; //	deviation of slopes across speices 
//	real<lower=0> sigma_a_spp; // (SD) variation of intercept among species; (sd of random effect)
	real<lower=0> sigma_b_spp; // (SD) variation of slope among species; (sd of random effect)
	
	real a_spp[Nspp]; //intercept for species within studies
}


transformed parameters{
	//Varying intercepts
	//real a_study[Nstudy]; //intercept for studies
	//real a_spp[Nspp]; //intercept for species within studies
	
	//Varying slopes
	real b_yr_study[Nstudy]; //slope of temperature (or year) at study level (random effect)(Dans b_warm_sp)	
	real b_yr_spp[Nspp]; //slope of temperature (or year) for species within studies
	
	//Individual mean
	real ypred[N];

	//Top Level/Level 3/Study level. Random intercepts (a) and slopes
	for (k in 1:Nstudy){
	//	a_study[k]=a_0+dev_a_study[k];
		b_yr_study[k]=b_yr_0+dev_b_yr_study[k]; //Random slopes (e.g. temperature change);
	}
	
	//Level 2: spp level. Random intercepts and slopes
	for (j in 1:Nspp){
	//	a_spp[j]=a_study[studyid[j]]+dev_a_spp[j];
		b_yr_spp[j]=b_yr_study[studyid[j]]+dev_b_yr_spp[j]; //Random slopes (e.g. temperature change);
	}
	
	//Level 1, Individual mean
	for (i in 1:N){
		ypred[i]=a_spp[species[i]]+b_yr_spp[species[i]]*year[i];
	}
}

model{
	//Prior part of Bayesian inference
	//Random effects distribution; Distribution of the varying slopes
	
//	dev_a_study~normal(0, sigma_a_study); // to allow pooling of intercepts
	dev_b_yr_study~normal(0, sigma_b_yr_study); // to allow pooling of slopes # or try b_yr_0
	
//	dev_a_spp~normal(0, sigma_a_spp); // to allow pooling of intercepts
	dev_b_yr_spp~normal(0, sigma_b_spp); // to allow pooling of slopes # or try b_yr_0

//	sigma_a_study~normal(0,10);
	sigma_b_yr_study~normal(0,10);
//	sigma_a_spp~normal(0,10);
	sigma_b_spp~normal(0,10);

	//Likelihood part of Bayesian inference
		y~normal(ypred, sigma_y);
}
