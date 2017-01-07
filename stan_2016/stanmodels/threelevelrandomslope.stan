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

data{
//counters
	int<lower=0> N; 		//LEVEL 1: number of observations (2000 iterations * number of spp)
	int<lower=0> Nspp; 		//LEVEL 2: number of species (i.e. grouping factor)
	int<lower=0> Nstudy; 	//LEVEL 3: number of studies (i.e. grouping factor)
//	int<lower=0> p;	 //Number of fixed effect parameters + intercept
	
//Group ids
	int<lower=1> species[N];			//spp identity
	int<lower=1> studyid[Nspp]; 	//vector of unique studyids for each species	
	
// predictors
//vector[N] temp; //temp change of interaction; Continuous
//real desMat[N,p]; //Design matrix
	vector[N] year; 	//year of data point

// response
	real y[N]; 		//mean synch change for each interaction; Continuous
}


parameters{
	//Fixed effects
	//Level 1:
	real a;  // population/overall intercept
	real mu_b;  // population/overall slope- here for temperature effect (need one for each predictor)
	real<lower=0> sigma_y; 		// measurement error, noise etc. (Level 1) population sd
	
	
	//LEVEL 2: interaction or species level
	real b_spp[Nspp]; 		//the slope of temperature change for each species (random effect)
	real<lower=0> sigma_b_spp; // variation of slope among species; (sd of random effect)
	
	//LEVEL 3: study level
	real b_study[Nstudy];	//the slope of temperature change for each study (random effect)
	real<lower=0> sigma_b_study; /// variation of slope among studies; 	
}


transformed parameters{
	//Varying slopes
	real beta_b_spp[Nspp]; //slope of temperature change at interaction level
	real beta_b_study[Nstudy]; //slope of temperature change at study level (Dan's b_warm_sp)
	
	//Individual mean
	real ypred[N];


	//Varying SLOPES definition
	//Level 3: study level
	for (k in 1:Nstudy){
		beta_b_study[k]<-mu_b+b_study[k]; //Random slopes (temperature change);
	}
	
	//Level 2: interaction level
	for (j in 1:Nspp){
		beta_b_spp[j]<-beta_b_study[studyid[j]]+b_spp[j]; //Random slopes (temperature change);
	}
	
	//Level 1, Individual mean
	for (i in 1:N){
		ypred[i]<-a[species[i]]+beta_b_spp[species[i]]*year[i];
	}
}

model{
	//Prior part of Bayesian inference
	//Random effects distribution
	//a_study~normal(a_study, sigma_a_study); //intercept for indiv study (level 3)
	//a_spp~normal(a_spp, sigma_a_spp); //intercept for indivi interaction (level 2)
	
	b_study~normal(0, sigma_b_study); //slope for indiv study (level 3)
	b_spp~normal(0, sigma_b_spp); //slope for indiv interaction (level 2)
	
	//Likelihood part of Bayesian inference
		y~normal(ypred, sigma_y);
}
		
	
	