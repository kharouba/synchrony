//SIGMA= SD!!!
// sigma^2= variance
//level=intercept
//trend=slope
//ALPHA= SLOPE
//BETA= INTERCEPT
//_0=intercept
//_1=slope

// Three-level (two hierarchical groups) nested random intercept AND slope model

data{
//counters
	int<lower=0> N; 		//LEVEL 1: number of observations (2000 iterations * number of interactions)
	int<lower=0> Nspp; 		//LEVEL 2: number of interactions (i.e. grouping factor)
	int<lower=0> Nstudy; 	//LEVEL 3: number of studies (i.e. grouping factor)
	int<lower=0> p;	 //Number of fixed effect parameters + intercept
	
//Group ids
	int<lower=1> species[N];			//interaction identity
	int<lower=1> studyid[Nspp]; 	//vector of uNque studyids for each species	
	
// predictors
//vector[N] temp; //temp change of interaction; Continuous
	real desMat[N,p]; //Design matrix
	//vector[N] year; 	//year of data point

// response
	real y[N]; 		//mean synch change for each interaction; Continuous
}


parameters{
	//Fixed effects
	//real mu_a[p];  // population/overall intercept
	real mu_a;  // population/overall intercept
	real mu_b;  // population/overall slope- here for temperature effect (need one for each predictor)
	//real mu_b[p];
	real<lower=0> sigma_y; 		// measurement error, noise etc. (Level 1) population sd
	
	//LEVEL 2: interaction or species level
	real a_spp[Nspp]; 		//the intercept for each species (random effect)
	real<lower=0> sigma_a_spp; // variation of intercept among species; (sd of random effect)
	real b_spp[Nspp]; 		//the slope of temperature change for each species (random effect)
	real<lower=0> sigma_b_spp; // variation of slope among species; (sd of random effect)
	
	//LEVEL 3: study level
	real a_study[Nstudy];	// intercept for each study (random effect)
	real<lower=0> sigma_a_study; // variation of intercept among studies; 
	real b_study[Nstudy];	//the slope of temperature change for each study (random effect)
	real<lower=0> sigma_b_study; /// variation of slope among studies; 	
}


transformed parameters{
	//Varying intercepts
	real beta_a_spp[Nspp]; //intercept at interaction level
	real beta_a_study[Nstudy]; //intercept at study level
	
	//Varying slopes
	real beta_b_spp[Nspp]; //slope of temperature change at interaction level
	real beta_b_study[Nstudy]; //slope of temperature change at study level (Dan's b_warm_sp)
	
	//Individual mean
	real ypred[N];

	//Varying intercepts AND SLOPES definition
	//Level 3: study level
	for (k in 1:Nstudy){
		beta_a_study[k]<-mu_a+a_study[k]; //Random intercepts; 
		beta_b_study[k]<-mu_b+b_study[k]; //Random slopes (temperature change);
	}
	
	//Level 2: interaction level
	for (j in 1:Nspp){
		beta_a_spp[j]<-beta_a_study[studyid[j]]+a_spp[j]; //Random intercepts;
		beta_b_spp[j]<-beta_b_study[studyid[j]]+b_spp[j]; //Random slopes (temperature change);
	}
	
	//Level 1, Individual mean
	for (i in 1:N){
		ypred[i]<-beta_a_spp[species[i]]+beta_b_spp[species[i]]*desMat[i,p];
	}
}

model{
	//Prior part of Bayesian inference
	//Random effects distribution
	a_study~normal(a_study, sigma_a_study); //intercept for indiv study (level 3)
	a_spp~normal(a_spp, sigma_a_spp); //intercept for indivi interaction (level 2)
	
	b_study~normal(b_study, sigma_b_study); //slope for indiv study (level 3)
	b_spp~normal(b_spp, sigma_b_spp); //slope for indiv interaction (level 2)
	
	//Likelihood part of Bayesian inference
		y~normal(ypred, sigma_y);
}
		
	
	