// Three-level (two hierarchical groups) nested random intercept AND slope model
// Trying to add ecosystem into phenological model
// Copied from threelevelrandomeffects2.stan
// Lizzie got as far as I could on 23 Sept 2016 but cannot figure out p and desMat below!

data{
//counters
	int<lower=0> N; 		//LEVEL 1: number of observations (total years of data across species)
	int<lower=0> Nspp; 		//LEVEL 2: number of spcecies (i.e. grouping factor)
	int<lower=0> Neco; 	//LEVEL 3: species is aquatic or terrestrial (i.e. grouping factor)
	int<lower=0> p;	 //Number of fixed effect parameters + intercept
	
//Group ids
	int<lower=1> species[N];			//interaction identity
	int<lower=1> eco[Nspp]; 	//vector of 0 or 1 for aquatic or terrestrial for each species
	
// predictors
	real desMat[N,p]; //Design matrix

// response
	real y[N]; 		//mean synch change for each interaction; Continuous
}


parameters{
	//Fixed effects
	//real mu_a[p];  // population/overall intercept
	real mu_a;  // population/overall intercept
	real mu_b;  // population/overall slope- here for year effect (need one for each predictor)
	//real mu_b[p];
	real<lower=0> sigma_y; 		// measurement error, noise etc. (Level 1) population sd
	
	//LEVEL 2: interaction or species level
	real a_spp[Nspp]; 		//the intercept for each species (random effect)
	real<lower=0> sigma_a_spp; // variation of intercept among species; (sd of random effect)
	real b_spp[Nspp]; 		//the slope of phen change for each species (random effect)
	real<lower=0> sigma_b_spp; // variation of slope among species; (sd of random effect)
	
	//LEVEL 3: ecosystem level
	real a_eco[Neco];	// intercept for each ecosystem (random effect)
	real<lower=0> sigma_a_eco; // variation of intercept among ecosystems; 
	real b_eco[Neco];	//the slope of temperature change for each ecosystem (random effect)
	real<lower=0> sigma_b_eco; /// variation of slope among ecosystems; 	
}


transformed parameters{
	//Varying intercepts
	real beta_a_spp[Nspp]; //intercept at species level
	real beta_a_eco[Neco]; //intercept at ecosystem level
	
	//Varying slopes
	real beta_b_spp[Nspp]; //slope of phen change at species level
	real beta_b_eco[Neco]; //slope of phen change at ecosystem level (Dan's b_warm_sp)
	
	//Individual mean
	real ypred[N];

	//Varying intercepts AND SLOPES definition
	//Level 3: ecosystem level
	for (k in 1:Neco){
		beta_a_eco[k]=mu_a+a_eco[k]; //Random intercepts; 
		beta_b_eco[k]=mu_b+b_eco[k]; //Random slopes (phen change);
	}
	
	//Level 2: species level
	for (j in 1:Nspp){
		beta_a_spp[j]=beta_a_eco[eco[j]]+a_spp[j]; //Random intercepts;
		beta_b_spp[j]=beta_b_eco[eco[j]]+b_spp[j]; //Random slopes;
	}
	
	//Level 1, Individual mean
	for (i in 1:N){
		ypred[i]=beta_a_spp[species[i]]+beta_b_spp[species[i]]*desMat[i,p];
	}
}

model{
	//Prior part of Bayesian inference
	//Random effects distribution
	a_eco~normal(a_eco, sigma_a_eco); //intercept for ecosystem (level 3)
	a_spp~normal(a_spp, sigma_a_spp); //intercept for indivi interaction (level 2)
	
	b_eco~normal(b_eco, sigma_b_eco); //slope for indiv ecosystem (level 3)
	b_spp~normal(b_spp, sigma_b_spp); //slope for indiv interaction (level 2)
	
	//Likelihood part of Bayesian inference
		y~normal(ypred, sigma_y);
}