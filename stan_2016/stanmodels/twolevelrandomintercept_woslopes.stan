//March 17 2017
//Two-level (1 hierarchical grouping) random intercept model ==> intercepts vary and are pooled, slopes do NOT vary (and of course not pooled)

data{
	int<lower=0> N; 				//Level 1: Number of observations
	int<lower=0> Nspp; 				//Level 2: Number of species (Grouping factor)
	int species[N]; 	//species identity, coded as int
	
	//predictors
	vector[N] year; 	//year of data point
	
	//response
	real y[N]; 		//DOY of pheno event (OR temperature for temperature change model)

}

parameters{	
	// hyperparameters
	real mu_a;				//mean intercept across species (population)
	real mu_b; 		//overall slope
	
	real<lower=0> sigma_y; 	//measurement error, noise etc. (population standard deviation)
	
	real a[Nspp]; 		//the intercept for each species
	 
	real<lower=0> sigma_a;	//variation of intercept among species; [sd of random effects]
	
}

transformed parameters{
 //Individual mean
 real ypred[N];
 
//Individual mean
for (i in 1:N){
		ypred[i]<-a[species[i]]+mu_b*year[i];
	}
	
}

model{
	//Random effects distribution	
	a~normal(mu_a, sigma_a);
	y~normal(ypred, sigma_y);
}