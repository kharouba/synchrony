//Two-level (1 hierarchical grouping) random intercept AND slope model

data{
	int<lower=0> N; 				//Level 1: number of observations
	int<lower=0> Nspp; 				//Level 2: number of species (i.e. grouping factor)
	int species[N]; 	//species identity
	
	//predictors
	vector[N] year; 	//year of data point
	//int<lower=0> year; 
	
	//response
	real y[N]; 		//DOY of pheno event (OR temperature for temperature change model)
}

parameters{
	// hyperparameters
	real mu_a;            // mean intercept across species (population intercept)
  	real mu_b;                 // mean slope across species  (population slope)
  	real<lower=0> sigma_y; 		//measurement error, noise etc. [population sd]

  	real a[Nspp]; 		//the intercept for each species (random effect)
	real b[Nspp]; 		//the slope for each species (random effect)
  	real<lower=0> sigma_a;     // variation of intercept among species; (sd of random effect)
  	real<lower=0> sigma_b;     // variation of slope among species; (sd of random effect)	
}

transformed parameters{
 //Individual mean
 real ypred[N];
 

//Individual mean
for (i in 1:N){
		ypred[i]<-a[species[i]]+b[species[i]]*year[i];
	}
	
}


model{
	//Random effects distribution
	a~normal(mu_a, sigma_a);		
	b~normal(mu_b, sigma_b);
	
	y~normal(ypred, sigma_y);
}
	