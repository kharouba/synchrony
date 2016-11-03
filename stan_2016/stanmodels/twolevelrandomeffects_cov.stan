//Two-level (1 hierarchical grouping) random slope model

data{
	int<lower=0> N; 				//Level 1: Number of observations
	int<lower=0> Nspp; 				//Level 2: Number of species (Grouping factor)
	int species[N]; 	//species identity, coded as int
	
	//predictors
	int<lower=0> p; 	//number of fixed effect parameters
	real desMat[N, p]; //Design matrix
	
	//response
	real y[N]; 		//DOY of pheno event (OR temperature for temperature change model)

}

parameters{	
	// hyperparameters
	real mu_b;				//mean slope across species (population)
	real mu_a;				//mean intercept across species (population)
	real<lower=0> sigma_y; 	//measurement error, noise etc. (population standard deviation)
	
	real a[Nspp]; 		//the intercept for each species
	real b[Nspp]; 		//the slope for each species 
	real<lower=0> sigma_b;	//variation of slope among species; [sd of random effects]
	real<lower=0> sigma_a;     // variation of intercept among species; (sd of random effect)

}

transformed parameters{
 //Individual mean
 real ypred[N];
 

//Individual mean
for (i in 1:N){
		ypred[i]<-a[species[i]]+b[species[i]]*desMat[i,2]+b[species[i]]*desMat[i,3]+b[species[i]]*desMat[i,4];
	}
	
}

model{
	//Random effects distribution	
	a~normal(mu_a, sigma_a);
	b~normal(mu_b, sigma_b); #not drawing from normal distribution, it is calculating the join t distribution function (in log space); it is evaluating the normal prpopability distribution function of b given 0 and sigma
	
	y~normal(ypred, sigma_y);
}	
