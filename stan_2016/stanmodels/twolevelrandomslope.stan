//Two-level (1 hierarchical grouping) random slope model

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
	real mu_b;				//mean slope across species (population)
	real<lower=0> sigma_y; 	//measurement error, noise etc. (population standard deviation)
	
	real a[Nspp]; 		//the intercept for each species
	real b[Nspp]; 		//the slope for each species 
	#real<lower=0> sigma_a;	//variation of intercept among species; [sd of random effects]
	real<lower=0> sigma_b;	//variation of slope among species; [sd of random effects]

}

transformed parameters{
 //Varying intercepts
 real beta_a[Nspp];

 //Varying slopes
 real beta_0j[Nspp];
 
 
 //Individual mean
 real ypred[N];
 
 //Varying intercepts definition
 for (j in 1:Nspp){
 	beta_a[j]<-mu_a+a[j];
 	}

 //Varying slopes definition
 for (j in 1:Nspp){
 	beta_0j[j]<-mu_b+b[j];
 }

//Individual mean
for (i in 1:N){
		ypred[i]<-beta_a[species[i]]+beta_0j[species[i]]*year[i];
	}
	
}

model{
	//Random effects distribution	
	b~normal(mu_b, sigma_b); #not drawing from normal distribution, it is calculating the join t distribution function (in log space); it is evaluating the normal prpopability distribution function of b given 0 and sigma
	
	y~normal(ypred, sigma_y);
}	

	