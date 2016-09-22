//Simple linear regression model

data{
	int<lower=0> N; 				//Level 1: Number of observations
	
	//predictors
	vector[N] year; 	//year of data point
	
	//response
	real y[N]; 		//DOY of pheno event (OR temperature for temperature change model)

}

parameters{	
	real mu_a;				//mean intercept across species (population)
	real mu_b;				//mean slope across species (population)
	real<lower=0> sigma_y; 	//measurement error, noise etc. (population standard deviation)
}

transformed parameters{
 
 //Individual mean
 real ypred[N];
 

//Individual mean
for (i in 1:N){
		ypred[i]<-mu_a+mu_b*year[i];
	}
	
}

model{
	y~normal(ypred, sigma_y);
}	

	