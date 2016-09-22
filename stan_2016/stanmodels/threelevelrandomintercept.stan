//SIGMA= SD!!!
//ALPHA= SLOPE
//BETA= INTERCEPT
//_0=intercept
//_1=slope

// Three-level nested random intercept model

data{
//counters
	int<lower=0> N; 				//LEVEL 1: number of observations 
	int<lower=0> Nspp; 				//LEVEL 2: number of interactions (i.e. grouping factor)
	int<lower=0> Nstudy; 				//LEVEL 3: number of studies (i.e. grouping factor)
	int<lower=0> p;
	
	int species[N];			//interaction identity
	int studyid[N]; 	//study identity, coded as int

// predictors
	real year[N]; 
	
// Design matrix
//  real desMat[N,p];

// response
	real y[N]; 		//mean synch change for each interaction; Continuous
}


parameters{
	// hyperparameters
	real mu_a;				//mean intercept across species (population)
	real mu_b;				//mean slope across species (population)
	real<lower=0> sigma_y; 	//measurement error, noise etc. (population standard deviation)
		
	//LEVEL 2:
	real a_spp[Nspp]; 		//the intercept for each interaction
	real<lower=0> sigma_a_spp; // variation of intercept among interactions; 
	
	//LEVEL 3:
	real a_study[Nstudy];	// intercept for each study
	real<lower=0> sigma_a_study; // variation of intercept among studies; 
	}


transformed parameters{
	//Varying intercepts
	real beta_0spp[Nspp];
	real beta_0study[Nstudy];
	
		
	//Individual mean
	real ypred[N];

	//Varying intercepts definition
	//Level 3: 
	for (k in 1:Nstudy){
		beta_0study[k]<-mu_a+a_study[k]; 
	}
	
	//Level 2:
	for (j in 1:Nspp){
		beta_0spp[j]<-beta_0study[studyid[j]]+a_spp[j];
	}
	
	//Level 1, Individual mean
	for (i in 1:N){
		ypred[i]<-beta_0spp[species[i]]+mu_b*year[i];
		}
	}

model{
	//Prior part of Bayesian inference
	//Flat prior for mu
	
	//Random effects distribution
	a_study~normal(a_study, sigma_a_study); //intercept for study (level 3)
	a_spp~normal(a_spp, sigma_a_spp); //intercept for interaction (level 2)
	
	y~normal(ypred, sigma_y);
	
}
		
	
	