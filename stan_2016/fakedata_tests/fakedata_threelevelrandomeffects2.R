FOR RANDOM SLOPE and INTERCEPT MODEL- matches threelevelrandomeffects2.stan

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(rstan)
library(shinyStan)

#Needed to create parameters
Nstudy<- 10 # number of studies
Nspp<- 50 # number of species
Nsppstudy<- 5 # number of species per study

speciesid<-1:50
studyid<-rep(1:Nstudy, each=Nsppstudy) # Assign each species to a study

#####################
# True parameter values
mu_a<- 121.2615 # mean doy - for intercept 
mu_b<- -0.349733 # mean slope from lm fits (based on hinge data)
sigma_y <- 42.30465 # sd associated with response, doy

sigma_a_study<-3
a_study<-rnorm(Nstudy, mu_a, sigma_a_study) #different intercept for each study

sigma_b_study<-0.2
b_study<-rnorm(Nstudy, mu_b, sigma_b_study) #different slope for each study

#make variation in intercepts among species with regards to study
a_spp <- rep(NA, Nspp)
for(j in 1:Nstudy){
	for(k in 1:Nsppstudy){
	a_spp[speciesid==k]<-a_study[j]+rnorm(1, 10,3) 
	}
	Nsppstudy<-Nsppstudy+5
}

#make variation in slopes among species with regards to study
b_spp<-rep(NA, Nspp)
Nsppstudy<-5
for(j in 1:Nstudy){
	for(k in 1:Nsppstudy){	
	b_spp[speciesid==k]<-b_study[j]+rnorm(1, 0,1)
	}
	Nsppstudy<-Nsppstudy+5
}
sigma_a_spp<-2
sigma_b_spp<-0.3

# Simulate/create the data
year_0 <- 1981 #hinge year
n_spp_per_study<- round(runif(Nstudy, 5, 5)) # how many sp per study?
n_yrs_per_species <- round(runif(Nspp, 10, 10)) # how many years per sp.?

species <- rep(1:Nspp, n_yrs_per_species) #adds sppid
N <- length(species) #nrow of 'dataframe'


#####################
year <- rep(NA, N)
for (j in 1:Nspp){
  year[species==j] <- rev(2009 - 1:(n_yrs_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}


## Hinge model
for (j in 1:Nspp){
   w<-year[species==j]<=0 # or 1981
   x<-which(w==FALSE)
   k<-length(w)-length(x)
   if (k>=5){
  	d<-year[species==j]
  	d[1:k]<-0 # or 1981
   	year[species==j]<-d
   }
   }


# STAN model
ypred <- length(N)
for (i in 1:N){ 
	s <- species[i] #sppid for each row
   ypred[i] <- a_spp[species[s]] + b_spp[species[s]]*year[i];
}


y <- rnorm(N, ypred, sigma_y);
desMat <- model.matrix(object = ~ 1 + year)
p<-ncol(desMat)


fit_simple<-stan("threelevelrandomeffects2.stan", data=c("N","Nspp","Nstudy","species", "studyid","y","year","p","desMat"), iter=3000, chains=4)

print(fit_simple, pars=c("mu_a","mu_b","sigma_y","sigma_a_study","sigma_b_study","sigma_a_spp","sigma_b_spp"))

