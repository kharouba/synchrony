FOR RANDOM SLOPE and INTERCEPT MODEL- matches threelevelrandomeffects.stan

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(rstan)
library(shinyStan)

#Needed to create parameters
#Nspp<-50
Nstudy<- 10 # number of studies


# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
n_data_per_study<- round(runif(Nstudy, 2, 8)) # how many sp per study?
studyid<- rep(1:Nstudy, n_data_per_study) ## Create a vector of study IDs where j-th element gives studyid ID for spp ID j; length of species, every species gets studyid
Nspp <- length(studyid) # creates number of species based on data structure for studyid

n_data_per_species <- round(runif(Nspp, 5, 40)) # how many years per sp.?

species <- rep(1:Nspp, n_data_per_species) #adds sppid-HK added
N <- length(species) #nrow of 'dataframe'

#####################
# True parameter values
mu_a<-1; mu_a[1:N] <- 121.2615 # mean doy - for intercept 
mu_b<-1; mu_b[1:N]<- -0.349733 # mean slope from lm fits (based on hinge data)
sigma_y <- 42.30465 # sd associated with response, doy

sigma_a_spp<-5 #sd of intercepts
a_spp<-rnorm(Nspp, mu_a, sigma_a_spp)

sigma_b_spp<-0.737724 #- sd of mean slopes
b_spp<-rnorm(Nspp, mu_b, sigma_b_spp); #generate slopes for each species

sigma_a_study<-3
a_study<-rnorm(Nstudy, mu_a, sigma_a_study)

sigma_b_study<-0.2
b_study<-rnorm(Nstudy, mu_b, sigma_b_study)

#####################
year <- rep(NA, N)
for (j in 1:Nspp){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}
ypred <- length(N) # Lizzie added
for (i in 1:N){ # actual STAN model
	s <- species[i] #sppid for each row
   ypred[i] <- a_spp[species[s]] + b_spp[species[s]]*year[i]; #mean? prediction is a function of vairance associated with species, fits slope with species random slope model, n loop, create data 
}


y <- rnorm(N, ypred, sigma_y);
desMat <- model.matrix(object = ~ 1 + year)
p<-ncol(desMat)


fit_simple<-stan("threelevelrandomeffects.stan", data=c("N","Nspp","Nstudy","species", "studyid","year","p","desMat"), iter=2000, chains=4)

print(fit_simple)
#print(fit_simple, pars=c("mu_b","a","b","sigma_y"))
