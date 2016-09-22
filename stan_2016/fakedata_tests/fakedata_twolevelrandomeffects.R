FOR RANDOM SLOPE and INTERCEPT MODEL- matches twolevelrandomeffects.stan

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(rstan)
library(shinyStan)


Nspp <- 50 # number of species (needed to generate parameters)

# True parameter values
mu_a <- 121.2615 # mean doy - for intercept 
mu_b<- -0.349733 # mean slope from lm fits (based on hinge data)
sigma_y <- 41 # sd associated with response, doy # actual= 42.30465

sigma_a<-2 #sd of intercepts; actual= 5
a<-rnorm(Nspp, mu_a, sigma_a) #generate intercepts for each species

sigma_b<-0.499 #- sd of mean slopes; actual = 0.737724
b<-rnorm(Nspp, mu_b, sigma_b); #generate slopes for each species


# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
n_data_per_species <- round(runif(Nspp, 10, 10)) # how many years per sp.?
species <- rep(1:Nspp, n_data_per_species) #adds sppid-HK added
N <- length(species) #nrow of 'dataframe'
year <- rep(NA, N)
for (j in 1:Nspp){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
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


ypred <- length(N) # Lizzie added
for (i in 1:N){ # actual STAN model
	s <- species[i] #sppid for each row
   ypred[i] <- a[species[s]] + b[species[s]]*year[i]; #mean? prediction is a function of vairance associated with species, fits slope with species random slope model, n loop, create data 
}
y <- rnorm(N, ypred, sigma_y);


fit_simple<-stan("twolevelrandomeffects.stan", data=c("N","y","Nspp","species","year"), iter=3000, chains=4)

print(fit_simple, pars=c("mu_a","mu_b","sigma_a","sigma_b","sigma_y"))

