For two level random intercepts mODEL- matches twolevelrandomsintercepts_woslopes.stan
Varying and pooled intercepts ONLY, (only 1 slope estimated)
==> USED FOR COVARIATE MODEL (PHENO CHANGE~ TEMP CHANGE)

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(rstan)
library(shinyStan)
set_cppo("fast")  # for best running speed


Nspp <- 50 # number of species

# True parameter values
mu_a <- 0.71 # mean pheno change - for intercept (mean(abs(tdata2$pheno.change)))
mu_b<- -1.47 # slope from complete pooling 
sigma_y <- 0.55 # sd associated with response, sd(abs(tdata2$pheno.change))

sigma_a<- 5 #- sd of mean intercept
a<-rnorm(Nspp, mu_a, sigma_a); #generate intercepts for each species
sigma_b<-0.015 #sd(abs(tdata2$temp.change))

# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
species <- rep(1:Nspp, 1000) 
N <- length(species) #nrow of 'dataframe'
year<-rep(NA, N)
year<-rnorm(N, mean(abs(tdata2$temp.change)), sd(abs(tdata2$temp.change))); #generate data

ypred <- length(N) 
for (i in 1:N){ # actual STAN model
	s <- species[i] #sppid for each row
   ypred[i] <- a[species[s]] + mu_b*year[i]; #model
   }
y <- rnorm(N, ypred, sigma_y);



fit_simple<-stan("stanmodels/twolevelrandomintercept_woslopes.stan", data=c("N","y","Nspp","species","year"), iter=2000, chains=4)

print(fit_simple, pars=c("mu_a","mu_b","sigma_a","sigma_y"))

