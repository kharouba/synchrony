For two level RANDOM Intercepts MODEL- matches twolevelrandomsintercepts.stan
N.B true mu_b out of bounds- value is higher than range

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(rstan)
library(shinyStan)
set_cppo("fast")  # for best running speed


Nspp <- 50 # number of species

# True parameter values
mu_a <- 121.2615 # mean doy - for intercept 
mu_b<--0.349733 # mean slope from lm fits (based on hinge data)
sigma_y <- 42.30465 # sd associated with response, doy

sigma_a<- 5 #- sd of mean intercept
a<-rnorm(Nspp, mu_a, sigma_a); #generate intercepts for each species

# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
n_data_per_species <- round(runif(Nspp, 5, 40)) # how many years per sp.?
#n_data_per_species <- round(runif(Nspp, 5, 40)) # structure based on dataset; how many years per sp.?
species <- rep(1:Nspp, n_data_per_species) 
N <- length(species) #nrow of 'dataframe'
year <- rep(NA, N)
for (j in 1:Nspp){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0
}
ypred <- length(N) # Lizzie added
for (i in 1:N){ # actual STAN model
	s <- species[i] #sppid for each row
   ypred[i] <- a[species[s]] + mu_b*year[i]; #model
   }
y <- rnorm(N, ypred, sigma_y);



fit_simple<-stan("twolevelrandomintercept.stan", data=c("N","y","Nspp","species","year"), iter=2000, chains=4)

print(fit_simple, pars=c("mu_a","mu_b","sigma_a","sigma_y"))

