#FOR two level RANDOM SLOPE MODEL- matches twolevelrandomslope.stan


rm(list=ls()) 
library(rstan)
library(shinyStan)
set_cppo("fast")  # for best running speed


Nspp <- 50 # number of species

# True parameter values
mu_a<- 121.2615 # mean doy - for intercept 
sigma_a<-5
mu_b<- -0.349733 # mean slope from lm fits (based on hinge data)
sigma_b<-0.1 #- sd of mean slopes actual=0.737724
sigma_y <- 42.30465 # sd associated with response, doy

a<-rnorm(Nspp, mu_a, sigma_a);
b<-rnorm(Nspp, mu_b, sigma_b); #generate slopes for each species

# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
#n_data_per_species <- round(runif(Nspp, 10, 10)) #balanced; how many years per sp.?
n_data_per_species <- round(runif(Nspp, 5, 35)) #based on dataset; how many years per sp.?
species <- rep(1:Nspp, n_data_per_species) #adds sppid
N <- length(species) #nrow of dataframe
year <- rep(NA, N)
for (j in 1:Nspp){
  #year[species==j] <- 1:(n_data_per_species[j])
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1981, sort in descending order, series of years for each spp
}
ypred <- length(N) 
for (i in 1:N){ 
	s <- species[i] #sppid for each row
   ypred[i] <- a[species[s]] + b[species[s]]*year[i]; #model
   }
y <- rnorm(N, ypred, sigma_y);

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



fit_simple<-stan("twolevelrandomslope.stan", data=c("N","y","Nspp","species","year"), iter=4000, chains=4)

print(fit_simple, pars=c("mu_b","sigma_b","sigma_y"))
