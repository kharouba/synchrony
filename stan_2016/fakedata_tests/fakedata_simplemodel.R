For SIMPLE LINEAR REGRESSION MODEL- matches simplemodel.stan


rm(list=ls()) 
library(rstan)
library(shinyStan)
set_cppo("fast")  # for best running speed


Nspp <- 50 # number of species

# True parameter values
mu_a <- 121.2615 # mean doy - for intercept 
mu_b <- -0.349733 # mean slope from lm fits (based on hinge data)
sigma_y <- 42.30465 # sd associated with response, doy


# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
n_data_per_species <- round(runif(Nspp, 10, 10)) # how many years per sp.?
species <- rep(1:Nspp, n_data_per_species) #adds sppid-HK added
N <- length(species) #nrow of 'dataframe'
year <- rep(NA, N)
for (j in 1:Nspp){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}
ypred <- length(N) # Lizzie added

ypred<- mu_a + mu_b*year; #model

y <- rnorm(N, ypred, sigma_y);


fit_simple<-stan("simplemodel.stan", data=c("N","y","year"), iter=2000, chains=4)

print(fit_simple, pars=c("mu_a","mu_b","sigma_y"))

goo <- extract(fit_simple) 
hist(goo$mu_b) # examp
hist(goo$mu_b); abline(v = mean(mu_b), col = "blue", lwd = 2) 
