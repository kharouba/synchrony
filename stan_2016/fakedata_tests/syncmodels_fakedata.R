FOR RANDOM SLOPE MODEL- matches synchrony_apr_nocov

rm(list=ls()) 


J <- 91 # number of species (needed to create parameters)

# True parameter values
a<-1 #need to create vector first
a[1:J] <- 121.2615 # mean doy - for intercept

mu_b<--0.349733 # mean slope from lm fits (based on hinge data)
sigma_y <- 42.30465 # sd associated with response, doy

sigma_b<-0.737724 #sd of mean slopes
b<-rnorm(J, mu_b, sigma_b); #generate slopes for each species


# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
n_data_per_species <- round(runif(J, 5, 40)) # how many years per sp.?
species <- rep(1:J, n_data_per_species) #adds sppid-HK added
N <- length(species) #nrow of 'dataframe'
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}
ypred <- length(N) # Lizzie added
for (i in 1:N){ # actual STAN model
	s <- species[i] #sppid for each row
   ypred[i] <- a[species[s]] + b[species[s]]*year[i]; #mean? prediction is a function of vairance associated with species, fits slope with species random slope model, n loop, create data 
}
y <- rnorm(N, ypred, sigma_y);
nVars<-1 

fit_simple<-stan("synchrony_apr_nocov.stan", data=c("N","y","J","species","year","nVars"), iter=2000, chains=4)

print(fit_simple)
#print(fit_simple, pars=c("mu_b","a","b","sigma_y"))
