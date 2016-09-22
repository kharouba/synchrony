FOR RANDOM SLOPE MODEL- matches threelevelrandomslope.stan

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(rstan)
library(shinyStan)

#Needed to create parameters
Nspp<-50
#Nk<- 10 # number of studies

#####################
# True parameter values
mu_a<-1; mu_a[1:Nspp] <- 121.2615 # mean doy - for intercept 
mu_b<- -0.349733 # mean slope from lm fits (based on hinge data)
sigma_y <- 42.30465 # sd associated with response, doy
sigma_b<-0.737724 #- sd of mean slopes

b<-rnorm(Nj, mu_b, sigma_b); #generate slopes for each species

# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
n_data_per_study<- round(runif(Nk, 4, 4)) # how many sp per study?
#n_data_per_study<- round(runif(Nk, 2, 8)) # how many sp per study?
studyid<- rep(1:Nk, n_data_per_study) ## Create a vector of study IDs where j-th element gives studyid ID for spp ID j; length of species, every species gets studyid
Nj <- length(studyid) # creates number of species based on data structure for studyid

#n_data_per_species <- round(runif(Nj, 5, 40)) # how many years per sp.?
n_data_per_species <- round(runif(Nj, 10, 10)) # how many years per sp.?

species <- rep(1:Nj, n_data_per_species) #adds sppid-HK added
Ni <- length(species) #nrow of 'dataframe'



#####################
year <- rep(NA, Ni)
for (j in 1:Nj){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}
ypred <- length(Ni) # Lizzie added
for (i in 1:Ni){ # actual STAN model
	s <- species[i] #sppid for each row
   ypred[i] <- mu_a + b[species[s]]*year[i]; #mean? prediction is a function of vairance associated with species, fits slope with species random slope model, n loop, create data 
}


y <- rnorm(Ni, ypred, sigma_y);
desMat <- model.matrix(object = ~ 1 + year)
p<-ncol(desMat)


fit_simple<-stan("threelevelrandomslope.stan", data=c("Ni","Nj","Nk","species", "studyid","year","p","desMat"), iter=2000, chains=4)

print(fit_simple)
#print(fit_simple, pars=c("mu_b","a","b","sigma_y"))
