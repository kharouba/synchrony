FOR RANDOM INTERCEPT MODEL- matches threelevelrandomintercept.stan

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(rstan)
library(shinyStan)

#Needed to create parameters
Nstudy<- 10 # number of studies

# Start to simulate/create the data
year_0 <- 1981 
n_data_per_study<- round(runif(Nstudy, 4, 4)) # Balanced; how many sp per study?
#n_data_per_study<- round(runif(Nk, 2, 8)) # how many sp per study?
studyid<- rep(1:Nstudy, n_data_per_study) ## Create a vector of study IDs where j-th element gives studyid ID for spp ID j; length of species, every species gets studyid
Nspp <- length(studyid) # creates number of species based on data structure for studyid

#n_data_per_species <- round(runif(Nj, 10, 10)) # how many years per sp.?
n_data_per_species <- round(runif(Nspp, 5, 40)) # how many years per sp.?

species <- rep(1:Nspp, n_data_per_species) #adds sppid-HK added
N <- length(species) #nrow of 'dataframe'


#####################
# True parameter values
mu_a<-121.2615 # mean doy - for intercept 
mu_b<- -0.349733 # mean slope from lm fits (based on hinge data)
sigma_y <- 42.30465 # sd associated with response, doy

sigma_a<-5 #sd of intercepts
a<-rnorm(Nspp, mu_a, sigma_a) #generate slopes for each species


#####################
#Continue to simulate/create the data
year <- rep(NA, N)
for (j in 1:Nspp){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}
ypred <- length(N) # Lizzie added
for (i in 1:N){ # actual STAN model
	s <- species[i] #sppid for each row
   ypred[i] <- a[species[s]] + mu_b*year[i];
}
plot(year, ypred)

y <- rnorm(N, ypred, sigma_y);
p<-2

desMat <- model.matrix(object = ~ 1 + year)

fit_simple<-stan("threelevelrandomintercept.stan", data=c("N","Nspp","Nstudy","species", "studyid","desMat","p"), iter=2000, chains=4)

print(fit_simple, pars=c("mu_a","mu_b","sigma_a", "sigma_y"))
