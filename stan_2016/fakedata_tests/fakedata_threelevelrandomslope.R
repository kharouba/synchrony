FOR RANDOM SLOPE MODEL- matches threelevelrandomslope2.stan

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(rstan)
library(shinystan)
library(rstanarm)

#Needed to create parameters
Nspp<- 50 #37 #(Old Nj)
Nstudy <- 10 #13 # number of studies (old Nk)
Nsppstudy<- 5 # number of species per study

#####################
# True parameter values
mu_a<- 121.2615 # mean doy - for intercept 
mu_b<- -0.349733 # mean slope from lm fits (based on hinge data) (SAME AS B_YR_0)
sigma_y <- 42.30465 # sd associated with response, doy

sigma_a_study<-3 #OR 5
sigma_b_study<-0.1 #- sd of mean slopes actual=0.737724

a_study<-rnorm(Nstudy, mu_a, sigma_a_study); #generate intercepts for each STUDY
b_study<-rnorm(Nstudy, mu_b, sigma_b_study); #generate slopes for each STUDY

#THEN make variation in intercepts among species with regards to study
a_spp <- rep(NA, Nspp)
for(j in 1:Nstudy){
	for(k in 1:Nsppstudy){
	a_spp[speciesid==k]<-a_study[j]+rnorm(1, 10,3) 
	}
	Nsppstudy<-Nsppstudy+5
}

#make variation in slopes among species with regards to study
b_spp<-rep(NA, Nspp)
for(j in 1:Nstudy){
	for(k in 1:Nsppstudy){	
	b_spp[speciesid==k]<-b_study[j]+rnorm(1, 0,1)
	}
	Nsppstudy<-Nsppstudy+5
}
sigma_a_spp<-2
sigma_b_spp<-0.3


# Simulate/create the data
year_0 <- 1981 # small numbers (like 0) are better than bigger numbers (like 1976)
n_spp_per_study<- round(runif(Nstudy, 5, 5)) # how many sp per study?
#n_spp_per_study<- round(runif(Nstudy, 2, 4)) # UNBALANCED: how many sp per study?

n_yrs_per_species <- round(runif(Nspp, 10, 10)) # how many years per sp.?
#n_yrs_per_species <- round(runif(Nspp, 5, 40)) # UNBALANCED: how many years per sp.?

studyid<- rep(1:Nstudy, n_spp_per_study) ## Create a vector of study IDs where j-th element gives studyid ID for spp ID j; length of species, every SPECIES gets studyid- No, see below!!!
#studyid<-c(1,1,1,1,2,2,2,3,3,3,3,4,4,4,5,5,5,6,6,7,8,8,9,9,9,10,10,10,10,11,11,11,12,12,13,13,13)
#Nspp <- length(studyid) # creates number of species based on data structure for studyid


species <- rep(1:Nspp, n_yrs_per_species) #adds sppid-HK added
N <- length(species) #nrow of 'dataframe'

uni<-as.data.frame(studyid)
uni$species<-unique(species)
uni2<-as.data.frame(species)
uni3<-merge(uni2, uni, by=c("species"))
studyid<-uni3$studyid #needs to be same length as number of lowest level/observations


#####################
year <- rep(NA, N)
for (j in 1:Nspp){
  year[species==j] <- rev(2009 - 1:(n_yrs_per_species[j])) - year_0 #assign 'new' year for each year/row for each species; from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}

# STAN model
ypred <- length(N) # Lizzie added
for (i in 1:N){ # actual STAN model
	s <- species[i] #sppid for each row
   ypred[i] <- a_spp[species[s]] + b_spp[species[s]]*year[i]; #mean? prediction is a function of vairance associated with species, fits slope with species random slope model, n loop, create data 
}
y <- rnorm(N, ypred, sigma_y);

#desMat <- model.matrix(object = ~ 1 + year)
#p<-ncol(desMat)


fit_simple<-stan("stanmodels/threelevelrandomslope2.stan", data=c("N","Nspp","Nstudy","species", "studyid","y","year"), iter=3000, chains=4)
fit_simple<-stan("stanmodels/threelevelrandomeffects3.stan", data=c("N","Nspp","Nstudy","species", "studyid","y","year"), iter=3000, chains=4)

#CHECK WHETHER PROBLEM IS THAT I'M CENTERING RANDOM EFFECTS DISTRIBUTIONS

print(fit_simple)
print(fit_simple, pars=c("b_yr_0","a_0","b_yr_spp","sigma_y"))

simple<-stan_lmer(y~year+(0+year|studyid)+(0+year|species), iter=3000)
simple<-stan_lmer(y~year+(0+year|studyid)+(0+year|species), iter=3000, prior_covariance=decov(scale=10))
simple<-stan_lmer(y~year+(0+year|studyid)+(0+year|species), iter=3000, prior_intercept=normal(0, scale=sigma), prior_covariance=decov(scale=10))

!!! NOTES
ranef(simple)
posterior_interval(womensrole_bglm_1, prob = 0.95, pars = "education")
launch_shinystan(womensrole_bglm_1)
preds <- posterior_predict(fit)
apply(preds, 2 quantile, prob = vector_of_probabilities))
prior_summary(simple)

SPECIFYING RANDOM EFFECTS
The intercept is implicitly added to the random effect parts. So (cov1 + cov2 | grouping1 )  is indeed (1 + cov1 + cov2 | grouping1 ). If I would like to get rid of the intercept I would have to write (0 + cov1 + cov2 | grouping1).


PRIORS
prior location= prior mean
prior specifies the distribution(s) for all slope parameters (coefficients) that do NOT vary by group
prior_intercept specifies the distribution for the model's GLOBAL intercept
prior_covariance specifies the ENTIRE structure for all parameters that vary by group because they have zero mean
decov priors= random effects= variance parameters on the diagnoal of the covariance matrix
there will be 1 scale parameter for each grouping (e.g. 2 because For (cov1 | grouping)+ (cov1 | grouping), there will be 2 scale parameter)




