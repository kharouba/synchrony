#STAN example - plots in sites

# FROM LIZZE AND MEGAN IN HAWAII- edited to work with synchrony data

library(rstan)
library(shinystan)
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")

nreps = 20  #replicate observations per plot
Nstudy = 10     #total number of sites
nindiv <- 8
Nspp = nindiv * Nstudy  #total number of plots 
N = nreps*nindiv*Nstudy 

year <- rgamma(N,5,2)

species <- rep(c(1:Nspp), each=nreps)   #list of plot numbers for each observation (length N)
studyid <- rep(c(1:Nstudy), each=nindiv)  #list of site numbers for each plot (length J)

#mu_a and sig_a are the mean and variance of the intercept across sites
mu_a <- 5
sig_a <- 2
#mu_b and sig_b are the mean and variance of the slope across sites
mu_b <- 4
sig_b <- 1

# To get top level (study or site)
#draw S=10 site means for slope and intercept
a_study <- rnorm(Nstudy,mu_a,sig_a)
b_study <- rnorm(Nstudy,mu_b,sig_b)
# draw S=10 site sds for slope and intercept (variability of plots within sites)
# sig_a_site <- rlnorm(S,0,1)
# sig_b_site <- rlnorm(S,0,1)
#OR draw a single sd for all sites
sig_a_study <- rlnorm(1,0,1)
sig_b_study <- rlnorm(1,0,1)

#for each plot, draw the slope and intercept from the appropriate site mean and sd
a_spp <- rep(0,Nspp)
b_spp <- rep(0,Nspp)
# for (j in 1:J){
#   a_plot[j] <- rnorm(1,a_site[sitenum[j]], sig_a_site[sitenum[j]]);
#   b_plot[j] <- rnorm(1,b_site[sitenum[j]], sig_b_site[sitenum[j]]);
# }
# Alternatively, assume same within-site variance for all sites
for (j in 1:Nspp){
  a_spp[j] <- rnorm(1,a_study[studyid[j]], sig_a_study[1]);
  b_spp[j] <- rnorm(1,b_study[studyid[j]], sig_b_study[1]);
}

#draw three observations from each plot with sd = sig_y
sig_y <- 1
y <- rep(0,N)
ypred <- rep(0,N)
for (n in 1:N){
  ypred[n] <- a_spp[species[n]] + b_spp[species[n]]*year[n]
  y[n] <- a_spp[species[n]] + b_spp[species[n]]*year[n] + rnorm(1,0,sig_y)
}

#plot data
allsites <- rep(c(1:10),each=12) #for plotting, list of site numbers for each obs
plot(y~x,col=allsites,type='n')
text(x,y,plotnum, cex=0.5, col=allsites)

#call stan model

# Random slopes
fitme<-stan("stanmodels/threelevelrandomslope3.stan", data=c("N","Nspp","Nstudy","y","species","studyid","year"), iter=4000, chains=4,control=list(adapt_delta = 0.9, stepsize = 0.5))

sig_a_study = 0.98
sig_b_study = 0.20
mu_b =4
sig_b <- 1
sig_y<-1


# Random slopes AND intercepts

fitme<-stan("stanmodels/threelevelrandomeffects4.stan", data=c("N","Nspp","Nstudy","y","species","studyid","year"), iter=4000, chains=4,control=list(adapt_delta = 0.9, stepsize = 0.5))

mu_a<-5
sig_a <- 2
mu_b<-4
sig_b <- 1
sig_y =1 

