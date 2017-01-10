#STAN example - plots in sites

# FROM LIZZE AND MEGAN IN HAWAII

library(rstan)
library(shinystan)

nreps = 20  #replicate observations per plot
S = 10     #total number of sites
nplots <- 8
J = nplots * S  #total number of plots 
N = nreps*nplots*S 

x <- rgamma(N,5,2)

plotnum <- rep(c(1:J), each=nreps)   #list of plot numbers for each observation (length N)
sitenum <- rep(c(1:S), each=nplots)  #list of site numbers for each plot (length J)

#mu_a and sig_a are the mean and variance of the intercept across sites
mu_a <- 5
sig_a <- 2
#mu_b and sig_b are the mean and variance of the slope across sites
mu_b <- 4
sig_b <- 1
#draw S=10 site means for slope and intercept
a_site <- rnorm(S,mu_a,sig_a)
b_site <- rnorm(S,mu_b,sig_b)
# draw S=10 site sds for slope and intercept (variability of plots within sites)
# sig_a_site <- rlnorm(S,0,1)
# sig_b_site <- rlnorm(S,0,1)
#OR draw a single sd for all sites
sig_a_site <- rlnorm(1,0,1)
sig_b_site <- rlnorm(1,0,1)

#for each plot, draw the slope and intercept from the appropriate site mean and sd
a_plot <- rep(0,J)
b_plot <- rep(0,J)
# for (j in 1:J){
#   a_plot[j] <- rnorm(1,a_site[sitenum[j]], sig_a_site[sitenum[j]]);
#   b_plot[j] <- rnorm(1,b_site[sitenum[j]], sig_b_site[sitenum[j]]);
# }
# Alternatively, assume same within-site variance for all sites
for (j in 1:J){
  a_plot[j] <- rnorm(1,a_site[sitenum[j]], sig_a_site[1]);
  b_plot[j] <- rnorm(1,b_site[sitenum[j]], sig_b_site[1]);
}

#draw three observations from each plot with sd = sig_y
sig_y <- 1
y <- rep(0,N)
yhat <- rep(0,N)
for (n in 1:N){
  yhat[n] <- a_plot[plotnum[n]] + b_plot[plotnum[n]]*x[n]
  y[n] <- a_plot[plotnum[n]] + b_plot[plotnum[n]]*x[n] + rnorm(1,0,sig_y)
}

#plot data
allsites <- rep(c(1:10),each=12) #for plotting, list of site numbers for each obs
plot(y~x,col=allsites,type='n')
text(x,y,plotnum, cex=0.5, col=allsites)

#call stan model
dat <- list(N=N,S=S,J=J,plotnum=plotnum, sitenum=sitenum,y=y,x=x) 
fitme <- stan("stanmodels/threelevel_plotsinsites.stan", data=c("N","J","S","plotnum", "sitenum","y","x"), iter=4000, chains=4, control=list(adapt_delta = 0.9, stepsize = 0.5))

dat <- list(N=N,S=S,J=J,plotnum=plotnum, sitenum=sitenum,y=y,x=x) 
fitme <- stan("stanmodels/threelevelrandomeffects4.stan", data=c("N","J","S","plotnum", "sitenum","y","x"), iter=4000, chains=4, control=list(adapt_delta = 0.9, stepsize = 0.5))

mu_a
mu_b
sig_a
sig_b
sig_y
a_plot
b_plot
a_site
b_site
sig_a_site
sig_b_site

launch_shinystan(fitme)
