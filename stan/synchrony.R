# Simulate fake data

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")


# RANDOM SLOPE AND INTERCEPT MODEL #
# Create the species-level parameters
J <- 50 # number of species 
type <- rep(c(1,2), c(25,25)) # consumer or resource
level <- c(125, 125) # mean
sigma_level <- c(25, 25) # variance
species_level <- rnorm(J, 0, 1) #species error
trend <- c(-.1, -.2) # mean slope?
sigma_trend <- c(.2, .2)  # variance around slope?
species_trend <- rnorm(J, 0, 1)
year_0 <- 1976
sigma_y <- 5


# Create the data
n_data_per_species <- round(runif(J, 5, 40))
species <- rep(1:J, n_data_per_species)
N <- length(species)
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0
}
ypred <- length(N) # Lizzie added
for (n in 1:N){
  s <- species[n]
  t <- type[s]
  ypred[n] <- level[t] + sigma_level[t]*species_level[s] +
    (trend[t] + sigma_trend[t]*species_trend[s])*year[n] # Lizzie changed year[i] to year[n]
}
y <- rnorm(N, ypred, sigma_y);

# Plot the data

pdf("raw1.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw data (types 1 and 2 are blue and red)")
for (j in 1:J)
  lines(year[species==j], y[species==j], col=colors[type[j]])
dev.off()

# Fit the model

library("rstan")
source("stan.R")

fit_simple <- stan("synchrony_simple.stan", data=c("N","y","J","species","year","type"), iter=1000, chains=4)
print(fit_simple)

# Look at the results
library("shinyStan")
launch_shinystan(fit_simple)

fit1 <- stan("synchrony1.stan", data=c("N","y","J","species","year","type"), iter=1000, chains=4)
print(fit1)
launch_shinystan(fit1)
