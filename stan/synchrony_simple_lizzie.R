#actual data
N <- nrow(longpheno)
y <- longpheno$measurement
J <- length(unique(longpheno$uniquespp))
species <- as.numeric(as.factor(longpheno$uniquespp))
year <- longpheno$yr1976
type <- as.numeric(as.factor(longpheno$variable))


# Simulate fake data

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")

# RANDOM INTERCEPT MODEL #
# Create the species-level parameters
J <- 50 # number of species
type <- rep(c(1,2), c(25,25)) # consumer or resource (1 or 2)
level <- c(125, 125) # mean doy for C and R
sigma_level <- c(25, 25) # variance for C and R
species_level <- rnorm(J, 0, 1) # error for each species, random draw

## no trend in simple
## so commenting out below
# trend <- c(-.1, -.2)
# sigma_trend <- c(.2, .2)
# species_trend <- rnorm(J, 0, 1)

year_0 <- 1976 # small numbers (like 0) are better than bigger numbers (like 1976)
sigma_y <- 5 # variance associated with response, doy

# Create the data
n_data_per_species <- round(runif(J, 5, 40)) # how many years per sp.?
species <- rep(1:J, n_data_per_species) #adds sppid-HK added
N <- length(species) #nrow of 'dataframe'
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0 #from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}
ypred <- length(N) # Lizzie added
for (n in 1:N){ # actual STAN model
  s <- species[n] #sppid for each row
  t <- type[s] # assigns whether consumer or resource for each species
  ypred[n] <- level[t] + sigma_level[t]*species_level[s] #mean? prediction is a function of whether you are Cor R and variance associated with Cor R and vairance associated with species, fits slope with species random intercept model, n loop, create data 
}
y <- rnorm(N, ypred, sigma_y); 

# Plot the data

pdf("raw1simple.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw data (types 1 and 2 are blue and red)")
for (j in 1:J)
  lines(year[species==j], y[species==j], col=colors[type[j]])
dev.off()

# Fit the model

library("rstan")

fit_simple <- stan("synchrony_simple_lizzie.stan", data=c("N","y","J","species","year","type"), iter=100, chains=4)
print(fit_simple)


# Look at the results
library("shinyStan")
launch_shinystan(fit_simple)

