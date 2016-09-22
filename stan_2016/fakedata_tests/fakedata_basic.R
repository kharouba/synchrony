# RANDOM SLOPE AND INTERCEPT MODEL #
> # Create the species-level parameters

J <- 50 # number of species 
a <- rnorm(J, 120, 41) #species error, intercept
b <- rnorm(J, -0.35, 0.499) # n=87 species
year_0 <- 1981 # hinge
sigma_y <- 42 #

# RANDOM SLOPE AND INTERCEPT MODEL #
# Create the species-level parameters

# Create the data
n_data_per_species <- round(runif(J, 5, 37)) #how many years per spp
species <- rep(1:J, n_data_per_species) #nrow of dataframe
N <- length(species)
year <- rep(NA, N)
for (j in 1:J){
 year[species==j] <- rev(2013 - 1:(n_data_per_species[j]))-year_0 
}

### HINGE MODEL ###
for (j in 1:J){
   w<-year[species==j]<=0 # or 1981
   x<-which(w==FALSE)
   k<-length(w)-length(x)
   if (k>=5){
  	d<-year[species==j]
  	d[1:k]<-0 # or 1981
   	year[species==j]<-d
   }
   }

 ypred <- length(N) # Lizzie added
for (n in 1:N){ #actual STAN model
  s <- species[n] # sppid for each row
  ypred[n] <- a[s] + b[s]*year[n]
   }
 y <- rnorm(N, ypred, sigma_y);
 
 
 # check the model works and returns the slopes it was given
# first, add in two things you need to run the model
 nVars <-1
 Imat <- diag(1, nVars)
 # then run it!

fit.hinge <- stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan/synchrony1_notype_wcovar.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)
