# Simulate fake data

rm(list=ls())
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")

# RANDOM SLOPE AND INTERCEPT MODEL #

# Create fake data from Stan model
goo <- extract(fit.hinge) #based on hinge model from STAN

goo <- extract(fit.notype) #regular, no hinge

sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
sigma_a <- mean(goo$sigma_a) # from stan output 
sigma_b <- mean(goo$sigma_b) # from stan output 
mu_b <- mean(goo$mu_b) # from stan output 
mu_a <- mean(goo$mu_a) # from stan output
sd_b<-sd(goo$mu_b)
sd_a <-sd(goo$mu_a)


J<-87
a <-rnorm(J, mu_a, sigma_a) # sd_a or sigma_a ?sd of
b <- rnorm(J, mu_b, sigma_b) #n=87 species, sd should be 0.43
year_0<-1981
sigma_y #overall measurement error should be 11.1575

# Create the data
n_data_per_species <- round(runif(J, 5, 37)) #how many years per spp #37
species <- rep(1:J, n_data_per_species) #nrow of dataframe
N <- length(species)
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2013 - 1:(n_data_per_species[j]))-year_0 #from first year of study, number of years that differ from 1976, rev:like sort in descending order-HK added, series of years for each spp
}


## NO HINGE ##
ypred <- length(N) # Lizzie added
for (n in 1:N){ #actual STAN model
  s <- species[n] #sppid for each row
  ypred[n] <- a[s] + b[s]*year[n]
   }
y <- rnorm(N, ypred, sigma_y);

nVars<-1
Imat <- diag(1, nVars)
fit.nohinge <- stan("synchrony1_notype_hk.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)

zzno<-summary(fit.nohinge)
intercepts_nohin <- as.vector(zzno[[1]][,1])[1:87]
spptrends_nohin<- as.vector(zzno[[1]][,1])[88:174]
sppSE_nohin<- as.vector(zzno[[1]][,2])[88:174]

#compare b (fake data) created from original stan model to new b from new stan model created from fake data
min<-min(b); max<-max(b)
plot(spptrends_nohin~b, ylim=c(min, max)); abline(0,1)


pdf("raw1_nohinge.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw data"); abline(b=0)
for (j in 1:J)
  lines(year[species==j], y[species==j])
dev.off()

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
  s <- species[n] #sppid for each row
  ypred[n] <- a[s] + b[s]*year[n]
   }
y <- rnorm(N, ypred, sigma_y);

nVars<-1
Imat <- diag(1, nVars)
fit.fakehinge <- stan("synchrony1_notype_hk.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)

zzfk<-summary(fit.fakehinge)
intercepts_hin <- as.vector(zzfk[[1]][,1])[1:87]
spptrends_hin<- as.vector(zzfk[[1]][,1])[88:174]
sppSE_hin<- as.vector(zzfk[[1]][,2])[88:174]

#compare b (fake data) created from original stan model to new b from new stan model created from fake data
min<-min(b); max<-max(b)
plot(spptrends_hin~b, ylim=c(min, max)); abline(0,1)

plot(b~as.vector(zzfk[[1]][,1])[88:174])


# Plot the data
pdf("raw1_hinge.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(x=year, y=y, type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw data", xlim=c(-17,31), ylim=c(39,209));
for (j in 1:J)
  lines(year[species==j], y[species==j])
dev.off()

#notes
# Create the species-level parameters
J <- 87 # number of species 
#level <- 120 # mean phenovalue and SE is 39
a <- rnorm(J, 120, 41) #species error, intercept
#trend <- rep(-0.3,J) # mean slope?; range is -1.75 to 1.21 chosen range based on slopes from synchrony_notype_hk
#sigma_trend <- rep(0.005,J)  # variance around slope? range is 0.0019-0.0096 chosen range from SE around slope from synchrony_notype_hk
b <- rnorm(J, -0.32, 0.499) #n=87 species
year_0 <- 1981 #nohinge
sigma_y<- 5 #??

