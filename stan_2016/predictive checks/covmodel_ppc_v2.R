# Posterior Predictive checks for covariate model (e.g. pheno change~ temp change)
#  varying and pooled intercepts, only single slope estimated

#complete pooling
comppool<-lm(abs(pheno.change)~abs(temp.change), data=tdata2)

# No pooling:
uniquespp<-unique(tdata2$species)
#uniquespp<-unique(rawlong.tot$newid)
slopefits<-rep(NA< length(uniquespp))
varfits<-rep(NA< length(uniquespp))
intfits<-rep(NA< length(uniquespp))
for(eachsp in 1:length(uniquespp)){
	lmhere<-lm(abs(pheno.change)~abs(temp.change), data=subset(tdata2, species==uniquespp[eachsp]))
	slopefits[eachsp]<-coef(lmhere)[2]
	varfits[eachsp]<-(summary(lmhere)$sigma)**2
	intfits[eachsp]<-coef(lmhere)[1]
}
dr<-data.frame(cbind(uniquespp, slopefits, varfits, intfits))
dr$slopefits<-as.numeric(dr$slopefits)
dr$intfits<-as.numeric(dr$intfits)
dr$varfits<-as.numeric(dr$varfits)

##For covariate model:
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=range(y), type="n", xlab="abs(temp change)", ylab="abs(pheno change)", bty="l", main="Raw phenology data")
for (j in 1:Nspp){
  lines(year[species==j], y[species==j])
}
# Now grab the stan output 
goo <- extract(cov.model) 

#extract the means for now:
mu_b <- mean(goo$mu_b) # from stan output (what if mu_b<- -1.47)
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
mu_a <- mean(goo$mu_a) # from stan output 
sigma_a <- mean(goo$sigma_a) # from stan output 

# Random intercepts only model:
a <- rnorm(Nspp, mean=mu_a, sd=sigma_a) # this is one way to create fake data from the Stan output to use in the PP check


## Create the PP data using new a nd b for each of 71 species
#This is just creating basically one new set of data, or one predictive check
# Varying intercepts NOT slopes:
ypred <- length(N) # Lizzie added
for (n in 1:N){
     s <- species[n]
    ypred[n] <- a[species[s]] + mu_b*year[n] # one way to handle the non-pooled intercepts, there may be other ways
}

y <- rnorm(N, ypred, sigma_y)
#Plot the data and see what the raw data predicted from the model looks like
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=range(y), type="n", xlab="abs(Temp change)", ylab="abs(Pheno change)",   bty="l", main="Data from posterior means")
for (j in 1:Nspp){
   lines(year[species==j], y[species==j])
}

