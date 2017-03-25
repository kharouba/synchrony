# Posterior Predictive checks for covariate model (e.g. pheno change~ temp change)
#  varying and pooled intercepts, only single slope estimated

#complete pooling

comppool<-lm(pheno.change~temp.change, data=tdata2)
with(tdata2, plot(pheno.change~temp.change)); with(tdata2, abline(lm(pheno.change~temp.change))

# No pooling:
uniquespp<-unique(tdata2$species)
#uniquespp<-unique(rawlong.tot$newid)
slopefits<-rep(1< length(uniquespp))
varfits<-rep(1< length(uniquespp))
intfits<-rep(1< length(uniquespp))
for(eachsp in 1:length(uniquespp)){
	lmhere<-lm(pheno.change~temp.change, data=subset(tdata2, species==uniquespp[eachsp]))
	slopefits[eachsp]<-coef(lmhere)[2]
	varfits[eachsp]<-(summary(lmhere)$sigma)**2
	intfits[eachsp]<-coef(lmhere)[1]
}
dr<-data.frame(cbind(uniquespp, slopefits, varfits, intfits))
dr$slopefits<-as.numeric(levels(dr$slopefits))[dr$slopefits]
dr$intfits<-as.numeric(levels(dr$intfits))[dr$intfits]
dr$varfits<-as.numeric(levels(dr$varfits))[dr$varfits] #factor to numeric

##For covariate model:

option 1:
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=range(y), type="n", xlab="Temp change", ylab="Pheno change", bty="l", main="No pooling- linear models")
for (j in 1:Nspp){
  abline(lm((y[species==j])~year[species==j]))
}

option 2:
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=range(y), type="n", xlab="Temp change", ylab="Pheno change", bty="l", main="No pooling intercept, complete pooling slope")
for (j in 1:Nspp){
	m1<-lm((y[species==j])~year[species==j])
  abline(a=coef(m1)[1], b=coef(comppool)[2])
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
plot(range(year), range(y), ylim=range(y), type="n", xlab="Temp change", ylab="Pheno change",   bty="l", main="Data from posterior means")
for (j in 1:Nspp){
   abline(a=(y[species==j]), b=mu_b)
}


### posterior slope from Stan vs histogram of no pooling slopes, and complete pooling
hist(dr$slopefits, main="Dist'n of no pooling slopes w posterior slope (blue) and complete pooling slope (red)")
abline(v = mu_b, col = "blue", lwd = 2) # mean of lm fits
abline(v=coef(comppool)[2], col="red", lwd=2)

##posterior intercepts from stan vs. histogram of no pooling slopes
hist(dr$intfits, col=rgb(1,0,0,0.5), main="Dist'n of no pool ints (pink) w post ints (blue), mu_a (blue), complete pool int (red)")
hist(a, xlim=range(dr$intfits), col=rgb(0,0,1,0.5), add=TRUE)
abline(v = mu_a, col = "blue", lwd = 2) # mean of lm fits
abline(v=coef(comppool)[1], col="red", lwd=2)

# okay, but predictive checks are about much more than ONE simulation, one draw ...
# Create the data using new a and b for each of 71 species, 100 times
# RANDOM intercept model
ypred<-length(N)
y.sd100 <- matrix(0, ncol=100, nrow=Nspp)
for (i in 1:100){
     for (n in 1:N){
         s <- species[n]
         ypred[n] <- a[s] + mu_b*year[n] # I think a[s] would also work 
     }
   y <- rnorm(N, ypred, sigma_y)
   y.df <- as.data.frame(cbind(y, species))
   y.sd <- aggregate(y.df["y"], y.df["species"], FUN=sd)
   y.sd100[,i] <- y.sd[,2] 
}
hist(rowMeans(y.sd100))

## looking at slopes
a100 <- matrix(0, ncol=100, nrow=Nspp)
for (j in 1:100){
     a100[,j] <- rnorm(Nspp, mean=mu_a, sd=sigma_a)
}

hist(rowMeans(a100), col=rgb(1,0,0,0.5), xlim=range(dr$intfits), main="Mean(a) from 100 random draws based on Stan model",  ylab="mean(a from one draw)") 
abline(v = mean(dr$intfits), col = "blue", lwd = 2) # mean of lm fits
hist(dr$intfits, col=rgb(0,0,1,0.5), add=TRUE)

#100 random draws instead of means
hist(a100, col=rgb(1,0,0,0.5), main="100 random draws based on Stan model",  ylab="") 
abline(v = mean(dr$intfits), col = "blue", lwd = 2) # mean of lm fits
hist(dr$intfits, col=rgb(0,0,1,0.5), xlim=range(a100), add=TRUE)

