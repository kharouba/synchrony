# Posterior Predictive checks for Synchrony model

# get estimates from no pooling (fixed effects) and complete pooling
# Complete pooling:
comppool<-lm(phenovalue~year, data=rawlong.nodups)
comppool.wtype<-lm(phenovalue~year*spp-1, data=rawlong.nodups)
#nopool<-lm(phenovalue~year*species, data=rawlong.nodups)
# No pooling:
#rawlong.tot$yr1981 <- rawlong.tot$newyear-1981
uniquespp<-unique(rawlong.tot2$species)
#uniquespp<-unique(rawlong.tot$newid)
slopefits<-rep(NA< length(uniquespp))
varfits<-rep(NA< length(uniquespp))
intfits<-rep(NA< length(uniquespp))
for(eachsp in 1:length(uniquespp)){
	lmhere<-lm(phenovalue~yr1981, data=subset(rawlong.tot2, species==uniquespp[eachsp]))
	slopefits[eachsp]<-coef(lmhere)[2]
	varfits[eachsp]<-(summary(lmhere)$sigma)**2
	intfits[eachsp]<-coef(lmhere)[1]
}
dr<-data.frame(cbind(uniquespp, slopefits, varfits, intfits))
dr$slopefits<-as.numeric(dr$slopefits)
dr$intfits<-as.numeric(dr$intfits)
dr$varfits<-as.numeric(dr$varfits)
#write.csv(dr, "lmfits.csv")

#First, plot the real data used in the model
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=c(0, 350), type="n", xlab="Year", ylab="Phenovalue", bty="l", main="Raw phenology data")
for (j in 1:Nspp){
  lines(year[species==j], y[species==j])
}


# Now grab the stan output 
goo <- extract(sync.model) 
#hist(goo$mu_b) # example!

#extract the means for now:
#mu_a <- mean(goo$mu_a)
mu_b <- mean(goo$mu_b) # from stan output 
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
sigma_b <- mean(goo$sigma_b) # from stan output
#sigma_a <- mean(goo$sigma_a) # from stan output 

# Random slopes only model:
a <- colMeans(goo$a) # this may not be ideal but deals with the non-pooling
#a <- rnorm(Nspp, mean=mu_a, sd=sigma_a)
b <- rnorm(Nspp, mean=mu_b, sd=sigma_b) # this is one way to create fake data from the Stan output to use in the PP check


## Create the PP data using new a nd b for each of 71 species
#This is just creating basically one new set of data, or one predictive check
# RANDOM SLOPE MODEL:
ypred <- length(N) # Lizzie added
for (n in 1:N){
     s <- species[n]
    ypred[n] <- a[s] + b[s]*year[n] # one way to handle the non-pooled intercepts, there may be other ways
}

y <- rnorm(N, ypred, sigma_y)


#Plot the data and see what the raw data predicted from the model looks like
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=c(0, 350), type="n", xlab="Year", ylab="Phenovalue",   bty="l", main="Data from posterior means")
for (j in 1:Nspp){
   lines(year[species==j], y[species==j])
}

 # compare a few things on this single new dataset
par(mfrow=c(1,2))
hist(b, main="slopes (b) from the stan model with mean from the raw data in blue")
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # less negative, slopes are generally pooled towards center which makes sense
hist(dr$varfits, main="sigma y (b) from the raw data with sigma_y from the model in blue")
abline(v=sigma_y, col = "blue", lwd = 2) 

# okay, but predictive checks are about much more than ONE simulation, one draw ...
# Create the data using new a and b for each of 71 species, 100 times
# RANDOM slope MODEL
ypred<-length(N)
y.sd100 <- matrix(0, ncol=100, nrow=Nspp)
for (i in 1:100){
     for (n in 1:N){
         s <- species[n]
         ypred[n] <- a[s] + b[s]*year[n] # I think a[s] would also work 
     }
   y <- rnorm(N, ypred, sigma_y)
   y.df <- as.data.frame(cbind(y, species))
   y.sd <- aggregate(y.df["y"], y.df["species"], FUN=sd)
   y.sd100[,i] <- y.sd[,2] 
}
hist(rowMeans(y.sd100))

# and here's the real data
real.sd <- aggregate(rawlong.tot["phenovalue"], rawlong.tot[c("studyid", "species")],  FUN=sd)
hist(rowMeans(y.sd100), col=rgb(1,0,0,0.5), main="Mean(SD of y) from 100 simulated datasets based on Stan model", ylab="mean(SD of one sim)")
abline(v = mean(real.sd$phenovalue, na.rm=TRUE), col = "blue", lwd = 2)

## looking at slopes
b100 <- matrix(0, ncol=100, nrow=Nspp)
for (j in 1:100){
     b100[,j] <- rnorm(Nspp, mean=mu_b, sd=sigma_b)
}

#RED= STAN MODEL
#BLUE= RAW DATA
hist(rowMeans(b100), col=rgb(1,0,0,0.5), xlim=range(dr$slopefits), main="Mean(b) from 100 random draws based on Stan model",  ylab="mean(b from one draw)") 
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # mean of lm fits
hist(dr$slopefits, col=rgb(0,0,1,0.5), add=TRUE)

#100 random draws instead of means
hist(b100, col=rgb(1,0,0,0.5), xlim=c(-2.7, 2.8), main="100 random draws based on Stan model",  ylab="") 
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # mean of lm fits
hist(dr$slopefits, col=rgb(0,0,1,0.5), xlim=c(-2.7, 2.8), add=TRUE)

#1 random draws instead of means
hist(b100[,2], col=rgb(1,0,0,0.5), xlim=c(-2.7, 2.8), main="1 random draw based on Stan model",  ylab="") 
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # mean of lm fits
hist(dr$slopefits, col=rgb(0,0,1,0.5), xlim=c(-2.7, 2.8), add=TRUE)


summ_studyspp <- subset(specieschar.formodel, select=c("studyid", "species")); summ_studyspp<-unique(summ_studyspp)
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- goo$b[i,]
    it1000[,(i-3000)] <- goo$b[i,]
}
summ_studyspp$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP

names(dr)[1]<-"species"
asdf<-unique(merge(dr[,c("species","slopefits","varfits")], summ_studyspp))
asdf2<-merge(asdf, specieschar.hin)

par(mfrow=c(1,2)); with(asdf2, plot(slopefits~phenovalue, ylim=c(-3, 3),xlab="number of years per spp", main="Slope coefficients from (a) simple lm (b) stan models")); abline(h=0); with(asdf2, plot(stanfit~phenovalue, ylim=c(-3,3) , xlab="number of years per spp")); abline(h=0); # adding a horizontal line which is mean from complete pooling doesn't make sense since temperatures are in really different ranges (water vs air etc.), 

#Correlation across species
with(asdf2, plot(slopefits~stanfit)); with(asdf2, abline(lm(slopefits~stanfit)))

#Relationship with sd
names(real.sd)[3]<-"sd"
asdf3<-merge(asdf2, real.sd, by=c("studyid","species"))
par(mfrow=c(1,2)); with(asdf3, plot(slopefits~sd, ylim=c(-5.2, 10.75))); abline(h=0); with(asdf3, plot(stanfit~sd, ylim=c(-5.2, 10.75))); abline(h=0);

