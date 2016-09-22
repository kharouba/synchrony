### Posterior predictive checks for temperature model##
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
###
# get estimates from no pooling (fixed effects) and complete pooling
uniquespp<-unique(dataset$datasetid)
slopefits<-rep(NA< length(uniquespp))
varfits<-rep(NA< length(uniquespp))
intfits<-rep(NA< length(uniquespp))
#clim3$z<-with(clim3, (envvalue-mean(envvalue))/sd(envvalue)) # normal standardization
# No pooling:
for(eachsp in 1:length(uniquespp)){
	#lmhere<-lm(envvalue~newyear, data=subset(clim3, species==uniquespp[eachsp]))
	lmhere<-lm(envvalue~yr1981, data=subset(dataset, datasetid==uniquespp[eachsp]))
	slopefits[eachsp]<-coef(lmhere)[2]
	varfits[eachsp]<-(summary(lmhere)$sigma)**2
	intfits[eachsp]<-coef(lmhere)[1]
}
dr<-data.frame(cbind(uniquespp, slopefits, varfits, intfits))
#dr$slopefits<-as.numeric(dr$slopefits)

#write.csv(dr, "lmfits_temp.csv")

temptrends <- read.csv("lmfits_temp.csv", header=TRUE)
#temptrends <- read.csv("lmfits_temp_center.csv", header=TRUE)

# load temp.model from tempmodels.R

#First, plot the real data used in the model
#y2<-subset(y, ) HAS 99999 IN ENVALUE
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Temperature", bty="l", main="Raw real temp data")
for (j in 1:J){
  lines(year[species==j], y[species==j])
}
#hist(y, xlab="Temperature (C)", main="Real temp data, cleaned (I assume)")

# Now grab the stan output 
goo <- extract(temp.model) 
#hist(goo$mu_b) # example!

#extract the means for now:
#mu_a <- mean(goo$mu_a)
mu_b <- mean(goo$mu_b) # from stan output 
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
#sigma_a <- mean(goo$sigma_a) # from stan output
sigma_b <- mean(goo$sigma_b) # from stan output 

# If random effects model:
#a <- rnorm(J, mean=mu_a, sd=sigma_a) # one way to create fake data from the Stan output to use in the PP check
#b <- rnorm(J, mean=mu_b, sd=sigma_b) # one way to create fake data from the Stan output to use in the PP check

# If random SLOPE model:
a <- colMeans(goo$a) # this may not be ideal but deals with the non-pooling
b <- rnorm(Nspp, mean=mu_b, sd=sigma_b) # one way to create fake data from the Stan output to use in the PP check

## Create the PP data using new a nd b for each of 71 species
#This is just creating basically one new set of data, or one predictive check
#
ypred <- length(N) # Lizzie added
for (n in 1:N){
     s <- species[n]
    ypred[n] <- a[s] + b[s]*year[n] # one way to handle the non-pooled intercepts, there may be other ways
}

y <- rnorm(N, ypred, sigma_y)

#Plot the data and see what the raw data predicted from the model looks like

#pdf("onepredictivecheck.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Temperature",   bty="l", main="Data from posterior means")
for (j in 1:J){
   lines(year[species==j], y[species==j])
   }

hist(y, xlab="Temperature", main="Data from posterior means") 
 dev.off()
 # Hmm, the model makes more data than we seem to have ... but I assume that is because of species rep issues? YES BECAUSE OF REPEATING SPECIES IN RAW DATA, IN STAN WHITE NOISE ADDED AND NOW SLOPES ESTIMATED SO NOW LINES ARE SEPARATED, RAW DATA SLOPES SIT ON EACH OTHER SO LOOKS LIKE LESS DATA

 
 # compare a few things on this single new dataset
par(mfrow=c(2,2))

#for random effects model look at (a)
hist(a, main="intercepts (a) from the stan model with mean from the raw data in blue")
abline(v = mean(temptrends$intfits), col = "blue", lwd = 2) 
#NOW RESOLVED- X WAS YEAR AND SHOULD BE YR1981 # not even on the graph! ** Are the YEARS centered differently between the model and the raw data?

hist(b, main="slopes (b) from the stan model with mean from the raw data in blue")
abline(v = mean(temptrends$slopefits), col = "blue", lwd = 2) # less negative, slopes are generally pooled towards center which makes sense
hist(temptrends$varfits, main="sigma y (b) from the raw data with sigma_y from the model in blue")
abline(v=sigma_y, col = "blue", lwd = 2) 

# okay, but predictive checks are about much more than ONE simulation, one draw ...
# Create the data using new a and b for each of 71 species, 100 times
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
real.sd <- aggregate(dataset["envvalue"], dataset[c("studyid", "species")],  FUN=sd)
hist(rowMeans(y.sd100), col=rgb(1,0,0,0.5), main="Mean(SD of y) from 100 simulated datasets based on Stan model", ylab="mean(SD of one sim)")
abline(v = mean(real.sd$envvalue, na.rm=TRUE), col = "blue", lwd = 2)

## looking at slopes
b100 <- matrix(0, ncol=100, nrow=J)
for (j in 1:100){
     b100[,j] <- rnorm(J, mean=mu_b, sd=sigma_b) ##generate slopes for each species based on stan model
}

#RED= STAN MODEL
#BLUE= RAW DATA
CHANGE XLIM EACH TIME
# Based on mean of random draws
hist(temptrends$slopefits, col=rgb(0,0,1,0.5), main="Mean(b) from 100 random draws based on Stan model",  ylab="mean(b from one draw)") #add=TRUE
hist(rowMeans(b100), col=rgb(1,0,0,0.5), xlim=range(temptrends$slopefits), add=TRUE) 
abline(v = mean(temptrends$slopefits), col = "blue", lwd = 2) # mean of lm fits
# Some very strong pooling of the slopes it looks like, both the negative fits and the positive high fits
# Is this what we want? Might be good to talk to Stan group about this

#Based on 1 random draws 
hist(temptrends$slopefits, col=rgb(0,0,1,0.5), main="Single random draw based on Stan model",  ylab="") #add=TRUE
hist(b100[,1], col=rgb(1,0,0,0.5), xlim=range(temptrends$slopefits), add=TRUE) 
abline(v = mean(temptrends$slopefits), col = "blue", lwd = 2) # mean of lm fits
# Some very strong pooling of the slopes it looks like, both the nega

#Mean of random draws vs. all random draws
hist(b100, col=rgb(1,0,0,0.5), xlim=range(temptrends$slopefits), main="100 random draws based on Stan model",  ylab="")
hist(temptrends$slopefits, col=rgb(0,0,1,0.5), xlim=range(temptrends$slopefits), add=TRUE) 
abline(v = mean(temptrends$slopefits), col = "blue", lwd = 2)

# Check that the pooling appears correct (e.g., are the shortest time series being pooled the most?)
specieschar.formodel <- aggregate(dataset["envvalue"], dataset[c("datasetid")], FUN=length) #number of years
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("datasetid"))
#look to see how stan is pooling, look at correlation between stan slope and lm slope
temp.change <- matrix(0, ncol=3000, nrow=Nspp) #2000 iterations for 21 interactions;
for (i in 3000:6000){ # 2000 iterations?
    specieschar.formodel.sm$model <- 5.1
    specieschar.formodel.sm$model<-goo$b[i,]
    temp.change[,(i-3000)] <- goo$b[i,]
    }
specieschar.formodel.sm$stanfit <- rowMeans(temp.change, na.rm=TRUE) #mean across iterations for EACH SPP; STAN SLOPE ESTIMATE PER SPECIES
names(temptrends)[2]<-"datasetid"
asdf<-unique(merge(temptrends[,c("datasetid","slopefits","varfits")], specieschar.formodel.sm))
asdf2<-merge(asdf, specieschar.formodel)

m1<-lm(envvalue~yr1981, data=dataset) #complete pooling

#to match fig 12.3 in Gelman and  Hill; No pooling (a) vs. partial pooling or multilevel model (b); Slope coefficients vs. number of obs:
par(mfrow=c(1,2)); with(asdf2, plot(slopefits~envvalue, ylim=range(temptrends$slopefits), xlab="number of years per spp", main="Slope coefficients from (a) simple lm (b) stan models")); abline(h=0); with(asdf2, plot(stanfit~envvalue, ylim=range(temptrends$slopefits), xlab="number of years per spp")); abline(h=0); # adding a horizontal line which is mean from complete pooling doesn't make sense since temperatures are in really different ranges (water vs air etc.), 

#Correlation across species
with(asdf2, plot(slopefits~stanfit)); with(asdf2, abline(lm(slopefits~stanfit)))

#look at those species with strong positive trend (warming) and those with negative trend (cooling)
high<-unique(subset(temptrends, slopefits>=0.08)); names(high)[2]<-"species"
new<-merge(high[,c("species","slopefits")], dataset, by="species")

low<-unique(subset(temptrends, slopefits<=0)); names(low)[2]<-"species"
new<-merge(low[,c("species","slopefits")], dataset, by="species")

ggplot(new, aes(x=yr1981, y=envvalue))+facet_wrap(~species, nrow=5)+  geom_point(size=1)+theme_bw()+geom_smooth(method="lm", se=FALSE)


#Relationship with sd
names(real.sd)[3]<-"sd"
asdf3<-merge(asdf2, real.sd, by=c("studyid","species"))
par(mfrow=c(1,2)); with(asdf3, plot(abs(slopefits)~sd, ylim=range(abs(temptrends$slopefits)))); abline(h=0); with(asdf3, plot(stanfit~sd, ylim=range(abs(temptrends$slopefits)))); abline(h=0);



###
### end posterior predictive checks
###

