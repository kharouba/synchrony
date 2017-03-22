# Posterior Predictive checks for covariate model (e.g. pheno change~ temp change)

# skipping ahead for a sec- covariate model- varying and pooled intercepts, only single slope estimated
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



## First for TEMP MODEL:
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=range(y), type="n", xlab="Year", ylab="Phenovalue", bty="l", main="Raw phenology data")
for (j in 1:Nspp){
  lines(year[species==j], y[species==j])
}

# Now grab the stan output 
goo <- extract(temp.model) 

#extract the means for now:
mu_b <- mean(goo$mu_b) # from stan output 
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
sigma_b <- mean(goo$sigma_b) # from stan output
mu_a <- mean(goo$mu_a) # from stan output 
sigma_a <- mean(goo$sigma_a) # from stan output 

# Random slopes only model:
#a <- colMeans(goo$a) # this may not be ideal but deals with the non-pooling
a <- rnorm(Nspp, mean=mu_a, sd=sigma_a) # this is one way to create fake data from the Stan output to use in the PP check
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
plot(range(year), range(y), ylim=range(y), type="n", xlab="Year", ylab="Phenovalue",   bty="l", main="Data from posterior means")
for (j in 1:Nspp){
   lines(year[species==j], y[species==j])
}


### Now for two level cov model
#First, plot the real data used in the model
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=c(range(y)), type="n", xlab="Year", ylab="Phenovalue", bty="l", main="Raw phenology data")
for (j in 1:Nspp){
  lines(year[species==j], y[species==j])
}

# Now grab the stan output 
goo <- extract(cov.model) 

#extract the means for now:
mu_a <- mean(goo$mu_a)
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
sigma_a <- mean(goo$sigma_a) # from stan output

# Random slopes only model:
a <- rnorm(Nspp, mean=mu_a, sd=sigma_a)
b <- colMeans(goo$b)# this is one way to create fake data from the Stan output to use in the PP check

ypred <- length(N) # Lizzie added
for (n in 1:N){
     s <- species[n]
    ypred[n] <- a[s] + b[s]*year[n] # one way to handle the non-pooled intercepts, there may be other ways
}

y <- rnorm(N, ypred, sigma_y)

#Plot the data and see what the raw data predicted from the model looks like
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=c(0, 4), type="n", xlab="Year", ylab="Phenovalue",   bty="l", main="Data from posterior means")
for (j in 1:Nspp){
   lines(year[species==j], y[species==j])
}

#### First for threelevelrandomslope.stan (for first model e.g. temp change)
# get estimates from no pooling (fixed effects) 
# No pooling:
uniquespp<-unique(clim3$species)
slopefits<-rep(NA< length(uniquespp))
varfits<-rep(NA< length(uniquespp))
intfits<-rep(NA< length(uniquespp))
for(eachsp in 1:length(uniquespp)){
	#lmhere<-lm(envvalue~yr1981, data=subset(clim3, species==uniquespp[eachsp])) #for temp change
	lmhere<-lm(phenovalue~yr1981, data=subset(clim3, species==uniquespp[eachsp])) #for pheno change
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
plot(range(year), range(y), ylim=c(range(y)), type="n", xlab="Year", ylab="Phenovalue", bty="l", main="Raw phenology data")
for (j in 1:Nspp){
  lines(year[species==j], y[species==j])
}


# Now grab the stan output 
goo <- extract(temp.model) 
goo <- extract(pheno.model) 

#extract the means for now:
mu_b <- mean(goo$mu_b) # from stan output 
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)

sigma_b <- mean(goo$sigma_b) # from stan output

sigma_a_study <- mean(goo$sigma_a_study) # from stan output
sigma_b_study <- mean(goo$sigma_b_study) # from stan output


#Top level
#a_study<-colMeans(goo$a_study) # this may not be ideal but deals with the non-pooling
a_study<-rnorm(Nstudy, mean=colMeans(goo$a_study), sd=sigma_a_study) #ASK LIZZIE IF THIS IS OK
b_study<-rnorm(Nstudy, mean=mu_b, sd=sigma_b) # this is one way to create fake data from the Stan output to use in the PP check


#2nd level
a_spp <- rep(0,Nspp)
b_spp <- rep(0,Nspp)
for (j in 1:Nspp){
  a_spp[j] <- rnorm(1,a_study[studyid[j]], sigma_a_study[1]);
  b_spp[j] <- rnorm(1,b_study[studyid[j]], sigma_b_study[1]);
}


## Create the PP data using new a nd b for each of 71 species
#This is just creating basically one new set of data, or one predictive check
# RANDOM SLOPE MODEL:
y <- rep(0,N)
ypred <- rep(0,N)
for (n in 1:N){
  ypred[n] <- a_spp[species[n]] + b_spp[species[n]]*year[n]
  y[n] <- a_spp[species[n]] + b_spp[species[n]]*year[n] + rnorm(1,0,sigma_y)
}


#Plot the data and see what the raw data predicted from the model looks like
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=c(range(y)), type="n", xlab="Year", ylab="Phenovalue",   bty="l", main="Data from posterior means")
for (j in 1:Nspp){
   lines(year[species==j], y[species==j])
}

 # compare a few things on this single new dataset
par(mfrow=c(1,2))
hist(b_spp, main="slopes (b) from the stan model with mean from the raw data in blue")
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # less negative, slopes are generally pooled towards center which makes sense
hist(dr$varfits, main="sigma y (b) from the raw data with sigma_y from the model in blue")
abline(v=sigma_y, col = "blue", lwd = 2) 

# okay, but predictive checks are about much more than ONE simulation, one draw ...
# Create the data using new a and b for each of 71 species, 100 times
# RANDOM slope MODEL
ypred <- rep(0,N)
y.sd100 <- matrix(0, ncol=100, nrow=Nspp)
for (i in 1:100){
     for (n in 1:N){
  ypred[n] <- a_spp[species[n]] + b_spp[species[n]]*year[n]
  y[n] <- a_spp[species[n]] + b_spp[species[n]]*year[n] + rnorm(1,0,sigma_y)
     }
   y.df <- as.data.frame(cbind(y, species))
   y.sd <- aggregate(y.df["y"], y.df["species"], FUN=sd)
   y.sd100[,i] <- y.sd[,2] 
}
hist(rowMeans(y.sd100))

# and here's the real data CHANGE VARIABLE NAMES!
real.sd <- aggregate(clim3["envvalue"], clim3[c("studyid", "species")],  FUN=sd)
hist(rowMeans(y.sd100), col=rgb(1,0,0,0.5), main="Mean(SD of y) from 100 simulated datasets based on Stan model", ylab="mean(SD of one sim)")
abline(v = mean(real.sd$envvalue, na.rm=TRUE), col = "blue", lwd = 2)

## looking at slopes
b100 <- matrix(0, ncol=100, nrow=Nspp)
for (j in 1:100){
     b100[,j] <- rnorm(Nspp, mean=b_study, sd=sigma_b_study) ##generate slopes for each species based on stan model
}
 a_spp[j] <- rnorm(1,a_study[studyid[j]], sigma_a_study[1]);

#RED= STAN MODEL
#BLUE= RAW DATA
hist(rowMeans(b100), col=rgb(1,0,0,0.5), xlim=range(dr$slopefits), main="Mean(b) from 100 random draws based on Stan model",  ylab="mean(b from one draw)") 
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # mean of lm fits
hist(dr$slopefits, col=rgb(0,0,1,0.5), add=TRUE)

#100 random draws instead of means
hist(b100, col=rgb(1,0,0,0.5), xlim=range(dr$slopefits), main="100 random draws based on Stan model",  ylab="") 
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # mean of lm fits
hist(dr$slopefits, col=rgb(0,0,1,0.5), xlim=c(-2.7, 2.8), add=TRUE)

#1 random draws instead of means
hist(b100[,2], col=rgb(1,0,0,0.5), xlim=range(dr$slopefits), main="1 random draw based on Stan model",  ylab="") 
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # mean of lm fits
hist(dr$slopefits, col=rgb(0,0,1,0.5), xlim=c(-2.7, 2.8), add=TRUE)


sock<-unique(clim3[,c("studyid","species")])
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    sock$model <- goo$b_spp[i,]
    it1000[,(i-3000)] <- goo$b_spp[i,]
}
sock$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP

names(dr)[1]<-"species"
asdf<-unique(merge(dr[,c("species","slopefits","varfits")], sock))
#asdf2<-merge(asdf, specieschar.hin)


#Correlation across species
with(asdf, plot(slopefits~stanfit)); with(asdf, abline(lm(slopefits~stanfit)))

#extra
par(mfrow=c(1,2)); with(asdf2, plot(slopefits~phenovalue, ylim=c(-3, 3),xlab="number of years per spp", main="Slope coefficients from (a) simple lm (b) stan models")); abline(h=0); with(asdf2, plot(stanfit~phenovalue, ylim=c(-3,3) , xlab="number of years per spp")); abline(h=0); # adding a horizontal line which is mean from complete pooling doesn't make sense since temperatures are in really different ranges (water vs air etc.), 


#Relationship with sd
names(real.sd)[3]<-"sd"
asdf3<-merge(asdf2, real.sd, by=c("studyid","species"))
par(mfrow=c(1,2)); with(asdf3, plot(slopefits~sd, ylim=c(-5.2, 10.75))); abline(h=0); with(asdf3, plot(stanfit~sd, ylim=c(-5.2, 10.75))); abline(h=0);


----------------------------------------------------
### 	SECOND- FOR THREE LEVEL RANDOM EFFECTS (e.g. pheno change~temp change)

N<-nrow(tdata2)
y <- abs(tdata2$pheno.change) #absolute value of pheno change!
Nspp <- length(unique(tdata2$id)) #J
species <- as.numeric(as.factor(tdata2$id))
#studyid <- as.numeric(as.factor(tdata2$studyid))
year <- abs(tdata2$temp.change) #absolute value of temp change
Nstudy <- length(unique(tdata2$studyid)) #J
sock<-unique(tdata2[,c("studyid","species")])
studyid <- as.numeric(as.factor(sock$studyid))


# get estimates from no pooling (fixed effects) 
# No pooling:
uniquespp<-unique(tdata2$species)
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

#First, plot the real data used in the model
#sub<-subset(tdata2, iteration<=10)
#y <- abs(sub$pheno.change) #absolute value of pheno change!
#species <- as.numeric(as.factor(sub$id))
#studyid <- as.numeric(as.factor(sub$studyid))
#year <- abs(sub$temp.change) #absolute value of temp change

par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(x=range(year), y=range(y), ylim=c(range(y)), type="n", xlab="abs(Temp change)", ylab="abs(Pheno change)", bty="l", main="Raw data")
for (j in 1:Nspp){
  abline(lm(y[species==j]~year[species==j]))
}

ggplot(tdata2, aes(x=abs(temp.change), y=abs(pheno.change), colour=factor(species)))+geom_point()+geom_smooth(method="lm", se=FALSE)+theme(legend.position="none")
ggplot(tdata2, aes(x=abs(temp.change), y=abs(pheno.change)))+geom_smooth(method="lm", se=FALSE, aes(colour=factor(species)))+geom_smooth(method="lm", se=FALSE)+theme(legend.position="none")

# Now grab the stan output 
goo <- extract(cov.model) 

#extract the means for now:
mu_a <- mean(goo$mu_a) #global intercept
mu_b <- mean(goo$mu_b) # global slope 
sigma_y <- mean(goo$sigma_y) # obervation error
sigma_a <- mean(goo$sigma_a) # sd around global intercept
sigma_b <- mean(goo$sigma_b) # sd around global slope

sigma_a_study <- mean(goo$sigma_a_study) # sd across species intercepts within each study
sigma_b_study <- mean(goo$sigma_b_study) # sd across species slopes within each study


#Top level
#a_study<-colMeans(goo$a_study) # this may not be ideal but deals with the non-pooling
a_study <- rnorm(Nstudy, mean=mu_a, sd=sigma_a)
b_study <- rnorm(Nstudy, mean=mu_b, sd=sigma_b) # this is one way to create fake data from the Stan output to use in the PP check

# draw S=10 site sds for slope and intercept (variability of plots within sites)
#sig_a_study <- rnorm(Nstudy,mean=mean(a_study),sd=sigma_a_study)
#sig_b_study <- rnorm(Nstudy,mean=mean(b_study),sd=sigma_b_study)
#sig_a_study <- rnorm(Nstudy,0,1)
#sig_b_study <- rnorm(Nstudy,0,1)
#sig_a_study <- rlnorm(Nstudy,0,sd=sigma_a_study) #need log norm so that sd aren't negative
#sig_b_study <- rlnorm(Nstudy,0,sd=sigma_b_study) #need log norm so that sd aren't negative
sig_a_study <- sigma_a_study
sig_b_study <- sigma_b_study



#2nd level
a_spp <- rep(0,Nspp)
b_spp <- rep(0,Nspp)
for (j in 1:Nspp){
  a_spp[j] <- rnorm(1,a_study[studyid[j]], sig_a_study[1]);
  b_spp[j] <- rnorm(1,b_study[studyid[j]], sig_b_study[1]);
}

#This is just creating basically one new set of data, or one predictive check
# RANDOM EFFECTS MODEL:
year <- abs(tdata2$temp.change) #absolute v
species <- as.numeric(as.factor(tdata2$id))
y <- rep(0,N)
ypred <- rep(0,N)
for (n in 1:N){
  ypred[n] <- a_spp[species[n]] + b_spp[species[n]]*year[n]
  y[n] <- a_spp[species[n]] + b_spp[species[n]]*year[n] + rnorm(1,0,sigma_y)
}



#Plot the data and see what the raw data predicted from the model looks like
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), ylim=c(range(y)), type="n", xlab="abs(temp change)", ylab="abs(pheno change)",   bty="l", main="Data from posterior means")
for (j in 1:Nspp){
   abline(lm(y[species==j]~year[species==j]))
}

sock<-unique(tdata2[,c("studyid","species")])
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    sock$model <- goo$b_spp[i,]
    it1000[,(i-3000)] <- goo$b_spp[i,]
}
sock$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP

names(dr)[1]<-"species"
asdf<-unique(merge(dr[,c("species","slopefits","varfits")], sock))
#asdf2<-

#Correlation across species
with(asdf, plot(slopefits~stanfit)); with(asdf, abline(lm(slopefits~stanfit)));
with(asdf, plot(slopefits~stanfit, xlim=range(slopefits), ylim=range(slopefits))); with(asdf, abline(lm(slopefits~stanfit))); abline(b=1)

ggplot(asdf, aes(stanfit, slopefits))+geom_point()+geom_smooth(method="lm", se=FALSE)+geom_text(aes(label=species), hjust=0, vjust=0)


#pooling based on sd?
asdf$fitdiff<-with(asdf, abs(slopefits)-abs(stanfit))
pheno.sd <- aggregate(tdata2["pheno.change"], tdata2["species"], FUN=sd); names(pheno.sd)[2]<-"pheno.sd"
asdf2<-merge(asdf, pheno.sd, by=c("species"))

ggplot(asdf2, aes(x=pheno.sd, y=fitdiff))+geom_point()+geom_smooth(method="lm", se=FALSE)+geom_text(aes(label=species), hjust=0, vjust=0) #+theme(main="Stan fits based ")



