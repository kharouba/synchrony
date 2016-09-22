### Started 1 April 2015 ###
### April Fool's Day! ###

## Trying to run STAN with synchrony data ##

rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(ggplot2)
library(rstan)
library(shinyStan)
library(grid)
library(nlme)
library(dplyr)
library(reshape)
set_cppo("fast")  # for best running speed
source("/users/kharouba/google drive/UBC/multiplot.R")
#library(reshape)
# library(lme4)

#get data WITH HINGE
source("tempmodels_clean.R")
source("syncmodels_clean.R")

rawlong.tot <- arrange(rawlong.tot, studyid, species)
clim3 <- arrange(clim3, studyid, species)
newdata<-merge(rawlong.tot, unique(clim3[,c("studyid","species")]), by=c("studyid","species"))
newdata.formodel<-merge(specieschar.formodel, unique(clim3[,c("studyid","species")]), by=c("studyid","species"))

#For predictive checks
uniquespp<-unique(newdata.formodel$species)
lmfits<-rep(NA< length(uniquespp))
varfits<-rep(NA< length(uniquespp))
for(eachsp in 1:length(uniquespp)){
	lmhere<-lm(phenovalue~year, data=subset(rawlong.nodups, species==uniquespp[eachsp]))
	lmfits[eachsp]<-coef(lmhere)[2]
	varfits[eachsp]<-(summary(lmhere)$sigma)**2
}
dr<-data.frame(cbind(uniquespp, lmfits, varfits))



############
#RUN synch models first
newdata$yr1981 <- newdata$newyear-1981
N <- nrow(newdata)
y <- newdata$phenovalue
specieschar.hin<- aggregate(newdata["phenovalue"], newdata[c("studyid", "species", "int_type", "terrestrial", "spp")], FUN=length) #number of years per species
J <- nrow(specieschar.hin)
species <- as.numeric(as.factor(newdata$species))
year <- newdata$yr1981
nVars <-1

#add warmin?
sync.model<-stan("synchrony_apr_nocov.stan", data=c("N","J","y","species","year","nVars"), iter=2000,chains=4)
#sync.model<-stan("synchrony_apr_nocov.stan", data=c("N","J","y","species","year","nVars"), iter=3000, warmup=1000, thin=10, chains=4)

fh.sim <- extract(sync.model, inc_warmup=FALSE) #only extracts parameters
dim(fh.sim$b) # number of iterations*number of species 
names(fh.sim)
fh.sim$b[2000,]

specieschar.formodel.sm <- subset(newdata.formodel, select=c("studyid", "species"))
intid <- read.csv("input/raw_april.csv", header=TRUE)
lal<-unique(newdata[,c("intid","terrestrial")])
intid2<-merge(intid, lal, by=c("intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction","terrestrial"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]


it1000 <- matrix(0, ncol=2000, nrow=22) #2000 iterations for 53 interactions;
for (i in 2000:4000){ # 2000 iterations?
    specieschar.formodel.sm$model <- fh.sim$b[i,]
    andtheanswer <- merge(intid.nodups, specieschar.formodel.sm, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, specieschar.formodel.sm, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
    it1000[,(i-2000)] <- andtheanswer$model.x-andtheanswer$model.y #model.x=spp1
}
syncs<-it1000

###############
# Temp model next
clim3$yr1981 <- clim3$newyear-1981
N <- nrow(clim3)
y <- clim3$envvalue
specieschar.hin<- aggregate(clim3["envvalue"], clim3[c("studyid", "species")], FUN=length) #number of years per species
J <- nrow(specieschar.hin)
species <- as.numeric(as.factor(clim3$species))
year <- clim3$yr1981
nVars <-1
Imat <- diag(1, nVars)


#New model as of April 2016
#Random slopes only, no random intercepts, hinge, no covariate matrix:
temp.model<-stan("synchrony_apr_nocov.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000,  chains=4)
#temp.model<-stan("synchrony_apr_nocov.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, warmup, thin=10, chains=4)

te.sim <- extract(temp.model, inc_warmup=FALSE) #only extracts parameters
dim(te.sim$b) # number of iterations*number of species 
names(te.sim)
te.sim$b[2000,]

specieschar.formodel <- aggregate(clim3["envvalue"], clim3[c("studyid", "species")], FUN=length) #number of years
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species"))
rawlong <- read.csv("rawlong.csv", header=TRUE)
sol<-merge(specieschar.formodel.sm, unique(rawlong[,c("studyid","intid","species")]), by=c("studyid","species"))
intid <- read.csv("input/raw_april.csv", header=TRUE)
intid2<-merge(intid, unique(sol[,c("studyid","intid")]), by=c("studyid","intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
temp_int<-intid.nodups #temp change interactions


temp.change <- matrix(0, ncol=2000, nrow=22) #2000 iterations for 21 interactions;
for (i in 2000:4000){ # 2000 iterations?
    specieschar.formodel.sm$model <- 5.1
    specieschar.formodel.sm$model<-te.sim$b[i,]
    andtheanswer <- merge(temp_int, specieschar.formodel.sm, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer2 <- merge(andtheanswer, specieschar.formodel.sm, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
        andtheanswer2[2,"model.x"]<- andtheanswer2[2,"model.y"] #Engraulis japonicus because first day of temp X
        andtheanswer2[9,"model.x"]<-andtheanswer2[9,"model.y"] #Diatom2 spp. because first day of temp X
        andtheanswer2[10,"model.x"]<-andtheanswer2[10,"model.y"]#Diatom2 spp.] because first day of temp X
    andtheanswer2$meanmodel <- rowMeans(andtheanswer2[,c(6:7)])
    temp.change[,(i-2000)] <- andtheanswer2$meanmodel
}

#########################
### COVARIATES MODEL ###
news<-t(syncs)
mdata <- melt(news)
names(mdata)[1]<-"iteration"; names(mdata)[2]<-"intid"; names(mdata)[3]<-"sync_change"

sup<-t(temp.change)
ndata<-melt(sup)
names(ndata)[1]<-"iteration"; names(ndata)[2]<-"intid"; names(ndata)[3]<-"temp_change"

newframe<-merge(mdata, ndata, by=c("iteration","intid"))
andtheanswer2$newintid<-1:nrow(andtheanswer2)
names(andtheanswer2)[4]<-"oldintid"; names(andtheanswer2)[9]<-"intid"
newframe2<- merge(newframe, unique(andtheanswer2[,c("intid","oldintid")]))
names(newdata)[4]<-"oldintid"
some<-unique(newdata[,c("studyid","oldintid")]); names(some)[1]<-"oldstudyid"; some$row<-9999 
somes<-unique(some[,c("oldstudyid","row")]); somes$studyid<-1:nrow(somes)
names(andtheanswer2)[1]<-"oldstudyid"
somess<-merge(unique(andtheanswer2[,c("oldstudyid","oldintid","intid")]), somes)
newframe2<-merge(newframe, somess)


#data for stan model
Ni<- nrow(newframe2)
Nj<-length(unique(newframe2$intid))
Nk<-length(unique(newdata$studyid))

desMat <- model.matrix(object = ~ 1 + temp_change, data = newframe)
p<-ncol(desMat)
year<-newframe2$temp_change

intid <- as.integer(as.factor(newframe2$intid))
studyid <- as.integer(as.factor(newframe2$studyid))
studyLookup <- unique(newframe2[c("intid","studyid")])[,"studyid"] # Create a vector of school IDs where j-th element gives school ID for class ID j

y <- newframe$sync_change
#nVars<-1
#Imat <- diag(1, nVars)

# need covariance matrix?
	stan.model<-stan("threelevelrandomeffects.stan", data=c("Ni","Nj","Nk","intid","studyLookup","year","y"), iter=2000, warmup=1000, thin=10,chains=4)

# results
fh.sim <- extract(stan.model, inc_warmup=FALSE) #only extracts parameters
dim(fh.sim$b) # number of iterations*number of species 
names(fh.sim)

beta_1[2]=overall slope
u_1jk= slope for each interaction
beta_1jk= slope of temperature change at interaction level

st.err <- function(x) {
    sd(x)/sqrt(length(x))
     }

sums<-aggregate(newframe2["sync_change"], newframe2[c("intid")], FUN=mean);
sums_se<-aggregate(newframe2["sync_change"], newframe2[c("intid")], FUN=st.err)
names(sums_se)[2]<-"se_sync"
temps<-aggregate(newframe2["temp_change"], newframe2[c("intid")], FUN=mean);
temps_se<-aggregate(newframe2["temp_change"], newframe2[c("intid")], FUN=st.err)
names(temps_se)[2]<-"se_temp"
tot<-merge(sums, sums_se); tots<-merge(tot, temps); tots2<-merge(tots, temps_se)

ggplot(data=tots2, aes(x=temp_change, y=sync_change,
       colour=factor(intid))) +
       geom_point(size=3)+geom_errorbarh(aes(xmin=temp_change-se_temp, xmax=temp_change+se_temp), colour="black")+geom_errorbar(aes(ymin=sync_change-se_sync, ymax=sync_change+se_sync), colour="black")+geom_abline(intercept=-1.659, slope=-4.658)+ylim(-2,2)
        
       geom_smooth(method="lm", se=FALSE, size=1, 
       aes(fill = factor(species)))+theme_bw()