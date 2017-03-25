##### TO GET TEMP CHANGE PER INTERACTION ######
### note that as of 23 Sep 2016 this file uses the same models and data as tempmodels_daasets.R ###
### But the order that the data go into the models in each file is different... ###

rm(list=ls())
options(stringsAsFactors=FALSE)


#row.names=FALSE
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
# setwd("~/Documents/git/projects/trophsynch/synchrony/stan_2016")
library(ggplot2)
library(rstan)
library(shinyStan)
library(grid)
library(nlme)
library(gridExtra)
library(reshape)
library(plyr)
library(dplyr)
# set_cppo("fast") 

## HAS TEMP CHANGED?
## source the temp dataset cleaning code

source("source/clean_tempdatasets2.R")

# check, which datasetids have more than one species
ddply(dataset, c("datasetid"), summarise,
    nspp=n_distinct(species))

dataset <- dataset[with(dataset, order(species, year)),]
N<-nrow(dataset)
Nspp <- length(unique(dataset$species)) 
y <- dataset$envvalue
species <- as.numeric(as.factor(dataset$species))
year <- dataset$yr1981

temp.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(temp.model, pars = c("mu_b", "sigma_y", "a", "b"))

summ_studyspp <- subset(dataset, select=c("studyid","species", "datasetid")); 
summ_studyspp <- unique(summ_studyspp)
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]

goo <- extract(temp.model) 
summ_studyspp<-unique(dataset[,c("studyid","species")])
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- goo$b[i,]
   it1000[,(i-3000)] <- goo$b[i,]
}
tempchange<-it1000
ndata<-melt(it1000[,2001:3000])
names(ndata)[1]<-"id"; names(ndata)[2]<-"iteration"; names(ndata)[3]<-"temp.change"; 
#ndata$intid<-rep(intid.nodups$intid, 1000)


## Run sync models




specieschar.formodel <- aggregate(dataset["envvalue"], dataset[c("studyid","species",
    "datasetid")], FUN=length)
# number of years (which is correct, I double-checked below)
ddply(dataset, c("studyid", "species", "datasetid"), summarise,
    nspp=n_distinct(year))
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species", "datasetid"))





# Step 1 in the merge: Get the intid
rawlong <- read.csv("input/rawlong2.csv", header=TRUE)
# use for 1981 hinge model
source("input/datacleaningmore.R")
sol<-merge(specieschar.formodel.sm, unique(rawlong[,c("studyid","intid")]),
   by=c("studyid"), all.x=TRUE)
#sol<- sol[with(sol,order(species)),]

# Step 2 in the merge: Get the rest of the info
intid <- read.csv("input/raw_oct.csv", header=TRUE)
intid2 <- merge(intid, unique(sol[,c("studyid","intid", "datasetid")]), by=c("studyid","intid"),
    all.y=TRUE)
dim(unique(intid2[,c("studyid","intid", "datasetid")]))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction", "datasetid"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
temp_int <- intid.nodups #temp change interactions
temp_int <- temp_int[with(temp_int,order(intid)),]

unique(temp_int[,c("datasetid")])

# Increase iterations someday!
# Get means for all interactions, including ...
# Those interactions where each species had a slightly different climate dataset

# Compare this model to the one (that should be the same!) in tempmodels_datasets
summ_studyspp$model.mean <- colMeans(goo$b)
othermodel <- read.csv("output/compare.tempmodels.csv", header=TRUE)
compare.samemodel <- merge(summ_studyspp, othermodel, by="datasetid")
plot(data.stanfit~model.mean, data=compare.samemodel)
abline(0,1)

# Now we need to get the averages since some species 
temp.change <- matrix(0, ncol=1000, nrow=18) # 1000 iterations for 22 interactions;
allspecieswcoef <- matrix(0, ncol=1000, nrow=18)
for (i in 2000:3000){ # 2000 iterations?
    summ_studyspp$model.est <- goo$b[i,]
    # step 1: Get back the missing species
    allspecieswcoef <- merge(summ_studyspp, datasetid.trans, by=c("datasetid",
         "studyid"), all.y=TRUE, suffixes=c(".model", ".data"))
    allspecieswcoef <- subset(allspecieswcoef, select=c("datasetid", "model.est", "species.data"))
    names(allspecieswcoef)[names(allspecieswcoef)=="species.data"] <- "species"
    # step 2 merge in spp1 and spp2 
    andtheanswer <- merge(temp_int, allspecieswcoef, by.x=c("datasetid", "spp1"),
        by.y=c("datasetid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, allspecieswcoef, by.x=c("datasetid","spp2"),
        by.y=c("datasetid", "species"), all.x=TRUE, suffixes=c(".spp1",".spp2"))
    # step 3 now we need to aggregate by spp1 and spp2 for the cases where we have
    # different estimates for each spp in an interaction
    models.intid <- data.frame(intid=rep(andtheanswer$intid, 2), model.est=c(
        andtheanswer$model.est.spp1, andtheanswer$model.est.spp2))
    modelmeanz <- ddply(models.intid , c("intid"), summarise,
        model.mean=mean(model.est, na.rm=TRUE))
    temp.change[,(i-2000)] <- modelmeanz$model.mean
}

row.names(temp.change) <- modelmeanz$intid
temp.change <- as.data.frame(temp.change)
temp.change$intid <- modelmeanz$intid

intid.wtempchange <- data.frame(intid=modelmeanz$intid, stanfit=rowMeans(temp.change, na.rm=TRUE)) # mean across iterations for EACH dataset; STAN SLOPE ESTIMATE PER dataset

intid.wtempchange1K.long <- melt(temp.change, id="intid")

write.csv(intid.wtempchange1K.long, "output/temp.change.1K.csv")
write.csv(intid.wtempchange, "output/temp.change.meanover1K.csv", row.names=FALSE)
