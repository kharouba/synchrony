### Started 6 May 2015 ###

## Trying to run STAN on linear models with synchrony data ##

options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")
library(ggplot2)
library(rstan)
library(shinyStan)

rawlong <- read.csv("indiv_sppdata2.csv", header=TRUE)
rawlong$int_type[which(rawlong$int_type=="poll")] <- "pollination"
rawlong.X <- NULL


#fix for southern hemisphere
sub<-subset(rawlong, studyid=="HMK018")
sub$phenovalue2<-with(sub, phenovalue+182) #number of days different between Dec 22 (SH summer solstice) and June 21 (NH summer solstice), to put all dates on same 'calendar'
sub$phenovalue3<-with(sub, phenovalue2-365) #now same number of days before NH summer solstice as it was before SH summer solstice


##
# figure out how many unique species
# and alter the data so each unique species shows up once
# watch out: if you leave in int_type you still get duplicated species...
# so it's not below until after duplicate species are removed!
specieschar.wdups <- aggregate(rawlong["phenovalue"],
    rawlong[c("studyid", "species", "spp")], FUN=length)
specieschar <- aggregate(specieschar.wdups["phenovalue"],
    specieschar.wdups[c("studyid", "species")], FUN=length)
dupspp <- subset(specieschar, phenovalue>1)
specieschar.wdups[which(specieschar.wdups$species %in% dupspp$species),] #there are some species with mult relationships
# delete duplicate species as spp1 (generally the shorter timeseries)
rawlong.nodups <- rawlong[-(which(rawlong$species %in% dupspp$species &
    rawlong$spp=="spp1")),]
# and order it!
rawlong.nodups <- rawlong.nodups[with(rawlong.nodups, order(species, year)),]
specieschar.formodel <- aggregate(rawlong.nodups["phenovalue"],
    rawlong.nodups[c("studyid", "species", "int_type", "spp")], FUN=length)

rawlong.nodups$yr1976 <- rawlong.nodups$year-1976
N <- nrow(rawlong.nodups)
y <- rawlong.nodups$phenovalue
J <- nrow(specieschar.formodel)
species <- as.numeric(as.factor(rawlong.nodups$species))
year <- rawlong.nodups$yr1976

#Run the model: complete pooling:
fit1<-stan("synchrony1_lmfit.stan",data=c("N","y","year"), iter=100, chains=4)
print(fit1) #Rhat should be ~1

# SAME MODEL AS:
# get estimate from no pooling (fixed effects) and complete pooling
comppool <- lm(phenovalue~year, data=rawlong.nodups)

### Stan model for each species ###
#cheating
rawlong.nodups$sp <- as.numeric(as.factor(rawlong.nodups$species))
Bgroups<-unique(rawlong.nodups$sp); b<-Bgroups; b<-as.character(b)
best<-data.frame(array(0, c(length(b), 3)))
names(best)[1]<-"species"; names(best)[2]<-"int"; names(best)[3]<-"slope"
#for(i in 1:length(b)){
for(i in 1:5){
mat<-rawlong.nodups[rawlong.nodups$sp==b[i],]
N <- nrow(mat)
y <- mat$phenovalue
year <- mat$yr1976
best[i,1]<-b[i]
fit1<-stan("synchrony1_lmfit.stan",data=c("N","y","year"), iter=100, chains=4)
qq<-summary(fit1)
best[i,2]<-qq[[1]][1,1] #alpha/int
best[i,3]<-qq[[1]][2,1] #beta/slope
}

# not cheating
Ndata <- nrow(rawlong.nodups)
Jspp <- nrow(specieschar.formodel)
fit1<-stan("synchrony1_lmfitwspp.stan",data=c("Jspp","Ndata","y","year","species"), iter=100, chains=4)



