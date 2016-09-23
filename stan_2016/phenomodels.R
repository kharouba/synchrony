##### ESTIMATES SPECIES PHENOLOGICAL SHIFT ####

rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis")
# setwd("~/Documents/git/projects/trophsynch/synchrony")
library(ggplot2)
library(rstan)
library(shinyStan)
library(grid)
library(nlme)
library(dplyr)
set_cppo("fast")  # for best running speed
source("/users/kharouba/google drive/UBC/multiplot.R")
#library(reshape)
# library(lme4)

#get data
source("datacleaning.R")
# in the data file: spp1 = neg= negative species e.g. resource

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
# setwd("~/Documents/git/projects/trophsynch/synchrony/stan_2016")

rawlong <- read.csv("rawlong.csv", header=TRUE)
##
# figure out how many unique species
# and alter the data so each unique species shows up once
# watch out: if you leave in int_type you still get duplicated species...
# so it's not below until after duplicate species are removed!
specieschar.wdups <- aggregate(rawlong["phenovalue"],
    rawlong[c("studyid", "species", "spp")], FUN=length)
specieschar <- aggregate(specieschar.wdups["phenovalue"],
    specieschar.wdups[c("studyid", "species")], FUN=length)
#dupspp <- subset(specieschar, phenovalue>1)
#specieschar.wdups[which(specieschar.wdups$species %in% dupspp$species),] #there are some species with mult relationships
# delete duplicate species as spp1 (generally the shorter timeseries)
#rawlong.nodups <- rawlong[-(which(rawlong$species %in% dupspp$species &
   # rawlong$spp=="spp1")),]
# and order it!
rawlong.nodups<-rawlong
rawlong.nodups <- rawlong.nodups[with(rawlong.nodups, order(species, year)),]
specieschar.formodel <- aggregate(rawlong.nodups["phenovalue"],
    rawlong.nodups[c("studyid", "species", "intid", "terrestrial","spp")], FUN=length) 
 specieschar.formodel <- specieschar.formodel[with(specieschar.formodel, order(species)),]   

#number of years per species
    
##

## Hinge models

##
# Heather's code to add in data formatted for hinge model
hinge <- subset(rawlong.nodups, intid=="170" | intid=="171" | intid=="177" |
    intid=="178" | intid=="179" | intid=="180" |intid=="181" | intid=="189" |
    intid=="191"|intid=="193" |intid=="194" | intid=="195" |intid=="196"|
    intid=="201" |intid=="207" |intid=="208")

hinge_non <- subset(rawlong.nodups, intid!="170" & intid!="171" & intid!="177"
     & intid!="178" & intid!="179" & intid!="180" & intid!="181" & intid!="189"
     & intid!="191" & intid!="193" & intid!="194" & intid!="195" & intid!="196"
     & intid!="201" & intid!="207" & intid!="208")

hinge_non$newyear<-hinge_non$year
hinge_pre<-subset(hinge, year<=1981); hinge_pre$newyear<-1981
hinge_post<-subset(hinge, year>1981); hinge_post$newyear<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)
rawlong.tot<-rbind(hinge_non, hinges)


# prep the data to fit the model including:
# aggregate to get species level characteristics
# subset down to the phenovalues

rawlong.tot <- arrange(rawlong.tot, species)
rawlong.tot$yr1981 <- rawlong.tot$newyear-1981
rawlong.tot <- rawlong.tot[with(rawlong.tot, order(species)),]
N <- nrow(rawlong.tot)
y <- rawlong.tot$phenovalue
specieschar.hin<- aggregate(rawlong.tot["phenovalue"], rawlong.tot[c("studyid","intid","species", "int_type", "terrestrial", "spp")], FUN=length) #number of years per species
specieschar.hin <- specieschar.hin[with(specieschar.hin, order(species)),]
specieschar.hin2 <- unique(specieschar.hin$species)
Nspp <- length(specieschar.hin2)
species <- as.numeric(as.factor(rawlong.tot$species))
year <- rawlong.tot$yr1981
#nVars <-1
# get ecosystem for each species
specieschar.hin3 <- subset(specieschar.hin, select=c("species", "terrestrial"))
specieschar.hin4 <- specieschar.hin3[!duplicated(specieschar.hin3),]
eco <- as.numeric(as.factor(specieschar.hin4$terrestrial))
Neco <- length(unique(specieschar.hin4$terrestrial))


#New model as of June 2016
#Random slopes only, no random intercepts, hinge, no covariate matrix:
#sync.model<-stan("synchrony1_notype_randslops_wcovar.stan", data=c("N","Nspp","y","species","year"), iter=2000, warmup=1000, thin=10, chains=4)
sync.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

## ADD ECOSYSTEM here
# current error: Error : variable does not exist; processing stage=data initialization; variable name=p; base type=int
# sync.model.eco <-stan("stanmodels/threelevelrandomslope_eco.stan", data=c("N","Nspp","y","species","year", "eco", "Neco"), iter=3000, chains=4)


#Match up interacting species and look at differences

fh.sim <- extract(sync.model)# inc_warmup=FALSE) #only extracts parameters, does not include warmup
dim(fh.sim$b) # number of iterations*number of species 
names(fh.sim)
#sd is not considered parameter (only rows in output are parameters)
# here's one iteration
fh.sim$b[2000,]

specieschar.formodel <- aggregate(rawlong.nodups["phenovalue"],
    rawlong.nodups[c("studyid", "species", "intid", "terrestrial","spp")], FUN=length) 
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species"))
specieschar.formodel.sm  <- specieschar.formodel.sm [with(specieschar.formodel.sm , order(species)),]
intid <- read.csv("input/raw_april.csv", header=TRUE)
lal<-unique(rawlong.tot[,c("intid","terrestrial")])
intid2<-merge(intid, lal, by=c("intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction","terrestrial"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
sync_int<-intid.nodups #synchrony change interactions


summ_studyspp <- subset(specieschar.formodel, select=c("studyid", "species")); summ_studyspp<-unique(summ_studyspp)
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]

## RESULTS !!	
# For spp' phenological change
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- fh.sim$b[i,]
    it1000[,(i-3000)] <- fh.sim$b[i,]
}
summ_studyspp$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP
mean(summ_studyspp$stanfit)

it1K <- matrix(0, ncol=1000, nrow=Nspp)
for (i in 3000:4000){ 
    summ_studyspp$model <- fh.sim$b[i,]
    it1K[,(i-3000)] <- fh.sim$b[i,]
}

it1K.out <- as.data.frame(it1K)
it1K.out$studyid <- summ_studyspp$studyid
it1K.out$species <- summ_studyspp$species

library(reshape)
it1K.out.long <- melt(it1K.out, id=c("studyid", "species"))

write.csv(it1K.out.long , "output/pheno.change.1K.csv")

#computation of the standard error of the mean
sem<-sd(summ_studyspp$stanfit)/sqrt(length(summ_studyspp$stanfit)); sem
#95% confidence intervals of the mean
c(mean(summ_studyspp$stanfit)-2*sem,mean(summ_studyspp$stanfit)+2*sem)

#aquatic vs. terrestrial
lol<-unique(specieschar.formodel[,c("studyid","species","terrestrial")])
summ_studyspp<-merge(summ_studyspp, lol, by=c("studyid", "species"))
data<-subset(summ_studyspp, terrestrial=="aquatic"); mean(data$stanfit); sd(data$stanfit)
data<-subset(summ_studyspp, terrestrial=="terrestrial"); mean(data$stanfit); sd(data$stanfit)

sem<-sd(data$stanfit)/sqrt(length(data$stanfit)); sem
#95% confidence intervals of the mean
c(mean(data$stanfit)-2*sem,mean(data$stanfit)+2*sem)
