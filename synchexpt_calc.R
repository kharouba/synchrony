 ### Started 23 January 2014 ###
### By Lizzie & Heather ###

## Cleaning (perhaps) the experimental data #
## And -- more importantly -- calculating a standardized fitness metric ##

## Last updated: 23 June 2014 (whoo, that is a long time since January) #
## see SynconryLizzieNotebooksJan2014_forcode. pdf for what I was trying to do #
## on the train to Avignon, on y danse, on y danse ... ##

## f(x)s adapted from timingmetafuncs.R (Lizzie's old code) ##

# basic housekeeping #
rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/")

library(reshape)

manip <- read.csv("manip.csv", header=TRUE)
intxn <- read.csv("intxn.csv", header=TRUE)

intxn$X <- NULL
intxnlong <- melt(intxn, id=c("studyid", "interaction", "intid"))
names(intxnlong)[names(intxnlong)=="variable"] <- "spp1or2"
names(intxnlong)[names(intxnlong)=="value"] <- "value"

# basic look-see of the data
unique(paste(manip$tempman, manip$phenoman)) # check out SET006 and see what is up there
unique(manip$secondfactor) # shaded or unshaded only

manip$latbi <- paste(manip$genus, manip$species)
manip$latbi2 <- paste(manip$genus2, manip$species2)

# fix one error in SET006
manip$tempman[manip$factormanipulated=="consumer" & manip$studyid=="SET006"] <- "NA"
manip$phenoman[manip$factormanipulated=="temperature" & manip$studyid=="SET006"] <- "NA"

manip<- subset(manip, secondfactor=="unshaded"|is.na(secondfactor)==TRUE)

# exclude studies with species data in genus2, species2
manip1spp <- subset(manip, is.na(genus2)==TRUE)

# exclude rows where effect of manipulation was on same species (synchrony with other species not considered)
manip<-subset(manip, factormanipulated!="same species")

#change treatment number so that mult treatments will be counted
manip$treatment[manip$treatment=="2" & manip$studyid=="KJB001"] <- "0"

manip$treatment[manip$treatment=="0" & manip$studyid=="HMK002"] <- "4"
manip$treatment[manip$treatment=="1" & manip$studyid=="HMK002"] <- "5"

# on 23 June 2014 I (Lizzie) added JMA002 to this exclude list #
# JMA002 is a lab study so it has no year so it does not work #
# (just add in a year and it would be fine!) #
#manip$year[manip$studyid=="JMA002"] <- "1992"

#choose 1 phenophase /study/interaction
manip2<-subset(manip, studyid!="HMK014")
su<-subset(manip, studyid=="HMK014" & phenophase=="peak density")
manip3<-rbind(manip2, su)

manip4<-subset(manip3, studyid!="HMK045")
su2<-subset(manip3, studyid=="HMK045" & phenophase!="95th percentile emergence")
manip5<- rbind(manip4, su2)

manip6<-subset(manip5, studyid!="HMK044")
su3<-subset(manip5, studyid=="HMK044" & phenophase!="end of budburst" & phenophase!="start of shoot elongation" & phenophase!="end of shoot elongation")
manip7<- rbind(manip6, su3)

manip8<-subset(manip7, studyid!="HMK041")
su4<-subset(manip7, studyid=="HMK041" & phenophase!="second peak")
manip9<- rbind(manip8, su4)

manip<-manip9

## choose 1 fitness estimate /species; randomly choose or go with positive interactor (ie. JEH001)
#manip2<-subset(manip, studyid!="HMK002")
#su<-subset(manip, studyid=="HMK002" & fitnesstype!="3")
#manip3<-rbind(manip2, su)
manip2<-subset(manip, studyid!="HMK002")
su<-subset(manip, studyid=="HMK002" & fitnesstype!="1")
manip3<-rbind(manip2, su)

#manip4<-subset(manip3, studyid!="JEH001")
#su2<-subset(manip3, studyid=="JEH001" & fitnesscomp!="host development time (ffinal instar to pupation)") #choose positive spp over negative spp
#manip5<- rbind(manip4, su2)

manip6<-subset(manip3, studyid!="JMA002")
su3<-subset(manip3, studyid=="JMA002" & fitnesstype=="1")
manip7<- rbind(manip6, su3)

manip8<-subset(manip7, studyid!="HMK044")
su4<-subset(manip7, studyid=="HMK044" & fitnesscomp=="larval survival")
#su4<-subset(manip7, studyid=="HMK044" & fitnesscomp=="pupal survival")
#su4<-subset(manip7, studyid=="HMK044" & fitnesscomp=="female adult dry mass")
su5<-subset(manip7, studyid=="HMK044" & fitnesscomp=="a")
manip9<- rbind(manip8, su4, su5)

manip10<-subset(manip9, studyid!="KJB001")
su6<-subset(manip9, studyid=="KJB001" & fitnesscomp=="fecundity")
manip11<- rbind(manip10, su6)

manip12<-subset(manip11, studyid!="HMK045")
su7<-subset(manip9, studyid=="HMK045" & fitnesscomp=="survivorship")
manip13<- rbind(manip12, su7)

manip<-manip13

# The following functions require:
# data: a dataframe object wth the following columns:
#   studyid
#   latbi1 and latbi2
#   phenophase
#   treatment
#   fitnesscomp
#   year -- MUST have year as a number (not NA as in JMA002) otherwise it will error
#   site
# You then have to tell it these things:
#   data: the dataframe
#   responsecol (fitness, phenovalue etc.)
#   treatcon: the treatment to be considered the control (0 or 1)
#   the 'heated treatment' is treated as >treatcon (and currently the code only allows for one)

#For fitness= heat-control/control
stdeffect.synch <- function(data, responsecol, treatcon){
sensexpt <- data.frame(cbind(studyid=character(0), site=character(0),
     year=numeric(0), respchange=numeric(0),
     uniqueintxn=character(0), heatlevel=numeric(0), control=numeric(0)))
print("Alert! This code does not work across different fitness types or phenophases currently. It just takes row 1. Beware.") 
for (uniquestudyid in seq_along(unique(data$studyid))){
    subsettostudyid <- subset(data, studyid==unique(data$studyid)[uniquestudyid])
    for (uniquesite in seq_along(unique(subsettostudyid$site))){
        subsettosite <- subset(subsettostudyid, site==unique(subsettostudyid$site)[uniquesite])
        for (uniqueyr in seq_along(unique(subsettosite$year))){
            subsettoyr <- subset(subsettosite, year==unique(subsettosite$year)[uniqueyr])     
            subsettoyr$uniqueintxn <- paste(subsettoyr$latbi, subsettoyr$latbi2,
                subsettoyr$phenophase, subsettoyr$fitnesscomp)
    control <- subset(subsettoyr, treatment==treatcon)
    heat1 <- subset(subsettoyr, treatment>treatcon)
    spp <- control[!duplicated(control$uniqueintxn),]
      for (speciesy in 1:nrow(spp)){
             controlsp <- subset(control, uniqueintxn==
                 spp$uniqueintxn[speciesy])
             heatsp <- subset(heat1, uniqueintxn==
                 spp$uniqueintxn[speciesy])
             for (treater in 1:length(unique(heat1$treatment))){
                 heatadd <- subset(heatsp, treatment==
                    unique(heat1$treatment)[treater])
                 sensadd <- data.frame(cbind(studyid=unique(data$studyid)[uniquestudyid],
                     site=unique(subsettoyr$site), # should be one only
                     year=unique(subsettoyr$year), # should be one only                      
                     respchange=((heatadd[[responsecol]][1]- 
                     controlsp[[responsecol]][1])/controlsp[[responsecol]][1]), uniqueintxn=
                     spp$uniqueintxn[speciesy], heatlevel=
                     unique(heat1$treatment)[treater]), control=
                     controlsp[[responsecol]][1])
             sensexpt <- rbind(sensexpt, sensadd)
               }
        }
  }
  }
  }
return(sensexpt)
}

# rm studies with no 0 control:
manip.sm <- subset(manip, studyid != "HMK003" & studyid != "HMK005"
    & studyid != "KJB001" & studyid!="HMK002" & studyid != "HMK044")
# note that HMK has a 0 control but not fitness data for the control


# now run stuff
# output: one df for pheno manip (response is fitness), one df for C manip (response is phenology or fitness -- just one study does both) 
man.fit <- stdeffect.synch(manip.sm, "fitnessvalue", 0)


# HMK002 has no fitness data for the control treatement so ...
# we assign 1 as 0 since that is closes to ambient
HMK002 <- subset(manip, studyid=="HMK002")
phenman.fit2 <- stdeffect.synch(HMK002, "fitnessvalue", 2)

KJB001 <- subset(manip, studyid=="KJB001")
phenman.fit3 <- stdeffect.synch(KJB001, "fitnessvalue", 0)

HMK003 <- subset(manip, studyid=="HMK003")
phenman.fit4 <- stdeffect.synch(HMK003, "fitnessvalue", 1)

HMK005 <- subset(manip, studyid=="HMK005")
phenman.fit5 <- stdeffect.synch(HMK005, "fitnessvalue", 1)

HMK044 <- subset(manip, studyid=="HMK044")
phenman.fit6 <- stdeffect.synch(HMK044, "fitnessvalue", 1)

man.fit2<-rbind(man.fit, phenman.fit2, phenman.fit3, phenman.fit4, phenman.fit5, phenman.fit6)
man.fit<-man.fit2
names(man.fit)[4]<-"fitchange"

#for Pheno (heat-control)
stdeffect.synch <- function(data, responsecol, treatcon){
sensexpt <- data.frame(cbind(studyid=character(0), site=character(0),
     year=numeric(0), respchange=numeric(0),
     uniqueintxn=character(0), heatlevel=numeric(0), control=numeric(0)))
print("Alert! This code does not work across different fitness types or phenophases currently. It just takes row 1. Beware.") 
for (uniquestudyid in seq_along(unique(data$studyid))){
    subsettostudyid <- subset(data, studyid==unique(data$studyid)[uniquestudyid])
    for (uniquesite in seq_along(unique(subsettostudyid$site))){
        subsettosite <- subset(subsettostudyid, site==unique(subsettostudyid$site)[uniquesite])
        for (uniqueyr in seq_along(unique(subsettosite$year))){
            subsettoyr <- subset(subsettosite, year==unique(subsettosite$year)[uniqueyr])     
            subsettoyr$uniqueintxn <- paste(subsettoyr$latbi, subsettoyr$latbi2,
                subsettoyr$phenophase, subsettoyr$fitnesscomp)
    control <- subset(subsettoyr, treatment==treatcon)
    heat1 <- subset(subsettoyr, treatment>treatcon)
    spp <- control[!duplicated(control$uniqueintxn),]
      for (speciesy in 1:nrow(spp)){
             controlsp <- subset(control, uniqueintxn==
                 spp$uniqueintxn[speciesy])
             heatsp <- subset(heat1, uniqueintxn==
                 spp$uniqueintxn[speciesy])
             for (treater in 1:length(unique(heat1$treatment))){
                 heatadd <- subset(heatsp, treatment==
                    unique(heat1$treatment)[treater])
                 sensadd <- data.frame(cbind(studyid=unique(data$studyid)[uniquestudyid],
                     site=unique(subsettoyr$site), # should be one only
                     year=unique(subsettoyr$year), # should be one only                      
                     respchange=(heatadd[[responsecol]][1]- 
                     controlsp[[responsecol]][1]), uniqueintxn=
                     spp$uniqueintxn[speciesy], heatlevel=
                     unique(heat1$treatment)[treater]), control=
                     controlsp[[responsecol]][1])
                     sensexpt <- rbind(sensexpt, sensadd)
               }
        }
  }
  }
  }
return(sensexpt)
}

# now run stuff
# output: one df for pheno manip (response is fitness), one df for C manip (response is phenology or fitness -- just one study does both) 
#man.phen <- stdeffect.synch(manip.sm, "phenovalue", 0) #just doy for each species
man.phen <- stdeffect.synch(manip.sm, "treatchangeday", 0) #difference in doy between species 1 and 2 for each treatment

# HMK002 has no fitness data for the control treatement so ...
# we assign 1 as 0 since that is closes to ambient
HMK002 <- subset(manip, studyid=="HMK002")
phenman.phen2 <- stdeffect.synch(HMK002, "treatchangeday", 2)

KJB001 <- subset(manip, studyid=="KJB001")
phenman.phen3 <- stdeffect.synch(KJB001, "treatchangeday", 0)

HMK003 <- subset(manip, studyid=="HMK003")
phenman.phen4 <- stdeffect.synch(HMK003, "treatchangeday", 1)

HMK005 <- subset(manip, studyid=="HMK005")
phenman.phen5 <- stdeffect.synch(HMK005, "treatchangeday", 1)

HMK044 <- subset(manip, studyid=="HMK044")
phenman.phen6 <- stdeffect.synch(HMK044, "treatchangeday", 1)

man.phen2<-rbind(man.phen, phenman.phen2, phenman.phen3, phenman.phen4, phenman.phen5, phenman.phen6)
man.phen<-man.phen2
names(man.phen)[4]<-"phenochange"

manoh<-merge(man.fit[,1:6], man.phen, by=c("studyid","site","year","uniqueintxn","heatlevel"))
manoh$fitchange<-as.numeric(manoh$fitchange)
manoh$phenochange<-as.numeric(manoh$phenochange)

# temp change
man.temp <- stdeffect.synch(manip.sm, "treatchangeC", 0) #difference in doy between species 1 and 2 for each treatment

# HMK002 has no fitness data for the control treatement so ...
# we assign 1 as 0 since that is closes to ambient
HMK002 <- subset(manip, studyid=="HMK002")
man.temp2 <- stdeffect.synch(HMK002, "treatchangeC", 2)

KJB001 <- subset(manip, studyid=="KJB001")
man.temp3 <- stdeffect.synch(KJB001, "treatchangeC", 0)

HMK003 <- subset(manip, studyid=="HMK003")
man.temp4 <- stdeffect.synch(HMK003, "treatchangeC", 1)

HMK005 <- subset(manip, studyid=="HMK005")
man.temp5 <- stdeffect.synch(HMK005, "treatchangeC", 1)

HMK044 <- subset(manip, studyid=="HMK044")
man.temp6 <- stdeffect.synch(HMK044, "treatchangeC", 1)

man.temp7<-rbind(man.temp, man.temp2, man.temp3, man.temp4, man.temp5, man.temp6)
man.temp<-man.temp7
names(man.temp)[4]<-"tempchange"

manoh2<-merge(manoh, man.temp[,c(1:6)], by=c("studyid","site","year","uniqueintxn","heatlevel"))
manoh2$tempchange<-as.numeric(manoh2$tempchange)
manoh<-manoh2

manip$uniqueintxn <- with(manip, paste(latbi, latbi2, phenophase, fitnesscomp))
manoh2<-merge(manoh, manip[,c("studyid","site","year","fitnesstype","uniqueintxn")], by=c("studyid","site","year","uniqueintxn"))

## make a smaller dataframe of the unique treatments info
## add site, year, species
#manip.infoall <- subset(manip, select=c("studyid", "site","year", "tempman",
    #"labfield", "phenophase","fitnesscomp","treatment", "treatchangeC", "treatchangeday"))
#manip.info <- manip.infoall[!duplicated(manip.infoall),]

