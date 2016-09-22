##### TO GET TEMP CHANGE PER INTERACTION ######

rm(list=ls())
#row.names=FALSE
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(ggplot2)
library(rstan)
library(shinyStan)
library(grid)
library(nlme)
library(gridExtra)
library(reshape)
set_cppo("fast") 

## HAS TEMP CHANGED?
clim<-read.csv("input/climate3.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
clim$name_env<-with(clim, paste(species, phenophase, sep="_"))
clim<-subset(clim, phenophase!="start" & extra!="10" & envunits=="C" & studyid!="HMK028" & studyid!="HMK029" & studyid!="HMK037") #029,037 because nutrients>phenology; 028- because temperature sum
clim$envvalue<-as.numeric(clim$envvalue)
clim$envfactor[clim$envfactor=="temperaure"] <- "temperature"

# Mean temp across sites for those studies with multiple sites
sites<-subset(clim, studyid=="HMK018" | studyid=="HMK019" | studyid=="HMK023" & site!="tomakomai")
nosites<-subset(clim, site=="tomakomai" | studyid!="HMK018" & studyid!="HMK019" & studyid!="HMK023")
nosites<-nosites[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")]

new<-with(sites, aggregate(envvalue, by=list(studyid, year, species), FUN=mean, na.rm=T)) # across sites
names(new)[1]<-"studyid"; names(new)[2]<-"year"; names(new)[3]<-"species"; names(new)[4]<-"envvalue"
new<-new[order(new$studyid),]
#'all' species: "HMK018" "HMK023" "HMK028" "HMK029" "HMK036" "HMK042" "HMK043" >> now fixed in climate3.csv

sites2<-merge(sites[,c("studyid","envfactor","envunits","envtype","year","species")], new, by=c("studyid","year","species"))
sites3<-unique(sites2[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")])

clim2<-rbind(nosites, sites3) # all years with data
#merge with spp data so calculating env change only over the years of interaction
total3<-read.csv("input/raw_april.csv", header=TRUE, na.strings="NA", as.is=TRUE)
total3<-na.omit(total3)
clim3<-merge(clim2, total3[,c("studyid","year")], by=c("studyid","year"))
clim3<-unique(clim3[,c("studyid","year","envfactor","envunits","envtype","species","envvalue")])
clim3<-na.omit(clim3)

#Hinge (11 spp)
hinge<-subset(clim3, species=="Acrocephalus arundinaceus" | species=="Acrocephalus scirpaceus" | species=="Copepod1 spp." | species=="Copepod2 spp." | species=="Cyclops vicinus"  | species=="Daphnia3a spp." | species=="Daphnia3b spp."| species=="Diatom3 spp." | species=="Perca fluviatillis" | species=="Phytoplankton1 spp." | species=="Pleurobrachia pileus" | species=="Pleurobrachia_a pileus")

#Non-hinge (11 spp)
hinge_non<-subset(clim3, species!="Acrocephalus arundinaceus" & species!="Acrocephalus scirpaceus"  & species!="Copepod1 spp." & species!="Copepod2 spp." & species!="Cyclops vicinus"  & species!="Daphnia3a spp." & species!="Daphnia3b spp." & species!="Diatom3 spp." & species!="Perca fluviatillis" & species!="Phytoplankton1 spp." & species!="Pleurobrachia pileus" & species!="Pleurobrachia_a pileus")

hinge_non$newyear<-hinge_non$year
hinge_pre<-subset(hinge, year<=1981); hinge_pre$newyear<-1981
hinge_post<-subset(hinge, year>1981); hinge_post$newyear<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)

clim4<-rbind(hinge_non, hinges);
clim3<-clim4
clim3$yr1981 <- clim3$newyear-1981

clim3 <- clim3[with(clim3, order(species, year)),]
N<-nrow(clim3)
y <- clim3$envvalue
specieschar.hin<- aggregate(clim3["envvalue"], clim3[c("species")], FUN=length)
specieschar.hin <- specieschar.hin[with(specieschar.hin, order(species)),]
Nspp <- nrow(specieschar.hin) #J
species <- as.numeric(as.factor(clim3$species))
year <- clim3$yr1981

temp.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

WHY AM I GETTING TWO DIFFERENT SLOPE ESTIMATES FROM SAME DATA, E.G .HMK042?????
ASSUME FOR NOW ONE ESTIMATE OF TEMP CHANGE PER STUDY!


!!!!!!! Mix up happens between species > interactions, beween summ_studyspp, temp_int

summ_studyspp <- subset(clim3, select=c("studyid","species")); 
summ_studyspp<-unique(summ_studyspp)
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]

goo <- extract(temp.model) 

specieschar.formodel <- aggregate(clim3["envvalue"], clim3[c("studyid","species")], FUN=length) #number of years
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid","species"))
rawlong <- read.csv("rawlong.csv", header=TRUE)
sol<-merge(specieschar.formodel.sm, unique(rawlong[,c("studyid","intid","species")]), by=c("studyid","species"))
#sol<- sol[with(sol,order(species)),]
intid <- read.csv("input/raw_april.csv", header=TRUE)
intid2<-merge(intid, unique(sol[,c("studyid","intid")]), by=c("studyid","intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
temp_int<-intid.nodups #temp change interactions
temp_int<- temp_int[with(temp_int,order(intid)),]

temp.change <- matrix(0, ncol=3000, nrow=22) #2000 iterations for 22 interactions;
for (i in 3000:6000){ # 2000 iterations?
    summ_studyspp$model <- goo$b[i,]
    andtheanswer <- merge(temp_int, summ_studyspp, by.x=c("studyid", "spp1"), by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, summ_studyspp, by.x=c("studyid","spp2"), by.y=c("studyid", "species"), all.x=TRUE)
       andtheanswer[andtheanswer$spp1=="Engraulis japonicus","model.x"]<- andtheanswer[andtheanswer$spp1=="Engraulis japonicus","model.y"] #Engraulis japonicus because first day of temp X
       andtheanswer[andtheanswer$spp1=="Diatom2a","model.x"]<-andtheanswer[andtheanswer$spp1=="Diatom2a","model.y"] #Diatom2 spp. because first day of temp X
       andtheanswer[andtheanswer$spp1=="Diatom2b","model.x"]<-andtheanswer[andtheanswer$spp1=="Diatom2b","model.y"]#Diatom2 spp.] because first day of temp X
    andtheanswer$meanmodel <- rowMeans(andtheanswer[,c(6:7)]) # take mean because want avg temp change PER interaction
    temp.change[,(i-3000)] <- andtheanswer$meanmodel
}

intid.nodups$stanfit <- rowMeans(temp.change, na.rm=TRUE) #mean across iterations for EACH dataset; STAN SLOPE ESTIMATE PER dataset
