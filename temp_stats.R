rm(list=ls())
#row.names=FALSE
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/2013_2015/input/")

## HAS TEMP CHANGED?
clim<-read.csv("climate3.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
clim$name_env<-with(clim, paste(species, phenophase, sep="_"))
clim<-subset(clim, phenophase!="start" & extra!="10")
clim$envvalue<-as.numeric(clim$envvalue)
clim$envfactor[clim$envfactor=="temperaure"] <- "temperature"

# Mean temp across sites for those studies with multiple sites
sites<-subset(clim, studyid=="HMK018" | studyid=="HMK019" | studyid=="HMK023" & site!="tomakomai")
nosites<-subset(clim, site=="tomakomai" | studyid!="HMK018" & studyid!="HMK019" & studyid!="HMK023")
nosites<-nosites[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")]

new<-with(sites, aggregate(envvalue, by=list(studyid, year, species), FUN=mean, na.rm=T)) # across sites
names(new)[1]<-"studyid"; names(new)[2]<-"year"; names(new)[3]<-"species"; names(new)[4]<-"envvalue"
new<-new[order(new$studyid),]
sites2<-merge(sites[,c("studyid","envfactor","envunits","envtype","year","species")], new, by=c("studyid","year","species"))
sites3<-unique(sites2[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")])

clim2<-rbind(nosites, sites3) # all years with data
#merge with spp data so calculating env change only over the years of interaction
total3<-read.csv("int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
total3<-na.omit(total3)
clim3<-merge(clim2, total3[,c("studyid","year")], by=c("studyid","year"))
clim3<-unique(clim3[,c("studyid","year","envfactor","envunits","envtype","species","envvalue")])
clemo<-unique(clim3[,c("studyid","year")])
mi<-aggregate(clemo["year"], by=clemo[c("studyid")], FUN=length); std(mi$year) #number of years
mimin<-aggregate(clemo["year"], by=clemo[c("studyid")], FUN=min) #number of years
mimax<-aggregate(clemo["year"], by=clemo[c("studyid")], FUN=max)
mimax$range<-mimax$year-mimin$year; mean(mimax$range); std(mimax$range)
std <- function(x) sd(x)/sqrt(length(x))


# ONLY FOR YEARS WITH SPECIES DATA
m1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype=="air" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-year); anova(m3,m1)
#water
m1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype=="water" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-year); anova(m3,m1)

m1<-lme(envvalue~envtype+year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype!="ground" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & envunits!="doy"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-envtype); anova(m3,m1)