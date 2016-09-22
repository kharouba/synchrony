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

#Hinge
hinge<-subset(clim3, species=="Acrocephalus arundinaceus" | species=="Acrocephalus scirpaceus" | species=="Copepod1 spp." | species=="Copepod2 spp." | species=="Cyclops vicinus"  | species=="Daphnia3 spp." | species=="Diatom3 spp." | species=="Perca fluviatillis" | species=="Phytoplankton1 spp." | species=="Pleurobrachia pileus" | species=="Pleurobrachia_a pileus")

hinge_non<-subset(clim3, species!="Acrocephalus arundinaceus" & species!="Acrocephalus scirpaceus"  & species!="Copepod1 spp." & species!="Copepod2 spp." & species!="Cyclops vicinus"  & species!="Daphnia3 spp." & species!="Diatom3 spp." & species!="Perca fluviatillis" & species!="Phytoplankton1 spp." & species!="Pleurobrachia pileus" & species!="Pleurobrachia_a pileus")

hinge_non$newyear<-hinge_non$year
hinge_pre<-subset(hinge, year<=1981); hinge_pre$newyear<-1981
hinge_post<-subset(hinge, year>1981); hinge_post$newyear<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)

clim4<-rbind(hinge_non, hinges);
clim3<-clim4