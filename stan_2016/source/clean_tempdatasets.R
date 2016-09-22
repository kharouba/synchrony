## taken from tempmodels_datasets.R ##

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


################################
### to get absolute value of temperature change (i.e. based on unique datasets) ###
#Stan model
#HMK011- 2 unique
#HMK019- 2 unique
#HMK038 - 2 uniqu3
#HMK051- 3
#HMK023-  1 unique; HMK031-1; HMK034- 1, HMK043=1; #HMK016 - 1 #HMK018 - 1; #HMK052- 1; #HMk036- 1; #HMK042-1

clim3 <- clim3[with(clim3, order(studyid, year)),]
sub1<-subset(clim3, studyid=="HMK011" & species=="Epirrita autumnata")
sub1$datasetid<-with(sub1, paste(studyid,"_",1))
sub2<-subset(clim3, studyid=="HMK011" & species=="Poecile montanus")
sub2$datasetid<-with(sub2, paste(studyid,"_",2))
sub3<-subset(clim3, studyid=="HMK019" & species=="Phytoplankton1 spp.")
sub3$datasetid<-with(sub3, paste(studyid,"_",3))
# for below EMW manually checked that all 3 remaining species (Perca fluviatillis, Daphnia3b spp., Daphnia3a spp., ) have same n years; then just picked one to subset on
sub4<-subset(clim3, studyid=="HMK019" & species=="Perca fluviatillis")
sub4$datasetid<-with(sub4, paste(studyid,"_",4)) 
sub5<-subset(clim3, studyid=="HMK038" & species=="Glis glis")
sub5$datasetid<-with(sub5, paste(studyid,"_",5))
# for below EMW manually checked that all 4 remaining species have same n years
sub6<-subset(clim3, studyid=="HMK038" & species=="Sitta europaea")
sub6$datasetid<-with(sub6, paste(studyid,"_",6))
sub7<-subset(clim3, studyid=="HMK051" & species=="Caterpillar4 spp.")
sub7$datasetid<-with(sub7, paste(studyid,"_",7))
sub8<-subset(clim3, studyid=="HMK051" & species=="Parus6 major")
sub8$datasetid<-with(sub8, paste(studyid,"_",8))
sub9<-subset(clim3, studyid=="HMK051" & species=="Parus3 caeruleus")
sub9$datasetid<-with(sub9, paste(studyid,"_",9))

# For all of the below I looked at all species, checked the n yrs were the same and...
# picked a species to subset on
sub10<-subset(clim3, studyid=="HMK023" & species=="Bombus spp.")
sub10$datasetid<-with(sub10, paste(studyid,"_",10))
sub11<-subset(clim3, studyid=="HMK031" & species=="Thermocyclops oithonoides")
sub11$datasetid<-with(sub11, paste(studyid,"_",11))
sub12<-subset(clim3, studyid=="HMK034" & species=="Cyclops vicinus")
sub12$datasetid<-with(sub12, paste(studyid,"_",12))
sub13<-subset(clim3, studyid=="HMK043" & species=="Beroe gracilis")
sub13$datasetid<-with(sub13, paste(studyid,"_",13))
sub14<-subset(clim3, studyid=="HMK016" & species=="Cerorhinca monocerata")
sub14$datasetid<-with(sub14, paste(studyid,"_",14))
sub15<-subset(clim3, studyid=="HMK018" & species=="Pygoscelis_b antarcticus")
sub15$datasetid<-with(sub15, paste(studyid,"_",15))
sub16<-subset(clim3, studyid=="HMK052" & species=="Pica pica")
sub16$datasetid<-with(sub16, paste(studyid,"_",16))
sub17<-subset(clim3, studyid=="HMK036" & species=="Copepod1 spp.")
sub17$datasetid<-with(sub17, paste(studyid,"_",17))
sub18<-subset(clim3, studyid=="HMK042" & species=="Acrocephalus arundinaceus")
sub18$datasetid<-with(sub18, paste(studyid,"_",18))

dataset<-rbind(sub1, sub2, sub3, sub4, sub5, sub6, sub7, sub8, sub9, sub10,
    sub11, sub12, sub13, sub14, sub15, sub16, sub17, sub18)

# Now! make a list of what datasetid each species belongs to, without dropping any species
# (note that we did drop species to make dataset)
# a translator so to speak!

datasetid.trans.full <- merge(dataset, clim3, by=c("species", "studyid"), all.y=TRUE, all.x=TRUE)
datasetid.trans <- subset(datasetid.trans.full, select=c("species", "studyid", "datasetid"))
datasetid.trans <- datasetid.trans[!duplicated(datasetid.trans),]

# get a list of the species x datasetids that we don't have above!
lookuptable <- c("Corydalis ambigua"="HMK023 _ 10", "Diatom2b spp."="HMK031 _ 11",
    "Diatom2a spp."="HMK031 _ 11", "Daphnia1 spp."="HMK031 _ 11",
    "Phytoplankton2 spp."="HMK034 _ 12", "Pleurobrachia pileus"="HMK043 _ 13",
    "Copepod2 spp."="HMK043 _ 13", "Pleurobrachia_a pileus"="HMK043 _ 13",
    "Engraulis japonicus"="HMK016 _ 14", "Pygoscelis adeliae"="HMK018 _ 15",
    "Pygoscelis_a antarcticus"="HMK018 _ 15", "Pygoscelis papua"="HMK018 _ 15",
    "Clamator glandarius"="HMK052 _ 16", "Diatom3 spp."="HMK036 _ 17",
    "Acrocephalus scirpaceus"="HMK042 _ 18", "Daphnia3b spp."="HMK019 _ 4",
    "Daphnia3a spp."="HMK019 _ 4", "Ficedula2 albicollis"="HMK038 _ 6", 
    "Parus caeruleus"="HMK038 _ 6", "Parus3 major"="HMK038 _ 6")

# now make a dataframe (not pretty, but works)
mergemein <- data.frame(species=row.names(as.data.frame(lookuptable)),
    datasetid.more=as.data.frame(lookuptable)$lookuptable)

datasetid.trans <- merge(datasetid.trans, mergemein, by=c("species"),
     all.y=TRUE, all.x=TRUE)

datasetid.trans$datasetid[is.na(datasetid.trans$datasetid)==TRUE] <-
    datasetid.trans$datasetid.more[is.na(datasetid.trans$datasetid)==TRUE]

datasetid.trans$datasetid.more <- NULL
