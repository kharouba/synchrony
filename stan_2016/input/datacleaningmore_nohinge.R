## Extracted from syncmodels.R ##
## on 27 Feb 2017 ##

#new changes to data (Oct 2016): No longer unique: Caterpillar 2a, 2b, 2c, 2d; Daphnia 3a, 3b; Paurs 2a, 2b caeruleus; Pyg. antarcticus a and b; Caterpillar4 now 4a and 4b
rawlong<-subset(rawlong, species!="asdf1" & species!="asdf2")
#asdf1= Diatom2b, asdf2= Thermocyclops oithonoides but first stage, in analysis now, last phenophase


#fix for Diatom4 spp. Dec 2016
sub<-subset(rawlong, species=="Diatom4 spp." & intid=="194");
sub[,c("species")]<-"Diatom4a spp."
rawlong<-rbind(rawlong, sub)
sub<-subset(rawlong, species=="Diatom4 spp." & intid=="195");
sub[,c("species")]<-"Diatom4b spp."
rawlong<-rbind(rawlong, sub)
sub<-subset(rawlong, species=="Diatom4 spp." & intid=="196");
sub[,c("species")]<-"Diatom4c spp."
rawlong<-rbind(rawlong, sub)
rawlong<-subset(rawlong, species!="Diatom4 spp.")




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
#rawlong.nodups$newid<-with(rawlong.nodups, paste(intid,"_",species))
#rawlong.nodups <- rawlong.tot[with(rawlong.nodups, order(newid)),]
#specieschar.formodel <- aggregate(rawlong.nodups["phenovalue"],
  #  rawlong.nodups[c("studyid", "species", "intid", "terrestrial","spp")], FUN=length) 
 #specieschar.formodel <- specieschar.formodel[with(specieschar.formodel, order(species)),]   

#number of years per species
    

#NO HINGE
rawlong.tot <- arrange(rawlong.tot, species)
rawlong.tot$yr1981 <- rawlong.tot$year-1981
rawlong.tot <- rawlong.tot[with(rawlong.tot, order(species)),]
rawlong.tot2<-unique(rawlong.tot[,c("studyid","species","phenovalue","year","yr1981")]) #CLEAN UP so only unique values across repeating species within studoes
N <- nrow(rawlong.tot2)
y <- rawlong.tot2$phenovalue
Nspp <- length(unique(rawlong.tot2$species)) #newid is character !
Nstudy<-length(unique(rawlong.tot2$studyid))
species <- as.numeric(as.factor(rawlong.tot2$species))
sock<-unique(rawlong.tot2[,c("studyid","species")])
studyid <- as.numeric(as.factor(sock$studyid))
year <- rawlong.tot2$yr1981

