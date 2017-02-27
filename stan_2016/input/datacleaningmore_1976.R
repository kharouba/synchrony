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
    
##

## Hinge models

#1976
# 3 groups of species: those with enough data <1976, those with some data but not enough <1976, and those with no data <1976
pre<-subset(rawlong.tot, year<=1976)
pre$count<-1
sss<- aggregate(pre["count"], pre[c("studyid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas

hinge<-merge(rawlong.tot, sss2, by=c("studyid", "species"))
hinge_pre<-subset(hinge, year<=1976); hinge_pre$newyear_76<-1976
hinge_post<-subset(hinge, year>1976); hinge_post$newyear_76<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)

hinge_non<-subset(sss, count<5) #datasets with NOT enough data pre climate change (n=21)
hinge_non$speciesid<-1:nrow(hinge_non) #number datas
hinge_non2<-merge(rawlong.tot, hinge_non, by=c("studyid", "species"))
hinge_non2$newyear_76<-hinge_non2$year
sums<-rbind(hinges, hinge_non2)
#names(sums)[6]<-"count"

letstry<-anti_join(rawlong.tot, sums)
letstry$newyear_76<-letstry$year

letstry2<-rbind(letstry,sums[,c("X","studyid","intid","repeatcode_study","figcode_study","int_type","spp","year","species","general","repeatcode","figcode","phenofreq","short_site","terrestrial","phenovalue","newyear","yr1981","newyear_76")])
rawlong.tot<-letstry2

#1976
rawlong.tot <- arrange(rawlong.tot, species)
rawlong.tot$yr1976 <- rawlong.tot$newyear-1976
rawlong.tot <- rawlong.tot[with(rawlong.tot, order(species)),]
rawlong.tot2<-unique(rawlong.tot[,c("studyid","species","phenovalue","yr1976","year")]) #CLEAN UP so only unique values across repeating species within studoes
N <- nrow(rawlong.tot2)
y <- rawlong.tot2$phenovalue
Nspp <- length(unique(rawlong.tot2$species)) #newid is character !
Nstudy<-length(unique(rawlong.tot2$studyid))
species <- as.numeric(as.factor(rawlong.tot2$species))
sock<-unique(rawlong.tot2[,c("studyid","species")])
studyid <- as.numeric(as.factor(sock$studyid))
year <- rawlong.tot2$yr1976
