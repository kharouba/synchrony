rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")
library(ggplot2)
library(rstan)
library(shinyStan)
# library(lme4)

# in the data file: spp1 = neg= negative species e.g. resource 
raw <- read.csv("raw_aug.csv", header=TRUE)
rawlong <- read.csv("indiv_sppdata_aug.csv", header=TRUE)
rawlong <- rawlong[with(rawlong, order(intid)),]
rawlong$int_type[which(rawlong$int_type=="poll")] <- "pollination"
rawlong.X <- NULL
change<-subset(rawlong, species=="Parus2 caeruleus" & intid=="216")
change$species[which(change$species=="Parus2 caeruleus")] <- "Parus2a caeruleus"

change2<-subset(rawlong, species=="Parus2 caeruleus" & intid=="220")
change2$species[which(change2$species=="Parus2 caeruleus")] <- "Parus2b caeruleus"

ori<-subset(rawlong, species!="Parus2 caeruleus")
ori2<-rbind(ori, change, change2)
rawlong<-ori2

change<-subset(rawlong, species=="Parus4 major" & intid=="221")
change$species[which(change$species=="Parus4 major")] <- "Parus4a major"

change2<-subset(rawlong, species=="Parus4 major" & intid=="217")
change2$species[which(change2$species=="Parus4 major")] <- "Parus4b major"

ori<-subset(rawlong, species!="Parus4 major")
ori2<-rbind(ori, change, change2)
rawlong<-ori2

#notes
#caterpillar spp.
#int177: 1961-2004
#int178: 1962-2006
#int179: 1961-2007
#for int209- Pleurobrachia pileus: 1976-2004 vs int208- 1976-2003
#Ficedula hypoleuca

# changed spp name because slightly different time series than other interaction in ms, now matches its interacting spp.
rawlong$rowid<-1:nrow(rawlong)
change<-subset(rawlong, species=="Caterpillar spp." & intid=="177")
change$species[which(change$species=="Caterpillar spp.")] <- "Caterpillar_a spp."
change$rowid<-2236:(2236+nrow(change)-1)
ori<-subset(rawlong, rowid<479 | rowid>514)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar spp." & intid=="178")
change$species[which(change$species=="Caterpillar spp.")] <- "Caterpillar_b spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<551 | rowid>584)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar spp." & intid=="179")
change$species[which(change$species=="Caterpillar spp.")] <- "Caterpillar_c spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<653 | rowid>686)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Pleurobrachia pileus" & intid=="209")
change$species[which(change$species=="Pleurobrachia pileus")] <- "Pleurobrachia_a pileus"
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1623 | rowid>1649)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Ficedula hypoleuca" & intid=="214")
change$species[which(change$species=="Ficedula hypoleuca")] <- "Ficedula_a hypoleuca"
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1677 | rowid>1692)
ori2<-rbind(ori, change)
rawlong<-ori2

#not different length of time series but easier to fix interactions HERE
#for int215- Parus ater
change<-subset(rawlong, species=="Parus ater" & intid=="215")
change$species[which(change$species=="Parus ater")] <- "Parus_a ater"
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1709 | rowid>1724)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar2 spp." & intid=="218")
change$species[which(change$species=="Caterpillar2 spp.")] <- "Caterpillar2a spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1773 | rowid>1791)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar2 spp." & intid=="219")
change$species[which(change$species=="Caterpillar2 spp.")] <- "Caterpillar2b spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1811 | rowid>1830)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar2 spp." & intid=="220")
change$species[which(change$species=="Caterpillar2 spp.")] <- "Caterpillar2c spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1851 | rowid>1870)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar2 spp." & intid=="221")
change$species[which(change$species=="Caterpillar2 spp.")] <- "Caterpillar2d spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1871 | rowid>1890)
ori2<-rbind(ori, change)
rawlong<-ori2

#check interactions to see if there are still 2 spp per interaction
#rename Parus4 major of int221

sub<-subset(rawlong, studyid=="HMK018")
sub$phenovalue2<-with(sub, phenovalue+182) #number of days different between Dec 22 (SH summer solstice) and June 21 (NH summer solstice), to put all dates on same 'calendar'
sub$phenovalue3<-with(sub, phenovalue2-365) #now same number of days before NH summer solstice as it was before SH summer solstice

rawlong2<-subset(rawlong, studyid!="HMK018")
rawlong2$phenovalue2<-rawlong2$phenovalue
rawlong2$phenovalue3<-rawlong2$phenovalue
rawlong3<-rbind(rawlong2, sub)
rawlong<-rawlong3
rawlong<-rawlong[,c("studyid","intid","int_type","spp","year","species","phenofreq","short_site","terrestrial","phenovalue3")]
names(rawlong)[10]<-"phenovalue"
