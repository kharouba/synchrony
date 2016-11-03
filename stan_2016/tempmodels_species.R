# To match individual species doy~year with temp~year

rm(list=ls()) 

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")


#step 1- sync climate and pheno data
clim<-read.csv("input/climate4.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
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

#load rawlong.tot from syncmodels.

clim3<-merge(clim2, unique(rawlong.tot[,c("studyid","species","year","phenovalue","terrestrial")]), by=c("studyid","species","year"))


#STEP 2- hinge model

#Hinge (11 spp)
hinge<-subset(clim3, species=="Acrocephalus arundinaceus" | species=="Acrocephalus scirpaceus" | species=="Copepod1 spp." | species=="Copepod2 spp." | species=="Cyclops vicinus"  | species=="Daphnia3 spp." | species=="Diatom3 spp." | species=="Perca fluviatillis" | species=="Phytoplankton1 spp." | species=="Pleurobrachia pileus" | species=="Pleurobrachia_a pileus")

#Non-hinge (11 spp)
hinge_non<-subset(clim3, species!="Acrocephalus arundinaceus" & species!="Acrocephalus scirpaceus"  & species!="Copepod1 spp." & species!="Copepod2 spp." & species!="Cyclops vicinus"  & species!="Daphnia3 spp." & species!="Diatom3 spp." & species!="Perca fluviatillis" & species!="Phytoplankton1 spp." & species!="Pleurobrachia pileus" & species!="Pleurobrachia_a pileus")

hinge_non$newyear<-hinge_non$year
hinge_pre<-subset(hinge, year<=1981); hinge_pre$newyear<-1981
hinge_post<-subset(hinge, year>1981); hinge_post$newyear<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)

clim4<-rbind(hinge_non, hinges);
clim3<-clim4
clim3$yr1981 <- clim3$newyear-1981


# STEP 3- Run stan on temp first
clim3 <- clim3[with(clim3, order(species, year)),]
clim3 <- na.omit(clim3)
N<-nrow(clim3)
y <- clim3$envvalue
Nspp <- length(unique(clim3$species)) #J
species <- as.numeric(as.factor(clim3$species))
year <- clim3$yr1981

temp.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

goo <- extract(temp.model)
asdf<-summary(temp.model)
print(temp.model, pars=c("mu_a", "mu_b", "sigma_y", "sigma_a","sigma_b"))

asdf<-summary(temp.model, pars="b")
#to get median coefficients from SUMMARY
median<-asdf[[1]][1:37]; #new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) 
min<-asdf[[1]][112:148]; e<-data.frame(min=unlist(min)) #-1.129 to -1.022
max<-asdf[[1]][260:296]; f<-data.frame(max=unlist(max)) #-0.3077 to 0.14
d$min<-e$min; d$max<-f$max;

summ_studyspp<-d;  names(summ_studyspp)[1]<-"tempchange"; names(summ_studyspp)[3]<-"temp_min"; names(summ_studyspp)[4]<-"temp_max"

summ_studyspp<-unique(clim3[,c("studyid","species")])

it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- goo$b[i,]
    it1000[,(i-3000)] <- goo$b[i,]
}
summ_studyspp$tempchange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP

# STEP 4- run stan on doy~year

N<-nrow(clim3)
y <- clim3$phenovalue
Nspp <- length(unique(clim3$species)) #J
species <- as.numeric(as.factor(clim3$species))
year <- clim3$yr1981

pheno.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

asdf<-summary(pheno.model, pars="b")
#to get median coefficients from SUMMARY
median<-asdf[[1]][1:37]; #new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) 
min<-asdf[[1]][112:148]; e<-data.frame(min=unlist(min)) #-1.129 to -1.022
max<-asdf[[1]][260:296]; f<-data.frame(max=unlist(max)) #-0.3077 to 0.14
d$min<-e$min; d$max<-f$max;

solo<-d; names(solo)[1]<-"phenochange"; names(solo)[3]<-"pheno_min"; names(solo)[4]<-"pheno_max"

data<-merge(summ_studyspp, solo, by=c("grp"))

faa <- extract(pheno.model)
solo<-unique(clim3[,c("studyid","species")])

it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    solo$model <- faa$b[i,]
    it1000[,(i-3000)] <- faa$b[i,]
}
solo$phenochange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP

data<-merge(summ_studyspp, solo[,c("studyid","species","phenochange")], by=c("studyid","species"))

*** check T oithonoides, stan not matching paper (stan does equal lm though)

ggplot(data, aes(y=phenochange,x=tempchange))+geom_errorbar(aes(ymin=pheno_min, ymax=pheno_max,width=.0025), colour="red")+geom_errorbarh(aes(xmin=temp_min, xmax=temp_max, height = .1), width=.0025,colour="red")+geom_point(size=2)
