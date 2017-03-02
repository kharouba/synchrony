### UPDATES FROM FEB 2017 WHEN CHANGES MADE TO MODEL TYPE 
## COVARIATE MODELS 
# To match individual species doy~year with temp~year

rm(list=ls()) 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")

# libraries
library(ggplot2)
library(rstan)
library(shinystan)
library(grid)
library(nlme)
library(dplyr)
library(ggrepel)
library(reshape)
set_cppo("fast")  # for best running speed
source("/users/kharouba/google drive/UBC/multiplot.R")

#load rawlong.tot from syncmodels.
#get data
rawlong <- read.csv("input/rawlong2.csv", header=TRUE)
# use for 1981 hinge model
source("input/datacleaningmore.R")

#step 1- sync climate and pheno data
clim<-read.csv("input/climate4.csv", header=TRUE, na.strings="<NA>", as.is=TRUE) #updated it to include raw data for HMK043 and NOT interannual difference
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


!!!!! Option: Choose 1 species/climate per study
Bgroups<-unique(clim3$studyid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0,c(nrow(clim3),ncol(clim3))))
rowcount<-1
for(i in 1:length(b)){
int_long<-subset(clim3, studyid==b[i])
Cgroups<-unique(int_long$species); c<-Cgroups; c<-as.character(c)
x<-c(1:length(c))
me<-sample(x, size=1)

spp<-subset(int_long, species==c[me])
asdf<-rowcount+(nrow(spp)-1)

new[rowcount:asdf,]<-spp
rowcount<-rowcount+nrow(spp)
asdf<-rowcount+(nrow(spp)-1)
}
names(new)<-names(clim3)
new<-subset(new, studyid!=0)
clim3<-new

# STEP 3- Run stan on temp first
clim3 <- clim3[with(clim3, order(species, year)),]
clim3 <- na.omit(clim3)
N<-nrow(clim3)
y <- clim3$envvalue
Nspp <- length(unique(clim3$species)) #J
Nstudy<-length(unique(clim3$studyid))
species <- as.numeric(as.factor(clim3$species))
#studyid <- as.numeric(as.factor(clim3$studyid))
sock<-unique(clim3[,c("studyid","species")])
studyid <- as.numeric(as.factor(sock$studyid))
year <- clim3$yr1981

#Feb 2017: don't have enough repeating species ACROSS studies to esimate error, therefore no study as grouping
temp.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(temp.model, pars = c("mu_b", "sigma_y", "a", "b"))

#to calculate temp sensitivity
y<- clim3$phenovalue
year<- clim3$envvalue

#Feb 2017: don't have enough repeating species ACROSS studies to esimate error, therefore no study as grouping
tempsens.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(tempsens.model, pars = c("mu_b", "sigma_y", "a", "b"))


uni<-unique(clim3[c("studyid","species")])

asdf<-summary(temp.model, pars="b_spp")
#to get median coefficients from SUMMARY
median<-asdf[[1]][1:37]; #new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) 
min<-asdf[[1]][112:148]; e<-data.frame(min=unlist(min)) #-1.129 to -1.022
max<-asdf[[1]][260:296]; f<-data.frame(max=unlist(max)) #-0.3077 to 0.14
d$min<-e$min; d$max<-f$max;

summ_studyspp<-d;  names(summ_studyspp)[1]<-"tempchange"; names(summ_studyspp)[3]<-"temp_min"; names(summ_studyspp)[4]<-"temp_max"

#dataset length
clim3$count<-1
tool<-with(clim3, aggregate(count, by=list(studyid, species), FUN=sum, na.rm=T)) 
names(tool)[1]<-"studyid"; names(tool)[2]<-"species"; names(tool)[3]<-"length"

d2<-cbind(d, tool)
m1<-lm(y~length, d2); summary(m1)

##alt
# for interpretation:
# for mu_b
fas<-extract(temp.model)
it1000 <- matrix(0, ncol=3000, nrow=1)
for (i in 3000:6000){ # 3000 iterations?
        it1000[,(i-3000)] <- fas$mu_b[i]
}
mean(rowMeans(it1000, na.rm=TRUE))
sem<-sd(it1000)/sqrt(length(it1000)); sem
#95% confidence intervals of the mean
c(mean(it1000)-2*sem,mean(it1000)+2*sem)


#for covariate model:
goo <- extract(temp.model)
summ_studyspp<-unique(clim3[,c("studyid","species")])
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- goo$b[i,]
   it1000[,(i-3000)] <- goo$b[i,]
}
tempchange<-it1000
ndata<-melt(it1000[,2001:3000])
names(ndata)[1]<-"id"; names(ndata)[2]<-"iteration"; names(ndata)[3]<-"temp.change"; 
#ndata$intid<-rep(intid.nodups$intid, 1000)


summ_studyspp$tempchange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP
mean(summ_studyspp$tempchange)
#computation of the standard error of the mean
sem<-sd(summ_studyspp$tempchange)/sqrt(length(summ_studyspp$tempchange)); sem
#95% confidence intervals of the mean
c(mean(summ_studyspp$tempchange)-2*sem,mean(summ_studyspp$tempchange)+2*sem)

#temp sens figure
asdf<-summary(temp.model, pars=c("a_spp","b_spp"))
lab<-as.data.frame(asdf[[1]][1:37]); names(lab)[1]<-"int"
lab$slope<-asdf[[1]][38:74]
lab$id<-unique(species)

me<-summary(temp.model, pars=c("a_study"))
mean(me[[1]][1:13]) #rough guess on mu_a

ggplot(clim3, aes(y=phenovalue,x=envvalue, colour=factor(species)))+geom_point()+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("Phenology (doy)")+theme(legend.position="none")+xlab(expression(paste("Temperature ",degree,"C")))+geom_abline(aes(intercept=int, slope=slope, colour=factor(id)), data=lab)+geom_abline(slope=-3.19, intercept=132, size=1.5)


#temp sens figure
asdf<-summary(temp.model, pars=c("a_spp","b_spp"))
lab<-as.data.frame(asdf[[1]][1:37]); names(lab)[1]<-"int"
lab$slope<-asdf[[1]][38:74]
lab$id<-unique(species)

me<-summary(temp.model, pars=c("a_study"))
mean(me[[1]][1:13]) #rough guess on mu_a

ggplot(clim3, aes(y=phenovalue,x=envvalue, colour=factor(species)))+geom_point()+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("Phenology (doy)")+theme(legend.position="none")+xlab(expression(paste("Temperature ",degree,"C")))+geom_abline(aes(intercept=int, slope=slope, colour=factor(id)), data=lab)+geom_abline(slope=-3.19, intercept=132, size=1.5)

# STEP 4- run stan on doy~year

N<-nrow(clim3)
y <- clim3$phenovalue
Nspp <- length(unique(clim3$species)) #J
species <- as.numeric(as.factor(clim3$species))
year <- clim3$yr1981

#Feb 2017: don't have enough repeating species ACROSS studies to esimate error, therefore no study as group
pheno.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(pheno.model, pars = c("mu_b", "sigma_y", "a", "b"))


asdf<-summary(pheno.model, pars="b_spp")
#to get median coefficients from SUMMARY
median<-asdf[[1]][1:37]; #new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) 
min<-asdf[[1]][112:148]; e<-data.frame(min=unlist(min)) #-1.129 to -1.022
max<-asdf[[1]][260:296]; f<-data.frame(max=unlist(max)) #-0.3077 to 0.14
d$min<-e$min; d$max<-f$max;

#dataset length
clim3$count<-1
tool<-with(clim3, aggregate(count, by=list(studyid, species), FUN=sum, na.rm=T)) 
names(tool)[1]<-"studyid"; names(tool)[2]<-"species"; names(tool)[3]<-"length"

d2<-cbind(d, tool)
m1<-lm(y~length, d2); summary(m1)


solo<-d; names(solo)[1]<-"phenochange"; names(solo)[3]<-"pheno_min"; names(solo)[4]<-"pheno_max"

data<-merge(summ_studyspp, solo, by=c("grp"))
data<-cbind(data, uni)

#non stan approach
m1<-lme(abs(phenochange)~abs(tempchange), random=~1|studyid, data=data); summary(m1)

#alt
faa <- extract(pheno.model)
solo<-unique(clim3[,c("studyid","species")])

it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    solo$model <- faa$b[i,]
    it1000[,(i-3000)] <- faa$b[i,]
}
solo$phenochange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP
mean(rowMeans(it1000, na.rm=TRUE))
data<-merge(summ_studyspp, solo[,c("studyid","species","phenochange")],by=c("studyid","species"))


phenochange<-it1000
mdata<-melt(it1000[,2001:3000])
#mdata<-melt(it1000[,2001:2100]) #test
names(mdata)[1]<-"id"; names(mdata)[2]<-"iteration"; names(mdata)[3]<-"pheno.change"; 
tdata<-merge(ndata, mdata, by=c("id","iteration"))

solo$id<-1:nrow(solo);
tdata2<-merge(tdata, solo[,c("studyid","species","id")], by=c("id"))
N<-nrow(tdata2)
y <- abs(tdata2$pheno.change) #absolute value of pheno change!
Nspp <- length(unique(tdata2$id)) #J
species <- as.numeric(as.factor(tdata2$id))
#studyid <- as.numeric(as.factor(tdata2$studyid))
year <- abs(tdata2$temp.change) #absolute value of temp change

# Feb 2017 decisions: (1) Not enough variation within species or studies therefore do not pool slopes; (2) Not enough repeating species ACROSS studies to justify including study as grouping, therefore two level random intercept model.
cov.model<-stan("stanmodels/twolevelrandomintercept2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(pheno.model, pars = c("mu_b", "sigma_y", "a", "b"))


print(cov.model)

What's important for posterior expectations is the effective
sample size (n_eff in Stan's output), because that's what
determines the MCMC standard error on the estimate.'
The simulation with higher effective sample size has a lower standard error of the mean and more stable estimates. 


#non stan approach
m1<-lme(abs(pheno.change)~abs(temp.change), random=~1|studyid/species, data=tdata2); summary(m1)

fun<-extract(cov.model)
sis<-as.data.frame(unique(tdata2[,c("id")]))

# for mu_b
it1000 <- matrix(0, ncol=3000, nrow=1)
for (i in 3000:6000){ # 3000 iterations?
    #sis$model <- fun$mu_b[i]
    it1000[,(i-3000)] <- fun$mu_b[i]
}
mean(rowMeans(it1000, na.rm=TRUE))
sem<-sd(it1000)/sqrt(length(it1000)); sem
#95% confidence intervals of the mean
c(mean(it1000)-2*sem,mean(it1000)+2*sem)

0.058, 0.054-0.062

# for mu_a
it2000 <- matrix(0, ncol=3000, nrow=1)
for (i in 3000:6000){ # 3000 iterations?
    #sis$model <- fun$mu_a[i]
    it2000[,(i-3000)] <- fun$mu_a[i]
}
mean(rowMeans(it2000, na.rm=TRUE))


#Figure-Relationship between temperature change and phenological change show 95%CI around (slope) estimates of phenological and temperature (i.e. b_spp), black line is fit with mu_a and mu_b
#load data from above
ggplot(data, aes(y=abs(phenochange),x=abs(tempchange)))+geom_errorbar(aes(ymin=abs(pheno_min), ymax=abs(pheno_max)), colour="grey")+geom_errorbarh(aes(xmin=temp_min, xmax=temp_max, height = .1), colour="grey")+geom_point(size=2)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab(expression(paste("Temperature change (",degree,"C/year)")))+geom_abline(slope=0.058, intercept=0.71, size=1.5)

#for display
study <- aggregate(tdata2["pheno.change"], tdata2["studyid"], FUN=mean); names(study)[2]<-"pheno"
study2 <- aggregate(tdata2["temp.change"], tdata2["studyid"], FUN=mean); names(study2)[2]<-"temp"
studs<-merge(study,study2, by=c("studyid"))
#by study
asdf<-summary(cov.model, pars=c("a_study","b_study"))
lab<-as.data.frame(asdf[[1]][1:13]); names(lab)[1]<-"int"
lab$slope<-asdf[[1]][14:26]
lab$id<-unique(tdata2$studyid)
ggplot(subset(tdata2, iteration==1), aes(y=abs(pheno.change),x=abs(temp.change), colour=factor(studyid)))+geom_point(size=2)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")+geom_abline(aes(intercept=int, slope=slope, colour=factor(id)), data=lab)+geom_abline(slope=0.07, intercept=0.71, size=2)

+geom_abline(intercept=mean(rowMeans(it2000, na.rm=TRUE)), slope=mean(rowMeans(it1000, na.rm=TRUE)), colour="red")

#one study
ggplot(subset(tdata2, iteration==1 & studyid=="HMK038"), aes(y=abs(pheno.change),x=abs(temp.change), colour=factor(studyid)))+geom_point(size=2)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")+geom_smooth(method="lm", se=FALSE)+geom_abline(slope=0.11, intercept=0.687)

#for b_spp
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    sis$model <- fun$b_spp[i,]
    it1000[,(i-3000)] <- fun$b_spp[i,]
}
slopes<-melt(it1000[,2001:3000]) #for figure
names(slopes)[1]<-"id"; names(slopes)[2]<-"iteration"; names(slopes)[3]<-"slope"; #for figure

sis$pheno.temp <- rowMeans(it1000, na.rm=TRUE); mean(rowMeans(it1000, na.rm=TRUE)) #mean across iterations for EACH SPP
#computation of the standard error of the mean
sem<-sd(sis$pheno.temp)/sqrt(length(sis$pheno.temp)); sem
#95% confidence intervals of the mean
c(mean(sis$pheno.temp)-2*sem,mean(sis$pheno.temp)+2*sem)

new_y<-extract(cov.model, pars="ypred")

*** check T oithonoides, [stan not matching paper (stan does equal lm though)]- CHECKED- DOES MATCH PAPER BECAUSE CHOSE LAST PHENOPHASE OF SEASON (NOT START) AND LAST STAGE HAS GOTTEN LATER- SHOULD FIX FOR FINAL PAPER
** UPDATE- CHANGED OITHONOIDES SO THAT NOW FIRST PHENOPHASE (DEC 2016)

*** SPP30, 31 ARE OUTLIERS: Pleurobrachia pileus and Pleurobrachia_a pileus from HMK043

ggplot(data, aes(y=abs(phenochange),x=tempchange))+geom_errorbar(aes(ymin=abs(pheno_min), ymax=abs(pheno_max)), colour="red")+geom_errorbarh(aes(xmin=temp_min, xmax=temp_max, height = .1), colour="red")+geom_point(size=2)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")

OR

ggplot(data, aes(y=abs(phenochange),x=tempchange))+geom_point(size=3)+theme(text = element_text(size=20), axis.text.x=element_text(size=20), axis.title.y=element_text(size=20, angle=90))+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")

yep<-unique(clim3[,c("studyid", "species", "terrestrial")])
yep2<-merge(data, yep, by=c("studyid","species"))

ggplot(subset(yep2, terrestrial=="terrestrial"), aes(y=abs(phenochange),x=tempchange))+geom_point(size=3)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")

