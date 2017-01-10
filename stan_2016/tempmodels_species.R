# To match individual species doy~year with temp~year

rm(list=ls()) 

library(reshape)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")

#load rawlong.tot from syncmodels.

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
studyid <- as.numeric(as.factor(clim3$studyid))
#sock<-unique(clim3[,c("studyid","species")])
#studyid <- as.numeric(as.factor(sock$studyid))
year <- clim3$yr1981

#to calculate temp sensitivity
#y<- clim3$phenovalue
#year<- clim3$envvalue

#temp.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=4000, chains=4)
temp.model<-stan("stanmodels/threelevelrandomslope3.stan", data=c("N","Nspp","Nstudy","y","species","studyid","year"), iter=4000, chains=4)


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

#alt
goo <- extract(temp.model)
summ_studyspp<-unique(clim3[,c("studyid","species")])
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- goo$b_spp[i,]
   it1000[,(i-3000)] <- goo$b_spp[i,]
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


# STEP 4- run stan on doy~year

N<-nrow(clim3)
y <- clim3$phenovalue
Nspp <- length(unique(clim3$species)) #J
species <- as.numeric(as.factor(clim3$species))
year <- clim3$yr1981

#pheno.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)
pheno.model<-stan("stanmodels/threelevelrandomslope3.stan", data=c("N","Nspp","Nstudy","y","species","studyid","year"), iter=4000, chains=4)


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
    solo$model <- faa$b_spp[i,]
    it1000[,(i-3000)] <- faa$b_spp[i,]
}
solo$phenochange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP
mean(rowMeans(it1000, na.rm=TRUE))
data<-merge(summ_studyspp, solo[,c("studyid","species","phenochange")],by=c("studyid","species"))


phenochange<-it1000
mdata<-melt(it1000[,2001:3000])
mdata<-it1000[,2001:2100]
names(mdata)[1]<-"id"; names(mdata)[2]<-"iteration"; names(mdata)[3]<-"pheno.change"; 
tdata<-merge(ndata, mdata, by=c("id","iteration"))

solo$id<-1:nrow(solo);
tdata2<-merge(tdata, solo[,c("studyid","species","id")], by=c("id"))
N<-nrow(tdata2)
y <- abs(tdata2$pheno.change) #absolute value of pheno change!
Nspp <- length(unique(tdata2$id)) #J
species <- as.numeric(as.factor(tdata2$id))
studyid <- as.numeric(as.factor(tdata2$studyid))
year <- abs(tdata2$temp.change) #absolute value of temp change


#cov.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)
cov.model<-stan("stanmodels/threelevelrandomeffects4.stan", data=c("N","Nspp","Nstudy","y","species","studyid","year"), iter=4000, chains=4)

N<-nrow(tdata2)
y <- abs(tdata2$pheno.change) #absolute value of pheno change!
J <- length(unique(tdata2$id)) #J
S <- length(unique(tdata2$studyid))
plotnum <- as.numeric(as.factor(tdata2$id))
sock<-unique(tdata2[,c("studyid","species")])
sitenum <- as.numeric(as.factor(sock$studyid))
#sitenum <- as.numeric(as.factor(tdata2$studyid))
x <- abs(tdata2$temp.change) #absolute value of temp change


#dat <- list(N=N,S=S,J=J,plotnum=plotnum, sitenum=sitenum,y=y,x=x) 
fitme <- stan("stanmodels/threelevel_plotsinsites.stan", data=c("N","J","S","plotnum", "sitenum","y","x"), iter=3000, chains=4, control=list(adapt_delta = 0.9, stepsize = 0.5))


fun<-extract(cov.model)
sis<-as.data.frame(unique(tdata[,c("id")]))

it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    sis$model <- fun$b[i,]
    it1000[,(i-3000)] <- fun$b[i,]
}
mean(sis$model)
#computation of the standard error of the mean
sem<-sd(sis$model)/sqrt(length(sis$model)); sem
#95% confidence intervals of the mean
c(mean(sis$model)-2*sem,mean(sis$model)+2*sem)



*** check T oithonoides, [stan not matching paper (stan does equal lm though)]- CHECKED- DOES MATCH PAPER BECAUSE CHOSE LAST PHENOPHASE OF SEASON (NOT START) AND LAST STAGE HAS GOTTEN LATER- SHOULD FIX FOR FINAL PAPER
** UPDATE- CHANGED OITHONOIDES SO THAT NOW FIRST PHENOPHASE (DEC 2016)

*** SPP30, 31 ARE OUTLIERS: Pleurobrachia pileus and Pleurobrachia_a pileus from HMK043

ggplot(data, aes(y=abs(phenochange),x=tempchange))+geom_errorbar(aes(ymin=abs(pheno_min), ymax=abs(pheno_max)), colour="red")+geom_errorbarh(aes(xmin=temp_min, xmax=temp_max, height = .1), colour="red")+geom_point(size=2)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")

OR

ggplot(data, aes(y=abs(phenochange),x=tempchange))+geom_point(size=3)+theme(text = element_text(size=20), axis.text.x=element_text(size=20), axis.title.y=element_text(size=20, angle=90))+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")

yep<-unique(clim3[,c("studyid", "species", "terrestrial")])
yep2<-merge(data, yep, by=c("studyid","species"))

ggplot(subset(yep2, terrestrial=="terrestrial"), aes(y=abs(phenochange),x=tempchange))+geom_point(size=3)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")

