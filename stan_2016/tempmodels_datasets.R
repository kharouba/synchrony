### to get absolute value of temperature change (i.e. based on unique datasets) ###

rm(list=ls())
#row.names=FALSE
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
# setwd("~/Documents/git/projects/trophsynch/synchrony/stan_2016")
library(ggplot2)
library(rstan)
library(shinystan)
library(grid)
library(nlme)
library(gridExtra)
library(plyr)
library(dplyr)
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

dataset<-rbind(sub1, sub2, sub3, sub4, sub5, sub6, sub7, sub8, sub9, sub10, sub11, sub12, sub13, sub14, sub15, sub16, sub17, sub18)

# check, which datasetids have more than one species
ddply(dataset, c("datasetid"), summarise,
    nspp=n_distinct(species))

#dataset<-unique(clim3[,c("studyid","year","envfactor","envunits","envtype","envvalue","newyear","yr1981")])
dataset <- dataset[with(dataset, order(datasetid)),]
N<-nrow(dataset)
y <- dataset$envvalue
specieschar.hin<- aggregate(dataset["envvalue"], dataset[c("datasetid")], FUN=length) 
specieschar.hin <- specieschar.hin[with(specieschar.hin, order(datasetid)),]
Nspp <- nrow(specieschar.hin) #J
#J <- nrow(specieschar.hin)
species <- as.numeric(as.factor(dataset$datasetid))
year <- dataset$yr1981

temp.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)


goo <- extract(temp.model)


#data <- subset(specieschar.hin, select=c("studyid"))
data <- unique(dataset[,c("studyid","datasetid","envtype")])
#look to see how stan is pooling, look at correlation between stan slope and lm slope
temp.change <- matrix(0, ncol=3000, nrow=Nspp) #2000 iterations for 21 interactions;
for (i in 3000:6000){ # 2000 iterations?
    #data$indiv_model <- 5.1
    #data$indiv_model<-goo$b[i,]
    temp.change[,(i-3000)] <- goo$b[i,]
    }
 
#temps<-t(temp.change)
data$stanfit <- rowMeans(temp.change, na.rm=TRUE) #mean across iterations for EACH dataset; STAN SLOPE ESTIMATE PER dataset
mean(data$stanfit)

sem<-sd(data$stanfit)/sqrt(length(data$stanfit)); sem
#95% confidence intervals of the mean
c(mean(data$stanfit)-2*sem,mean(data$stanfit)+2*sem)

# lizzie checks against some linear models
idshere <- unique(dataset$datasetid)
lmfits <- data.frame(datasetid=NA, lmfit=NA)

for (uniqdataset in seq_along(idshere)){
    subby <- subset(dataset, datasetid==idshere[uniqdataset])
    modelhere <- lm(envvalue~yr1981, data=subby)
    lmfits[uniqdataset,] <- data.frame(datasetit=idshere[uniqdataset], lmfit=coef(modelhere)[2])
}

compare.models <- cbind(lmfits, data$stanfit, idshere)

plot(lmfit~data$stanfit, data=compare.models)
abline(0,1)
abline(lm(lmfit~data$stanfit, data=compare.models))

# looking at few datasets

# start with one long-term one
# this looks pretty clear and strong and the agreement between stan and lm is good
subby <- subset(dataset, datasetid=="HMK019 _ 3")
plot(envvalue~yr1981, data=subby)
abline(lm(envvalue~yr1981, data=subby))

# one more long term one
# yeah, this model is totally working!
subby <- subset(dataset, datasetid=="HMK031 _ 11")
plot(envvalue~yr1981, data=subby)
abline(lm(envvalue~yr1981, data=subby))

# now look at one where a big increase get decreased by stan
# okay, I could see that getting smaller
subby <- subset(dataset, datasetid=="HMK011 _ 1")
plot(envvalue~yr1981, data=subby)
abline(lm(envvalue~yr1981, data=subby))

# now look at one where decrease gets pulled to the mean (positive change) by stan
# right
subby <- subset(dataset, datasetid=="HMK052 _ 16")
plot(envvalue~yr1981, data=subby)
abline(lm(envvalue~yr1981, data=subby))


#by group
sub<-subset(data, envtype=="air"); mean(sub$stanfit)
sub<-subset(data, envtype=="water"); mean(sub$stanfit)
sem<-with(sub, sd(stanfit)/sqrt(length(stanfit)));
#95% confidence intervals of the mean
with(sub, c(mean(stanfit)-2*sem,mean(stanfit)+2*sem))

news<-t(temp.change)
mdata <- melt(news)
names(mdata)[1]<-"iteration"; names(mdata)[2]<-"studyid"

#New model as of June 2016
#Random slopes only, no random intercepts, hinge, no covariate matrix:
temp.model<-stan("twolevelrandomslope.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

#New model as of April 2016
#Random slopes only, no random intercepts, hinge, no covariate matrix:
#temp.model<-stan("synchrony_apr_nocov.stan", data=c("N","J","y","species","year"), iter=2000, chains=4)
#temp.model<-stan("twolevelrandomeffects.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)



#RESULTS!!!
plot(temp.model, pars=c("b"))
new_y<-extract(temp.model, pars="ypred")
pred<-apply(new_y[[1]],2, quantile, probs=c(0.025, 0.975))
with(dataset, plot(env)


par(mfrow=(c(1,2)))
#95% posterior interval:
print(temp.model, "b", probs=c(.025,.975))


##

#########
#Plot datasets
asdf<-summary(temp.model, pars="b")
#to get median coefficients from SUMMARY
median<-asdf[[1]][1:13]; new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) #-0.715645239
min<-asdf[[1]][40:52]; e<-data.frame(min=unlist(min)) #-1.1222299
max<-asdf[[1]][92:104]; f<-data.frame(max=unlist(max)) #-0.316333047
d$min<-e$min; d$max<-f$max;
details<-unique(clim3[,c("studyid","envtype")])
details <- details[with(details, order(studyid)),]
d$envtype<-details$envtype

ggplot(d, aes(x=factor(grp), y=y, ymin=min, ymax=max, colour=envtype))+geom_pointrange(aes(width=.05, colour=factor(envtype)))+theme_bw()+geom_point(size=3, aes(colour=factor(envtype)))+xlab("datasets")+ylab("Temperature change (C/year)")+geom_hline(xintercept=0, linetype="longdash", colour="grey")+theme(legend.position="false")+facet_wrap(~envtype, scale="free")+ylim(-0.1, 0.2)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Exploratory figures
ggplot(subset(clim3, envtype=="air"), aes(x=year, y=envvalue))+
geom_point(aes(colour=factor(species)))+xlim(1969, 2013)+geom_vline(xintercept=1981, colour="grey", linetype = "longdash")+theme_bw()#+geom_smooth(method="lm", se=FALSE, aes(colour = factor(species)))

#by species
ggplot(subset(clim3, envtype=="air"), aes(x=year, y=envvalue))+
geom_point()+xlim(1969, 2013)+geom_vline(xintercept=1981, colour="grey", linetype = "longdash")+facet_wrap(~species)+theme_bw()

ggplot(subset(clim3, envtype=="water"), aes(x=year, y=envvalue))+
geom_point()+xlim(1969, 2013)+geom_vline(xintercept=1981, colour="grey", linetype = "longdash")+facet_wrap(~species)+theme_bw()

#by unique weather station
#air- 
ggplot(subset(clim3, envtype=="air"), aes(x=year, y=envvalue))+
geom_point(aes(colour=factor(studyid)))+facet_wrap(~sitecode)
#water- by unique weather station
ggplot(subset(clim3, envtype=="water"), aes(x=year, y=envvalue))+
geom_point(aes(colour=factor(species)))+facet_wrap(~sitecode)
#geom_path(data=yo, aes(x=SiteTypeBinary, y=log(std), group=factor(sppid), colour=factor(sppid)), size=1)+
#geom_point(data=yo, aes(x=SiteTypeBinary, y=log(std), colour=factor(sppid)))

potential hinge
air:
Acrocephalus arundinaceus (207)
Acrocephalus scirpaceus (207)
water:
Copepod1 spp. (193)
Copepod2 spp. (208)
Cyclops vicinus (191)
Daphnia3 spp.
Diatom3 spp.
Perca fluviatillis
Phytoplankton1 spp.
Pleurobrachia pileus
Pleurobrachia_a pileus
