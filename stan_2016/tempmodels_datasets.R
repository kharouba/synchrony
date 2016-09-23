### to get absolute value of temperature change (i.e. based on unique datasets) ###
### Note that as of 23 Sep 2016 this file is super similar to tempmodels_interactions.R ###
### Uniquely, however, this file looks at the model output a little to make sure  #
## ...that the models appear to be working. ###
### But the order that the data go into the models in each file is different... ###
### In tempmodels_interactions.R I check that we're getting similar answers between the two code files though ###

rm(list=ls())
options(stringsAsFactors=FALSE)

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
# set_cppo("fast") 

## source the temp dataset cleaning code

source("source/clean_tempdatasets.R")

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

write.csv(compare.models, "output/compare.tempmodels.csv", row.names=FALSE)

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

## OLDER code below, Lizzie did not work on or check this ... 


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
