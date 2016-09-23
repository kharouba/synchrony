rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(ggplot2)
library(rstan)
library(shinyStan)
library(grid)
library(nlme)
library(reshape)
library(dplyr)
set_cppo("fast")  # for best running speed
source("/users/kharouba/google drive/UBC/multiplot.R")
#library(reshape)
# library(lme4)

#Step 1: Synchrony data, use all 3000 iterations, no means taken per interaction
#use synch.model loaded from syncmodels.R

### need to match with the right intid!!!

ndata<-melt(synch.model[,2001:3000])
names(ndata)[1]<-"id"; names(ndata)[2]<-"iteration"; names(ndata)[3]<-"sync.change"; 
ndata$intid<-rep(intid.nodups$intid, 1000) # 1,19,144,145 x3000

interact <- read.csv("input/raw_april.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue) #spp1-spp2
#positive phenodiff= consumer emerges BEFORE resource
#negative phenodiff= resource emerges BEFORE consumer
season <- aggregate(interact["newphenodiff"], interact[c("studyid", "intid")], FUN=mean)

sync<-merge(ndata, season, by=c("intid"))
#increasing synchrony (i.e. species are getting closer together)
inc1<-subset(sync, newphenodiff>0 & sync.change<0)
inc2<-subset(sync, newphenodiff<0 & sync.change>0)
inc<-rbind(inc1, inc2)
inc$sync.change<-abs(inc$sync.change)
# decreasing synchrony (i.e. speceis are getting farther apart)
dec1<-subset(sync, newphenodiff>0 & sync.change>0)
dec2<-subset(sync, newphenodiff<0 & sync.change<0)
dec<-rbind(dec1, dec2)
dec$sync.change<-abs(dec$sync.change)
dec$sync.change<--dec$sync.change
tog<-rbind(dec, inc)

#to calculate mismatch
interact$count<-1
length<-aggregate(interact["count"], interact[c("studyid", "intid")], FUN=sum)
names(length)[3]<-"total"
neg<-subset(interact, newphenodiff<0)
length_neg<-aggregate(neg["count"], neg[c("studyid", "intid")], FUN=sum)
close<-merge(length, length_neg, by=c("studyid","intid"))
all<-merge(sync, close, by=c("studyid","intid"))
all$count_diff<-with(all, count-total)

#fix mismatchs
sub<-subset(tog, intid=="1" | intid=="175" | intid=="235")
sub$sync.change<--sub$sync.change
sub2<-subset(tog, intid!="1" & intid!="175" & intid!="235")
tog<-rbind(sub, sub2)

tog  <- tog[with(tog , order(intid, iteration)),]

write.csv(tog, "/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/output/sync.change.1K.csv")

tog <- read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/output/sync.change.1K.csv", header=TRUE)


# Step 2- Load temp change data- from tempmodels_interactions.R
mdata <- read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/output/temp.change.1K.csv", header=TRUE)

#try<-as.data.frame(temp.change)
#try$is<-row.names(temp.change)
#goober<-melt(try, id="is")

#mdata<-melt(temp.change[,2:1001])
names(mdata)[1]<-"id"; names(mdata)[3]<-"iteration"; names(mdata)[4]<-"temp.change";
 #mdata$intid<-rep(intid.nodups$intid, 1000) # 1,19,144,145 x3000
 mdata<- mdata[with(mdata, order(intid)),]

# Step 3- Reduce synchrony interactions to only those with climate data
all<-merge(mdata[,2:4], tog, by=c("intid"))
all<- all[with(all , order(intid, iteration)),]

new<-inner_join(mdata, tog, by = "intid")

# Model
N <- nrow(all)
y <- all$sync.change
year<-all$temp.change
Nspp <- length(unique(all$intid))
species <- as.numeric(as.factor(all$intid))

cov.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

