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

interact <- read.csv("input/raw_oct.csv", header=TRUE)
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

syncdat.full <- read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/output/sync.change.1K.csv", header=TRUE)

# Step 2- Load temp change data- from tempmodels_interactions.R

tempdat.full <- read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/output/temp.change.1K.csv", header=TRUE, col.names=c("rowhere","intid", "iteration", "temp.change"))

syncdat <- subset(syncdat.full, select=c("intid", "sync.change", "iteration"))
tempdat <- subset(tempdat.full, select=c("intid", "temp.change", "iteration"))

tempdat$iteration <- as.numeric(as.factor(tempdat$iteration))
mode(syncdat$iteration)


# Step 3- Reduce synchrony interactions to only those with climate data
tryme <- inner_join(syncdat, tempdat, by=c("intid", "iteration"))


# Model
N <- nrow(tryme)
y <- tryme$sync.change
year<-tryme$temp.change
Nspp <- length(unique(tryme$intid))
species <- as.numeric(as.factor(tryme$intid))

cov.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomeffects.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

launch_shinystan(cov.model)

goo <- extract(cov.model)# dim(fh.sim$b) # number of iterations*number of species 

data<-unique(tryme[,c("intid")]); data<-as.data.frame(data)
data$stan<-colMeans(goo$b)

print(cov.model, pars=c("mu_a", "mu_b", "sigma_y", "sigma_a","sigma_b"))


#Figure
data1<-aggregate(tryme["sync.change"], list(tryme[,c("intid")]), FUN=mean)
data2<-aggregate(tryme["temp.change"], list(tryme[,c("intid")]), FUN=mean)
data3<-merge(data1, data2, by=c("Group.1"))

ggplot(data3, aes(y=abs(sync.change),x=temp.change))+geom_point(size=3)+theme(text = element_text(size=20), axis.text.x=element_text(size=20), axis.title.y=element_text(size=20, angle=90))+ylab("abs(Phenological change (days/yr))")+theme_bw()+xlab("Temperature change (C/yr)")
