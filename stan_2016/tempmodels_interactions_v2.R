##### TO GET TEMP CHANGE PER INTERACTION ######
### note that as of 23 Sep 2016 this file uses the same models and data as tempmodels_daasets.R ###
### But the order that the data go into the models in each file is different... ###

rm(list=ls())
options(stringsAsFactors=FALSE)


#row.names=FALSE
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
# setwd("~/Documents/git/projects/trophsynch/synchrony/stan_2016")
library(ggplot2)
library(rstan)
library(shinyStan)
library(grid)
library(nlme)
library(gridExtra)
library(reshape)
library(plyr)
library(dplyr)
# set_cppo("fast") 

#Load sync data first
#get data
rawlong <- read.csv("input/rawlong2.csv", header=TRUE)

# use for 1981 hinge model
source("input/datacleaningmore.R")

#Load temp data next
## source the temp dataset cleaning code
source("source/clean_tempdatasets2.R")


#MERGE!

merges<-merge(rawlong.tot2, clim3, by=c("studyid","species","year","yr1981"))
merges <- merges[with(merges, order(species, year)),]
merges2<-merge(rawlong.tot, clim3, by=c("studyid","species","year","newyear","yr1981")) #to get intid
merges2 <- merges2[with(merges2, order(species, year)),]

#Temp model#

# check, which datasetids have more than one species
#ddply(dataset, c("datasetid"), summarise,
#    nspp=n_distinct(species))

## FOR REPORTING TEMP CHANGE IN THE PAPER: only on unique datasets (n=18)
temp.id<-unique(datasetid.trans2[,1:3])
temp.id2<-merge(merges[,c("studyid","species","year","yr1981","envvalue","envtype")], temp.id)
temp.id3<-unique(temp.id2[,c("studyid","datasetid","year","yr1981","envvalue","envtype")])
temp.id3<-temp.id3[order(temp.id3$datasetid),]

N<-nrow(temp.id3)
Nspp <- length(unique(temp.id3$datasetid)) 
y <- temp.id3$envvalue
species <- as.numeric(as.factor(temp.id3$datasetid))
year <- temp.id3$yr1981

temp.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(temp.model, pars = c("mu_b", "sigma_y", "a", "b"))

temp.est <- extract(temp.model)# 
it1000 <- matrix(0, ncol=3000, nrow=1) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    it1000[,(i-3000)] <- temp.est$mu_b[i]
}
mean(rowMeans(it1000, na.rm=TRUE))



## OTHERWISE-- CLEAN UP TEMP DATA- get one dataset per interaction ##

temp.id<-datasetid.trans2[,1:3]
temp.id2<-merge(unique(merges2[,c("studyid","species","intid")]), temp.id)
temp.id2<-temp.id2[order(temp.id2$intid),]

id<-as.numeric(as.factor(merges$species))
merges$id<-id
spp.id<-merge(unique(merges[,c("studyid","species","id")]), temp.id2)
spp.id2<-merge(spp.id, unique(rawlong.tot[,c("studyid","species","spp")]))
spp.id2$row<-1:nrow(spp.id2)
#assign one role for competitions:
spp.id3<-subset(spp.id2, row!=7 & row!=10 & row!=14 & row!=15)


#mult temp change per interaction
#mult<-subset(spp.id, intid=="155" | intid=="170" | intid=="197" | intid=="198" | intid=="199" | intid=="200" | intid=="236" | intid=="237")
#single<-subset(spp.id, intid!="155" & intid!="170" & intid!="197" & intid!="198" & intid!="199" & intid!="200" & intid!="236" & intid!="237")
#choose random species per interaction (matters for where mult temp change per interaction but easy to do it where same temp.change repeated 2x)
uni_int<- data.frame(array(0, c(nrow=length(unique(spp.id3$intid)), ncol=5))) 
#b<-as.character(unique(mult$intid)); Bgroups<-length(b)
b<-as.character(unique(spp.id3$intid)); Bgroups<-length(b)
count<-1
for(i in 1:length(b)){
	int<-subset(spp.id3, intid==b[i])
	#int$row<-1:2 #for random
	#sub<-subset(int, row==sample(1:2, 1)) #for random
	#sub<-subset(int, spp=="spp1") #choose resource
	sub<-subset(int, spp=="spp2") #choose consumer
	uni_int[count,]<-sub[,1:5]
	count<-count+1
}
names(uni_int)<-names(spp.id)


merges_again<-merge(merges, uni_int[,c("studyid","species","id")]) #22 interactions, 18 species 
merges_again<-merges_again[order(merges_again$species),]

#estimate temp change on single species per interaction
N<-nrow(merges_again)
Nspp <- length(unique(merges_again$species)) 
y <- merges_again$envvalue
species <- as.numeric(as.factor(merges_again$species))
year <- merges_again$yr1981

temp.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(temp.model, pars = c("mu_b", "sigma_y", "a", "b"))

#for covariate model: #assign temp.change back to interaction
#load intid.nodups2
goo <- extract(temp.model)
summ_studyspp<-unique(merges_again[,c("studyid","species")])
it1000 <- matrix(0, ncol=3000, nrow=22)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- goo$b[i,]
    andtheanswer <- merge(intid.nodups2, summ_studyspp, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, summ_studyspp, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
        for(j in 1:nrow(andtheanswer)){
andtheanswer$final[j]<-with(andtheanswer, mean(c(model.x[j],model.y[j]), na.rm=TRUE))
    }
   it1000[,(i-3000)] <- andtheanswer$final
}

ndata<-melt(it1000[,2001:3000])
names(ndata)[1]<-"id"; names(ndata)[2]<-"iteration"; names(ndata)[3]<-"temp.change"; 
#ndata$intid<-rep(intid.nodups$intid, 1000)


#PHENO MODEL# - fit on all of 37 species because need both estimates of pheno change to get sync change
N<-nrow(merges)
y <- merges$phenovalue
Nspp <- length(unique(merges$species)) #J
species <- as.numeric(as.factor(merges$species))
year <- merges$yr1981

#Feb 2017: don't have enough repeating species ACROSS studies to esimate error, therefore no study as group
pheno.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(pheno.model, pars = c("mu_b", "sigma_y", "a", "b"))


### ESTIMATE SYNC CHANGE ###

fh.sim <- extract(pheno.model)# inc_warmup=FALSE) #only extracts parameters, does not include warmup

specieschar.formodel <- aggregate(merges2["phenovalue"], merges2[c("studyid", "species", "intid", "terrestrial","spp")], FUN=length) 
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species","intid"))
specieschar.formodel.sm  <- specieschar.formodel.sm [with(specieschar.formodel.sm , order(intid)),]
intid_full <- read.csv("input/raw_oct.csv", header=TRUE)
#fix for diatom4 spp. (Dec 2016)
sub<-subset(intid_full, intid=="194");
sub[,c("spp1")]<-"Diatom4a spp."
intid_full<-rbind(intid_full, sub)
sub<-subset(intid_full, intid=="195");
sub[,c("spp1")]<-"Diatom4b spp."
intid_full<-rbind(intid_full, sub)
sub<-subset(intid_full, intid=="196");
sub[,c("spp1")]<-"Diatom4c spp."
intid_full<-rbind(intid_full, sub)
intid_full<-subset(intid_full, spp1!="Diatom4 spp.")
intid_full<-subset(intid_full, spp1!="asdf1" & spp2!="asdf2")

lal<-unique(rawlong.tot[,c("intid","terrestrial")])
intid2<-merge(intid_full, lal, by=c("intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction","terrestrial"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
intid.nodups2<-merge(unique(merges2[,c("studyid","intid")]), intid.nodups)

summ_studyspp <- subset(specieschar.formodel, select=c("studyid", "species")); summ_studyspp<-unique(summ_studyspp)
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]

#summ_studyspp2<-merge(unique(merges2[,c("studyid","species")]), summ_studyspp)

#For interactions AND synchrony change- KEEP ALL INTERACTIONS
it1000 <- matrix(0, ncol=3000, nrow=length(unique(intid.nodups2$intid))) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    summ_studyspp$model <- fh.sim$b[i,]
    andtheanswer <- merge(intid.nodups2, summ_studyspp, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, summ_studyspp, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
    it1000[,(i-3000)] <- andtheanswer$model.x-andtheanswer$model.y #model.x=spp1
}
synch.model<-it1000

mdata<-melt(it1000[,2001:3000])
#mdata<-melt(it1000[,2001:2100]) #test
names(mdata)[1]<-"id"; names(mdata)[2]<-"iteration"; names(mdata)[3]<-"sync.change"; 
mdata$intid<-unique(andtheanswer$intid)

# For interpretation -- KEEP ALL INTERACTIONS
#Negative difference= decreasing synchrony ()
#Positive difference= increasing synchrony ()
# spp1 = negative species in trophic interaction (RESOURCE) and spp2= positive species (CONSUMER)
#spp1-spp2=#resource-consumer

## Figure out seasonal order for each year across interactions (positive phenodiff= consumer emerges before resource; negative phenodiff= resource emerges before consumer)
intid_full$newphenodiff<- with(intid_full, neg_phenovalue-pos_phenovalue) #spp1-spp2
season <- aggregate(intid_full["newphenodiff"], intid_full[c("studyid", "intid")], FUN=mean) ## Take average seasonal order for each interaction

mdata2<-merge(mdata, season)

#Group all scenarios where there has been an increase synchrony (i.e. species are getting closer together): e.g. one is where consumer emerges before resource and consumer has had no pheno shift or delayed and resource has advanced  second is where resource emerges
inc1<-subset(mdata2, newphenodiff>0 & sync.change<0)
inc2<-subset(mdata2, newphenodiff<0 & sync.change>0)
inc<-rbind(inc1, inc2)
inc$sync.change<-abs(inc$sync.change) # give all these scenarios the same sign since we know that they are all increasing
# # Group all scenarios where there has been a decrease in synchrony (i.e. speceis are getting farther apart)- e.g. one is where consumer emerges before resource and consumer has had a larger advance than resource
dec1<-subset(mdata2, newphenodiff>0 & sync.change>0)
dec2<-subset(mdata2, newphenodiff<0 & sync.change<0)
dec<-rbind(dec1, dec2)
dec$sync.change<-abs(dec$sync.change) # give all these scenarios the same sign since we know that they are all decreasing
dec$sync.change<--dec$sync.change #add negative value to distinguish from increasing synchrony changes
tog<-rbind(dec, inc)

#to calculate mismatch
intid_full$count<-1
length<-aggregate(intid_full["count"], intid_full[c("studyid", "intid")], FUN=sum)
names(length)[3]<-"total"
neg<-subset(intid_full, newphenodiff<0)
length_neg<-aggregate(neg["count"], neg[c("studyid", "intid")], FUN=sum)
close<-merge(length, length_neg, by=c("studyid","intid"))

all<-merge(mdata2, close, by=c("studyid","intid"))
all$count_diff<-with(all, count-total)

#fix mismatchs
sub<-subset(tog, intid=="175")
sub$sync.change<--sub$sync.change
sub2<-subset(tog, intid!="175")
tog<-rbind(sub, sub2); nrow(tog) #### FINAL SYNCHRONY DATASET
tog<-tog[order(tog$id),]

## NOW BUILD COVARIATE MODEL
tdata<-merge(ndata, tog, by=c("id","iteration"))
#solo$id<-1:nrow(solo);
#tdata2<-merge(tdata, solo[,c("studyid","species","id")], by=c("id"))
N<-nrow(tdata)
Nspp <- length(unique(tdata$id)) #J
species <- as.numeric(as.factor(tdata$id))
y <- tdata$sync.change 
year <- tdata$temp.change #absolute value of temp change

# March 2017 decision- not enough variation within species to assign each species its own slope (i.e. no varying slopes) 
cov.model<-stan("stanmodels/twolevelrandomintercept_woslopes.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(cov.model, pars = c("mu_a", "mu_b","sigma_y", "a"))


