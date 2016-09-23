### Started 1 April 2015 ###
### April Fool's Day! ###

## Run STAN with synchrony data ##

rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
library(ggplot2)
library(rstan)
library(shinyStan)
library(grid)
library(nlme)
library(dplyr)
set_cppo("fast")  # for best running speed
source("/users/kharouba/google drive/UBC/multiplot.R")
#library(reshape)
# library(lme4)

#get data
source("datacleaning.R")
# in the data file: spp1 = neg= negative species e.g. resource 

rawlong <- read.csv("rawlong.csv", header=TRUE)
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

##
# Heather's code to add in data formatted for hinge model
hinge <- subset(rawlong.nodups, intid=="170" | intid=="171" | intid=="177" |
    intid=="178" | intid=="179" | intid=="180" |intid=="181" | intid=="189" |
    intid=="191"|intid=="193" |intid=="194" | intid=="195" |intid=="196"|
    intid=="201" |intid=="207" |intid=="208")

hinge_non <- subset(rawlong.nodups, intid!="170" & intid!="171" & intid!="177"
     & intid!="178" & intid!="179" & intid!="180" & intid!="181" & intid!="189"
     & intid!="191" & intid!="193" & intid!="194" & intid!="195" & intid!="196"
     & intid!="201" & intid!="207" & intid!="208")

hinge_non$newyear<-hinge_non$year
hinge_pre<-subset(hinge, year<=1981); hinge_pre$newyear<-1981
hinge_post<-subset(hinge, year>1981); hinge_post$newyear<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)
rawlong.tot<-rbind(hinge_non, hinges)


# prep the data to fit the model including:
# aggregate to get species level characteristics
# subset down to the phenovalues
# Run stan on 108 unique intid-species interactions (n=108) (NOT unique species (n=91))
rawlong.tot <- arrange(rawlong.tot, species)
rawlong.tot$yr1981 <- rawlong.tot$newyear-1981
rawlong.tot <- rawlong.tot[with(rawlong.tot, order(species)),]
#rawlong.tot$newid<-with(rawlong, paste(intid,"_",species)); 
#rawlong.tot <- rawlong.tot[with(rawlong.tot, order(species)),]
N <- nrow(rawlong.tot)
y <- rawlong.tot$phenovalue
Nspp <- length(unique(rawlong.tot$species)) #newid is character !
species <- as.numeric(as.factor(rawlong.tot$species))
year <- rawlong.tot$yr1981

#specieschar.hin<- aggregate(rawlong.tot["phenovalue"], rawlong.tot[c("studyid","intid","species", "int_type", "terrestrial", "spp")], FUN=length) #number of years per species
#specieschar.hin <- specieschar.hin[with(specieschar.hin, order(intid)),]
#specieschar.hin2<-unique(specieschar.hin$species)


#New model as of June 2016
#Random slopes only, no random intercepts, hinge, no covariate matrix:
#sync.model<-stan("synchrony1_notype_randslops_wcovar.stan", data=c("N","Nspp","y","species","year"), iter=2000, warmup=1000, thin=10, chains=4)
sync.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

#Match up interacting species and look at differences

fh.sim <- extract(sync.model)# inc_warmup=FALSE) #only extracts parameters, does not include warmup
dim(fh.sim$b) # number of iterations*number of species 
names(fh.sim)
#sd is not considered parameter (only rows in output are parameters)
# here's one iteration
fh.sim$b[2000,]

specieschar.formodel <- aggregate(rawlong.nodups["phenovalue"], rawlong.nodups[c("studyid", "species", "intid", "terrestrial","spp")], FUN=length) 
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species","intid"))
specieschar.formodel.sm  <- specieschar.formodel.sm [with(specieschar.formodel.sm , order(intid)),]
intid <- read.csv("input/raw_april.csv", header=TRUE)
lal<-unique(rawlong.tot[,c("intid","terrestrial")])
intid2<-merge(intid, lal, by=c("intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction","terrestrial"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
#sync_int<-intid.nodups #synchrony change interactions


summ_studyspp <- subset(specieschar.formodel, select=c("studyid", "species")); summ_studyspp<-unique(summ_studyspp)
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]

#For interactions AND synchrony change
it1000 <- matrix(0, ncol=3000, nrow=length(unique(intid.sm$intid))) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    summ_studyspp$model <- fh.sim$b[i,]
    andtheanswer <- merge(intid.nodups, summ_studyspp, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, summ_studyspp, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
    it1000[,(i-3000)] <- andtheanswer$model.x-andtheanswer$model.y #model.x=spp1
}
synch.model<-it1000

#it1000 <- matrix(0, ncol=3000, nrow=Nspp) #2000 iterations for 53 interactions;
#for (i in 3000:6000){ # 2000 iterations?
 #   summ_studyspp$model <- fh.sim$b[i,]
  #  it1000[,(i-3000)] <- summ_studyspp$model
#}
#synch.model<-it1000
#dr$stan<-rowMeans(it1000, na.rm=TRUE)



# FOR INTERPRETATION
#Negative difference= decreasing synchrony ()
#Positive difference= increasing synchrony ()
# spp1 = negative species in trophic interaction and spp2= positive species

#spp1-spp2=#resource-consumer
meanchange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH INTXN
andtheanswer$meanchange<-meanchange

interact <- read.csv("input/raw_april.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue) #spp1-spp2
#positive phenodiff= consumer emerges BEFORE resource
#negative phenodiff= resource emerges BEFORE consumer
season <- aggregate(interact["newphenodiff"], interact[c("studyid", "intid")], FUN=mean)

andtheanswer2<-merge(andtheanswer, season, by=c("studyid","intid"))

#increasing synchrony (i.e. species are getting closer together)
inc1<-subset(andtheanswer2, newphenodiff>0 & meanchange<0)
inc2<-subset(andtheanswer2, newphenodiff<0 & meanchange>0)
inc<-rbind(inc1, inc2)
inc$meanchange<-abs(inc$meanchange)
# decreasing synchrony (i.e. speceis are getting farther apart)
dec1<-subset(andtheanswer2, newphenodiff>0 & meanchange>0)
dec2<-subset(andtheanswer2, newphenodiff<0 & meanchange<0)
dec<-rbind(dec1, dec2)
dec$meanchange<-abs(dec$meanchange)
dec$meanchange<--dec$meanchange
tog<-rbind(dec, inc)


#to calculate mismatch
interact$count<-1
length<-aggregate(interact["count"], interact[c("studyid", "intid")], FUN=sum)
names(length)[3]<-"total"
neg<-subset(interact, newphenodiff<0)
length_neg<-aggregate(neg["count"], neg[c("studyid", "intid")], FUN=sum)
close<-merge(length, length_neg, by=c("studyid","intid"))
all<-merge(andtheanswer2, close, by=c("studyid","intid"))
all$count_diff<-with(all, count-total)

#fix mismatchs
sub<-subset(tog, intid=="1" | intid=="175" | intid=="235")
sub$meanchange<--sub$meanchange
sub2<-subset(tog, intid!="1" & intid!="175" & intid!="235")
tog<-rbind(sub, sub2)

#KEY RESULT:
#Direction and magnitude
mean(tog$meanchange, na.rm=TRUE) #mean difference across interactions; they drift apart by half a day a decade
#computation of the standard error of the mean
sem<-sd(tog$meanchange)/sqrt(length(tog$meanchange)); sem
#95% confidence intervals of the mean
c(mean(tog$meanchange)-2*sem,mean(tog$meanchange)+2*sem)

#Magnitude only
mean(abs(tog$meanchange), na.rm=TRUE) #mean difference across interactions; they drift apart by half a day a decade
sem<-sd(tog$meanchange)/sqrt(length(tog$meanchange)); sem
#95% confidence intervals of the mean
c(mean(tog$meanchange)-2*sem,mean(tog$meanchange)+2*sem)

#by group
sub<-subset(tog, terrestrial=="terrestrial"); mean(sub$meanchange); sd(sub$meanchange)
sub<-subset(tog, terrestrial=="aquatic"); mean(sub$meanchange); mean(sub$meanchange); sd(sub$meanchange)
sem<-with(sub, sd(meanchange)/sqrt(length(meanchange))); sem
#95% confidence intervals of the mean
with(sub, c(mean(meanchange)-2*sem,mean(meanchange)+2*sem))

library(grid)
text_high <- textGrob("Smaller interval", gp=gpar(fontsize=13, fontface="bold"))
text_low <- textGrob("Larger interval", gp=gpar(fontsize=13, fontface="bold"))
ggplot(tog, aes(x=meanchange, fill=terrestrial))+geom_histogram(binwidth=.5, alpha=.5, position="identity", colour="black")+theme_bw()+geom_vline(xintercept=0, linetype="dashed",size=1)+annotation_custom(text_high, xmin=2, xmax=2, ymin=-0.5, ymax=-0.5)+annotation_custom(text_low, xmin=-2, xmax=-2, ymin=-0.5, ymax=-0.5)

#Link to individual species
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- fh.sim$b[i,]
    it1000[,(i-3000)] <- summ_studyspp$model
}
summ_studyspp$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP
#mean(summ_studyspp$stanfit)


indiv_intxn <- merge(tog, summ_studyspp[,c("studyid","species","stanfit")], by.x=c("studyid", "spp1"), by.y=c("studyid", "species"), all.x=TRUE)
indiv_intxn <- merge(indiv_intxn, summ_studyspp[,c("studyid","species","stanfit")], by.x=c("studyid", "spp2"), by.y=c("studyid", "species"), all.x=TRUE)
#indiv_intxn$direction<-ifelse(indiv_intxn$meanchange>0, 1, 0)
#with(indiv_intxn, cor(stanfit.x, stanfit.y))


ggplot(indiv_intxn, aes(x=stanfit.x, y=stanfit.y))+geom_point(aes(size=2))+geom_vline(xintercept=0, linetype="dashed")+geom_hline(yintercept=0, linetype="dashed")+theme_bw()+theme(legend.position="false")+xlab("resource")+ylab("consumer")

#annotate(geom="text", x=2, y=-1, label="smaller interval", colour="black")
#theme(plot.margin = unit(c(1,1,2,1), "lines"))


#######   Magnitude of change: Phenodiff- year 1 difference ########
interact <- read.csv("input/raw_april.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue) #spp
sp<-merge(rawlong.tot, interact[,c("studyid","intid","year","newphenodiff")], by=c("studyid","intid","year"))
baseline<-aggregate(sp["newphenodiff"], sp[c("studyid", "intid")], FUN=length)


sp<-sp[order(sp$intid,sp$year),]; sp<-na.omit(sp)
best<-data.frame(array(0, c(nrow(sp), 4)))
names(best)[1]<-"intid"; names(best)[2]<-"year"; names(best)[3]<-"phenodiff_base"; names(best)[4]<-"base"
Bgroups<-unique(sp$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-sp[sp$intid==b[i],]
asdf<-rowcount+(nrow(spp)-1)
best[rowcount:asdf,1]<-b[i]
best[rowcount:asdf,2]<-spp[,c("year")]
best[rowcount:asdf,3]<- spp$newphenodiff-spp[1,c("newphenodiff")]
best[rowcount:asdf,4]<- spp[1,c("newphenodiff")]
rowcount<-rowcount+nrow(spp)
}
best2<-unique(best)
sp2<-merge(sp, best2, by=c("intid","year"))
sp<-sp2


sp <- sp[with(sp, order(intid, species)),]
sp2<-unique(sp[,c("studyid","intid","yr1981","phenodiff_base")])
N <- nrow(sp2)
y <- abs(sp2$phenodiff_base)
#specieschar.hin<- aggregate(sp2["phenodiff_base"], sp2[c("studyid", "intid")], FUN=length) #number of years per species
#specieschar.hin <- specieschar.hin[with(specieschar.hin, order(intid)),]
Nspp <- length(unique(sp2$intid))
species <- as.numeric(as.factor(sp2$intid))
year <- sp2$yr1981

baseline.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)


goo <- extract(baseline.model)# dim(fh.sim$b) # number of iterations*number of species 


specieschar.formodel <- aggregate(sp2["phenodiff_base"], sp2[c("studyid", "intid")], FUN=length) 
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "intid"))
specieschar.formodel.sm  <- specieschar.formodel.sm [with(specieschar.formodel.sm , order(intid)),]
intid <- read.csv("input/raw_april.csv", header=TRUE)
lal<-unique(rawlong.tot[,c("intid","terrestrial")])
intid2<-merge(intid, lal, by=c("intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction","terrestrial"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
sync_int<-intid.nodups #synchrony change interactions


summ_studyspp <- subset(specieschar.formodel, select=c("studyid", "intid")); summ_studyspp<-unique(summ_studyspp)
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(intid)),]

#For interactions AND synchrony change
it1000 <- matrix(0, ncol=3000, nrow=length(unique(intid.sm$intid))) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    it1000[,(i-3000)] <- goo$b[i,]
       #it1000[,(i-3000)] <- andtheanswer$model.x-andtheanswer$model.y #model.x=spp1
}
synch.model<-it1000
summ_studyspp$meanslope<-abs(rowMeans(synch.model)) # take absolute value because we want to know how much relative timing has changed since year 1

*** to compare stan fits to linear fits : scatterplot with 1:1 line (abline(0,1))

#Results
median(summ_studyspp$meanslope)
sem<-sd(summ_studyspp$meanslope)/sqrt(length(summ_studyspp$meanslope)); sem

because hist is non-normal use: Bootstrap CI calculations
boot.ci(boot(summ_studyspp$meanslope,function(x,i) median(x[i]), R=1000))
OR quantiles of binomial distribution:
sort(summ_studyspp$meanslope)[qbinom(c(.025,.975), length(summ_studyspp$meanslope), 0.5)]

#95% confidence intervals of the mean
c(mean(summ_studyspp$meanslope)-2*sem,mean(summ_studyspp$meanslope)+2*sem)


print(baseline.model,"b",probs=c(.025,.975))

##### FIGURES #####
#plot all species
asdf<-summary(sync.model, pars="b")
#to get median coefficients from SUMMARY
median<-asdf[[1]][1:91]; new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) #-0.715645239
min<-asdf[[1]][274:364]; e<-data.frame(min=unlist(min)) #-1.1222299
max<-asdf[[1]][638:728]; f<-data.frame(max=unlist(max)) #-0.316333047
d$min<-e$min; d$max<-f$max;
details<-unique(specieschar.formodel[,c("studyid","species","terrestrial")])
details <- details[with(details, order(species)),]
d$terrestrial<-details$terrestrial; d$species<-details$species

#asdf[[1]][182:273] #to get sd
#min<-asdf[[1]][274:364]; min<-as.data.frame(min) #to get 2.5%
#asdf[[1]][365:455] #to get 25 (-0.658)
#asdf[[1]][456:546] # to get 50
#asdf[[1]][547:637] # to get 75
#max<-asdf[[1]][638:728]; max<-as.data.frame(max) # to get 97.5

d$group <- factor(d$terrestrial, levels = d$terrestrial[order(d$terrestrial)])

ggplot(d, aes(x=factor(species), y=y, ymin=min, ymax=max, colour=terrestrial))+geom_pointrange(aes(width=.05, colour=factor(terrestrial)))+theme_bw()+geom_point(size=3, aes(colour=factor(terrestrial)))+xlab("species")+ylab("Phenological change (doy/year)")+geom_hline(xintercept=0, linetype="longdash", colour="grey")+theme(legend.position="false")+coord_flip()

OR
ggplot(d, aes(x=factor(species), y=y, ymin=min, ymax=max, colour=terrestrial))+geom_pointrange(aes(width=.05, colour=factor(terrestrial)))+theme_bw()+geom_point(size=3, aes(colour=factor(terrestrial)))+xlab("species")+ylab("Phenological change (doy/year)")+geom_hline(xintercept=0, linetype="longdash", colour="grey")+theme(legend.position="false")+facet_wrap(~terrestrial, scale="free")+theme(axis.text.x = element_text(angle = 45, hjust = 1))

______________________________
EXTRA!!!!
#to get variance (sd) estimate of slope
sd_it1000<-matrix(0,ncol=2000, nrow=89)
for (i in 2000:4000){
sd_it1000[,(i-2000)]<-fh.sim$b[i,]
}
specieschar.formodel.sm$sd<-apply(sd_it1000, 1, sd)

###
par(mfrow=(c(1,2)))
# get the mean slope shifts
mean(colMeans(fh.sim$b))*10 #mean across iterations for EACH SPECIES; they shift about 4 days/decade
hist(colMeans(fh.sim$b)*10, main="", xlab="change in phenology (days/decade)", breaks=20)
hist(meanchange, main="", xlab="change in synchrony (days/decade)", breaks=20)



#95% posterior interval:
print(sync.model, "b", probs=c(.025,.975))
fh.sim

andtheanswer <- cbind(andtheanswer, meanchange)
m1<-lme(meanchange~1, random=~1|studyid, na.action=na.omit, method="ML", andtheanswer); summary(m1)


ggplot(andtheanswer, aes(x=meanchange, fill=interaction)) +
    geom_histogram(binwidth=0.5, alpha=0.5, position="identity")

ggplot(andtheanswer, aes(x=meanchange, colour=interaction)) + geom_density() # I don't like geom_density but my default histogram is awful

# histograms
# Interaction type
preds <- subset(andtheanswer, interaction=="predation")
comps <- subset(andtheanswer, interaction=="competition")
polln <- subset(andtheanswer, interaction=="pollination")
herbiv <- subset(andtheanswer, interaction=="herbivory")
paras <- subset(andtheanswer, interaction=="parasitism")
mutul <- subset(andtheanswer, interaction=="mutuliasm")

par(mfrow=c(2,3))
xlim <- c(-17,17)
hist(preds$meanchange, main="predators", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(comps$meanchange, main="competitors", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(polln$meanchange, main="pollination", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(herbiv$meanchange, main="herbivory", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(paras$meanchange, main="parasitism", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(mutul$meanchange, main="mutualism", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)

ggplot(andtheanswer, aes(x=meanchange, fill=interaction)) + 
    geom_histogram(data = , fill = "red", alpha = 0.2) + 
    geom_histogram(data = subset(andtheanswer, interaction==""), fill = "blue", alpha = 0.2) + 
    geom_histogram(data = subset(andtheanswer, interaction==""), fill = "purple", alpha = 0.2) + 
    geom_histogram(data = subset(andtheanswer, interaction=="herbivory"), fill = "green", alpha = 0.2) +
    geom_histogram(data = subset(andtheanswer, interaction=="pollination"), fill = "orange", alpha = 0.2) + 
    geom_histogram(data = subset(andtheanswer, interaction=="parasitism"), fill = "yellow", alpha = 0.2) 

#Ecosystem
aqua <- subset(andtheanswer, terrestrial=="aquatic")
terra <- subset(andtheanswer, terrestrial=="terrestrial")
par(mfrow=c(1,2))
xlim<-c(-16,20)
hist(aqua$meanchange, main="aquatic", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)
hist(terra$meanchange, main="terrestrial", xlab="change in phenology (days/decade)", breaks=20, xlim=xlim)

# EXPLORATORY FIGURES:
#Repeating species
so<-subset(rawlong.tot, spp=="spp1"); spp1<-unique(so[,c("studyid","intid","figcode")]); names(spp1)[3]<-"figcode_1"
sp<-subset(rawlong.tot, spp=="spp2"); spp2<-unique(sp[,c("studyid","intid","figcode")]); names(spp2)[3]<-"figcode_2"
tots<-merge(spp1,spp2, by=c("studyid","intid"))
comp<-merge(andtheanswer, tots, by=c("studyid","intid"))

re_uni<-ggplot(comp, aes(x=model.x, fill=figcode_1)) + 
     geom_histogram(data =subset(comp, figcode_1=="unique") , fill = "grey", col="black", alpha=0.2)+theme_bw()+ggtitle("Resource_unique")+xlim(-1.5,1.5)+ylim(0,5)
re<-ggplot(comp, aes(x=model.x, fill=figcode_1))+
     geom_histogram(data = subset(comp, figcode_1!="unique"))+theme_bw()+ggtitle("Resource_repeat")+theme(legend.position="none")+xlim(-1.5,1.5)+ylim(0,5)
co_uni<-ggplot(comp, aes(x=model.x, fill=figcode_2)) + 
     geom_histogram(data =subset(comp, figcode_2=="unique") , fill = "grey", col="black", alpha=0.2)+theme_bw()+ggtitle("Consumer_unique")+theme(legend.position="none")+xlim(-1.5,1.5)
co<-ggplot(comp, aes(x=model.x, fill=figcode_2)) + 
     geom_histogram(data = subset(comp, figcode_2!="unique"))+theme_bw()+ggtitle("Consumer_repeat")+xlim(-1.5,1.5)+theme(legend.position="none")+ylim(0,6)
multiplot(re_uni, co_uni, re, co, cols=2)

#Repeating interactions
so<-subset(rawlong.tot, spp=="spp1"); spp1<-unique(so[,c("studyid","intid","figcode_study")]); 
comp<-merge(unique(andtheanswer[,c("studyid","intid","meanchange")]), spp1, by=c("studyid","intid"))


re_uni<-ggplot(comp, aes(x=meanchange, fill=figcode_study)) + 
     geom_histogram(data =subset(comp, figcode_study=="unique") , fill = "grey", col="black", alpha=0.2)+theme_bw()+ggtitle("Studies_unique")+xlim(-12,8)+ylim(0,4)
re<-ggplot(comp, aes(x=meanchange, fill=figcode_study))+
     geom_histogram(data = subset(comp, figcode_study!="unique"))+theme_bw()+ggtitle("Studies_repeat")+xlim(-12,9)+ylim(0,4)+theme(legend.position="none")
multiplot(re_uni, re, cols=1)


## OLD MODELS ###
#New model as Apri. 2016 here!
#Random slopes only, no type (aka Margaret Kosmala's model) and with hinge
fit.hinge.rs<-stan("synchrony_apr_nocov.stan", data=c("N","J","y","species","year"), iter=2000, chains=4)

#New model as of Dec 2015 here!
#Random slopes only, no type (aka Margaret Kosmala's model) and with hinge
Imat <- diag(1, nVars)
fit.hinge.rs<-stan("synchrony1_notype_randslops_wcovar.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)

launch_shinystan(fit.hinge.rs)

#  Random slopes, random intercepts, no type (aka Margaret Kosmala's model) and with hinge
fit.hinge<-stan("synchrony1_notype_wcovar.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)



#
andtheanswer <- cbind(meanchange, andtheanswer)
bigchanges <- subset(andtheanswer, abs(meanchange)>5)
smallerchanges <- subset(andtheanswer, abs(meanchange)<5)
andtheanswer.sm <- subset(andtheanswer, is.nan(meanchange)==FALSE)