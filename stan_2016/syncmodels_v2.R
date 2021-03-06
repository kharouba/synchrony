### Started 1 April 2015 ###
### April Fool's Day! ###

## Run STAN with synchrony data ##

rm(list=ls()) 
options(stringsAsFactors = FALSE)


# Set working directory: 
if(length(grep("Lizzie", getwd())>0)) {    setwd("~/Documents/git/projects/trophsynch/synchrony/stan_2016") 
} else 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")


# libraries
library(ggplot2)
library(rstan)
library(shinystan)
library(grid)
library(nlme)
library(dplyr)
library(ggrepel)
set_cppo("fast")  # for best running speed
source("..//multiplot.R")
#library(reshape)
# library(lme4)

#get data
rawlong <- read.csv("input/rawlong2.csv", header=TRUE)

# use for 1981 hinge model
source("input/datacleaningmore.R")


#New model as of June 2016
#Random slopes only, no random intercepts, hinge, no covariate matrix:
#sync.model<-stan("synchrony1_notype_randslops_wcovar.stan", data=c("N","Nspp","y","species","year"), iter=2000, warmup=1000, thin=10, chains=4)

sync.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(sync.model, pars = c("mu_b", "sigma_y", "a", "b"))


#just for pheno change ie. >1981
aftercc<-subset(rawlong.tot2, year>1981)
N <- nrow(aftercc)
y <- aftercc$phenovalue
Nspp <- length(unique(aftercc$species)) #newid is character !
species <- as.numeric(as.factor(aftercc$species))
year <- aftercc$yr1981
pheno.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)
print(pheno.model, pars = c("mu_b", "sigma_y", "a", "b"))

#For pheno change
#Population-level slope estimates in text are from “print(stan)”
fh.sim <- extract(pheno.model)# 
it1000 <- matrix(0, ncol=3000, nrow=1) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    it1000[,(i-3000)] <- fh.sim$mu_b[i]
}
mean(rowMeans(it1000, na.rm=TRUE))
sem<-sd(it1000)/sqrt(length(it1000)); sem
#95% confidence intervals of the mean
c(mean(it1000)-2*sem,mean(it1000)+2*sem)

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

alls <- matrix(0, ncol=4, nrow=54);
for(i in 1:nrow(synch.model)){
sem<-sd(synch.model[i,])/sqrt(length(synch.model[i,])); sem
#95% confidence intervals of the mean
ci<-c(mean(synch.model[i,])-2*sem,mean(synch.model[i,])+2*sem)
alls[i,2]<-ci[1]
alls[i,3]<-ci[2]
alls[i,4]<-mean(synch.model[i,])
}


# FOR INTERPRETATION
#Negative difference= decreasing synchrony ()
#Positive difference= increasing synchrony ()
# spp1 = negative species in trophic interaction (RESOURCE) and spp2= positive species (CONSUMER)

#spp1-spp2=#resource-consumer
meanchange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH INTXN
andtheanswer$meanchange<-meanchange

## Figure out seasonal order for each year across interactions (positive phenodiff= consumer emerges before resource; negative phenodiff= resource emerges before consumer)
intid_full$newphenodiff<- with(intid_full, neg_phenovalue-pos_phenovalue) #spp1-spp2
season <- aggregate(intid_full["newphenodiff"], intid_full[c("studyid", "intid")], FUN=mean) ## Take average seasonal order for each interaction

andtheanswer2<-merge(andtheanswer, season, by=c("studyid","intid"))

#Group all scenarios where there has been an increase synchrony (i.e. species are getting closer together): e.g. one is where consumer emerges before resource and consumer has had no pheno shift or delayed and resource has advanced  second is where resource emerges
inc1<-subset(andtheanswer2, newphenodiff>0 & meanchange<0)
inc2<-subset(andtheanswer2, newphenodiff<0 & meanchange>0)
inc<-rbind(inc1, inc2)
inc$meanchange<-abs(inc$meanchange) # give all these scenarios the same sign since we know that they are all increasing
# # Group all scenarios where there has been a decrease in synchrony (i.e. speceis are getting farther apart)- e.g. one is where consumer emerges before resource and consumer has had a larger advance than resource
dec1<-subset(andtheanswer2, newphenodiff>0 & meanchange>0)
dec2<-subset(andtheanswer2, newphenodiff<0 & meanchange<0)
dec<-rbind(dec1, dec2)
dec$meanchange<-abs(dec$meanchange) # give all these scenarios the same sign since we know that they are all decreasing
dec$meanchange<--dec$meanchange #add negative value to distinguish from increasing synchrony changes
tog<-rbind(dec, inc)


#to calculate mismatch
intid_full$count<-1
length<-aggregate(intid_full["count"], intid_full[c("studyid", "intid")], FUN=sum)
names(length)[3]<-"total"
neg<-subset(intid_full, newphenodiff<0)
length_neg<-aggregate(neg["count"], neg[c("studyid", "intid")], FUN=sum)
close<-merge(length, length_neg, by=c("studyid","intid"))
all<-merge(andtheanswer2, close, by=c("studyid","intid"))
all$count_diff<-with(all, count-total)

#fix mismatchs
sub<-subset(tog, intid=="1" | intid=="175" | intid=="235")
sub$meanchange<--sub$meanchange
sub2<-subset(tog, intid!="1" & intid!="175" & intid!="235")
tog<-rbind(sub, sub2); nrow(tog) #### FINAL SYNCHRONY DATASET

#### KEY RESULT:
#Direction and magnitude
mean(tog$meanchange, na.rm=TRUE) #mean difference across interactions; they drift apart by half a day a decade
!number decreasing: 
sub<-subset(tog, meanchange<0); nrow(sub)
!number increasing:
sub<-subset(tog, meanchange>0); nrow(sub)
#computation of the standard error of the mean
sem<-sd(tog$meanchange)/sqrt(length(tog$meanchange)); sem
#95% confidence intervals of the mean
c(mean(tog$meanchange)-2*sem,mean(tog$meanchange)+2*sem)

#by group- with direction and magnitude
sub<-subset(tog, terrestrial=="terrestrial"); mean(sub$meanchange); sd(sub$meanchange)
sem<-sd(sub$meanchange)/sqrt(length(sub$meanchange)); sem
#95% confidence intervals of the mean
c(mean(sub$meanchange)-2*sem,mean(sub$meanchange)+2*sem)
sub<-subset(tog, terrestrial=="aquatic"); mean(sub$meanchange); sd(sub$meanchange)
sem<-sd(sub$meanchange)/sqrt(length(sub$meanchange)); sem
#95% confidence intervals of the mean
c(mean(sub$meanchange)-2*sem,mean(sub$meanchange)+2*sem)


#Magnitude only
mean(abs(tog$meanchange), na.rm=TRUE) #mean difference across interactions; they drift apart by half a day a decade
sem<-sd(abs(tog$meanchange)/sqrt(length(abs(tog$meanchange)))); sem
#95% confidence intervals of the mean
c(mean(abs(tog$meanchange))-2*sem,mean(abs(tog$meanchange))+2*sem)

#by group- magnitude
sub<-subset(tog, terrestrial=="terrestrial"); mean(abs(sub$meanchange)); sd(abs(sub$meanchange))
sub<-subset(tog, terrestrial=="aquatic"); mean(abs(sub$meanchange)); sd(abs(sub$meanchange))
sem<-with(sub, sd(abs(meanchange))/sqrt(length(abs(meanchange)))); sem
#95% confidence intervals of the mean
with(sub, c(mean(meanchange)-2*sem,mean(meanchange)+2*sem))


############################
####Covariates
rawlong.tot$count<-1
startdate<-aggregate(rawlong.tot["year"], rawlong.tot[c("studyid", "intid")], FUN=min); names(startdate)[3]<-"minyear"
ts_length<-aggregate(rawlong.tot["count"], rawlong.tot[c("studyid", "intid","spp")], FUN=sum); names(ts_length)[4]<-"length"
cov<-merge(startdate, unique(ts_length[,c("studyid","intid","length")]), by=c("studyid","intid"))
covs<-merge(cov, unique(rawlong.tot[,c("studyid","intid","phenofreq")]), by=c("studyid","intid"))
tog2<-merge(tog, covs, by=c("studyid","intid"))
tog2$minyr1981<-tog2$minyear-1981

#stan
simple linear models
x<-with(tog2, minyear-1951)
x<-tog2$minyear
x<-tog2$length
x<-tog2$phenofreq
x<-ifelse(tog2$phenofreq == "daily", 1, ifelse(tog2$phenofreq == "weekly", 2, 3))
#x1<-tog2$length
y<-tog2$meanchange
N<-nrow(tog2)

stan.lm<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/simplelinearmodel.stan", data=c("N","x","y"), iter=2000, chains=4)
print(stan.lm)
lol<-extract(stan.lm)
mean(lol$beta)

m1<-lm(meanchange~minyear, data=tog2); summary(m1)
with(tog2, plot(abs(meanchange)~length)); 

ggplot(tog2, aes(x=minyear, y=meanchange))+geom_point()+geom_abline(slope=0.01, intercept=-23.27)

# multiple factor model
N <- nrow(tog2)
y <- abs(tog2$meanchange)
Nspp <- length(unique(tog2$studyid)) #newid is character !
species <- as.numeric(as.factor(tog2$studyid))
p<-4 #number of predictors +1
desMat <- model.matrix(object = ~ 1 + minyr1981 + phenofreq + length, data = tog2)

cov.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomeffects_cov.stan", data=c("N","Nspp","y","species","p", "desMat"), iter=3000, chains=4)

OR
#single factor model
#year<-tog2$phenofreq #HOW TO CODE CATEGORICAL??? I think int is ok
m<-desMat[,3]; year<-as.vector(t(m)); 
year<-tog2$length
year<-as.numeric(tog2$minyr1981)

>1981
tog3<-subset(tog2, minyear>1981)
N <- nrow(tog3)
y <- abs(tog3$meanchange)
Nspp <- length(unique(tog3$studyid)) #newid is character !
species <- as.numeric(as.factor(tog3$studyid))
year<-as.numeric(tog3$minyear)

cov.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomeffects.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

print(cov.model, pars=c("mu_a", "mu_b", "sigma_y", "sigma_a","sigma_b"))

---
#FIGURES-v2 (days/decade)
# All groups together
library(grid)
#text_high <- textGrob("Closer together", gp=gpar(fontsize=12, fontface="bold"))
#text_low <- textGrob("Further apart", gp=gpar(fontsize=12, fontface="bold"))
text_high <- textGrob("Closer together", gp=gpar(fontsize=10))
text_low <- textGrob("Further apart", gp=gpar(fontsize=10))
null<-read.csv("nullmodel_sync.csv", header=TRUE, na.strings="NA", as.is=TRUE) #values for null model are from here
tog$meanperdec<-with(tog, meanchange*10)
tog$null<-nots$meanchange
tog$nullperdec<-with(tog, null*10)
#magnitude
a<-ggplot(tog, aes(x=abs(meanperdec)))+geom_histogram(aes(x=abs(nullperdec), colour="red"), binwidth=0.2)+geom_histogram(binwidth=2, alpha=.5, position="identity", colour="black")+theme_bw()+geom_vline(xintercept=4.4, linetype="solid",size=1)+geom_vline(xintercept=0.23, linetype=2,size=0.75)+theme(legend.position="none", axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/decade")+ annotation_custom(grob = textGrob(label = "a)", hjust = 0, gp = gpar(cex = 1.5)), ymin = 21, ymax = 21, xmin = -1.5, xmax = -1.5)


c<-ggplot(tog, aes(x=meanperdec))+geom_histogram(aes(x=nullperdec, colour="red"), binwidth=0.15)+geom_histogram(binwidth=2.25, alpha=.5, position="identity", colour="black")+theme_bw()+geom_vline(xintercept=-0.51, linetype="solid",size=1)+geom_vline(xintercept=0.048, linetype=2,size=0.75)+annotation_custom(text_high, xmin=10.2, xmax=13, ymin=-0.28, ymax=-0.5)+annotation_custom(text_low, xmin=-13, xmax=-13, ymin=-0.23, ymax=-0.5)+theme(legend.position="none", axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/decade")+ annotation_custom(grob = textGrob(label = "c)", hjust = 0, gp = gpar(cex = 1.5)), ymin = 17, ymax = 17, xmin = -17.5, xmax = -17.5)
#all five iterations:
+geom_vline(xintercept=-0.084, linetype=2,size=0.05)
geom_vline(xintercept=0.073, linetype=2,size=0.05)
geom_vline(xintercept=-0.058, linetype=2,size=0.05)
+geom_vline(xintercept=0.094, linetype=2,size=0.05)


# Figure - RElationship between individual species change and synchrony change. Showing stan slope estimates with 95% credible intervals.
#Link to individual species
tog2<-tog[,c(1:6,9:10)]
summ_studyspp2<-summ_studyspp[,1:2]
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp2$model <- fh.sim$b[i,]
    it1000[,(i-3000)] <- summ_studyspp2$model
}
summ_studyspp2$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP
#mean(summ_studyspp$stanfit)

asdf<-summary(sync.model, pars="b")
#to get median coefficients from SUMMARY
median<-asdf[[1]][1:88]; #new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) 
min<-asdf[[1]][265:352]; e<-data.frame(min=unlist(min)) #-1.129 to -1.00
max<-asdf[[1]][617:704]; f<-data.frame(max=unlist(max)) #-0.3109 to 0.125
d$min<-e$min; d$max<-f$max;


summ_studyspp2$min<-d$min; summ_studyspp2$max<-d$max; summ_studyspp2$median<-d$y


indiv_intxn <- merge(tog, summ_studyspp[,c("studyid","species","stanfit","min","max")], by.x=c("studyid", "spp1"), by.y=c("studyid", "species"), all.x=TRUE)
indiv_intxn <- merge(indiv_intxn, summ_studyspp[,c("studyid","species","stanfit","min","max")], by.x=c("studyid", "spp2"), by.y=c("studyid", "species"), all.x=TRUE)
#indiv_intxn$direction<-ifelse(indiv_intxn$meanchange>0, 1, 0)
#with(indiv_intxn, cor(stanfit.x, stanfit.y))


tryagain<-merge(summ_studyspp2, unique(rawlong.tot[,c("intid","species","spp")]), by="species")
tryagain2<-merge(tryagain, tog[,c("studyid","intid","meanchange")], by=c("studyid","intid"))
uni<-as.data.frame(unique(tryagain2$species))
uni$label<-1:nrow(uni); names(uni)[1]<-"species"
tryagain3<-merge(tryagain2, uni, by=c("species"))

#sub1$label<-letters

#NEWEST :#ggtitle("Stan slope estimates with 95% credible intervals")
text_small <- textGrob("Smaller synchrony change", rot=90, gp=gpar(fontsize=10, fontface="bold"))
text_large <- textGrob("Larger synchrony change", rot=90, gp=gpar(fontsize=10, fontface="bold"))
text_left <- textGrob("Advancement", rot=0, gp=gpar(fontsize=8, fontface="bold"))
text_right <- textGrob("Delay", rot=0, gp=gpar(fontsize=8, fontface="bold"))
tryagain3$stanfit_dec<-with(tryagain3, stanfit*10)
tryagain3$min_dec<-with(tryagain3, min*10)
tryagain3$max_dec<-with(tryagain3, max*10)
b<-ggplot(tryagain3, aes(x=factor(reorder(intid, abs(meanchange))), y=stanfit_dec, label = label))+geom_errorbar(aes(ymin=min_dec, ymax=max_dec, linetype=factor(spp)), width=.0025, colour="black")+geom_hline(yintercept=0, linetype="dashed")+geom_point(size=4, aes(order=abs(meanchange), colour=factor(spp), shape=factor(spp)))+theme_bw()+theme(legend.position=c(.83, .1))+xlab("interactions")+ylab("phenological change (days/decade)")+coord_flip()+scale_linetype(name="Species role", labels=c("Resource","Consumer"))+scale_shape(name="Species role", labels=c("Resource", "Consumer"))+scale_color_discrete(name="Species role", labels=c("Resource", "Consumer"))+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=8), axis.title.y=element_text(size=15, angle=90))+geom_text(size=3, hjust=2)+annotation_custom(text_small, xmin=9, xmax=9, ymin=1, ymax=47)+annotation_custom(text_large, xmin=47, xmax=47, ymin=1, ymax=47)+annotation_custom(text_left, xmin=1, xmax=1, ymin=-26, ymax=-20)+annotation_custom(text_right, xmin=1, xmax=1, ymin=1, ymax=40)+annotation_custom(grob = textGrob(label = "b)", hjust = 0, gp = gpar(cex = 1.5)), ymin = -4, ymax = -50, xmin = 51, xmax = 54)

multiplot(a,c,cols=1)
multiplot(b,cols=1)

#Figures- V1
# All groups together
library(grid)
text_high <- textGrob("Closer together", gp=gpar(fontsize=13, fontface="bold"))
text_low <- textGrob("Further apart", gp=gpar(fontsize=13, fontface="bold"))
null<-read.csv("nullmodel_sync.csv", header=TRUE, na.strings="NA", as.is=TRUE) #values for null model are from here
tog$meanperdec<-with(tog, meanchange*10)
a<-ggplot(tog, aes(x=meanchange))+geom_histogram(binwidth=.5, alpha=.5, position="identity", colour="black")+theme_bw()+geom_vline(xintercept=-0.052, linetype="solid",size=1)+geom_vline(xintercept=0.0073, linetype=2,size=0.05)+geom_vline(xintercept=-0.0084, linetype=2,size=0.05)+geom_vline(xintercept=-0.0058, linetype=2,size=0.05)+geom_vline(xintercept=0.0094, linetype=2,size=0.05)+geom_vline(xintercept=0.0023, linetype=2,size=0.05)+annotation_custom(text_high, xmin=1.2, xmax=1.2, ymin=-0.5, ymax=-0.5)+annotation_custom(text_low, xmin=-1.2, xmax=-1.2, ymin=-0.5, ymax=-0.5)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/year")+ annotation_custom(grob = textGrob(label = "a", hjust = 0, gp = gpar(cex = 1.5)), ymin = 15, ymax = 15, xmin = -2, xmax = -2)



#Magnitude
ggplot(tog, aes(x=abs(meanchange)))+geom_histogram(binwidth=.5, alpha=.5, position="identity", colour="black")+theme_bw()+geom_vline(xintercept=0.44, linetype="solid",size=1)+geom_vline(xintercept=0.042, linetype=2,size=1)+annotation_custom(text_high, xmin=1.5, xmax=1.5, ymin=-0.5, ymax=-0.5)+annotation_custom(text_low, xmin=-1.5, xmax=-1.5, ymin=-0.5, ymax=-0.5)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/year")


#Aquatic vs. terrestrial
text_high <- textGrob("Closer together", gp=gpar(fontsize=13, fontface="bold"))
text_low <- textGrob("Further apart", gp=gpar(fontsize=13, fontface="bold"))
ggplot(tog, aes(x=meanchange, fill=terrestrial))+geom_histogram(binwidth=.5, alpha=.5, position="identity", colour="black")+theme_bw()+geom_vline(xintercept=0, linetype="dashed",size=1)+annotation_custom(text_high, xmin=1.5, xmax=1.5, ymin=-0.5, ymax=-0.5)+annotation_custom(text_low, xmin=-1.4, xmax=-1.6, ymin=-0.5, ymax=-0.5)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/year")+theme(legend.position="none")


##OLDER version
ggplot(tryagain, aes(x=factor(reorder(intid, abs(meanchange))), y=stanfit, label = label))+geom_errorbar(aes(ymin=min, ymax=max, linetype=factor(spp)), width=.0025, colour="black")+geom_hline(yintercept=0, linetype="dashed")+geom_point(size=4, aes(order=abs(meanchange), colour=abs(meanchange), shape=factor(spp)))+theme_bw()+theme()+xlab("interactions")+ylab("phenological change (days/year)")+coord_flip()+scale_colour_continuous(name="abs(synchrony change)")+scale_linetype(name="Species role", labels=c("Resource","Consumer"))+scale_shape(name="Species role", labels=c("Resource", "Consumer"))+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=10), axis.title.y=element_text(size=15, angle=90))+geom_text(size=3, hjust=2)
#geom_label_repel()
#+geom_text(hjust=3)

#Extra
#Just on pre_cc data 
rawlong.tot2<-unique(rawlong.tot[,c("studyid","species","phenovalue","year","yr1981")]) 
rawlong.tot2$count<-1
pre<-subset(rawlong.tot2, year<=1981)
sss<- aggregate(pre["count"], pre[c("studyid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas
pre_cc<-merge(rawlong.tot2, sss2, by=c("studyid", "species"))
pre_cc<-subset(pre_cc, year<=1981)
rawlong.tot<-pre_cc
rawlong.tot <- arrange(rawlong.tot, species)
rawlong.tot2<-unique(rawlong.tot[,c("studyid","species","phenovalue","yr1981")]) #CLEAN UP so only unique values across repeating species within studoes
N <- nrow(rawlong.tot2)
y <- rawlong.tot2$phenovalue
Nspp <- length(unique(rawlong.tot2$species)) #newid is character !



#Appendix
startdate<-aggregate(rawlong.tot["year"], rawlong.tot[c("studyid", "intid")], FUN=min)
tryagain2<-merge(tryagain, startdate, by=c("studyid","intid"))
tryagain2<-unique(tryagain2[,c("studyid","intid","meanchange","year")])
a<-ggplot(tryagain2, aes(x=year, y=abs(meanchange)))+geom_point(size=2)+theme_bw()
rawlong.tot$count<-1
length<-aggregate(rawlong.tot["count"], rawlong.tot[c("studyid", "intid","spp")], FUN=sum)
length2<-unique(length[,c("studyid","intid","count")])
tryagain2<-merge(tryagain2, length[1:54,], by=c("studyid","intid"))
b<-ggplot(tryagain2, aes(x=count, y=abs(meanchange)))+geom_point(size=2)+theme_bw()
multiplot(a,b, cols=2)


# Old Scatterplots
ggplot(indiv_intxn, aes(x=stanfit.x, y=stanfit.y))+geom_point(aes(size=2))+geom_vline(xintercept=0, linetype="dashed")+geom_hline(yintercept=0, linetype="dashed")+theme_bw()+theme(legend.position="false")+xlab("resource")+ylab("consumer")

#with labels:
indiv_intxn$pattern<-with(indiv_intxn, ifelse(meanchange<=0, "decreasing", "increasing"))
#indiv_intxn$pattern<-with(indiv_intxn, recode(meanchange,"<0='decreasing';>=0='increasing'"))

indiv_intxn$seasonal_order<-with(indiv_intxn, ifelse(newphenodiff<0, "consumer_first", "resource_first"))
#positive phenodiff= consumer emerges BEFORE resource
#negative phenodiff= resource emerges BEFORE consumer

ggplot(indiv_intxn, aes(x=stanfit.x, y=stanfit.y, label = intid))+geom_vline(xintercept=0, linetype="dashed")+geom_hline(yintercept=0, linetype="dashed")+theme_bw()+theme()+xlab("resource")+ylab("consumer")+geom_text(aes(colour=factor(pattern)), check_overlap=TRUE)+geom_abline(intercept=0, slope=1)

ggplot(subset(indiv_intxn, seasonal_order=="consumer_first"), aes(x=stanfit.x, y=stanfit.y, label = intid))+geom_vline(xintercept=0, linetype="dashed")+geom_hline(yintercept=0, linetype="dashed")+theme_bw()+theme()+xlab("resource")+ylab("consumer")+geom_text(aes(colour=factor(pattern)), check_overlap=TRUE)+geom_abline(intercept=0, slope=1)






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

################################################

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
