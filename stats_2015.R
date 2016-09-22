library(nlme)
library(car)
library(MASS)
library(MuMIn)

rm(list=ls())
row.names=FALSE
# run synchbuildintxn_hk to produce interactions
# THEN run datacleaning
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis")
#mat2<-merge(mat, stud, by="studyid")

total3<-read.csv("int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)

### to deal with pseudo-replication for INTERACTIONS
#remove Keratella cochlearis from either HMK029 or HMK037
#total<-subset(total3, spp2!="Keratella1 cochlearis") # from HMK029
total<-subset(total3, intid!="182") #see right above for reason
total<-subset(total, studyid!="HMK026")
total<-subset(total, studyid!="HMK049") #HMK049- remove because predicted relationship


total3<-total
total3$intxn<-with(total3, paste(spp1,spp2,interaction))
total3<-subset(total3, intxn!="Caterpillar2 spp. Parus4 major predation") #same interaction as HMK027, excluded because shorter time series

total3<-subset(total3, intxn!="Quercus3 robur Caterpillar2 spp. herbivory")

##number of interactions by study
#ubtloc <- unique(total3[,c("studyid","intid")])
#yo<-aggregate(ubtloc["intid"], by=ubtloc[c("studyid")], FUN=length)

##to see whether breakpt possible
bk<-subset(total3, intid=="170" | intid=="171" | intid=="177" | intid=="178" | intid=="179" | intid=="180" | intid=="189" | intid=="195" | intid=="196")
pre<-subset(bk, year<1976)
post<-subset(bk, year>1975)
#sub<-subset(yo3, length>4)
m2<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, method="ML", pre, na.action=na.omit);

ggplot(bk, aes(x=year, y=sqrt(abs(phenodiff))))+geom_point() + facet_wrap(~intid)+stat_smooth(method="lm", se=FALSE, formula = y ~ x + I(x^2))+theme_bw()

#number of years by INTERACTION (not study!)
total4<-na.omit(total3)
yo<-aggregate(total4["year"], by=total4[c("studyid","intid")], FUN=length)
names(yo)[names(yo)=="year"]<-"length"
yo2<-merge(total4,yo, by=c("studyid","intid"))

#first year of study WITH INTERACTION
yo<-aggregate(total4["year"], by=total4[c("studyid","intid")], FUN=min)
names(yo)[names(yo)=="year"]<-"minyear"
yo3<-merge(yo2,yo, by=c("studyid","intid"))

#re-calculate phenodiff from min year
yo3<-yo3[order(yo3$intid,yo3$year),]; yo3<-na.omit(yo3)
best<-data.frame(array(0, c(nrow(yo3), 4)))
names(best)[1]<-"intid"; names(best)[2]<-"year"; names(best)[3]<-"phenodiff_base"; names(best)[4]<-"base"
Bgroups<-unique(yo3$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-yo3[yo3$intid==b[i],]
asdf<-rowcount+(nrow(spp)-1)
best[rowcount:asdf,1]<-b[i]
best[rowcount:asdf,2]<-spp[,c("year")]
best[rowcount:asdf,3]<- with(spp, phenodiff-spp[1,c("phenodiff")])
best[rowcount:asdf,4]<- spp[1,c("phenodiff")]
rowcount<-rowcount+nrow(spp)
}
yo4<-merge(yo3, best, by=c("intid","year"))
yo3<-yo4

#re-calculate phenodiff from 3 years
yo3<-yo3[order(yo3$intid,yo3$year),]; yo3<-na.omit(yo3)
best<-data.frame(array(0, c(nrow(yo3), 3)))
names(best)[1]<-"intid"; names(best)[2]<-"year"; names(best)[3]<-"phenodiff_baseavg"
Bgroups<-unique(yo3$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-yo3[yo3$intid==b[i],]
spp2<-spp[1:3,]
yo<-aggregate(spp2["phenodiff"], by=spp2[c("studyid","intid")], FUN=mean)
asdf<-rowcount+(nrow(spp)-1)
best[rowcount:asdf,1]<-b[i]
best[rowcount:asdf,2]<-spp[,c("year")]
best[rowcount:asdf,3]<- with(spp, phenodiff-yo[1,3])
rowcount<-rowcount+nrow(spp)
}
yo4<-merge(yo3, best, by=c("intid","year"))
yo3<-yo4

#add mismatch interactions- 
#more conservative- if lm fits of species crosses
yo3$mismatch<-0
yo3$mismatch[yo3$intid=="1"] <- 1
yo4<-subset(yo3, intid!="1")
int1<-subset(yo3, intid=="1"); min(int1$phenodiff)
int1$phenodiff<-with(int1, phenodiff+abs(min(int1$phenodiff)))
yo5<-rbind(int1, yo4)
yo3<-yo5

yo3$mismatch[yo3$intid=="192"] <- 1
yo4<-subset(yo3, intid!="192")
int192<-subset(yo3, intid=="192"); min(int192$phenodiff)
int192$phenodiff<-with(int192, phenodiff+abs(min(int192$phenodiff)))
yo5<-rbind(int192, yo4)
yo3<-yo5
yo3$mismatch[yo3$intid=="175"] <- 1
yo3$mismatch[yo3$intid=="232"] <- 1
yo3$mismatch[yo3$intid=="235"] <- 1

#redefined as any time series with one year minimum above <>0
yo3$mismatch<-0
yo3$mismatch[yo3$intid=="1"] <- 1
yo3$mismatch[yo3$intid=="19"] <- 1
yo3$mismatch[yo3$intid=="144"] <- 1
yo3$mismatch[yo3$intid=="146"] <- 1
yo3$mismatch[yo3$intid=="165"] <- 1
yo3$mismatch[yo3$intid=="168"] <- 1
yo3$mismatch[yo3$intid=="171"] <- 1
yo3$mismatch[yo3$intid=="175"] <- 1
yo3$mismatch[yo3$intid=="178"] <- 1
yo3$mismatch[yo3$intid=="179"] <- 1
yo3$mismatch[yo3$intid=="182"] <- 1
yo3$mismatch[yo3$intid=="183"] <- 1
yo3$mismatch[yo3$intid=="184"] <- 1
yo3$mismatch[yo3$intid=="186"] <- 1
yo3$mismatch[yo3$intid=="189"] <- 1
yo3$mismatch[yo3$intid=="190"] <- 1
yo3$mismatch[yo3$intid=="192"] <- 1
yo3$mismatch[yo3$intid=="195"] <- 1
yo3$mismatch[yo3$intid=="196"] <- 1
yo3$mismatch[yo3$intid=="197"] <- 1
yo3$mismatch[yo3$intid=="201"] <- 1
yo3$mismatch[yo3$intid=="207"] <- 1
yo3$mismatch[yo3$intid=="208"] <- 1
yo3$mismatch[yo3$intid=="209"] <- 1
yo3$mismatch[yo3$intid=="218"] <- 1
yo3$mismatch[yo3$intid=="219"] <- 1
yo3$mismatch[yo3$intid=="220"] <- 1
yo3$mismatch[yo3$intid=="221"] <- 1
yo3$mismatch[yo3$intid=="222"] <- 1
yo3$mismatch<-as.factor(yo3$mismatch)

####################################################################################
##what is the LONG-TERM trend in the annual variation in interacting spp' phenology
# tells us about shift in synchrony NOT whether switch happened
#phenofreq should account for the fact that some studies only measured phenology weekly or monthly creating large differences
#minyear and length correlated (-0.71)
# ESA abstract based on studies up to and including HMK036
yo3$phenofreq<-as.factor(yo3$phenofreq)
yo3$taxatype<-as.factor(yo3$taxatype)
yo3$interaction[yo3$interaction=="poll"] <- "pollination"
yo3$interaction<-as.factor(yo3$interaction)

###change in absolute relative timing ###
m2<-lme(sqrt(abs(phenodiff))~phenofreq+minyear+year, random=~1|studyid/intid, method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); plot(m2)


m1<-lme(sqrt(abs(phenodiff))~phenofreq+minyear+year, random=~1|studyid/intid, method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")); summary(m1); # mismatch!="1
allmodels<-dredge(m1, rank="AIC"); allmodels # all model combinations
suballmodels<-subset(allmodels, delta<2); suballmodels # only those models with <2 deltaAIC
m1avg<-model.avg(suballmodels)
m2<-get.models(allmodels, 1)[[1]]; summary(m2) # best model
m4<-lme(sqrt(abs(phenodiff))~phenofreq+year, random=~1|studyid/intid, method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4) # mismatch!="1
m2<-update(m4, ~.-year); anova(m4,m2)

###change in absolute relative timing WHILE ACCOUNTING FOR FIRST YEAR DIFFERENCE ###
m4<-lme(sqrt(abs(phenodiff))~abs(base)+year, random=~1|studyid/intid, method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4) #
OR
yo4<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
yo4$base_c<-scale(yo4$base, center=TRUE, scale=FALSE) #center variable by subracting mena of all data points from each individual data point
yo4$year_c<-scale(yo4$year, center=TRUE, scale=FALSE)
!!!! scale(A, center=TRUE, scale=TRUE) # to get zscores in one line (data point-mean/sd) 
m4<-lme(sqrt(abs(phenodiff))~abs(base_c)*year_c, random=~1|studyid/intid, method="ML", yo4, na.action=na.omit); summary(m4) #

###change in absolute relative timing RELATIVE TO FIRST YEAR ###
#incl. covariates first because variance structure affected by model AND can't assess importance of covariates if model doesn't meet assumptions
m1<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year, random=~1|studyid/intid, method="REML", data=subset(yo3, length>4), na.action=na.omit); summary(m1)
m2<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year, random=~1|studyid/intid, weights=varExp(form=~year), method="REML", data=subset(yo3, length>4), na.action=na.omit); 
m3<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year, random=~1|studyid/intid, weights=varExp(form=~year, fixed=0.5), method="REML", data=subset(yo3, length>4), na.action=na.omit); 
m4<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year, random=~1|studyid/intid, weights=varPower(form=~year), method="REML", data=subset(yo3, length>4), na.action=na.omit);
m5<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year, random=~1|studyid/intid, weights=varPower(form=~year, fixed=0.5), method="REML", data=subset(yo3, length>4), na.action=na.omit);
m6<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year, random=~1|studyid/intid, weights=varFixed(~year), method="REML", data=subset(yo3, length>4), na.action=na.omit);
anova(m1,m4, m5, m6)

m4<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")); summary(m4)
allmodels<-dredge(m4, rank="AIC"); allmodels # all model combinations
suballmodels<-subset(allmodels, delta<2); suballmodels 
m2<-get.models(allmodels, 1)[[1]]; summary(m2) # best model
avg<-model.avg(suballmodels); summary(avg) 
m4<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")); summary(m4)
m2<-update(m4, ~.-year); anova(m4,m2)

#type of taxa
m4<-lme(sqrt(abs(phenodiff))~phenofreq+taxatype*year, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4) # mismatch!="1
m4<-lme(sqrt(abs(phenodiff))~phenofreq+year, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & taxatype=="spp"), na.action=na.omit); summary(m4) # mismatch!="1
m2<-update(m4, ~.-year); anova(m4,m2)

m4<-lme(sqrt(abs(phenodiff_base))~taxatype+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")); summary(m4)

#without mismatch studies- 
yo4<-subset(yo3, intid!="187" & intid!="188")
m4<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, method="ML", data=subset(yo4, length>4 & mismatch!="1"), na.action=na.omit); summary(m4)

m4<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, method="ML", data=subset(yo4, length>4 & mismatch!="1"), na.action=na.omit); summary(m4)

##just temp sensitivity studies
mom2<-unique(mom[,c("studyid","intid","year")])
yess<-merge(yo3, mom2, by=c("studyid","intid","year"))
m4<-lme(abs(phenodiff)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yess, length>4 & mismatch!="1" & studyid!="SET005"), na.action=na.omit); summary(m4)

m4<-lme(abs(phenodiff_base)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-phenofreq); anova(m4,m2) # SIG SO KEEP!!

#first year and length are correlated!
m4<-lme(abs(phenodiff_base)~length+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-length); anova(m4,m2) # LENGTH NS and quadratic term

m4<-lme(abs(phenodiff_base)~minyear+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-minyear); anova(m4,m2) #NS so don't keep!

#check effect of excluding weekly and monthly surveys DON"T INCLUDE PHENOFREQ in model bc then only level of factor
m4<-lme(abs(phenodiff)~minyear+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & phenofreq=="daily"), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-minyear); anova(m4,m2) #now Sig so keep!

#check effect of HMK031 which has really temporally segragated events
m4<-lme(abs(phenodiff)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-minyear); anova(m4,m2) #now Sig so keep!

m4<-lme(abs(phenodiff)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & studyid!="HMK031"), na.action=na.omit); summary(m4)

# Temporal autocorrelation
m4<-lme(abs(phenodiff)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), data=subset(yo3, length>4), na.action=na.omit)
m5<-lme(abs(phenodiff)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), correlation=corAR1(), data=subset(yo3, length>4), na.action=na.omit)
AIC(m4, m5)

#Exclude studies < 10 years
m4<-lme(abs(phenodiff)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>=10), na.action=na.omit); summary(m4)

## 
phenodiff z score
newdata<-with(yo3, aggregate(phenodiff_base, by=list(studyid,intid), FUN=mean, na.rm=T))
names(newdata)[1]<-"studyid"; names(newdata)[2]<-"intid"; names(newdata)[3]<-"phenodiff_mean"

newdata2<-with(yo3, aggregate(phenodiff_base, by=list(studyid,intid), FUN=sd, na.rm=T))
names(newdata2)[1]<-"studyid"; names(newdata2)[2]<-"intid"; names(newdata2)[3]<-"phenodiff_sd"

newdata3<-merge(yo3, newdata, by=c("studyid","intid"))
newdata4<-merge(newdata3, newdata2, by=c("studyid","intid"))

newdata4$phenodiff_z<-with(newdata4, (phenodiff_base-phenodiff_mean)/phenodiff_sd)
newdata<-newdata4

m4<-lme(abs(phenodiff_z)~year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(newdata, length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)

# VARIANCE #
newdata<-with(yo3, aggregate(phenodiff, by=list(studyid,intid), FUN=mean, na.rm=T))
names(newdata)[1]<-"studyid"; names(newdata)[2]<-"intid"; names(newdata)[3]<-"phenodiff_mean"
newdata3<-merge(yo3, newdata, by=c("studyid","intid"))
newdata3$phenodiff_var<-with(newdata3, (phenodiff-phenodiff_mean))

m4<-lme(abs(phenodiff_var)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(newdata3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4) 

#phenodiff_base based on 3 year mean
m4<-lme(abs(phenodiff_baseavg)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)

#includ whether temp has changed (cat)
yo5$tempcat<-as.factor(yo5$tempcat)
m4<-lme(abs(phenodiff_base)~phenofreq+year+tempcat, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo5, length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)

m4<-lme(abs(phenodiff_base)~phenofreq+year+tempcat, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo5, tempcat=="nodata", length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)

####
#study characteristics
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(yo3, stud[,c("studyid","biome2","terrestrial", "terrestrial2", "big_interaction")], by="studyid")
yep$biome2<-as.factor(yep$biome2)
yep$terrestrial<-as.factor(yep$terrestrial)
yep$terrestrial2<-as.factor(yep$terrestrial2)

# NEGATIVE trophic interaction
yep2<-subset(yep, interaction=="predation" | interaction=="herbivory")
m4<-lme(sqrt(abs(phenodiff))~terrestrial+year*interaction, random=~1|studyid/intid, method="ML", data=subset(yep2, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4); AIC(m4)
m4<-lme(sqrt(abs(phenodiff))~terrestrial+year*interaction, random=~1|studyid/intid, method="ML", weights=varComb(varIdent(form=~1|interaction), varPower(form=~year)), data=subset(yep2, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); plot(m4); summary(m4)
m2<-update(m4, ~.-year:interaction); anova(m2,m4)

yep2<-subset(yep, interaction=="predation" | interaction=="herbivory")
m4<-lme(sqrt(abs(phenodiff))~terrestrial+year*interaction, random=~1|studyid/intid, method="ML", data=subset(yep2, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4); AIC(m4)

ggplot(data=subset(yep2, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(x=year, y=sqrt(abs(phenodiff)), colour=factor(interaction))) +  geom_point(aes(shape=factor(interaction))) + geom_smooth(method="lm", se=FALSE, aes(fill = factor(interaction)))+theme_bw()

m4<-lme(sqrt(abs(phenodiff))~terrestrial+year, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & interaction=="predation"), na.action=na.omit); plot(m4); summary(m4)
m2<-update(m4, ~.-year); anova(m2,m4)

m4<-lme(sqrt(abs(phenodiff))~terrestrial+year, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & interaction=="herbivory"), na.action=na.omit); plot(m4); summary(m4)
m2<-update(m4, ~.-year); anova(m2,m4)

#differences across biomes?
m1<-lme(sqrt(abs(phenodiff))~year*biome2, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-year*biome2); anova(m2,m1)

ggplot(data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(x=year, y=sqrt(abs(phenodiff)), colour=factor(biome2))) +  geom_point(aes(shape=factor(biome2))) + geom_smooth(method="lm", se=FALSE, aes(fill = factor(biome2)))+theme_bw()

#Aquatic vs terrestrial # All studies with non-daily surveys are aquatic!!
yep$terrestrial<-as.factor(yep$terrestrial)
m1<-lme(sqrt(abs(phenodiff_base))~year*terrestrial, random=~1|studyid/intid, method="REML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
m2<-lme(sqrt(abs(phenodiff_base))~year*terrestrial, random=~1|studyid/intid, weights=varPower(form=~year), method="REML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit);
m3<-lme(sqrt(abs(phenodiff_base))~year*terrestrial, random=~1|studyid/intid, weights=varIdent(form=~1|terrestrial), method="REML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit);
m4<-lme(sqrt(abs(phenodiff_base))~year*terrestrial, random=~1|studyid/intid, weights=varComb(varIdent(form=~1|terrestrial), varPower(form=~year)), method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)
anova(m1,m2,m3, m4)
m2<-update(m4, ~.-year:terrestrial); anova(m2,m4)

ggplot(data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(x=year, y=sqrt(abs(phenodiff_base)), colour=factor(terrestrial))) +  geom_point(shape=1) + geom_smooth(method="lm", se=FALSE, aes(fill = factor(terrestrial)))+theme_bw()

m4<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="aquatic"), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-year); anova(m4,m2)

m4<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="terrestrial"), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-year); anova(m4,m2)

#fw vs marine
m4<-lme(sqrt(abs(phenodiff))~year*terrestrial2, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="aquatic"), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-year); anova(m4,m2)

m4<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial2=="fresh"), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-year); anova(m4,m2)



#Broad class of interactions?
m1<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year+big_interaction, random=~1|studyid/intid, weights=varIdent(form=~1|big_interaction), method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
#m1<-lme(abs(phenodiff_base)~phenofreq+year+big_interaction, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yep, length>4), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-big_interaction); anova(m2,m1)

##### What is the interannual variability in pheno diff #####

m4<-lme(abs(phenodiff)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4), na.action=na.omit); summary(m4)

not<-subset(yo3, length>4); not<-na.omit(not)
not$res<-residuals(m4)
m5<-lme(res~year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=not); summary(m5)

#############################################################################
## Has rate of phenological change differed between interacting species ? ###
#############################################################################
so<-read.csv("spp_phenodata2.csv", header=TRUE, na.strings="NA", as.is=TRUE)
so<-so[order(so$studyid),]
### to deal with pseudo-replication for INTERACTIONS
so<-subset(so, intid!="182")
so<-subset(so, studyid!="HMK026")
so<-subset(so, studyid!="HMK049") #HMK049- remove because predicted relationship

#intxn="Caterpillar2 spp. Parus4 major predation") =intid==221
#intxn!="Quercus3 robur Caterpillar2 spp. herbivory") intid==222

# Only years with interaction #
# to incorporate number of years by study
total3<-read.csv("int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
total4<-na.omit(total3)
yo<-aggregate(total4["year"], by=total4[c("studyid","intid")], FUN=length)
names(yo)[names(yo)=="year"]<-"length"
so2<-merge(so,yo, by=c("studyid","intid"))
so3<-merge(so2,unique(total4[,c("studyid","intid","phenofreq")]), by=c("studyid","intid"))
so3$intxn<-with(so3, paste(species,int_type, spp))

ubtloc <- unique(total4[,c("studyid","intid","year")])
so3<-merge(so3, ubtloc, by=c("studyid","intid","year"))
so3<-so3[order(so3$studyid, so3$spp),]

####
#study characteristics
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(so3, stud[,c("studyid","short_site","biome","terrestrial","big_interaction")], by="studyid")
so3<-yep

m1<-lme(phenovalue~year*spp, random=~1|studyid/intid/species, method="ML", data=subset(so3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intid!="221" & intid!="222"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-year:spp); anova(m2,m1)

m1<-lme(phenovalue~year*spp, random=~1|studyid/intid/species, method="ML", data=subset(so3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intid!="221" & intid!="222" & terrestrial=="aquatic"), na.action=na.omit); summary(m1)

sup<-subset(so3, int_type=="predation" | int_type=="herbivory")
m1<-lme(phenovalue~year*spp*int_type, random=~1|studyid/intid/species, method="ML", data=subset(sup, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intid!="221" & intid!="222"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-year:spp:int_type); anova(m2,m1)

m1<-lme(phenovalue~year*spp, random=~1|studyid/intid/species, method="ML", data=subset(sup, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intid!="221" & intid!="222" & int_type=="predation"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-year:spp); anova(m2,m1)

m1<-lme(phenovalue~year*spp, random=~1|studyid/intid/species, method="ML", data=subset(sup, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intid!="221" & intid!="222" & int_type=="herbivory"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-year:spp); anova(m2,m1)

m1<-lme(phenovalue~year, random=~1|studyid/intid/species, method="ML", data=subset(sup, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intid!="221" & intid!="222" & int_type=="predation" & spp=="spp1"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-year:spp); anova(m2,m1)
# spp1 = negative species in trophic interaction and spp2= positive species

#calculate pheno shift across unique spp (some spp are in mult interacitons) i.e. remove intid as random effect + RENAME some
so3<-so3[order(so3$species),]
so3$species[so3$species=="Ficedula1 albicollis"] <- "Ficedula albicollis"
so3$species[so3$species=="Ficedula2 albicollis"] <- "Ficedula albicollis"
so3$species[so3$species=="Keratella2 cochlearis"]<-"Keratella cochlearis"
so3$species[so3$species=="Parus1 major"]<-"Parus major"
so3$species[so3$species=="Parus2 major"]<-"Parus major"
so3$species[so3$species=="Parus3 major"]<-"Parus major"
so3$species[so3$species=="Parus4 major"]<-"Parus major"
so3$species[so3$species=="Quercus1 robur"]<-"Quercus robur"
so3$species[so3$species=="Quercus2 robur"]<-"Quercus robur"
so3$species[so3$species=="Quercus3 robur"]<-"Quercus robur"
so3$species[so3$species=="Mnemiopsis2 leidyi"]<-"Mnemiopsis leidyi"
so3$species[so3$species=="Rissa1 tridactyla"]<-"Rissa tridactyla"
so3$species[so3$species=="Rissa2 tridactyla"]<-"Rissa tridactyla"
so3$species[so3$species=="Parus2 caeruleus"]<-"Parus caeruleus"
so3<-subset(so3, species!="Mnemiopsis1 leidyi")
so3<-so3[order(so3$species, so3$short_site, so3$year),]

#to get rid of species with multiple roles/study
so4<-subset(so3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
sure<-unique(so4[,c("studyid","year","species","phenovalue","phenofreq","short_site", "terrestrial")])

m1<-lme(phenovalue~year, random=~1|species/short_site, method="ML", sure, na.action=na.omit); summary(m1)
m2<-update(m1, ~.-year); anova(m2,m1)
m1<-lme(phenovalue~year*terrestrial, random=~1|species/short_site, method="ML", sure, na.action=na.omit); summary(m1)
 m2<-update(m1, ~.-year:terrestrial); anova(m2,m1)

m1<-lme(phenovalue~year, random=~1|species/short_site, method="ML", subset(sure, terrestrial=="aquatic"), na.action=na.omit); summary(m1) 
 m2<-update(m1, ~.-year); anova(m2,m1)

ggplot(sure, aes(x=year, y=phenovalue, colour=factor(terrestrial))) +  geom_point(aes(shape=factor(terrestrial))) + geom_smooth(method="lm", se=FALSE, aes(fill = factor(terrestrial)))+theme_bw() 

# Does pheno advance change by taxonomic group?
taxer <- read.csv("taxa.csv", header=TRUE)
taxer$latbi <- paste(taxer$genus, taxer$species)
taxer<-taxer[,c("studyid","taxa","latbi")]
names(taxer)[3]<-"species"
so4<-merge(so3, taxer, by=c("studyid","species"))

m1<-lme(phenovalue~year+spp+taxa, random=~1|studyid/intid/species, method="REML", data=subset(so4, length>4), na.action=na.omit); summary(m1)
m2<-lme(phenovalue~year+spp+taxa, random=~1|studyid/intid/species, method="ML", data=subset(so4, length>4), weights=varIdent(form=~1|spp), na.action=na.omit); summary(m2); anova(m1, m2)
m3<-update(m2, ~.-taxa); anova(m3,m2)

so5<-subset(so4, length>4)
sun<-unique(so5[,c("studyid","species","slope","taxa")])
m1<-lme(slope~taxa, random=~1|studyid, method="ML", data=subset(sun, species!="Thermocyclops oithonoides" & species!="Plant1 spp."), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-taxa); anova(m1,m2)

---------------------------------------------------------------


