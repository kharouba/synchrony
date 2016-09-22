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
total<-subset(total3, spp2!="Keratella1 cochlearis") # from HMK029
#remove HMK026
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
#total3<-subset(total3, year<1976)
#bk<-subset(total3, intid=="170" | intid=="171" | intid=="177" | intid=="178" | intid=="179" | intid=="180" | intid=="189" | intid=="195" | intid=="196")
#total3<-subset(total3, year>1975)
#sub<-subset(yo3, length>4)

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

#add mismatch interactions- redefined as any time series with one year minimum above <>0
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
yo3$interaction[yo3$interaction=="poll"] <- "pollination"
yo3$interaction<-as.factor(yo3$interaction)

###change in absolute relative timing ###
m1<-lme(abs(phenodiff)~phenofreq+minyear+year, random=~1|studyid/intid, method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4); qqnorm(m4) # mismatch!="1

m1<-lme(abs(phenodiff)~year, random=~1|studyid/intid, method="REML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
m2<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, weights=varExp(form=~year), method="REML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); 
m3<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, weights=varExp(form=~year, fixed=0.5), method="REML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); 
m4<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, weights=varPower(form=~year), method="REML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit);
m5<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, weights=varPower(form=~year, fixed=0.5), method="REML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit);
m6<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, weights=varFixed(~year), method="REML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit);
anova(m1,m4, m5, m6)

m4<-lme(abs(phenodiff)~phenofreq+minyear+year,  weights=varPower(form=~year), random=~1|studyid/intid, method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4); qqnorm(m4) # mismatch!="1
allmodels<-dredge(m4, rank="AICc"); allmodels # all model combinations
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
m1<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, method="REML", data=subset(yo3, length>4), na.action=na.omit); summary(m1)
m2<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, weights=varExp(form=~year), method="REML", data=subset(yo3, length>4), na.action=na.omit); 
m3<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, weights=varExp(form=~year, fixed=0.5), method="REML", data=subset(yo3, length>4), na.action=na.omit); 
m4<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, weights=varPower(form=~year), method="REML", data=subset(yo3, length>4), na.action=na.omit);
m5<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, weights=varPower(form=~year, fixed=0.5), method="REML", data=subset(yo3, length>4), na.action=na.omit);
m6<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, weights=varFixed(~year), method="REML", data=subset(yo3, length>4), na.action=na.omit);
anova(m1,m4, m5, m6)

m4<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-year); anova(m4,m2)

#without mismatch studies- 
yo4<-subset(yo3, intid!="187" & intid!="188")
m4<-lme(abs(phenodiff_base)~phenofreq+year, random=~1|studyid/intid, method="ML", data=subset(yo4, length>4 & mismatch!="1"), na.action=na.omit); summary(m4)

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

m4<-lme(abs(phenodiff_var)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(newdata3, length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4) 

#phenodiff_base based on 3 year mean
m4<-lme(abs(phenodiff_baseavg)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo3, length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)

#includ whether temp has changed (cat)
yo5$tempcat<-as.factor(yo5$tempcat)
m4<-lme(abs(phenodiff_base)~phenofreq+year+tempcat, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo5, length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)

m4<-lme(abs(phenodiff_base)~phenofreq+year+tempcat, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yo5, tempcat=="nodata", length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)

####
#study characteristics
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(yo3, stud[,c("studyid","biome","terrestrial","big_interaction")], by="studyid")

#differences across biomes?
T=tropical, S=subtropical/mediteranean, Te=temperatre, B=boreal, A=arctic/alpine, M=marine, f=freshwater
m1<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year+biome, random=~1|studyid/intid, method="ML", weights=varPower(form=~year), data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-biome); anova(m2,m1)

#Aquatic vs terrestrial # All studies with non-daily surveys are aquatic!!
m1<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year+terrestrial, random=~1|studyid/intid, method="REML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
m2<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year+terrestrial, random=~1|studyid/intid, weights=varPower(form=~year), method="REML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit);
m3<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year+terrestrial, random=~1|studyid/intid, weights=varIdent(form=~1|terrestrial), method="REML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit);
m4<-lme(sqrt(abs(phenodiff_base))~phenofreq+minyear+year+terrestrial, random=~1|studyid/intid, weights=varComb(varIdent(form=~1|terrestrial), varPower(form=~year)), method="ML", data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)
anova(m1,m2,m3, m4)
m2<-update(m4, ~.-terrestrial); anova(m2,m4)

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
#remove Keratella cochlearis from either HMK029 or HMK037
total<-subset(so, species!="Keratella1 cochlearis") # from HMK029
#remove HMK026
total<-subset(total, studyid!="HMK026")
total<-subset(total, studyid!="HMK049") #HMK049- remove because predicted relationship
so<-total


## Across all years with data ##
# to incorporate number of years by study
total3<-read.csv("int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
yo<-aggregate(total3["year"], by=total3[c("studyid","intid")], FUN=length)
names(yo)[names(yo)=="year"]<-"length"
so2<-merge(so,yo, by=c("studyid","intid"))
so3<-merge(so2,unique(total3[,c("studyid","intid","phenofreq")]), by=c("studyid","intid"))
so3$intxn<-with(so3, paste(species,int_type, spp))
#so3<-subset(so3, intxn!="Parus4 major predation spp2") #or HMK048same interaction as HMK027, excluded because shorter time series
#so3<-subset(so3, intxn!="Quercus3 robur herbivory spp1")


m1<-lme(phenovalue~year*spp, random=~1|studyid/intid/species, method="ML", data=subset(so3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intxn!="Parus4 major predation spp2" & intxn!="Quercus3 robur herbivory spp1"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-phenofreq); anova(m2,m1)
m2<-update(m1, ~.-year:spp); anova(m2,m1)
#m2<-lme(phenovalue~year*spp, random=~1|studyid/intid, method="ML",correlation=corAR1(), data=subset(so2, length>4), na.action=na.omit); summary(m2)
m1<-lme(phenovalue~year+spp, random=~1|studyid/intid/species, method="ML", data=subset(so3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)

#only years with interaction
total4<-na.omit(total3)
ubtloc <- unique(total4[,c("studyid","intid","year")])
so3<-merge(so3, ubtloc, by=c("studyid","intid","year"))
so3<-so3[order(so3$studyid, so3$spp),]

####
#study characteristics
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(so3, stud[,c("studyid","short_site","biome","terrestrial","big_interaction")], by="studyid")
so3<-yep

m1<-lme(phenovalue~year*spp, random=~1|studyid/intid/species, method="ML", data=subset(so3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intxn!="Parus4 major predation spp2" & intxn!="Quercus3 robur herbivory spp1"), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-year:spp); anova(m2,m1)
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
sure<-unique(so4[,c("studyid","year","species","phenovalue","phenofreq","short_site")])

m1<-lme(phenovalue~year, random=~1|species/short_site, method="ML", sure, na.action=na.omit); summary(m1)
m2<-lme(phenovalue~year, random=~1|studyid/species, method="REML", data=subset(sure, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), weights=varIdent(form=~1|spp), na.action=na.omit); m2<-update(m1, ~.-year); anova(m2,m1) 

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
--------------------------------------------
				 ENV DATA
--------------------------------------------
## HAS TEMP CHANGED?
clim<-read.csv("climate.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
clim<-subset(clim, phenophase!="start" & extra!="10")
clim$envvalue<-as.numeric(clim$envvalue)
clim$envfactor[clim$envfactor=="temperaure"] <- "temperature"
sites<-subset(clim, studyid=="HMK018" | studyid=="HMK019" | studyid=="HMK023" & site!="tomakomai")
nosites<-subset(clim, site=="tomakomai" | studyid!="HMK018" & studyid!="HMK019" & studyid!="HMK023")
nosites<-nosites[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")]

new<-with(sites, aggregate(envvalue, by=list(studyid, year, species), FUN=mean, na.rm=T)) # across sites
names(new)[1]<-"studyid"; names(new)[2]<-"year"; names(new)[3]<-"species"; names(new)[4]<-"envvalue"
new<-new[order(new$studyid),]
sites2<-merge(sites[,c("studyid","envfactor","envunits","envtype","year","species")], new, by=c("studyid","year","species"))
sites3<-unique(sites2[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")])

clim2<-rbind(nosites, sites3) # all years with data
#merge with spp data so calculating env change only over the years of interaction
total3<-read.csv("int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
total3<-na.omit(total3)
clim3<-merge(clim2, total3[,c("studyid","year")], by=c("studyid","year"))
clim3<-unique(clim3[,c("studyid","year","envfactor","envunits","envtype","species","envvalue")])
clemo<-unique(clim3[,c("studyid","year")])
mi<-aggregate(clemo["year"], by=clemo[c("studyid")], FUN=length); std(mi$year) #number of years
mimin<-aggregate(clemo["year"], by=clemo[c("studyid")], FUN=min) #number of years
mimax<-aggregate(clemo["year"], by=clemo[c("studyid")], FUN=max)
mimax$range<-mimax$year-mimin$year; mean(mimax$range); std(mimax$range)
std <- function(x) sd(x)/sqrt(length(x))


# ONLY FOR YEARS WITH SPECIES DATA
m1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype=="air" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-year); anova(m3,m1)
#water
m1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype=="water" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-year); anova(m3,m1)

m1<-lme(envvalue~envtype+year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype!="ground" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & envunits!="doy"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-envtype); anova(m3,m1)

################################################################
#### For interactions with same cue, does sensitivity differ?###
################################################################
#so<-read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/spp_phenodata2.csv", header=TRUE, na.strings="NA", as.is=TRUE)
#so<-so[order(so$studyid),]
#bo<-merge(lets[,c("studyid","year","envvalue","z")], so, by=c("studyid","year"));
#bo<-bo[order(bo$studyid, bo$studyid),]

#temp only
HMK029- stratification better predictor than temp # NOW INCLUDED_  nov 2014
HMK038- diff temp cues but z scores are only r=0.49 correlated
HMK035- no correlation in paper and temp measured as doy NOT C
HMK031- temp measured as doy
clim8<-subset(tot, envfactor=="temperature" | studyid=="HMK016" | studyid=="HMK034");
#bo<-merge(clim8[,c("studyid","year","envvalue","z")], so, by=c("studyid","year"));
clim8<-clim8[order(clim8$studyid, clim8$species),]
#bo2<-unique(bo[,c("studyid","year","envvalue","z","species","phenovalue","intid","spp")])
write.csv(clim8, "/users/kharouba/google drive/UBC/synchrony project/analysis/tempsens_nov.csv")

mom<-read.csv("tempsens_nov.csv", header=TRUE, na.strings="NA", as.is=TRUE)
mom<-na.omit(mom)
mom$spp<-as.factor(mom$spp)
#mom<-read.csv("tempsens.csv", header=TRUE, na.strings="NA", as.is=TRUE)
#mom[mom$intid=="202","intid"]<-"207"
#mom[mom$intid=="203","intid"]<-"208"
#mom[mom$intid=="204","intid"]<-"209"

Bgroups<-unique(mom$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0, c(length(b), 5)))
names(new)[1] = "intid"; names(new)[2] = "studyid"; names(new)[3]="spp1"; names(new)[4]<-"spp2"; names(new)[5] = "intp";                                                                                                                                                                                                                    
for(i in 1:length(b)) {
spp<-mom[mom$intid==b[i],]; spp<-na.omit(spp)
spp$intid<-as.character(spp$intid);
m1<-lm(phenovalue~z*species, data=spp, na.action=na.omit);
 new[i,1]<-spp$intid[1]
new[i,2]<-spp$studyid[1]
new[i,3]<-summary(m1)$coefficients[3]
new[i,5]<-anova(m1)$'Pr(>F)'[5]
}


m1<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
m2<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m2)
m3<-update(m2, ~.-z:spp); anova(m2,m3)
only include unique spp! i.e. remove intid as random effect
m2<-lme(phenovalue~z+spp, random=~1|studyid/species, method="ML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m2)
m3<-update(m2, ~.-z); anova(m2,m3)
#remove HMK038
m2<-lme(phenovalue~z+spp, random=~1|studyid/intid, method="REML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK038"), na.action=na.omit); summary(m2)

#Prep for overall model:
#Step 1 get temp sens per species
rowcount<-1
Bgroups<-unique(mom$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0, c(80, 6)))
names(new)[1] = "studyid"; names(new)[2]<-"species"; names(new)[3]<-"spp"; names(new)[4]<-"intid"; names(new)[5]<-"doybytemp"; names(new)[6]<-"rsquared"
for(i in 1:length(b)) { 
yo<-mom[mom$intid==b[i],]
Cgroups<-unique(yo$species); c<-Cgroups; c<-as.character(c)
for(j in 1:length(c)){
	spp<-yo[yo$species==c[j],]
new[rowcount+length(c),1]<-spp[1,c("studyid")]
new[rowcount+length(c),2]<-spp[1,c("species")]
new[rowcount+length(c),3]<-spp[1,c("spp")]
new[rowcount+length(c),4]<-spp[1,c("intid")]
m1<-with(spp, lm(phenovalue~z))
new[rowcount+length(c),5]<-summary(m1)$coefficients[2]
new[rowcount+length(c),6]<-summary(m1)$r.squared
rowcount<-rowcount+1
}
}
new<-subset(new, studyid!="0")

Bgroups<-unique(new$intid); b<-Bgroups; b<-as.character(b)
no<-data.frame(array(0, c(length(b), 6)))
names(no)[1] = "studyid"; names(no)[2]<-"intid"; names(no)[3]<-"doybytemp";  names(no)[4]<-"rsquarediff"; names(no)[5]<-"proprsquare"; names(no)[6]<-"medrsquare"
rowcount<-1
for(i in 1:length(b)) { 
yo<-new[new$intid==b[i],]
asdf<-(rowcount+nrow(yo))-1
no[rowcount:asdf,1]<-yo[1,c("studyid")]
no[rowcount:asdf,2]<-yo[1,c("intid")]
spp1<-subset(yo, spp=="1")
spp2<-subset(yo, spp=="2")
no[rowcount:asdf,3]<-spp1[1,5]-spp2[1,5]
no[rowcount:asdf,4]<-spp1[1,6]-spp2[1,6]
no[rowcount:asdf,6]<-median(c(spp1[1,6],spp2[1,6]))

if(spp1[1,6]>spp2[1,6]){
no[rowcount:asdf,5]<-spp2[1,6]/spp1[1,6]	
}

if(spp2[1,6]>spp1[1,6]){
no[rowcount:asdf,5]<-spp1[1,6]/spp2[1,6]	
}
rowcount<-rowcount+nrow(yo)
}
no<-unique(no[,c("studyid","intid","doybytemp","rsquarediff","proprsquare","medrsquare")])

yum<-merge(no, yo3, by=c("studyid","intid"))

NEED TO SQUARE ROOT
m4<-lme(abs(phenodiff_base)~phenofreq+year+medrsquare, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(yum, length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m4)

#pheno shift
rowcount<-1
Bgroups<-unique(mom$intid); b<-Bgroups; b<-as.character(b)
nup<-data.frame(array(0, c(80, 5)))
names(nup)[1] = "studyid"; names(nup)[2]<-"species"; names(nup)[3]<-"spp"; names(nup)[4]<-"intid"; names(nup)[5]<-"doybyyr"
for(i in 1:length(b)) { 
yo<-mom[mom$intid==b[i],]
Cgroups<-unique(yo$species); c<-Cgroups; c<-as.character(c)
for(j in 1:length(c)){
	spp<-yo[yo$species==c[j],]
nup[rowcount+length(c),1]<-spp[1,c("studyid")]
nup[rowcount+length(c),2]<-spp[1,c("species")]
nup[rowcount+length(c),3]<-spp[1,c("spp")]
nup[rowcount+length(c),4]<-spp[1,c("intid")]
m1<-with(spp, lm(phenovalue~year))
nup[rowcount+length(c),5]<-summary(m1)$coefficients[2]
rowcount<-rowcount+1
}
}
nup<-subset(nup, studyid!="0")
son<-merge(new, nup, by=c("studyid","species","spp","intid"))

#temp change BY SPECIES
rowcount<-1
Bgroups<-unique(mom$intid); b<-Bgroups; b<-as.character(b)
sup<-data.frame(array(0, c(80, 5)))
names(sup)[1] = "studyid"; names(sup)[2]<-"species"; names(sup)[3]<-"spp"; names(sup)[4]<-"intid"; names(sup)[5]<-"tempbyyr"
for(i in 1:length(b)) { 
yo<-mom[mom$intid==b[i],]
Cgroups<-unique(yo$species); c<-Cgroups; c<-as.character(c)
for(j in 1:length(c)){
	spp<-yo[yo$species==c[j],]
sup[rowcount+length(c),1]<-spp[1,c("studyid")]
sup[rowcount+length(c),2]<-spp[1,c("species")]
sup[rowcount+length(c),3]<-spp[1,c("spp")]
sup[rowcount+length(c),4]<-spp[1,c("intid")]
m1<-with(spp, lm(z~year))
sup[rowcount+length(c),5]<-summary(m1)$coefficients[2]
rowcount<-rowcount+1
}
}
sup<-subset(sup, studyid!="0")
son<-merge(son, sup, by=c("studyid","species","spp","intid"))

#Step 2 get synchrony change per interaction
Bgroups<-unique(yo3$intid); b<-Bgroups; b<-as.character(b)
yes<-data.frame(array(0, c(length(b), 4)))
names(yes)[1] = "studyid"; names(yes)[2]<-"intid"; names(yes)[3]<-"phenodiffbyyr"; names(yes)[4]<-"p"
rowcount<-1
for(i in 1:length(b)) { 
sub<-yo3[yo3$intid==b[i],]
yes[i,1]<-sub[1,c("studyid")]
yes[i,2]<-sub[1,c("intid")]
m1<-with(sub, lm(abs(phenodiff)~year))
yes[i,3]<-summary(m1)$coefficients[2]
yes[i,4]<-summary(m1)$coefficients[8]
rowcount<-rowcount+1
}

lala<-merge(yes, unique(yo3[,c("studyid","intid","phenofreq","length")]), by=c("studyid","intid"))
lala2<-subset(lala, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
sig<-subset(lala2, p<0.07)
pos<-subset(lala2, p<0.07 & phenodiffbyyr>0)

#get temp change per interaction
Bgroups<-unique(mom$intid); b<-Bgroups; b<-as.character(b)
soop<-data.frame(array(0, c(length(b), 4)))
names(soop)[1] = "studyid"; names(soop)[2]<-"intid"; names(soop)[3]<-"tempbyyr"; names(soop)[4]<-"tempp"
for(i in 1:length(b)) { 
sub<-mom[mom$intid==b[i],]
soop[i,1]<-sub[1,c("studyid")]
soop[i,2]<-sub[1,c("intid")]
#m1<-with(sub, lm(envvalue~year))
m1<-with(sub, lm(z~year))
soop[i,3]<-summary(m1)$coefficients[2]
soop[i,4]<-summary(m1)$coefficients[8]
rowcount<-rowcount+1
}

all<-merge(no, yes, by=c("studyid","intid"))
all2<-merge(all, soop, by=c("studyid","intid"))


# get z per year
#mom[mom$intid=="202","intid"]<-"207"
#mom[mom$intid=="203","intid"]<-"208"
#mom[mom$intid=="204","intid"]<-"209"
yo3$intid<-as.factor(yo3$intid)
#solo$intid<-as.factor(solo$intid)
solo<-merge(yo3, mom[,c("studyid","year","intid","envvalue","z")], by=c("studyid","year","intid"))
solo2<-unique(solo[,c("studyid","year","intid","interaction","phenodiff","phenodiff_base","length","phenofreq","envvalue","z")])
solo2$phenofreq<-as.factor(solo2$phenofreq)
m1<-lme(sqrt(abs(phenodiff))~phenofreq+abs(z), random=~1|studyid/intid, method="ML", data=subset(solo2, length>6 & studyid!="EMW001" & studyid!="HMK035" & studyid!="175")); summary(m1)
m3<-update(m1, ~.-abs(z)); anova(m1,m3)
#m1<-lme(sqrt(abs(phenodiff_base))~phenofreq+abs(z), random=~1|studyid/intid, data=subset(solo2, length>6 & studyid!="EMW001" & studyid!="HMK035" & studyid!="175")); summary(m1)
#just those with C
m1<-lme(sqrt(abs(phenodiff_base))~phenofreq+abs(envvalue), random=~1|studyid/intid, weights=varPower(form=~envvalue), data=subset(solo2, length>6 & studyid!="EMW001" & studyid!="HMK035" & studyid!="175" & studyid!="HMK028")); summary(m1)

#overall z change
solo2$z<-as.numeric(solo2$z)
m1<-lme(z~year, random=~1|studyid/intid,data=subset(solo2, length>6 & studyid!="EMW001" & studyid!="HMK035" & studyid!="175")); summary(m1)


ggplot(data=all2, aes(abs(tempbyyr), abs(phenodiffbyyr))) +  geom_point(size=3) +
geom_smooth(method="lm", se=FALSE)+theme_bw()

#tempchange
m1<-lme(sqrt(abs(phenodiffbyyr))~abs(tempbyyr), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1) #& studyid!="175"
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)
sig<-subset(all2, p<0.07)

#only sig
m1<-lme(sqrt(abs(phenodiffbyyr))~abs(tempbyyr), random=~1|studyid, data=subset(all2, p<0.07 & studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1) #& studyid!="175"
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)
sig<-subset(all2, p<0.07)

#total difference in temp sensitivity (magnitude + direction)
m1<-lme(sqrt(abs(phenodiffbyyr))~abs(doybytemp), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp)); anova(m1,m3)

m1<-lme(sqrt(abs(phenodiffbyyr))~abs(rsquarediff), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)

#can temp sens predict spp' pheno shift?
m1<-lme(sqrt(abs(doybyyr))~abs(doybytemp), random=~1|studyid, data=subset(son,studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp)); anova(m1,m3)

#can temp change predict spp' pheno shift?
m1<-lme(sqrt(abs(doybyyr))~abs(tempbyyr), random=~1|studyid, data=subset(son, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)

m1<-lme(rsquarediff~1, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-1); anova(m1,m3)

m1<-lme(proprsquare~1, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-1); anova(m1,m3)

m1<-lme(sqrt(abs(phenodiffbyyr))~rsquarediff, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-1); anova(m1,m3)


# temp sens vary by taxa?
taxer <- read.csv("taxa.csv", header=TRUE)
taxer$latbi <- paste(taxer$genus, taxer$species)
taxer<-taxer[,c("studyid","taxa","latbi")]
names(taxer)[3]<-"species"
so4<-merge(mom, taxer, by=c("studyid","species"))
m2<-lme(phenovalue~z+spp+taxa, random=~1|studyid/intid/species, method="REML", weights=varIdent(form=~1|spp), data=subset(so4, studyid!="HMK035" & studyid!="EMW001"), na.action=na.omit); summary(m2)
m3<-update(m2, ~.-taxa); anova(m3,m2)
m3<-lme(phenovalue~z+taxa+spp, random=~1|studyid/intid/species, method="REML", weights=varIdent(form=~1|taxa), data=subset(so4, studyid!="HMK035" & studyid!="EMW001"), na.action=na.omit); summary(m3)
m4<-lme(phenovalue~z+spp+taxa, random=~1|studyid/intid/species, method="REML", weights=varComb(varIdent(form=~1|taxa),varIdent(form=~1|spp)), data=subset(so4, studyid!="HMK035" & studyid!="EMW001"), na.action=na.omit);

#phenodiff changes for those int with temp sens?
mom2<-unique(mom[,c("studyid","intid")])
mom3<-merge(mom2,yo3, by=c("studyid","intid"))
m1<-lme(phenovalue~z+spp, random=~1|studyid/intid/species, method="REML", data=subset(mom, studyid!="HMK035" & studyid!="EMW001"), na.action=na.omit); summary(m1)
m4<-lme(abs(phenodiff)~phenofreq+year, random=~1|studyid/intid, weights=varPower(form=~year), method="ML", data=subset(mom3, length>4 & studyid!="HMK035" & studyid!="EMW001"), na.action=na.omit);

#### Do interannual deviations in synchro correlate with interannual dev in envionmental condtions? ####
mom<-read.csv("int_phenoenv.csv", header=TRUE, na.strings="NA", as.is=TRUE)
m1<-lme(abs(phenodiff)~z, random=~1|studyid/intid, method="ML", data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
m2<-lme(abs(phenodiff)~z+I(z^2), random=~1|studyid/intid, method="ML", data=mom, na.action=na.omit); summary(m2)
AIC(m1,m2)

m1<-lme(z_pheno~z, random=~1|studyid/intid, method="ML", data=mom, na.action=na.omit); summary(m1)

anova(m1,m2)

ggplot(data=mom, aes(z, abs(phenodiff),colour=factor(intid)))+
	geom_point(size=1)+geom_smooth(method="lm", formula=y~poly(x,2), se=FALSE, aes(fill=factor(intid)), size=1)+theme_bw()+theme(legend.position="none")
	
ggplot(data=mom, aes(z, abs(phenodiff)))+
	geom_point(size=1)+geom_smooth(method="lm", formula=y~poly(x,2), se=FALSE)+theme_bw()+theme(legend.position="none")
	
ggplot(data=mom, aes(z, z_pheno, colour=factor(intid)))+
	geom_point(size=1)+geom_smooth(method="lm", se=FALSE, size=1)+theme_bw()+theme(legend.position="none")


## OLD ##################################
#differences in study length
mat2<-subset(yo2, length>=4)
mat3<-subset(yo2, length>=10)
par(mfrow=c(2,2)); with(yo2, plot(phenodiff~year)); with(mat2, plot(phenodiff~year)); with(mat3, plot(phenodiff~year))


#differences across interaction type?
m1<-lme(abs(phenodiff)~minyear+interaction+year, random=~1|studyid, method="ML", data=subset(yo3, studyid!=c("HMK006","HMK007","HMK010","HMK031","HMK032","HMK036")), na.action=na.omit); summary(m1)
m2<-update(m1, ~.-interaction); anova(m2,m1)

yess<-subset(yep, s))


model<-lme(abs(phenodiff)~year+length, random=~1|studyid/intid, method="ML", data=subset(yo2, length>=4 & studyid!="SET005"), na.action=na.omit); summary(model)
model2<-update(model, ~.-length); anova(model2, model)

model<-lme(abs(phenodiff)~year+min_year, random=~1|studyid/intid, method="ML", data=subset(yeo, length>=4 & studyid!="SET005"), na.action=na.omit); summary(model)
model2<-update(model, ~.-min_year); anova(model2, model)

model<-lme(abs(phenodiff)~year+max_year, random=~1|studyid/intid, method="ML", data=subset(yeo, length>=4 & studyid!="SET005"), na.action=na.omit); summary(model)
model2<-update(model, ~.-max_year); anova(model2, model)

yeo$interaction<-as.factor(yeo$interaction)
model<-lme(abs(phenodiff)~interaction+year, random=~1|studyid/intid, method="ML", data=subset(yeo, length>=4 & studyid!="SET005"), na.action=na.omit); summary(model)
model2<-update(model, ~.-interaction); anova(model2, model)

model<-lme(abs(phenodiff)~year, random=~1|studyid/intid, method="ML", data=subset(yo2, length>=4 & studyid!="SET005"), na.action=na.omit); summary(model)
model2<-update(model, ~.-year); anova(model2, model)

anova(model)
Anova(model)

endstart<-read.csv("endstart.csv", header=TRUE, na.strings="NA", as.is=TRUE)
yeo<-merge(yo2, endstart[,c("studyid","intid","min_year","max_year")])

# how are individual species responding over time
all<-read.csv("spp_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
all<-all[order(all$studyid),]
total2<-read.csv("int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
yo<-aggregate(total2["year"], by=total2[c("studyid","intid")], FUN=length)
names(yo)[names(yo)=="year"]<-"length"
all2<-merge(all, yo, by=c("studyid","intid"))

model<-lme(phenovalue~year, random=~1|studyid/intid/species, data=all)
model<-lme(phenovalue~year, random=~1|studyid/intid/species, data=subset(all2, )

# what are slope diff for interacting spp
total2<-read.csv("int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
ubtloc <- unique(total2[,c("studyid","intid","neg_slope","pos_slope")])
ubtloc<-na.omit(ubtloc)
ubtloc$slope_diff<-with(ubtloc, neg_slope-pos_slope)
yo<-aggregate(total2["year"], by=total2[c("studyid","intid")], FUN=length)
names(yo)[names(yo)=="year"]<-"length"
ubt2<-merge(ubtloc, yo, by=c("studyid","intid"))
ubt3<-subset(ubt2, length>=10)

