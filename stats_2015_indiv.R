library(nlme)
library(lme4)
library(car)
library(MASS)
library(MuMIn)
library(mgcv)

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

#number of years by INTERACTION (not study!)
total4<-na.omit(total3)
yo<-aggregate(total4["year"], by=total4[c("studyid","intid")], FUN=length)
names(yo)[names(yo)=="year"]<-"length"
yo2<-merge(total4,yo, by=c("studyid","intid"))

#first year of study WITH INTERACTION
yo<-aggregate(total4["year"], by=total4[c("studyid","intid")], FUN=min)
names(yo)[names(yo)=="year"]<-"minyear"
yo3<-merge(yo2,yo, by=c("studyid","intid"))

yo3$phenofreq<-as.factor(yo3$phenofreq)
yo3$taxatype<-as.factor(yo3$taxatype)
yo3$interaction[yo3$interaction=="poll"] <- "pollination"
yo3$interaction<-as.factor(yo3$interaction)

# calculate indiv slopes
yo3<-yo3[order(yo3$intid,yo3$year),]; yo3<-na.omit(yo3)
mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
best<-data.frame(array(0, c(length(unique(mat$intid)), 4)))
names(best)[1]<-"intid"; names(best)[2]<-"coef"; names(best)[3]<-"se"; names(best)[4]<-"p"
Bgroups<-unique(mat$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-mat[mat$intid==b[i],]
best[i,1]<-b[i]
m1<-lm(phenodiff~year, spp)
best[i,1]<-b[i]
best[i,2]<-summary(m1)$coefficients[2]
best[i,3]<-summary(m1)$coefficients[4]
best[i,4]<-summary(m1)$coefficients[8]
}
inc<-subset(best, intid=="19" | intid=="144" | intid=="145" | intid=="146" | intid=="155" | intid=="165" | intid=="167"  | intid=="168"  | intid=="169"  | intid=="170"  | intid=="175"  | intid=="179" | intid=="180" | intid=="185" | intid=="193" | intid=="197" | intid=="198"  | intid=="199" | intid=="200" | intid=="207" | intid=="219" | intid=="220" | intid=="231" | intid=="232" | intid=="233" | intid=="234") # ie. towards 0
inc$coef<-abs(inc$coef)
dec<-subset(best, intid=="1" | intid=="171" | intid=="177"| intid=="178"| intid=="181" | intid=="183"| intid=="189"| intid=="190"| intid=="191"| intid=="192"| intid=="194"| intid=="195"| intid=="196"| intid=="201"| intid=="208"| intid=="209"| intid=="214"| intid=="215"| intid=="216"| intid=="217"| intid=="218"| intid=="235") #away from 0
dec$coef<-abs(dec$coef); dec$coef<--(dec$coef) # away from 0
tot<-rbind(inc, dec)

yolo<-merge(unique(mat[,c("studyid","intid","interaction","taxatype","phenofreq","intxn","length","minyear")]), best, by=c("intid"))

m1<-lme(coef~1, random=~1|studyid, method="ML", yolo, na.action=na.omit); summary(m1)
m2<-lme(coef~1, random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), yolo, na.action=na.omit);
AIC(m1,m2)
m4<-update(m2, ~.-1); anova(m2,m4)

#without non-linear intxns
m1<-lme(coef~1, random=~1|studyid, method="ML", subset(yolo, intid!="177" & intid!="197"), na.action=na.omit); summary(m1)
m2<-lme(coef~1, random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), subset(yolo, intid!="177" & intid!="197"), na.action=na.omit); summary(m2)
AIC(m1,m2)
m4<-update(m2, ~.-1); anova(m2,m4)

m2<-lme(coef~phenofreq+minyear+length, random=~1|studyid, weights=varFixed(~1/log(se)), method="ML", yolo, na.action=na.omit); summary(m2)
plot(m2)
m4<-update(m2, ~.-minyear); anova(m4,m2)
m1<-lme(coef~minyear, random=~1|studyid, method="ML", yolo, na.action=na.omit); summary(m1)
m2<-lme(coef~minyear, random=~1|studyid, method="REML", weights=varExp(form=~minyear), yolo, na.action=na.omit); 
m3<-lme(coef~minyear, random=~1|studyid, method="REML", weights=varExp(form=~minyear, fixed=0.5), yolo, na.action=na.omit); 
m4<-lme(coef~minyear, random=~1|studyid, method="REML", weights=varPower(form=~minyear), yolo, na.action=na.omit); 
m5<-lme(coef~minyear, random=~1|studyid, method="REML", weights=varPower(form=~minyear, fixed=0.5), yolo, na.action=na.omit); 
m6<-lme(coef~minyear, random=~1|studyid, weights=varFixed(~minyear), method="REML", yolo, na.action=na.omit);
AIC(m1,m4,m5, m6)

m1<-lme(coef~minyear, random=~1|studyid, method="ML", yolo, na.action=na.omit);
m2<-lme(coef~minyear, random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), yolo, na.action=na.omit);
m4<-update(m1, ~.-1); anova(m1,m4)

#taxatype
m1<-lme(coef~taxatype, random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), yep, na.action=na.omit); summary(m1)
m4<-lme(abs(coef)~taxatype*terrestrial, random=~1|studyid, method="ML", yep, weights=varFixed(~1/log(se)), na.action=na.omit); summary(m4)
m4<-lme(abs(coef)~terrestrial, random=~1|studyid, method="ML", subset(yep, taxatype=="spp"), weights=varFixed(~1/log(se)), na.action=na.omit); summary(m4) 
m4<-lme(abs(coef)~terrestrial, random=~1|studyid, method="ML", subset(yep, taxatype=="genus"), weights=varFixed(~1/log(se)), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-terrestrial); anova(m4,m2) 

m4<-lme(abs(coef)~length+taxatype*terrestrial, random=~1|studyid, method="ML", yep, weights=varFixed(~1/log(se)), na.action=na.omit); summary(m4)
m2<-update(m4, ~.-taxatype:terrestrial); anova(m4,m2) 

#study characteristics
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(yolo, stud[,c("studyid","biome2","terrestrial", "terrestrial2", "big_interaction")], by="studyid")
yep$biome2<-as.factor(yep$biome2)
yep$terrestrial<-as.factor(yep$terrestrial)
yep$terrestrial2<-as.factor(yep$terrestrial2)

m4<-lme(abs(coef)~terrestrial, random=~1|studyid, method="ML", yep, weights=varFixed(~1/log(se)), na.action=na.omit); summary(m4) 
m5<-lme(abs(coef)~terrestrial, random=~1|studyid, method="ML", yep, weights=varIdent(~terrestrial), na.action=na.omit); 
AIC(m4,m5)

m4<-lme(abs(coef)~terrestrial, random=~1|studyid, method="ML", yep, weights=varFixed(~1/log(se)), na.action=na.omit); summary(m4) 
m2<-update(m4, ~.-terrestrial); anova(m4,m2)

#without non-linear
m4<-lme(abs(coef)~terrestrial, random=~1|studyid, method="ML", subset(yep, intid!="177" & intid!="197"), weights=varFixed(~1/log(se)), na.action=na.omit); summary(m4) 
m2<-update(m4, ~.-terrestrial); anova(m4,m2)

pos<-subset(yep, coef>=0 & p<=0.05); pos$count<-1; pos$direc<-"posi"; with(pos, table(count, terrestrial))
neg<-subset(yep, coef<=0 & p<=0.05); neg$count<-1; neg$direc<-"nega"; with(neg, table(count, terrestrial)); 
ns<-subset(yep, p>0.05);
sig<-rbind(pos,neg); sig$direc<-as.factor(sig$direc); sig$terrestrial<-as.factor(sig$terrestrial)

m4<-glmer(direc~terrestrial+(1|studyid), family=binomial, sig, na.action=na.omit); summary(m4)
m2<-update(m4, ~.-terrestrial); anova(m4,m2)
library(blme)
m4<-bglmer(direc~terrestrial+(1|studyid), family=binomial, sig, cov.prior = NULL, fixef.prior = normal);

http://stackoverflow.com/questions/25985970/generalised-linear-mixed-model-error-binary-response
Ben Bolker suggests using bglmer from the blme package to impose zero-mean Normal priors on the fixed effects (a 4 \( \times \) 4 diagonal matrix with diagonal elements equal to 9, for variances of 9 or standard deviations of 3), or add a B element to the priors list for MCMCglmm to specify the same priors:



pos<-c(6,5)
neg<-c(2,6)
chisq.test(pos,neg)

###change in absolute relative timing WHILE ACCOUNTING FOR FIRST YEAR DIFFERENCE ###
mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
best<-data.frame(array(0, c(length(unique(mat$intid)), 3)))
names(best)[1]<-"intid"; names(best)[2]<-"studyid"; names(best)[3]<-"base"
Bgroups<-unique(mat$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-mat[mat$intid==b[i],]
spp2<-spp[1:3,]
best[i,1]<-b[i]
best[i,2]<-spp[1,c("studyid")]
#best[i,3]<-spp[1,c("phenodiff")] # first year
best[i,3]<-mean(spp2$phenodiff) #avg of 3 first years
}
best2<-merge(best, yolo, by=c("studyid","intid"))
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(best2, stud[,c("studyid","biome2","terrestrial", "terrestrial2", "big_interaction")], by="studyid")
yep$biome2<-as.factor(yep$biome2)
yep$terrestrial<-as.factor(yep$terrestrial)
yep$terrestrial2<-as.factor(yep$terrestrial2)


m1<-lme(abs(coef)~abs(base), random=~1|studyid, method="ML", best2, na.action=na.omit); summary(m1)
m2<-lme(abs(coef)~abs(base), random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), best2, na.action=na.omit); summary(m2)
m2<-lme(log(abs(coef))~abs(base), random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), best2, na.action=na.omit); summary(m2)
AIC(m1,m2)
m4<-update(m2, ~.-abs(base)); anova(m2,m4)

m2<-lme(log(abs(coef))~abs(base)*terrestrial, random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), yep, na.action=na.omit); summary(m2)
m4<-update(m2, ~.-abs(base):terrestrial); anova(m2,m4)

m2<-lme(log(abs(coef))~abs(base), random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), subset(yep, terrestrial=="terrestrial"), na.action=na.omit); summary(m2)
m4<-update(m2, ~.-abs(base)); anova(m2,m4)

## FROM BASELINE ##
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

# calculate indiv slopes
yo3<-yo3[order(yo3$intid,yo3$year),]; yo3<-na.omit(yo3)
mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
best<-data.frame(array(0, c(length(unique(mat$intid)), 4)))
names(best)[1]<-"intid"; names(best)[2]<-"coef"; names(best)[3]<-"se"; names(best)[4]<-"p"
Bgroups<-unique(mat$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-mat[mat$intid==b[i],]
best[i,1]<-b[i]
m1<-lm(abs(phenodiff_base)~year, spp)
best[i,1]<-b[i]
best[i,2]<-summary(m1)$coefficients[2]
best[i,3]<-summary(m1)$coefficients[4]
best[i,4]<-summary(m1)$coefficients[8]
}

yolo<-merge(unique(mat[,c("studyid","intid","interaction","taxatype","phenofreq","intxn","length","minyear")]), best, by=c("intid"))

m1<-lme(coef~1, random=~1|studyid, method="ML", yolo, na.action=na.omit); summary(m1)
m2<-lme(coef~1, random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), yolo, na.action=na.omit);summary(m2)
AIC(m1,m2)

m1<-lme(log(coef+1)~1, random=~1|studyid, method="ML", yolo, na.action=na.omit); summary(m1)
m2<-lme(log(coef+1)~1, random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), yolo, na.action=na.omit);summary(m2)
AIC(m1,m2)
m4<-update(m2, ~.-1); anova(m4,m2)

#without non-linear intxns
m1<-lme(log(coef+1)~1, random=~1|studyid, method="ML", subset(yolo, intid!="165" & intid!="177" & intid!="185" & intid!="197" & intid!="198" & intid!="199"), na.action=na.omit); summary(m1)
m2<-lme(log(coef+1)~1, random=~1|studyid, method="ML", weights=varFixed(~1/log(se)), subset(yolo, intid!="165" & intid!="177" & intid!="185" & intid!="197" & intid!="198" & intid!="199"), na.action=na.omit);summary(m2)
AIC(m1,m2)
m4<-update(m2, ~.-1); anova(m4,m2)

m2<-lme(coef~phenofreq+minyear+length, random=~1|studyid, weights=varFixed(~1/log(se)), method="ML", yolo, na.action=na.omit); summary(m2)

#study characteristics
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(yolo, stud[,c("studyid","biome2","terrestrial", "terrestrial2", "big_interaction")], by="studyid")
yep$biome2<-as.factor(yep$biome2)
yep$terrestrial<-as.factor(yep$terrestrial)
yep$terrestrial2<-as.factor(yep$terrestrial2)

m4<-lme(log(coef+1)~terrestrial, random=~1|studyid, method="ML", yep, weights=varFixed(~1/log(se)), na.action=na.omit); summary(m4) # mismatch!="1
m5<-lme(abs(coef)~terrestrial, random=~1|studyid, method="ML", yep, weights=varIdent(~terrestrial), na.action=na.omit); 



#######################
#      Non-linear     #
######################
# calculate indiv slopes
yo3<-yo3[order(yo3$intid,yo3$year),]; yo3<-na.omit(yo3)
mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
best<-data.frame(array(0, c(length(unique(mat$intid)), 5)))
names(best)[1]<-"intid"; names(best)[2]<-"coef"; names(best)[3]<-"se"; names(best)[4]<-"p"
Bgroups<-unique(mat$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-mat[mat$intid==b[i],]
if (nrow(spp)<10){
best[i,1]<-b[i]
m1<-lm(phenodiff~year, spp)
best[i,1]<-b[i]
best[i,2]<-summary(m1)$coefficients[2]
best[i,3]<-summary(m1)$coefficients[4]
best[i,4]<-summary(m1)$coefficients[8]
best[i,5]<-NA
}
if (nrow(spp)>=10){
best[i,1]<-b[i]
m1<-lm(phenodiff~year, spp)
#m2<-gam(phenodiff ~ s(year, k=-1), na.action=na.omit, data=spp)
m3<-lm(phenodiff ~ poly(year,2), na.action=na.omit, data=spp) #doesn't need to be >10 pts
best[i,1]<-b[i]
best[i,2]<-summary(m1)$coefficients[2]
best[i,3]<-summary(m1)$coefficients[4]
best[i,4]<-summary(m1)$coefficients[8]
best[i,5]<-summary(m3)$coefficients[12]
#best[i,5]<-AIC(m1)-AIC(m3)
}
}
nl<-subset(best, X5<=0.06); 
nl<-subset(best, X5>0); nl2<-subset(best, X5>=2);nl2
with(subset(yo3, intid=="165"), plot(phenodiff~year));
with(subset(yo3, intid=="165"), abline(gam(phenodiff~s(year, k=-1))))
ggplot(data=subset(yo3, intid=="197"), aes(year, phenodiff))+geom_point()+geom_smooth(method="lm", formula=y ~ poly(x, 2), se=FALSE,size=0.5, colour="black") + theme_bw()

ggplot(data=subset(yo3, intid=="177"), aes(year, phenodiff))+geom_point()+geom_smooth(method="gam", formula=y ~ s(x), se=FALSE,size=0.5, colour="black") + theme_bw()
	
#breakpt studies
bk<-subset(total3, intid=="170" | intid=="171" | intid=="177" | intid=="178" | intid=="179" | intid=="180" | intid=="189" | intid=="195" | intid=="196")
mat<-bk
best<-data.frame(array(0, c(length(unique(mat$intid)), 5)))
names(best)[1]<-"intid"; names(best)[2]<-"coef"; names(best)[3]<-"se"; names(best)[4]<-"p"
Bgroups<-unique(mat$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-mat[mat$intid==b[i],]
if (nrow(spp)<10){
best[i,1]<-b[i]
m1<-lm(phenodiff~year, spp)
best[i,1]<-b[i]
best[i,2]<-summary(m1)$coefficients[2]
best[i,3]<-summary(m1)$coefficients[4]
best[i,4]<-summary(m1)$coefficients[8]
best[i,5]<-NA
}
if (nrow(spp)>=10){
best[i,1]<-b[i]
m1<-lm(phenodiff~year, spp)
m2<-gam(phenodiff ~ s(year, k=-1), na.action=na.omit, data=spp)
m3<-lm(phenodiff ~ poly(year,2), na.action=na.omit, data=spp)
best[i,1]<-b[i]
best[i,2]<-summary(m1)$coefficients[2]
best[i,3]<-summary(m1)$coefficients[4]
best[i,4]<-summary(m1)$coefficients[8]
best[i,5]<-AIC(m1)-AIC(m3)
}
}
nl<-subset(best, X5>0); nl2<-subset(best, X5>=2); nl2

#baseline
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

yo3<-yo3[order(yo3$intid,yo3$year),]; yo3<-na.omit(yo3)
mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
best<-data.frame(array(0, c(length(unique(mat$intid)), 5)))
names(best)[1]<-"intid"; names(best)[2]<-"coef"; names(best)[3]<-"se"; names(best)[4]<-"p"
Bgroups<-unique(mat$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-mat[mat$intid==b[i],]
best[i,1]<-b[i]
m1<-lm(phenodiff_base~year, spp)
#m2<-gam(abs(phenodiff_base) ~ s(year, k=-1), na.action=na.omit, data=spp)
m3<-lm(phenodiff_base ~ poly(year,2), na.action=na.omit, data=spp)
best[i,1]<-b[i]
best[i,2]<-summary(m1)$coefficients[2]
best[i,3]<-summary(m1)$coefficients[4]
best[i,4]<-summary(m1)$coefficients[8]
best[i,5]<-AIC(m1)-AIC(m3)
}
nl<-subset(best, X5>70); nl2<-subset(best, X5>=2);nl2

ggplot(data=subset(yo3, intid=="165"), aes(year, abs(phenodiff_base)))+geom_point()+geom_smooth(method="lm", formula=y~poly(x,2),se=FALSE,size=0.5, colour="black") + theme_bw()


### ####
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
sure<-unique(so4[,c("studyid","intid","int_type","spp","year","species","phenovalue","phenofreq","short_site", "terrestrial")])
write.csv(sure, "/users/kharouba/google drive/UBC/synchrony project/analysis/stan/indiv_sppdata2.csv")


best<-data.frame(array(0, c(length(unique(sure$species)), 4)))
names(best)[1]<-"species"; names(best)[2]<-"coef"; names(best)[3]<-"se"; names(best)[4]<-"p"
Bgroups<-unique(sure$species); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-sure[sure$species==b[i],]
best[i,1]<-b[i]
m1<-lm(phenodiff~year, spp)
best[i,1]<-b[i]
best[i,2]<-summary(m1)$coefficients[2]
best[i,3]<-summary(m1)$coefficients[4]
best[i,4]<-summary(m1)$coefficients[8]
}

