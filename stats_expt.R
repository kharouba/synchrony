library(nlme)
library(car)

rm(list=ls())
row.names=FALSE

setwd("/Volumes/Music/UBC/synchrony project/analysis")
source("/users/kharouba/google drive/UBC/synchrony project/analysis/synchexpt_calc.R") # need to change all wd!!
#mat2<-merge(mat, stud, by="studyid")
std <- function(x) sd(x)/sqrt(length(x))

mat<-manoh[,1:7];mat<-na.omit(mat)
mat$fitnesstype<-c(1,1,1,1,3,3,3,3,3,3,1,1,1,1,1,1,1,1,2,1,1,3,3,2)
mat$fitnesstype<-as.factor(mat$fitnesstype)
mat<-subset(mat, site!="Mustaseika, Finland" & site!="farm" & uniqueintxn!="Asphondylia aucubae Aucuba japonica timing of oviposition oviposition success" & uniqueintxn!="Choristoneura fumiferana Picea glauca start of larval development larval survival" & uniqueintxn!="Malacosoma  californicum pluviale NA NA larval emergence larval weight (log)") #excluded non-linear
mat<-subset(mat, site!="Mustaseika, Finland" & site!="farm" & uniqueintxn!="Malacosoma  californicum pluviale NA NA larval emergence larval weight (log)") 
Asphondylia aucubae Aucuba japonica timing of oviposition oviposition success

#whats average phenochange that studies manipulate?
mat<-subset(mat, studyid!="JEH001" & studyid!="HMK047")
mat<-subset(mat, studyid!="JEH001" & studyid!="HMK047" & uniqueintxn!="Choristoneura fumiferana Picea glauca start of larval development larval survival")
mean(abs(mat$phenochange)); std(abs(mat$phenochange))

# phenochange should be abs because negative values for some studies means synchrony improved, others worse
m1<-lme(fitchange~abs(phenochange), random=~1|uniqueintxn, data=subset(mat, fitchange<0.8)); summary(m1)
m1<-lme(fitchange~fitnesstype+abs(phenochange), random=~1|uniqueintxn, method="ML",data=subset(mat, fitchange<0.8)); summary(m1)
m2<-update(m1,~.-fitnesstype); anova(m1,m2)

mu<-subset(manoh, tempchange!="NA" & phenochange!="NA")
mu$unique<-with(mu, paste(uniqueintxn,"_",year))
m1<-lme(phenochange~tempchange, random=~1|studyid, data=subset(mu, phenochange>-40)); summary(m1)

#how much synchrony changes by deg C
mu<-subset(manip, treatchangeC!="NA" & treatchangeday!="NA")
mu$unique<-with(mu, paste(latbi,"_",latbi2))
mu2<-subset(mu,  rowid!="90" & rowid!="91" & site!="Haugastol, Norway" & rowid!="100" & rowid!="103" & genus!="Alnus" & site!="farm" & treatment!="0" & rowid!="86" & rowid!="228" & rowid!="279") # get rid of mult years for HK017- and years for HMK001 AND data for JEH001 AND Hmk047 & HMK014- all randomly chosen
m1<-lme(abs(treatchangeday)~treatchangeC, random=~1|studyid, data=mu2);summary(m1)





