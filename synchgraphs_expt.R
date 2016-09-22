rm(list=ls())
options(stringsAsFactors=FALSE)
library(ggplot2)

source("/Volumes/MUSIC/UBC/synchrony project/analysis/synchexpt_calc.R")
setwd("/Volumes/Music/UBC/synchrony project/")


#Relationship between warming and change in doy

mu<-subset(manip, treatchangeC!="NA" & treatchangeday!="NA")
mu$unique<-with(mu, paste(latbi,"_",latbi2))
ggplot(data=mu, aes(treatchangeC, abs(treatchangeday),colour=factor(unique)))+ geom_point(size=2.5)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(unique)))+theme_bw()

mu<-subset(manoh, tempchange!="NA" & phenochange!="NA")
ggplot(data=subset(mu, phenochange>-40), aes(tempchange, phenochange))+ geom_point(size=2.5, aes(colour=factor(studyid)))+geom_smooth(method="lm", se=FALSE)+xlab("Temperature change (C)")+theme_bw()

############################################
#Relationship between doychange and fitness#
############################################
mat<-manoh[,1:7];mat<-na.omit(mat)
mat<-subset(mat, site!="Mustaseika, Finland")
#only pos interactions
ggplot(data=subset(mat, studyid!="JMA002"), aes(phenochange, fitchange))+ geom_point(size=2.5, aes(colour=factor(uniqueintxn)))+geom_smooth(method="lm",se=FALSE)+geom_line(y=0, linetype="dashed")+geom_vline(linetype="dashed")+theme_bw()+theme(legend.position="none")

mats<-manip[,c("studyid","site","treatchangeday","fitnessvalue","latbi","latbi2")]
mats$unique<-with(mats, paste(latbi,"_",latbi2,"_",site))
mats<-subset(mats, treatchangeday!="NA" & fitnessvalue!="NA")
ggplot(data=mats, aes(treatchangeday, fitnessvalue, colour=factor(unique)))+ geom_point(size=2.5, aes(colour=factor(unique)))+geom_smooth(method="lm",se=FALSE, aes(fill=factor(unique)))+theme_bw()

#reorganize data
names(man.fit)[6]<-"treatment"   
mo<-merge(manip[,c("studyid","site","year","treatment","uniqueintxn","treatchangeday")], man.fit, by=c("studyid","site","year","uniqueintxn","treatment"), all.x=TRUE)
mo<-mo[,1:7]
mo$fitchange[mo$treatment=="0"] <- 0
mo$fitchange[mo$treatment=="2" & mo$studyid=="HMK002"] <- 0.00
mo$fitchange[mo$treatment=="1" & mo$studyid=="HMK003"] <- 0.00
mo$fitchange[mo$treatment=="1" & mo$studyid=="HMK005"] <- 0.00
mo$fitchange[mo$treatment=="1" & mo$studyid=="HMK044"] <- 0.00
mo<-na.omit(mo)
mo$fitchange<-as.numeric(mo$fitchange)
mo<-subset(mo, studyid!="HMK001" & studyid!="HMK014" & studyid!="HMK017" & studyid!="HMK022" & studyid!="HMK041")
mo$rowid<-1:nrow(mo)
mo<-subset(mo, rowid!="1" & rowid!="2" & rowid!="37")
mo$treatchangeday[mo$treatment=="0" & mo$studyid=="HMK045"] <- 0.00
mo$treatchangeday[mo$treatment=="1" & mo$studyid=="HMK045"] <- 52.0 # treatment is delayed
mo$fitchange[mo$treatment=="1" & mo$uniqueintxn=="Cotesia bignellii Euphydryas aurinia parasitoid hatching parasitoid development time"]<-0.11707364 # development time works opp to fecundity in that less devt time is better so need to flip sign


mo$treatment[mo$treatment=="2" & mo$studyid=="HMK002"] <-"0"
mo$treatment[mo$treatment=="1" & mo$studyid=="HMK003"] <-"0"
mo$treatment[mo$treatment=="1" & mo$studyid=="HMK005"] <-"0"
mo$treatment[mo$treatment=="1" & mo$studyid=="HMK044"] <-"0"
mo$treatchangeday2<-0
Bgroups<-unique(mo$uniqueintxn); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-mo[mo$uniqueintxn==b[i],]
asdf<-rowcount+(nrow(spp)-1)
yes<-subset(spp,treatment=="0")
no<-yes[1, c("treatchangeday")]
mo[rowcount:asdf,c("treatchangeday2")]<- with(spp, treatchangeday-no)
rowcount<-rowcount+nrow(spp)
}

#only positive interactions
new<-subset(mo, studyid!="JMA002" & uniqueintxn!="Euphydryas aurinia Cotesia bignellii host availability host development time (ffinal instar to pupation)")
fm1<-lme(fitchange~treatchangeday, random=~1|uniqueintxn, data=new);
newdat<-expand.grid(new$treatchangeday); names(newdat)[1]<-"treatchangeday"
newdat$pred<-predict(fm1, newdat, level=0)

Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)
ggplot(data=newdat, aes(treatchangeday, y=pred))+ geom_point(data=new, aes(treatchangeday, fitchange, colour=factor(uniqueintxn)), size=2.5)+geom_smooth(data=new, method="lm",se=FALSE, aes(treatchangeday, fitchange,colour=factor(uniqueintxn)), size=1)+geom_vline(size=0.25,linetype="dashed")+geom_hline(size=0.25,linetype="dashed")+
geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.2, fill="blue")+
geom_line(size=2)+
theme_bw()+opts(legend.position="none",axis.title.x =theme_text(size=17), axis.text.x=theme_text(size=17), axis.text.y=theme_text(size=17), axis.title.y=theme_text(size=17, angle=90))+ylab("Performance change")+xlab("Relative timing")
#formula=y~poly(x,2)

mo$treatchangeday<-as.numeric(mo$treatchangeday)
m1<-lme(fitchange~treatchangeday, random=~1|uniqueintxn, data=subset(mo, studyid!="JMA002" & uniqueintxn!="Euphydryas aurinia Cotesia bignellii host availability host development time (ffinal instar to pupation)")); summary(m1) 

ggplot(data=subset(mo, studyid!="JMA002"), aes(abs(treatchangeday2), fitchange, colour=factor(uniqueintxn)))+ geom_point(size=2.5)+geom_smooth(method="lm",se=FALSE, aes(fill=factor(uniqueintxn)))+theme_bw()+theme(legend.position="none")+xlab("Days apart")

#
doPlotphenodiff <- function(sel_name) {
   subby <- mo[mo$uniqueintxn == sel_name,]
   ggobj <- ggplot(data=subby, aes(treatchangeday, fitchange)) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE)+theme_bw()+opts(legend.position="none")
   print(ggobj)
   ggsave(sprintf("graphs/fitbydoy_expt/%s.pdf", sel_name))
}
   
lapply(unique(mo$uniqueintxn), doPlotphenodiff)

#negative
mo$fitchange[mo$treatment=="1" & mo$uniqueintxn=="Euphydryas aurinia Cotesia bignellii host availability host development time (ffinal instar to pupation)"]<-0.14734104 # development time works opp to fecundity in that less devt time is better so need to flip sign
ggplot(data=subset(mo, studyid=="JMA002" | uniqueintxn=="Euphydryas aurinia Cotesia bignellii host availability host development time (ffinal instar to pupation)"), aes(treatchangeday, fitchange, colour=factor(uniqueintxn)))+ geom_point(size=2.5)+geom_smooth(method="lm",se=FALSE, aes(fill=factor(uniqueintxn)))+geom_vline(size=0.25,linetype="dashed")+geom_hline(size=0.25,linetype="dashed")+theme_bw()+theme(legend.position="none")+ylab("Performance change")

mat<-unique(manip[,c("studyid","year","site","genus","species","genus2","species2","treatment","treatchangeday","fitnesscomp")]); mat<-na.omit(mat)
mat$unique<-with(mat, paste(studyid,"_",year,"_",site,"_", genus,"_",species2))
ggplot(data=mat, aes(treatment, treatchangeday,colour=factor(unique)))+ geom_point(size=1)+geom_line(aes(group=unique))+theme_bw()+opts(legend.position="none")


dater<-unique(manip[,c("studyid","year","site","genus","species","genus2","species2","treatment","treatchangeday")]); dater<-na.omit(dater)
dater$unique<-with(dater, paste(studyid,"_",year,"_",site,"_",genus,"_",species2))
doPlotphenodiff <- function(sel_name) {
   subby <- dater[dater$unique == sel_name,]
   ggobj <- ggplot(data=subby, aes(treatment, treatchangeday,
       group=1)) +  geom_point(shape=1) +
       geom_line(group=1)+theme_bw()+opts(legend.position="none")
   print(ggobj)
   ggsave(sprintf("graphs/doybytreat/%s.pdf", sel_name))
}
   
lapply(unique(dater$unique), doPlotphenodiff)




dater<-unique(manip[,c("studyid","year","site","uni","treatchangeday","fitnessvalue")]); dater<-na.omit(dater)
doPlotphenodiff <- function(sel_name) {
   subby <- manoh[manoh$uniqueintxn == sel_name,]
   ggobj <- ggplot(data=subby, aes(phenochange, fitchange)) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE)+theme_bw()+opts(legend.position="none")
   print(ggobj)
   ggsave(sprintf("graphs/fitbyphenochange_expt/%s.pdf", sel_name))
}
   
lapply(unique(manoh$uniqueintxn), doPlotphenodiff)

#relationship between doy/c vs fitness
mo<-subset(manip, studyid=="JEH001")
a<-ggplot(data=mo, aes(treatchangeC, treatchangeday))+ geom_point(size=2.5)+geom_line()+theme_bw()+opts(legend.position="none")
(treatchangeday(treatment=1)-treatchangeday(treatment=0))/(treatchangeC(treatment=1)-treatchangeC(treatment=0))
mo$doyc<-c(0,0.5379612)
b<-ggplot(data=mo, aes(doyc,fitnessvalue))+ geom_point(size=2.5)+geom_line()+theme_bw()+opts(legend.position="none")