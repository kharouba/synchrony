library(nlme)
library(ggplot2)
rm(list=ls())
row.names=FALSE

setwd("/users/kharouba/google drive/UBC/synchrony project/")
#mat2<-merge(mat, stud, by="studyid")

dater<-read.csv("analysis/obs_fit.csv", header=TRUE, na.strings="NA", as.is=TRUE)
dater<-subset(dater, rowid!="364")
dater$spp1<-with(dater, paste(genus,species))
dater$spp2<-with(dater, paste(genus2, species2))
dater$intxn<-with(dater, paste(spp1,spp2))

pheno<-read.csv("analysis/int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)

#chose 1 site for HMK023 based on linear relationhip (i.e. no mismatch) AND with greatest number of years
HMK023<-subset(dater, studyid=="HMK023")
#newHMK023<-subset(HMK023, site!="nopporo" & site!="tomakomai")
newHMK023<-subset(HMK023, site=="nopporo")
datero<-subset(dater, studyid!="HMK023")
datero2<-rbind(datero,newHMK023)
dater<-datero2

HMK026<-subset(dater, studyid=="HMK026")
newHMK026<-subset(HMK026, fitnesscomp=="calf production") # can't use survival because mortality
datero<-subset(dater, studyid!="HMK026")
datero2<-rbind(datero,newHMK026)
dater<-datero2

HMK030<-subset(dater, studyid=="HMK030")
newHMK030<-subset(HMK030, fitnesscomp=="calf production") # can't use survival because mortality
datero<-subset(dater, studyid!="HMK030")
datero2<-rbind(datero,newHMK030)
dater<-datero2

## USE THIS SERIES OR SERIES BELOW BUT NOT BOTH!!
### DON'T USE RANDOM FOR FINAL ANALYSIS ##
#randomly choose fitness component when multiple/intxn
HMK016<-subset(dater, studyid=="HMK016")
newHMK016<-subset(HMK016, fitnesstype=="survival")
datero<-subset(dater, studyid!="HMK016")
datero2<-rbind(datero,newHMK016)
dater<-datero2

HMK011<-subset(dater, studyid=="HMK011")
newHMK011<-subset(HMK011, fitnesscomp=="number of chicks hatched")
datero<-subset(dater, studyid!="HMK011")
datero2<-rbind(datero,newHMK011)
dater<-datero2

HMK029<-subset(dater, studyid=="HMK029")
newHMK029<-subset(HMK029, fitnesstype=="reproduction")
datero<-subset(dater, studyid!="HMK029")
datero2<-rbind(datero,newHMK029)
dater<-datero2

###
#Take estimate of performance most likely to represent fitness, population>survival>repro>growth
HMK016<-subset(dater, studyid=="HMK016")
newHMK016<-subset(HMK016, fitnesstype=="survival")
datero<-subset(dater, studyid!="HMK016")
datero2<-rbind(datero,newHMK016)
dater<-datero2
 
HMK011<-subset(dater, studyid=="HMK011")
newHMK011<-subset(HMK011, fitnesscomp=="fledged chicks")
datero<-subset(dater, studyid!="HMK011")
datero2<-rbind(datero,newHMK011)
dater<-datero2

HMK029<-subset(dater, studyid=="HMK029")
newHMK029<-subset(HMK029, fitnesstype=="demographic") #negative trend over time
datero<-subset(dater, studyid!="HMK029")
datero2<-rbind(datero,newHMK029)
dater<-datero2

#need to take avg across sites for HMK035 THEN join with phenodiff BECAUSE no phenodifferences were averaged across sites IN paper
HMK035<-subset(dater, studyid=="HMK035")
newHMK035<-with(HMK035, aggregate(fitnessvalue, by=list(year, spp1), FUN=mean, na.rm=T)) #
names(newHMK035)[1]<-"year"; names(newHMK035)[2]<-"spp1"; names(newHMK035)[3]<-"fitnessvalue"
newHMK035$studyid<-"HMK035"

yo<-merge(newHMK035, HMK035[,c(1:12,14:17)], by=c("studyid","year","spp1"))
yo2<-unique(yo[,c("studyid","year","spp1","fitnessvalue","fitnesscomp","fitnessunit","fitnesstype","spp2","intxn")])
yo2<-subset(yo2, spp1!="Acartia tonsa")
pho<-subset(pheno, studyid=="HMK035")
pho2<-merge(yo2, pho[,c("studyid","year","phenodiff")], by=c("studyid","year"))

dater2<-subset(dater, studyid!="HMK035")
dater3<-unique(dater2[,c("studyid","year","spp1","fitnessvalue","fitnesscomp","fitnessunit","fitnesstype","spp2","intxn","phenodiff")])
dater4<-rbind(dater3,pho2)

dater<-dater4
#dater<-na.omit(dater)
dater$phenodiff<-as.numeric(dater$phenodiff)

#first year of study
yo<-aggregate(dater["year"], by=dater[c("studyid","intxn")], FUN=min)
names(yo)[names(yo)=="year"]<-"minyear"
dater2<-merge(dater,yo, by=c("studyid","intxn"))
dater<-dater2

#standardize fitness based on MEAN temporal trend
newdata<-with(dater, aggregate(fitnessvalue, by=list(studyid,intxn, fitnesscomp), FUN=mean, na.rm=T))
names(newdata)[1]<-"studyid"; names(newdata)[2]<-"intxn"; names(newdata)[3]<-"fitnesscomp"; names(newdata)[4]<-"meanfitness"

newdata2<-with(dater, aggregate(fitnessvalue, by=list(studyid,intxn, fitnesscomp), FUN=sd, na.rm=T))
names(newdata2)[1]<-"studyid"; names(newdata2)[2]<-"intxn"; names(newdata2)[3]<-"fitnesscomp"; names(newdata2)[4]<-"sdfitness"

newdata3<-merge(dater, newdata, by=c("studyid","intxn","fitnesscomp"))
newdata4<-merge(newdata3, newdata2, by=c("studyid","intxn","fitnesscomp"))

newdata4$fitz<-with(newdata4, (fitnessvalue-meanfitness)/sdfitness)
newdata<-newdata4

#Don't use#
#standardize PHENODIFF based on MEAN temporal trend
yum<-with(newdata, aggregate(abs(phenodiff), by=list(studyid,intxn, fitnesscomp), FUN=mean, na.rm=T))
names(yum)[1]<-"studyid"; names(yum)[2]<-"intxn"; names(yum)[3]<-"fitnesscomp"; names(yum)[4]<-"meanphenodiff"

yum2<-with(newdata, aggregate(abs(phenodiff), by=list(studyid,intxn, fitnesscomp), FUN=sd, na.rm=T))
names(yum2)[1]<-"studyid"; names(yum2)[2]<-"intxn"; names(yum2)[3]<-"fitnesscomp"; names(yum2)[4]<-"sdphenodiff"

yum3<-merge(yum, newdata, by=c("studyid","intxn","fitnesscomp"))
newdata4<-merge(yum3, yum2, by=c("studyid","intxn","fitnesscomp"))

newdata4$phenodiffz<-with(newdata4, (abs(phenodiff)-meanphenodiff)/sdphenodiff)
newdata<-newdata4


##extra
#re-calculate phenodiff from min year WILL LOSE HMK030 BECAUSE NO FIRST YEAR!!
dater<-dater[order(dater$intxn,dater$year),];
best<-data.frame(array(0, c(nrow(dater), 3)))
names(best)[1]<-"intxn"; names(best)[2]<-"year"; names(best)[3]<-"phenodiff_base"
Bgroups<-unique(dater$intxn); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-dater[dater$intxn==b[i],]
asdf<-rowcount+(nrow(spp)-1)
best[rowcount:asdf,1]<-b[i]
best[rowcount:asdf,2]<-spp[,c("year")]
best[rowcount:asdf,3]<- with(spp, phenodiff-spp[1,c("phenodiff")])
rowcount<-rowcount+nrow(spp)
}
dater2<-merge(dater, best, by=c("intxn","year"))
dater<-dater2

#standardize fitness based on YEAR 1
dater$fitz_yr1<-1
Bgroups<-unique(dater$intxn); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-dater[dater$intxn==b[i],]
asdf<-rowcount+(nrow(spp)-1)
#dater[rowcount:asdf,c("fitz_yr1")]<- with(spp, fitnessvalue-spp$fitnessvalue[spp$phenodiff_base==0])
dater[rowcount:asdf,c("fitz_yr1")]<- with(spp, (fitnessvalue-spp$fitnessvalue[spp$phenodiff_base==0])/spp$fitnessvalue[spp$phenodiff_base==0])
rowcount<-rowcount+nrow(spp)
}

########################
#### STATS #############
########################
## only positive interactions
#phenodiff has to be absolute because of differences across studies in how they measured synchrony/mismatch (esp. if negative mismatch), now will mean the same thing.
#SET005 removed because synchrony change has opposite meaning,
#HMK031 removed because unclear whether copepod interacts with diatom
#HMK026 and HMK030 have the same species BUT no years with data; choose Hmk030 because more data
#HMk039 removed becasue negative interaction
#HMK023 (spp1!="Corydalis ambigua") kept because positive role
mat<-newdata[,c("studyid","fitz","fitnesstype","phenodiff","intxn","spp1")]; mat<-na.omit(mat)
mat$fitnesstype<-as.factor(mat$fitnesstype)
m1<-lme(fitz~fitnesstype+abs(phenodiff), random=~1|intxn, method="ML", data=subset(mat, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")); summary(m1)
m2<-update(m1, ~.-fitnesstype); anova(m1,m2)
m1<-lme(fitz~abs(phenodiff), random=~1|studyid/intxn, method="ML", data=subset(mat, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")); summary(m1)
m2<-update(m1, ~.-abs(phenodiff)); anova(m1,m2)

m1<-lme(fitz~phenodiffz, random=~1|intxn, method="ML", data=subset(mat, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")); summary(m1)
m2<-update(m1, ~.-abs(phenodiff)); anova(m1,m2)

#sig slopes
mat2<-subset(mat, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")
Bgroups<-unique(mat2$intxn); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0, c(length(b), 4)))
names(new)[1] = "intxn"; names(new)[2] = "studyid"; names(new)[3]<-"coef"; names(new)[4] = "intp";                                                                                                                                                                                                                    
for(i in 1:length(b)) {
spp<-mat2[mat2$intxn==b[i],]; spp<-na.omit(spp)
m1<-lm(fitz~abs(phenodiff), data=spp, na.action=na.omit);
 new[i,1]<-spp$intxn[1]
new[i,2]<-spp$studyid[1]
new[i,3]<-summary(m1)$coefficients[2]
new[i,4]<-summary(m1)$coefficients[8]
}

#did these interactions change in synchrony?
mat<-newdata[,c("studyid","fitz","fitnesstype","phenodiff","intxn","spp1", "year")]; mat<-na.omit(mat)
mat$fitnesstype<-as.factor(mat$fitnesstype)
m1<-lme(sqrt(abs(phenodiff))~year, random=~1|intxn, method="ML", data=subset(mat, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")); summary(m1)
m2<-update(m1, ~.-year); anova(m1,m2)

# have these interactions had a decrease in performance OVER TIME?
m1<-lme(fitz~year, random=~1|intxn, method="ML", data=subset(mat, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")); summary(m1)
m2<-update(m1, ~.-year); anova(m1,m2)

#get performance~mismatch per interaction
matdata<-subset(newdata, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")
mat<-matdata[,c("studyid","fitz","fitnesstype","phenodiff","intxn","spp1", "year")]; mat<-na.omit(mat)

mat$fitnesstype<-as.factor(mat$fitnesstype)
Bgroups<-unique(mat$intxn); b<-Bgroups; b<-as.character(b)
soop<-data.frame(array(0, c(length(b), 4)))
names(soop)[1] = "studyid"; names(soop)[2]<-"intxn"; names(soop)[3]<-"fitzpheno"; names(soop)[4]<-"p1"
rowcount<-1
for(i in 1:length(b)) { 
sub<-mat[mat$intxn==b[i],]
soop[i,1]<-sub[1,c("studyid")]
soop[i,2]<-sub[1,c("intxn")]
m1<-with(sub, lm(fitz~abs(phenodiff)))
soop[i,3]<-summary(m1)$coefficients[2]
soop[i,4]<-summary(m1)$coefficients[8]
rowcount<-rowcount+1
}

#get synchrony change per interaction
yes<-data.frame(array(0, c(length(b), 4)))
names(yes)[1] = "studyid"; names(yes)[2]<-"intxn"; names(yes)[3]<-"phenoshift"; names(yes)[4]<-"p2"
rowcount<-1
for(i in 1:length(b)) { 
sub<-mat[mat$intxn==b[i],]
yes[i,1]<-sub[1,c("studyid")]
yes[i,2]<-sub[1,c("intxn")]
m1<-with(sub, lm(abs(phenodiff)~year))
yes[i,3]<-summary(m1)$coefficients[2]
yes[i,4]<-summary(m1)$coefficients[8]
rowcount<-rowcount+1
}

#get fitness change per interaction
sup<-data.frame(array(0, c(length(b), 4)))
names(sup)[1] = "studyid"; names(sup)[2]<-"intxn"; names(sup)[3]<-"fitshift"; names(sup)[4]<-"p3"
rowcount<-1
for(i in 1:length(b)) { 
sub<-mat[mat$intxn==b[i],]
sup[i,1]<-sub[1,c("studyid")]
sup[i,2]<-sub[1,c("intxn")]
m1<-with(sub, lm(fitz~year))
sup[i,3]<-summary(m1)$coefficients[2]
sup[i,4]<-summary(m1)$coefficients[8]
rowcount<-rowcount+1
}
son<-merge(soop, yes, by=c("studyid","intxn"))
son2<-merge(son, sup, by=c("studyid","intxn"))


m1<-lm(abs(fitshift)~abs(fitzpheno), subset(son2, abs(fitzpheno)<0.23)); summary(m1)
son2$fitshift<-son2$fitshift+1
m1<-glm(abs(fitshift+1)~abs(fitzpheno), family=Gamma, son2); summary(m1)

m1<-lm(abs(fitshift)~abs(phenoshift), son2); summary(m1)



stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
newdata2<-merge(newdata, stud[,c("studyid","biome","terrestrial","big_interaction")], by="studyid")

ggplot(data=subset(newdata2, terrestrial!="aquatic"), aes(abs(phenodiff), fitz,colour=factor(intxn)))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE)+theme_bw()+theme(legend.position="none")

m1<-lme(fitz~terrestrial+abs(phenodiff), random=~1|intxn, data=newdata2); summary(m1)


fm1<-lme(fitz~abs(phenodiff), random=~1|intxn, weights=varPower(form=~year), method="ML", data=subset(newdata, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Corydalis ambigua" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031"), na.action=na.omit);

#####################
### FIGURES ########
#####################
pdf("graphs/fitzbyphenodiff/pos_intxn_type/allspp.pdf")
mat2<-subset(mat, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")
fm1<-lme(fitz~abs(phenodiff), random=~1|intxn, method="ML", mat2);
newdat<-expand.grid(mat2$phenodiff); names(newdat)[1]<-"phenodiff"
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)
ggplot(newdat, aes(x=abs(phenodiff),y=pred))+ 
geom_point(data=mat2, size=2.5, aes(x=abs(phenodiff), y=fitz, colour=factor(intxn)))+
geom_smooth(data=mat2, method="lm", se=FALSE, size=1, aes(x=abs(phenodiff), y=fitz, colour=factor(intxn), fill=factor(intxn)))+
geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")+
geom_line(size=2)+theme_bw()+geom_hline(linetype="dashed", size=0.5)+opts(legend.position="none",axis.title.x =theme_text(size=17), axis.text.x=theme_text(size=17), axis.text.y=theme_text(size=17), axis.title.y=theme_text(size=17, angle=90))+ylab("Performance change")+xlab("Relative timing (days)")
dev.off()

#both axes standardized
mat2<-subset(mat, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031" & studyid!="HMK026")
fm1<-lme(fitz~phenodiffz, random=~1|intxn, method="ML", mat2);
newdat<-expand.grid(mat2$phenodiffz); names(newdat)[1]<-"phenodiffz"
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)
ggplot(newdat, aes(x=phenodiffz,y=pred))+ 
geom_point(data=mat2, size=2.5, aes(x=phenodiffz, y=fitz, colour=factor(intxn)))+
geom_smooth(data=mat2, method="lm", se=FALSE, size=1, aes(x=phenodiffz, y=fitz, colour=factor(intxn), fill=factor(intxn)))+
geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")+
geom_line(size=2)+theme_bw()+geom_hline(linetype="dashed", size=0.5)+opts(legend.position="none",axis.title.x =theme_text(size=17), axis.text.x=theme_text(size=17), axis.text.y=theme_text(size=17), axis.title.y=theme_text(size=17, angle=90))+ylab("Performance change")+xlab("Relative timing")+xlim(0,100)+ylim(-2.5,3)


m1<-lme(abs(phenodiff)~fitz, random=~1|intxn, data=subset(newdata, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Corydalis ambigua" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031")); summary(m1)

## Change in fitness over time
ggplot(data=subset(newdata, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Corydalis ambigua" & spp1!="Diatom3 spp." & studyid!="HMK039" & studyid!="HMK031"), aes(year, fitz,colour=factor(intxn)))+ geom_point(size=2)+geom_smooth(method="lm", se=FALSE, size=1)+theme_bw()+opts(legend.position="none")
#FOR DISPLAY only
ggplot(data=subset(dater, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Corydalis ambigua" & spp1!="Diatom3 spp." & studyid!="HMK039" & intxn!="Thermocyclops oithonoides NA NA" & fitz_yr1<10), aes(year, fitz_yr1,colour=factor(intxn)))+ geom_point(size=2)+geom_smooth(method="lm", se=FALSE, size=1)+theme_bw()+opts(legend.position="none")+ylab("Change in performance")


data=subset(newdata, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Corydalis ambigua" & spp1!="Diatom3 spp." & studyid!="HMK039")
doPlotphenodiff <- function(sel_name) {
   subby <- data[data$intxn == sel_name,]
   ggobj <- ggplot(data=subby, aes(phenodiff, fitz)) +  geom_point(size=3) +
       geom_smooth(method="lm", se=FALSE)+theme_bw()
   print(ggobj)
   ggsave(sprintf("graphs/fitzbyphenodiff/pos_intxn_type/%s.pdf", sel_name))
}
   
lapply(unique(data$intxn), doPlotphenodiff)

data=subset(dater, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Corydalis ambigua" & spp1!="Diatom3 spp." & studyid!="HMK039")
doPlotphenodiff <- function(sel_name) {
   subby <- data[data$intxn == sel_name,]
   ggobj <- ggplot(data=subby, aes(year, fitnessvalue)) +  geom_point(size=3) +
       geom_smooth(method="lm", se=FALSE)+theme_bw()
   print(ggobj)
   ggsave(sprintf("graphs/fitbyyr/%s.pdf", sel_name))
}
   
lapply(unique(data$intxn), doPlotphenodiff)

#using phenodiff_base
pdf("graphs/fitzbyphenodiff/pos_intxn_type/allspp_baseline.pdf")
ggplot(data=subset(dater, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & spp1!="Ammodytes marinus" & spp1!="Corydalis ambigua" & spp1!="Diatom3 spp." & studyid!="HMK039"), aes(abs(phenodiff_base), fitz_yr1,colour=factor(intxn)))+ geom_point(size=2)+geom_smooth(method="lm", se=FALSE, size=1)+theme_bw()+opts(legend.position="none")+geom_vline(linetype="dashed", size=0.25)+geom_hline(linetype="dashed", size=0.25)
dev.off()

## fitchange by phenochange
Bgroups<-unique(dater$intxn); b<-Bgroups; b<-as.character(b)
newest<-data.frame(array(0, c(length(b)), 6))
names(newest)[1]<-"studyid"; names(newest)[2]<-"intxn"; names(newest)[3]<-"fitnesscomp"; names(newest)[4]<-"fitnesstype"; names(newest)[5]<-"fitbyyr"; names(newest)[6]<-"phenodiffbyyr"
rowcount<-1
for(i in 1:length(b)){
spp<-dater[dater$intxn==b[i],]
asdf<-rowcount+(nrow(spp)-1)
newest[i,c("studyid")]<-spp[1,c("studyid")]
newest[i,c("intxn")]<-b[i]
newest[i,c("fitnesscomp")]<-spp[1,c("fitnesscomp")] 
newest[i,c("fitnesstype")]<-spp[1,c("fitnesstype")] 
m1<-with(spp, lm(fitnessvalue~year))
newest[i,c("fitbyyr")]<-summary(m1)$coefficients[2]
m2<-with(spp, lm(phenodiff~year))
newest[i,c("phenodiffbyyr")]<-summary(m2)$coefficients[2]
}
ggplot(data=subset(newest, studyid!="SET005" & intxn!="Acartia hudsonica NA NA" & intxn!="Corydalis ambigua NA NA" & intxn!="Diatom3 spp. NA NA" & studyid!="HMK039"), aes(phenodiffbyyr, fitbyyr))+ geom_point(size=2)+geom_smooth(method="lm", se=FALSE, size=1)+theme_bw()+opts(legend.position="none")+geom_vline(linetype="dashed", size=0.25)+geom_hline(linetype="dashed", size=0.25)


