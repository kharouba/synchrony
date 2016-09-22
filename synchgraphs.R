### Started 23 January 2014 ###
### By Lizzie & Heather ###

## Hopefully this will be useful plotting f(x)s for the trophic synchrony work ##
## Right now, just for observations ##

## Be sure to set the wd and add folders for:
# input, output, within output put 'graphs'#
# within graphs put two folders: doybyyr, and phenodiff ##

options(stringsAsFactors=FALSE)
library(ggplot2)

#setwd("/Volumes/Music/UBC/synchrony project/")
setwd("/users/kharouba/google drive/UBC/synchrony project/")
source("/users/kharouba/google drive/UBC/multiplot.R")

dater <- read.csv("analysis/spp_phenodata2.csv", header=TRUE)
daterdiff <- read.csv("analysis/int_phenodata.csv", header=TRUE)

######################
### QUALITY CONTROL ##
#####################
# basic idea here (showing all data, it's slow!)
ggplot(data=dater, aes(year, phenovalue,
       colour=factor(species))) +
       geom_point(shape=1) + facet_wrap(~intid) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(species)))

# make a f(x), which I adapted from one I found online
# and use lapply
doPlot <- function(sel_name) {
   subby <- dater[dater$intid == sel_name,]
   ggobj <- ggplot(data=subby, aes(year, phenovalue,
       colour=factor(species))) +
       geom_point(size=3) + facet_wrap(~intid) +
       geom_smooth(method="lm", se=FALSE, size=1, 
       aes(fill = factor(species)))+theme_bw()+
       theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("Day of year")+xlab("year")
   print(ggobj)
   ggsave(sprintf("/users/kharouba/google drive/UBC/synchrony project/graphs/2016/int%s.pdf", sel_name))
}
   
lapply(unique(dater$intid), doPlot)

#load yo3 from stats.R
mat<-subset(yo3, length>4)
doPlotphenodiff <- function(sel_name) {
   subby <- mat[mat$intid == sel_name,]
   ggobj <- ggplot(data=subby, aes(year, phenodiff_base),
       colour=factor(intid)) +  geom_point(size=4) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(intid)))+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("Phenodiff")+xlab("year")
   print(ggobj)
   ggsave(sprintf("graphs/phenodiff/int%s_base.pdf", sel_name))
}
   
lapply(unique(mat$intid), doPlotphenodiff)

doPlotphenodiff <- function(sel_name) {
   subby <- yo3[yo3$intid == sel_name,]
   ggobj <- ggplot(data=subby, aes(abs(z), abs(phenodiff_base),
       colour=factor(intid))) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(intid)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("graphs/phenodiffbyz/int%s_base.pdf", sel_name))
}
   
lapply(unique(yo3$intid), doPlotphenodiff)


data=subset(newdata3, length>4 & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188")
doPlotphenodiff <- function(sel_name) {
   subby <- data[data$intid == sel_name,]
   ggobj <- ggplot(data=subby, aes(year, phenodiff_var),
       colour=factor(intid)) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(intid)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("graphs/phenodiffvarbyyr/int%s.pdf", sel_name))
}
   
lapply(unique(data$intid), doPlotphenodiff)

#with absolute values
daterdiff<-yo3
daterdiff<-subset(yo3, length>4)
doPlotphenodiff <- function(sel_name) {
   subby <- daterdiff[daterdiff$intid == sel_name,]
   ggobj <- ggplot(data=subby, aes(year, abs(phenodiff),
       colour=factor(intid))) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(intid)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("graphs/phenodiff/%s_abs.pdf", sel_name))
}
   
lapply(unique(daterdiff$intid), doPlotphenodiff)

#phenodiff_Z
#load newdata from stats.R
daterdiff<-newdata
daterdiff<-subset(daterdiff, length>4)
doPlotphenodiff <- function(sel_name) {
   subby <- daterdiff[daterdiff$intid == sel_name,]
   ggobj <- ggplot(data=subby, aes(year, phenodiff_z,
       colour=factor(intid))) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(intid)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("graphs/phenodiff_z/int%s.pdf", sel_name))
}
   
lapply(unique(daterdiff$intid), doPlotphenodiff)

# phenodiff ~ env deviations
doPlotphenodiff <- function(sel_name) {
   subby <- tot[tot$studyid == sel_name,]
   ggobj <- ggplot(data=subby, aes(z, abs(phenodiff),
       colour=factor(intid))) +  geom_point(shape=1) +
       geom_smooth(method="lm", formula=y~poly(x,2), se=FALSE, 
       aes(fill = factor(intid)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("graphs/phenodiff/z scores/%s_abs.pdf", sel_name))
}
   
lapply(unique(tot$studyid), doPlotphenodiff)

mom<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/int_phenoenv.csv", header=TRUE, na.strings="NA", as.is=TRUE)
# phenodiff ~ env deviations
doPlotphenodiff <- function(sel_name) {
   subby <- mom[mom$studyid == sel_name,]
   ggobj <- ggplot(data=subby, aes(z, z_pheno,
       colour=factor(intid))) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(intid)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("graphs/phenodiff/z scores/%s.pdf", sel_name))
}
   
lapply(unique(mom$studyid), doPlotphenodiff)

#phenodiff by z
mom2<-merge(mom, yo3[,c("studyid","intid","year","phenodiff")], by=c("studyid","intid","year"))
mom3<-unique(mom2[,c("studyid","intid","year","z","phenodiff")])
ggplot(mom3, aes(z, abs(phenodiff),colour=factor(intid)))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(intid)))+theme_bw()+opts(legend.position="none")


## phenovalue~z I.E. TEMP SENSITIVITY
mom<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/tempsens.csv", header=TRUE, na.strings="NA", as.is=TRUE)
doPlotphenodiff <- function(sel_name) {
   subby <- mom[mom$intid == sel_name,]
   ggobj <- ggplot(data=subby, aes(z, phenovalue,
       colour=factor(spp))) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(spp)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("/Volumes/Music/UBC/synchrony project/graphs/doybyz/int%s.pdf", sel_name))
}
lapply(unique(mom$intid), doPlotphenodiff)

mom<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/tempsens.csv", header=TRUE, na.strings="NA", as.is=TRUE)
su<-subset(mom, studyid!="HMK035" & studyid!="EMW001")
doPlotphenodiff <- function(sel_name) {
   subby <- su[su$intid == sel_name,]
   ggobj <- ggplot(data=subby, aes(envvalue, phenovalue,
       colour=factor(spp))) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(spp)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("/Volumes/Music/UBC/synchrony project/graphs/doybytemp/int%s.pdf", sel_name))
}
lapply(unique(su$intid), doPlotphenodiff)

# tens sens by taxa
mom<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/tempsens.csv", header=TRUE, na.strings="NA", as.is=TRUE)
taxer <- read.csv("taxa.csv", header=TRUE)
taxer$latbi <- paste(taxer$genus, taxer$species)
taxer<-taxer[,c("studyid","taxa","latbi")]
names(taxer)[3]<-"species"
so4<-merge(mom, taxer, by=c("studyid","species"))
ggplot(data=subset(so4, studyid!="HMK035" & studyid!="EMW001"), aes(z, phenovalue,colour=factor(species))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(taxa)))+theme_bw()+theme(legend.position="none")

#temp change
clim3$stud_spp<-with(clim3, paste(studyid,"_",species))
doPlotphenodiff <- function(sel_name) {
   subby <- clim3[clim3$stud_spp == sel_name,]
   ggobj <- ggplot(data=subby, aes(year, envvalue,
       colour=factor(stud_spp))) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(stud_spp)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("/Volumes/Music/UBC/synchrony project/graphs/tempbyyr/%s.pdf", sel_name))
}
with(data=subset(clim3, envfactor=="temperature" & envtype=="air" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), lapply(unique(stud_spp), doPlotphenodiff))

#envvalue change
doPlotphenodiff <- function(sel_name) {
   subby <- mom[mom$studyid == sel_name,]
   ggobj <- ggplot(data=subby, aes(year, envvalue,
       colour=factor(stud_spp))) +  geom_point(shape=1) +
       geom_smooth(method="lm", se=FALSE, 
       aes(fill = factor(stud_spp)))+theme_bw()
   print(ggobj)
   ggsave(sprintf("/Volumes/Music/UBC/synchrony project/graphs/envbyyr/%s.pdf", sel_name))
}
with(data=subset(studyid!="EMW001" & studyid!="HMK035"), lapply(unique(studyid), doPlotphenodiff))

################################
###### FIGURES, FIGURES #####
##############################
#figure for Andrew gelman- spp pheno shifts by intxn
rowcount<-1
mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
Bgroups<-unique(mat$intid); b<-Bgroups; b<-as.character(b)
nup<-data.frame(array(0, c(80, 4)))
names(nup)[1] = "studyid"; names(nup)[2]<-"species"; names(nup)[3]<-"intid"; names(nup)[4]<-"doybyyr"; 
for(i in 1:length(b)) { 
yo<-mat[mat$intid==b[i],]
nup[rowcount,1]<-yo[1,c("studyid")]
nup[rowcount+1,1]<-yo[1,c("studyid")]
nup[rowcount,2]<-yo[1,c("spp1")]
nup[rowcount+1,2]<-yo[1,c("spp2")]
nup[rowcount,3]<-yo[1,c("intid")]
nup[rowcount+1,3]<-yo[1,c("intid")]
m1<-with(yo, lm(neg_phenovalue~year))
m2<-with(yo, lm(pos_phenovalue~year))
nup[rowcount,4]<-summary(m1)$coefficients[2]
nup[rowcount+1,4]<-summary(m2)$coefficients[2]
rowcount<-rowcount+2
}
nup$spp<-with(nup, c(1:2))
ggplot(data=nup, aes(x=spp, y=doybydec))+
geom_path(aes(group=factor(intid), colour=factor(intid)), size=1)+
geom_point(aes(colour=factor(intid), size=1))+
theme_bw()+theme(legend.position="none", panel.grid.major = element_blank(), axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("Phenological shift (days/decade)")+xlab("species role")

#Version2#
----------
#Figure 1#
-----------
mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")

a<-ggplot(mat, aes(x=year, y=phenodiff))+
geom_point(size=2, aes(colour=factor(intid)))+
geom_smooth(size=0.75, method="lm", se=FALSE, aes(colour=factor(intid), fill=factor(intid)))+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("relative timing (days)")+xlab("year")

#baseline synchrony
b<-ggplot(mat, aes(x=year, y=abs(phenodiff_base)))+
geom_point(size=2, aes(colour=factor(intid)))+
geom_smooth(size=0.75, method="lm", se=FALSE, aes(colour=factor(intid), fill=factor(intid)))+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=14, angle=90))+ylim(0,115)+ylab("timing relative to baseline (days)")
multiplot(a,b, cols=2)

#baseline
ggplot(best, aes(x=coef))+
geom_histogram(colour="black", fill="grey")+
geom_vline(xintercept=0, linetype="solid",size=1)+theme_bw()+opts(legend.position="none",axis.title.x =theme_text(size=15), axis.text.x=theme_text(size=15), axis.text.y=theme_text(size=15), axis.title.y=theme_text(size=14, angle=90))+ylab("Number of interactions")

#shifts
c<-ggplot(yep, aes(x=coef, fill=terrestrial))+geom_histogram(alpha=0.3,  colour="black", position="identity")+theme_bw()+geom_vline(xintercept=0, linetype="solid",size=1)+theme(legend.justification=c(1,0), legend.position=c(1,0.6), axis.title.x =element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=14, angle=90))+ylab("Number of interactions")+xlab("slope (days/year)")
multiplot(b,a,c, cols=2)

##Figure 2##
mom<-read.csv("tempsens_nov2.csv", header=TRUE, na.strings="NA", as.is=TRUE)
mom<-na.omit(mom)
mom$spp<-as.factor(mom$spp)
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(mom, stud[,c("studyid","short_site","biome","terrestrial","big_interaction")], by="studyid")
mom<-yep

fm1<-lme(phenovalue~z, random=~1|studyid/species, method="ML", data=mat, na.action=na.omit);
#fm1<-lme(phenovalue~z+spp, random=~1|studyid/speci es, method="ML", weights=varIdent(form=~1|spp), data=mat, na.action=na.omit);
#newdat <- with(mat,expand.grid(spp=levels(spp),z=seq(min(z),max(z),25)))
newdat<-expand.grid(mat$z); names(newdat)[1]<-"z"
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

a<-ggplot(newdat, aes(x=z, y=pred))+
geom_point(data=mat, size=2, aes(x=z, y=phenovalue, colour=factor(intid), shape=factor(terrestrial)))+
geom_smooth(data=mat, size=0.75, method="lm", se=FALSE, aes(x=z, y=phenovalue, colour=factor(intid), fill=factor(intid)))+geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")+
geom_line(size=2)+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("doy")+xlab("z")
multiplot(a,b, cols=2)


#Version1#
# Figure 1 = phenodiff~year- different estimates
load yo3 from stats.R
mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")

fm1<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid, method="ML", mat, na.action=na.omit);
newdat<-expand.grid(mat$year); names(newdat)[1]<-"year"
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

a<-ggplot(newdat, aes(x=year, y=pred))+
geom_point(data=mat, size=2, aes(x=year, y=sqrt(abs(phenodiff)), colour=factor(intid)))+
geom_smooth(data=mat, size=0.75, method="lm", se=FALSE, aes(x=year, y=sqrt(abs(phenodiff)), colour=factor(intid), fill=factor(intid)))+geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")+geom_line(size=2)+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("sqrt(timing)")+xlab("year")

#baseline synchrony
fm1<-lme(sqrt(abs(phenodiff_base))~year, random=~1|studyid/intid, method="ML", mat, na.action=na.omit);
newdat<-expand.grid(mat$year); names(newdat)[1]<-"year"
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

b<-ggplot(newdat, aes(x=year, y=pred))+
geom_point(data=mat, size=2, aes(x=year, y=sqrt(abs(phenodiff_base)), colour=factor(intid)))+
geom_smooth(data=mat, size=0.75, method="lm", se=FALSE, aes(x=year, y=sqrt(abs(phenodiff_base)), colour=factor(intid), fill=factor(intid)))+geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")+
geom_line(size=2)+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=14, angle=90))+ylab("sqrt(timing relative to baseline)")
multiplot(a,b, cols=2)

#aquatic vs. terrestrial
mat<-subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
te<-subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="terrestrial")
aq<-subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="aquatic")
fm1<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid,  method="ML",  na.action=na.omit, te); summary(m4)
fm2<-lme(sqrt(abs(phenodiff))~year, random=~1|studyid/intid,  method="ML",  na.action=na.omit, aq); summary(m4)

#newdat<-expand.grid(mat$year, mat$terrestrial); names(newdat)[1]<-"year"; names(newdat)[2]<-"terrestrial"
newdat<-expand.grid(te$year); names(newdat)[1]<-"year";
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
#newdat$comb<-interaction(newdat$year, newdat$terrestrial)
pd<-position_dodge(width=0.4)

newdat2<-expand.grid(aq$year); names(newdat2)[1]<-"year";
newdat2$pred<-predict(fm2, newdat2, level=0)
Designmat2 <- model.matrix(eval(eval(fm2$call$fixed)[-2]), newdat2[-3])
predvar2 <- diag(Designmat2 %*% fm2$varFix %*% t(Designmat2))
newdat2$SE <- sqrt(predvar2) 
newdat2$SE2 <- sqrt(predvar2+fm2$sigma^2)
pd<-position_dodge(width=0.4)

ggplot(data=te, aes(x=year, y=sqrt(abs(phenodiff))))+geom_point(data=te, size=2, aes(x=year, y=sqrt(abs(phenodiff))))+
geom_smooth(data=newdat, size=0.75, method="lm", se=FALSE, aes(x=year, y=pred, colour="orange"))+geom_ribbon(data=newdat,aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")

ggplot(data=aq, aes(x=year, y=sqrt(abs(phenodiff))))+geom_point(data=te, size=2, aes(x=year, y=sqrt(abs(phenodiff))))+
geom_smooth(data=newdat2, size=0.75, method="lm", se=FALSE, aes(x=year, y=pred, colour="blue"))+geom_ribbon(data=newdat2, aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")
+ylim(0,125)+geom_line(size=2)

c<-ggplot(data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(x=year, y=sqrt(abs(phenodiff)), colour=factor(terrestrial))) +  geom_point(size=2) + geom_smooth(method="lm", se=FALSE, aes(fill = factor(terrestrial)), size=1)+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("sqrt(timing)")+xlab("year")

OR
ggplot(data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(x=year, y=abs(phenodiff), shape=factor(terrestrial))) +  geom_point(data=subset(yep, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), size=2, aes(colour=factor(intid))) + geom_smooth(method="lm", se=FALSE, aes(fill = factor(terrestrial)), size=1)+theme_bw()+opts(legend.position="none",axis.title.x =theme_text(size=15), axis.text.x=theme_text(size=15), axis.text.y=theme_text(size=15), axis.title.y=theme_text(size=15, angle=90))+ylab("sqrt(timing)")+xlab("year")


yep2<-subset(yep, interaction=="predation" | interaction=="herbivory")
d<-ggplot(data=subset(yep2, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(x=year, y=sqrt(abs(phenodiff)), colour=factor(interaction))) +  geom_point(size=2, aes(shape=factor(terrestrial))) + geom_smooth(method="lm", se=FALSE, aes(fill = factor(interaction)), size=1)+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("sqrt(timing)")+xlab("year")

multiplot(b,a,c,d, cols=2)


#just those with temp sensitivities
mom2<-unique(mom[,c("studyid","intid","year")])
yess<-merge(yo3, mom2, by=c("studyid","intid","year"))
ggplot(data=subset(yess, length>4 & mismatch!="1" & studyid!="SET005"), aes(year, abs(phenodiff),colour=factor(intid)))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(intid)))+theme_bw()+theme(legend.position="none")


ggplot(newdat, aes(x=year, y=pred))+
geom_point(data=subset(yo3, length>4 & mismatch!="1" & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(x=year, y=abs(phenodiff), colour=factor(intid)))+
geom_smooth(data=subset(yo3, length>4 & mismatch!="1"  & studyid!="SET005" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), method="lm", se=FALSE, aes(x=year, y=abs(phenodiff), colour=factor(intid), fill=factor(intid)))+
geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")+
geom_line(size=2)+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("abs(Difference in doy)")+xlab("year")

----------------------------
# figure 2== phenovalue~year
----------------------------
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

mat<-subset(so3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188" & intid!="221" & intid!="222")

fm1<-lme(phenovalue~year, random=~1|studyid/species, method="REML", data=mat, na.action=na.omit); 
newdat<-expand.grid(mat$year); names(newdat)[1]<-"year"
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

b<-ggplot(newdat, aes(x=year, y=pred))+
geom_point(data=mat, aes(x=year, y=phenovalue, colour=factor(intid), shape=factor(terrestrial)))+
geom_smooth(data=mat, method="lm", se=FALSE, aes(x=year, y=phenovalue, colour=factor(intid), fill=factor(intid)))+
geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")+
geom_line(size=2)+theme_bw()+scale_x_continuous(breaks=seq(1950, 2012, 15))+ 
theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("Day of year")+xlab("year")

mom<-read.csv("tempsens_nov.csv", header=TRUE, na.strings="NA", as.is=TRUE)
mom<-na.omit(mom); mom$spp<-as.factor(mom$spp)
#study characteristics
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(mom, stud[,c("studyid","short_site","biome","terrestrial","big_interaction")], by="studyid")
mom<-yep
mat<-subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188")

fm1<-lme(phenovalue~z, random=~1|studyid/species, method="ML", data=mat, na.action=na.omit);
#fm1<-lme(phenovalue~z+spp, random=~1|studyid/speci es, method="ML", weights=varIdent(form=~1|spp), data=mat, na.action=na.omit);
#newdat <- with(mat,expand.grid(spp=levels(spp),z=seq(min(z),max(z),25)))
newdat<-expand.grid(mat$z); names(newdat)[1]<-"z"
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

a<-ggplot(newdat, aes(x=z, y=pred))+
geom_point(data=mat, size=2, aes(x=z, y=phenovalue, colour=factor(intid), shape=factor(terrestrial)))+
geom_smooth(data=mat, size=0.75, method="lm", se=FALSE, aes(x=z, y=phenovalue, colour=factor(intid), fill=factor(intid)))+geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.3, fill="blue")+
geom_line(size=2)+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("doy")+xlab("z")
multiplot(a,b, cols=2)


pdf("graphs/doybyyr/allspp.pdf")
ggplot(data=subset(so3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(year, phenovalue,colour=factor(species)))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(species)))+theme_bw()+theme(legend.position="none")
dev.off()

pdf("graphs/doybyyr/all.pdf")
ggplot(data=subset(so2, length>4), aes(year, phenovalue))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE)+theme_bw()+theme(legend.position="none")
dev.off()
---------------------------------
# Figure 3: temp change
--------------------------------
source("/Volumes/Music/UBC/multiplot.R")
clim<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/climate.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
clim<-subset(clim, phenophase!="start" & extra!="10")
clim$envvalue<-as.numeric(clim$envvalue)
clim$envfactor[clim$envfactor=="temperaure"] <- "temperature"
clim$stud_block<-with(clim, paste(studyid,"_",seasonalblock))
sites<-subset(clim, studyid=="HMK018" | studyid=="HMK019" | studyid=="HMK023" & site!="tomakomai")
nosites<-subset(clim, site=="tomakomai" | studyid!="HMK018" & studyid!="HMK019" & studyid!="HMK023" & stud_block!="HMK029 _ april" & stud_block!="HMK029 _ may" & stud_block!="HMK029 _ june")
nosites<-nosites[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")]

new<-with(sites, aggregate(envvalue, by=list(studyid, year, species), FUN=mean, na.rm=T)) # across sites
names(new)[1]<-"studyid"; names(new)[2]<-"year"; names(new)[3]<-"species"; names(new)[4]<-"envvalue"
new<-new[order(new$studyid),]
sites2<-merge(sites[,c("studyid","envfactor","envunits","envtype","year","species")], new, by=c("studyid","year","species"))
sites3<-unique(sites2[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")])

clim2<-rbind(nosites, sites3) # all years with data
#merge with spp data so calculating env change only over the years of interaction
total3<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
total3<-na.omit(total3)
clim3<-merge(clim2, total3[,c("studyid","year")], by=c("studyid","year"))
clim3<-unique(clim3[,c("studyid","year","envfactor","envunits","envtype","species","envvalue")])
clim3$stud_spp<-with(clim3, paste(studyid,"_",species))

###AIR###
fm1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype=="air" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), na.action=na.omit);
newdat<-expand.grid(clim3$year); names(newdat)[1]<-"year"
newdat$pred<-predict(fm1, newdat, level=0)

Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

air<-ggplot(newdat, aes(x=year, y=pred))+
geom_point(data=subset(clim3, envfactor=="temperature" & envtype=="air" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), aes(x=year, y=envvalue, colour=factor(stud_spp)))+
geom_smooth(data=subset(clim3, envfactor=="temperature" & envtype=="air" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), method="lm", se=FALSE, aes(x=year, y=envvalue, colour=factor(stud_spp), fill=factor(stud_spp)))+
geom_line(size=2)+geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.1, fill="blue")+theme_bw()+xlim(1970,2013)+theme(legend.position="none", axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab(expression(paste("Temperature (",degree,"C)")))+xlab("year")

##WATER##
fm1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype=="water" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis"), na.action=na.omit);
newdat<-expand.grid(clim3$year); names(newdat)[1]<-"year"
newdat$pred<-predict(fm1, newdat, level=0)

Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

water<-ggplot(newdat, aes(x=year, y=pred))+
geom_point(data=subset(clim3, envfactor=="temperature" & envtype=="water" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), aes(x=year, y=envvalue, colour=factor(stud_spp)))+
geom_smooth(data=subset(clim3, envfactor=="temperature" & envtype=="water" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), method="lm", se=FALSE, aes(x=year, y=envvalue, colour=factor(stud_spp), fill=factor(stud_spp)))+
geom_line(size=2)+geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.1, fill="blue")+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_blank())+xlab("year")
multiplot(air, water, cols=2)

##BOTH##
fm1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis"), na.action=na.omit);
newdat<-expand.grid(mom$z); names(newdat)[1]<-"z"
newdat$pred<-predict(fm1, newdat, level=0)
Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

ggplot(newdat, aes(x=year, y=pred))+
geom_point(data=subset(clim3, envfactor=="temperature" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), aes(x=year, y=envvalue, colour=factor(stud_spp)))+
geom_smooth(data=subset(clim3, envfactor=="temperature" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), method="lm", se=FALSE, aes(x=year, y=envvalue, colour=factor(stud_spp), fill=factor(stud_spp)))+
geom_line(size=2)+geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.1, fill="blue")+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab(expression(paste("Temperature (",degree,"C)")))+xlab("year")

pdf("graphs/tempbyyr.pdf")
ggplot(data=subset(clim2, envfactor=="temperature" & envtype=="air" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), aes(year, envvalue, colour=factor(stud_spp)))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(stud_spp)))+theme_bw()+theme(legend.position="none")
dev.off()

# Figure 4: temp sensitivity
mom<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/tempsens.csv", header=TRUE, na.strings="NA", as.is=TRUE)

fm1<-lme(phenovalue~z, random=~1|studyid/intid/species, method="REML", data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(fm1)
newdat<-expand.grid(mom$z); names(newdat)[1]<-"z"
newdat$pred<-predict(fm1, newdat, level=0)

Designmat <- model.matrix(eval(eval(fm1$call$fixed)[-2]), newdat[-3])
predvar <- diag(Designmat %*% fm1$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+fm1$sigma^2)
pd<-position_dodge(width=0.4)

ggplot(newdat, aes(x=z, y=pred))+
geom_line(size=2)+
geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.1, fill="blue")+theme_bw()+ylim(0,300)+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("Day of year")+xlab("z score")

ggplot(newdat, aes(x=z, y=pred))+
geom_point(data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), aes(x=z, y=phenovalue, colour=factor(species)))+
geom_smooth(data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), method="lm", se=FALSE, aes(x=z, y=phenovalue, colour=factor(species), fill=factor(species)))+
geom_line(size=2)+
geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.1, fill="blue")+theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("Day of year")+xlab("z score")

---
pdf("/Volumes/Music/UBC/synchrony project/graphs/doybyz/allspp.pdf")
ggplot(data=mom, aes(z, phenovalue, colour=factor(intid)))+ ylim(30,200)+geom_point(size=1)+geom_line(pred)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(intid)))+theme_bw()+theme(legend.position="none")
dev.off()


#Figure x- phenodiff ~ env deviations
mom<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/int_phenoenv.csv", header=TRUE, na.strings="NA", as.is=TRUE)
pdf("graphs/z_phenobyz.pdf")
ggplot(data=mom, aes(z_pheno, z, colour=factor(intid)))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(intid)))+theme_bw()+theme(legend.position="none")
dev.off()

pdf("graphs/z_phenobyz.pdf")
ggplot(data=mom, aes(z_pheno, z))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE)+theme_bw()+theme(legend.position="none")
dev.off()


pdf("graphs/phenodiff/allspp.pdf")
ggplot(data=subset(yo3, length>4), aes(year, abs(phenodiff_base),colour=factor(intid)))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(intid)))+theme_bw()+theme(legend.position="none")
dev.off()

pdf("graphs/phenodiff/allspp_terres.pdf")
ggplot(data=subset(yo2, length>4), aes(year, abs(phenodiff),colour=factor(terrestrial)))+ geom_point(size=1)+geom_smooth(method="lm", se=FALSE, aes(fill=factor(terrestrial)))+theme_bw()+theme(legend.position="none")
dev.off()


# MAP  #
library(rworldmap)
newmap<-getMap(resolution="coarse")
dad<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/studies_geog.csv", header=TRUE, na.strings="NA", as.is=TRUE)
coordinates(dad)<-c("x","y")
plot(newmap)
plot(dad, add=T, pch=20, col="red")

