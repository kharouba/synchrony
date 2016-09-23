library(ggplot2)
library(lattice)
 #run lizzie's code: synchbuildintxn_hk

## Lizzie made some updates to this code so it will run on the new file structure ##
## I updated it so it can files and had it now write to a new output folder ##
## but it still stalls out on all the figures and such so it really needs to be cleaned up ##

# rm(list=ls())

source("synchbuildintxn_hk.R")

row.names=FALSE

#studies
#stud<-read.csv("studies.csv", header=TRUE, na.strings="NA", as.is=TRUE)

#OBSERVATIONAL
obs<-read.csv("2013_2015/input/obs_start.csv", header=TRUE, na.strings="NA", as.is=TRUE)
#species across multiple studies- solution-add unique id
obs$spp<-with(obs, paste(genus, species, sep=" "))
obs$spp2<-with(obs, paste(genus_2nd, spp_2nd, sep=" "))
obs<-subset(obs, phenophase!="spring increase") #to choose one phenophase for HMK035, int192

obs$spp[obs$spp=="Hylurgops rugipennis pinnifex"] <- "Hylurgops rugipennis pinifex"
obs$spp[obs$spp=="Monochamus spp"] <- "Monochamus spp."
obs$spp[obs$spp=="Quercus spp"] <- "Quercus spp."
#obs2<-read.csv("rough_taxa.csv", header=TRUE, na.strings="NA", as.is=TRUE)
#obs2$spp<-with(obs2, paste(genus, species, sep=" "))


int<-intxn
int<-int[order(int$studyid),]
int$intid<-1:nrow(int)
# spp1 = negative species in trophic interaction and spp2= positive species

#remove one row from HMK004
#int<-int[c(1:6,8:nrow(int)),]
int$spp1[int$spp1=="Quercus pubescens"] <- "Quercus spp." #HMK004- does not specify oak spp in figure
int$spp1[int$spp1=="Quercus petreae"] <- "Quercus spp."
int<-subset(int, intid!="20") #& intid!="154" changed feb 2015
# & intid!="181" & intid!="182")

some<-subset(int, intid>=154 & intid<=236)
some$intid<-with(some, intid+1)

some2<-subset(int, intid<154)
int<-rbind(some2, some)

write.csv(int, "output/intxn.csv")



#roughintxn2<-unique(roughintxn[,c("studyid","neglatbi","poslatbi")])

#only studies that included doy as units (e.g HMK013 uses proportion) and multiple years
obs<-subset(obs, phenounits=="doy" & predictvalue!="NA")

#to eliminate non-obs studies
yes<-unique(obs[,c("studyid")]);yes<-as.data.frame(yes); names(yes)[1]<-"studyid"
int2<-merge(yes,int, by="studyid")
#int2<-subset(int2, neglatbi!="Rubus spp.") #& neglatbi!="Ammodytes marinus"

#Multiple sites
int2<-subset(int2,studyid!="EMW002")
new<-data.frame(array(0, c(nrow(obs), 5)))
names(new)[1] = "studyid"; names(new)[2] = "year"; names(new)[3] = "spp"; names(new)[4] = "pheno"; names(new)[5] = "phenovalue";
rowcount<-1
Bgroups<-unique(int2$studyid); b<-Bgroups; b<-as.character(b)
for(i in 1:length(b)) { 
yo<-obs[obs$studyid==b[i],]
q<-unique(yo$site)
if(length(q)==1){ #if there is only 1 phenophase/study
asdf<-rowcount+(nrow(yo)-1)
new[rowcount:asdf,1]<-yo[,1]
new[rowcount:asdf,2]<-yo[,c("predictvalue")]
new[rowcount:asdf,3]<-yo[,c("spp")]
new[rowcount:asdf,4]<-yo[,c("pheno")]
new[rowcount:asdf,5]<-yo[,c("phenovalue")]
rowcount<-rowcount+nrow(yo)
}
if(length(q)!=1){ #if there is > 1 site/study
lol<-with(yo, aggregate(phenovalue, by=list(spp,predictvalue,pheno), FUN=median, na.rm=T))
names(lol)[1]<-"spp"; names(lol)[2]<-"year"; names(lol)[3]<-"pheno"; names(lol)[4]<-"phenovalue"
asdf2<-rowcount+(nrow(lol)-1)
new[rowcount:asdf2,1]<-yo[1:nrow(lol),1]
new[rowcount:asdf2,2]<-lol[,c("year")]
new[rowcount:asdf2,3]<-lol[,c("spp")]
new[rowcount:asdf2,4]<-lol[,c("pheno")]
new[rowcount:asdf2,5]<-lol[,c("phenovalue")]
rowcount<-rowcount+nrow(lol)
}
}

sites<-new
sites<-subset(sites, studyid!="0")

#deal with pheno obs on multiple hosts
pop<-subset(obs, studyid=="EMW002")
lol<-with(pop, aggregate(phenovalue, by=list(spp,predictvalue,pheno), FUN=median, na.rm=T))
names(lol)[1]<-"spp"; names(lol)[2]<-"year"; names(lol)[3]<-"pheno"; names(lol)[4]<-"phenovalue"
lol2<-subset(lol, spp!="Acer saccharum" & spp!="Ostrya virginiana" & spp!="Rubus spp.")
lol3<-subset(lol2, pheno=="start")
lol4<-subset(lol, spp!="Barypeithes pellucidus" & spp!="Phyllobius oblongus" & spp!="Polydrusus sericeus" & spp!="Sciaphilus asperatus")
lol5<-rbind(lol3,lol4)
lol5$studyid<-"EMW002"
lol6<-lol5[,c("studyid","year","spp","pheno","phenovalue")]
sites2<-rbind(lol6, sites)

#deal with multiple obs per year but at same site -- problem getting indiv years from figure
yo<-with(sites2, aggregate(phenovalue, by=list(studyid,year,spp), FUN=median, na.rm=T))
names(yo)[names(yo)=="Group.1"]<-"studyid"
names(yo)[names(yo)=="Group.2"]<-"year"
names(yo)[names(yo)=="Group.3"]<-"spp"
names(yo)[names(yo)=="x"]<-"phenovalue"
yo<-yo[order(yo$studyid),]
sites2<-yo

#Organize data and run linear model for both species
Bgroups<-unique(int2$spp1); b<-Bgroups; b<-as.character(b) #negative or first species in pair
Cgroups<-unique(int2$spp2); c<-Cgroups; c<-as.character(c) #positive or second species in pair

neg<-data.frame(array(0, c(nrow(sites2), 7)))
names(neg)[1] = "studyid"; names(neg)[2] = "neglatbi"; names(neg)[3] = "year"; 
names(neg)[4] = "neg_phenovalue"; names(neg)[5]<-"neg_slope"; names(neg)[6] = "intid"; names(neg)[7] = "int_type"
rowcount<-1
for(i in 1:length(b)) {  
species<-sites2[sites2$spp==b[i],]
a<-int2[int2$spp1==b[i],c("intid")]
to<-int2[int2$spp1==b[i],c("interaction")]  
if(length(a)==1){ #if there is only 1 interaction/study
asdf<-rowcount+(nrow(species)-1)
neg[rowcount:asdf,1]<-species[1,c("studyid")]
neg[rowcount:asdf,2]<-b[i] #neg spp
neg[rowcount:asdf,3]<-species[,c("year")]
neg[rowcount:asdf,4]<-species[,c("phenovalue")]
lm<-with(species, lm(phenovalue~year))
neg[rowcount:asdf,5]<-summary(lm)$coefficients[2]
neg[rowcount:asdf,6]<-int2[int2$spp1==b[i],c("intid")]
neg[rowcount:asdf,7]<-int2[int2$spp1==b[i],c("interaction")]
rowcount<-rowcount+nrow(species)
}

if(length(a)>1){ #if multiple interactions per spp
asdf<-rowcount+(nrow(species)-1)
asdf2<-rowcount+(nrow(species)*length(a))-1
neg[rowcount:asdf2,1]<-species[1,c("studyid")]
neg[rowcount:asdf2,2]<-b[i] #neg spp
neg[rowcount:asdf2,3]<-species[,c("year")]
neg[rowcount:asdf2,4]<-species[,c("phenovalue")]
lm<-with(species, lm(phenovalue~year))
neg[rowcount:asdf2,5]<-summary(lm)$coefficients[2]
for(j in 1:length(a)){
neg[rowcount:asdf,6]<-a[j]
neg[rowcount:asdf,7]<-to[j]
rowcount<-rowcount+nrow(species)
asdf<-rowcount+(nrow(species)-1)
}
}
}

#ts<-with(spp, aggregate(late_spr_mean, by=list(station), FUN=mean))

pos<-data.frame(array(0, c(nrow(obs), 7)))
names(pos)[1] = "studyid"; names(pos)[2] = "poslatbi"; names(pos)[3] = "year"; 
names(pos)[4] = "pos_phenovalue"; names(pos)[5]<-"pos_slope"; names(pos)[6] = "intid"; names(pos)[7] = "int_type";
rowcount<-1
for(i in 1:length(c)) {  
species<-sites2[sites2$spp==c[i],]
a<-int2[int2$spp2==c[i],c("intid")]
to<-int2[int2$spp2==c[i],c("interaction")]  
if(length(a)==1){
asdf<-rowcount+(nrow(species)-1)
pos[rowcount:asdf,1]<-species[1,c("studyid")]
pos[rowcount:asdf,2]<-c[i] #pos spp
pos[rowcount:asdf,3]<-species[,c("year")]
pos[rowcount:asdf,4]<-species[,c("phenovalue")]
lm<-with(species, lm(phenovalue~year))
pos[rowcount:asdf,5]<-summary(lm)$coefficients[2]
pos[rowcount:asdf,6]<-int2[int2$spp2==c[i],c("intid")]
pos[rowcount:asdf,7]<-int2[int2$spp2==c[i],c("interaction")]
rowcount<-rowcount+nrow(species)
}
if(length(a)>1){ #if multiple interactions per spp
asdf<-rowcount+(nrow(species)-1)
asdf2<-rowcount+(nrow(species)*length(a))-1
pos[rowcount:asdf2,1]<-species[1,c("studyid")]
pos[rowcount:asdf2,2]<-c[i] #pos spp
pos[rowcount:asdf2,3]<-species[,c("year")]
pos[rowcount:asdf2,4]<-species[,c("phenovalue")]
lm<-with(species, lm(phenovalue~year))
pos[rowcount:asdf2,5]<-summary(lm)$coefficients[2]
for(j in 1:length(a)){
pos[rowcount:asdf,6]<-a[j]
pos[rowcount:asdf,7]<-to[j]
rowcount<-rowcount+nrow(species)
asdf<-rowcount+(nrow(species)-1)
}
}
}


#raw data synchronized by year and interaction
total<-merge(int2,neg[,c("year","neg_phenovalue","neg_slope","intid")],by=c("intid"), all.x=T)
total2<-merge(total,pos[,c("year","pos_phenovalue","pos_slope","intid")],by=c("intid","year"), all.x=T)

total2$phenodiff<-with(total2, neg_phenovalue-pos_phenovalue)

blerk<-unique(obs[,c("studyid","phenofreq")])
total3<-merge(total2, blerk, by="studyid")

write.csv(total3, "output/int_phenodata.csv")

neg$spp<-"spp1"
pos$spp<-"spp2"
names(neg)[names(neg)=="neglatbi"]<-"species"
names(neg)[names(neg)=="neg_phenovalue"]<-"phenovalue"
names(neg)[names(neg)=="neg_slope"]<-"slope"

names(pos)[names(pos)=="poslatbi"]<-"species"
names(pos)[names(pos)=="pos_phenovalue"]<-"phenovalue"
names(pos)[names(pos)=="pos_slope"]<-"slope"

all<-rbind(neg,pos)
all<-subset(all, studyid!="0")
write.csv(all, "output/spp_phenodata2.csv")

##Figure 1- raw data doy~year
ggplot(all, aes(x=year, y=phenovalue)) + 
geom_point(aes(color=species, size=1))+
geom_smooth(aes(group=species),method="lm", size=1, colour="black", se=FALSE)

##Figure 2- yearly phenodiff~year
library(ggplot2)
mat<-total2[,c("year","studyid","phenodiff","intid")]
mat<-na.omit(mat)
#all
ggplot(mat, aes(x=year, y=phenodiff)) + 
geom_point(aes(color=studyid, size=1))+
geom_hline(yintercept=0, linetype="dashed")+
#geom_smooth(method="lm", size=1, colour="black", se=FALSE)+
theme_bw()+xlab("Year") + ylab("Difference in doy")+
# opts(axis.title.x = theme_text(size=15), axis.text.x=theme_text(size=15), axis.text.y=theme_text(size=15), axis.title.y=theme_text(size=15, angle=90))
#by study
ggplot(mat, aes(x=year, y=phenodiff)) + 
geom_point(aes(color=studyid, size=1))+facet_grid(~studyid)+
geom_hline(yintercept=0, linetype="dashed")+
#geom_smooth(method="lm", size=1, colour="black", se=FALSE)+
theme_bw()+xlab("Year") + ylab("Difference in doy")+
# opts(axis.title.x = theme_text(size=15), axis.text.x=theme_text(size=15), axis.text.y=theme_text(size=15), axis.title.y=theme_text(size=15, angle=90))

xyplot(abs(phenodiff)~year|studyid, groups=intid, data=mat)

#long
# mat2


######
#to get start and end difference in pheno
new<-data.frame(array(0, c(nrow(total2), 10)))
names(new)[1]<-"studyid"; names(new)[2]<-"intid"; names(new)[3]<-"spp1"; names(new)[4]<-"spp1_slope"; names(new)[5] = "spp2"; names(new)[6]<-"spp2_slope"; names(new)[7] = "min_year"; names(new)[8] = "max_year"; names(new)[9] = "start_phenodiff"; names(new)[10] = "end_phenodiff";
#Bgroups<-unique(total2$studyid); b<-Bgroups; b<-as.character(b)
Bgroups<-unique(total2$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)) {  
study<-total2[total2$intid==b[i],]
study<-na.omit(study)
s<-with(study,min(year)) #start
s2<-subset(study, year==s) #start
s2$start_phenodiff<-with(s2, neg_phenovalue-pos_phenovalue) #start
e<-with(study, max(year)) #end
e2<-subset(study, year==e)
e2$end_phenodiff<-with(e2, neg_phenovalue-pos_phenovalue)
asdf<-rowcount+(nrow(e2)-1)
new[rowcount:asdf,1]<-s2[,c("studyid")]
new[rowcount:asdf,2]<-b[i]
new[rowcount:asdf,3]<-s2[,c("spp1")]
new[rowcount:asdf,4]<-s2[,c("neg_slope")]
new[rowcount:asdf,5]<-s2[,c("spp2")]
new[rowcount:asdf,6]<-s2[,c("pos_slope")]
new[rowcount:asdf,7]<-s
new[rowcount:asdf,8]<-e
new[rowcount:asdf,9]<-s2[,c("start_phenodiff")]
new[rowcount:asdf,10]<-e2[,c("end_phenodiff")]
rowcount<-rowcount+nrow(e2)
}
new$slopediff<-with(new, spp1_slope-spp2_slope)
new$endstart<-with(new, end_phenodiff-start_phenodiff)
new<-new[order(new$start_phenodiff, new$end_phenodiff),]
new<-subset(new, studyid!="0")
write.csv(new, "output/endstart.csv")

plot(slopediff~endstart, ylim=c(-5,5), data=subset(new, studyid!="HMK006"))
xyplot(slopediff~endstart|studyid, groups=intid, data=subset(new, studyid!="HMK006"))

yo<-aggregate(total2["year"], by=total2[c("studyid","intid")], FUN=length)
names(yo)[names(yo)=="year"]<-"length"
new2<-merge(new, yo, by=c("studyid","intid"))
new3<-subset(new2, length>=4)

##Figure 3-  

#EXPERIMENTS
man<-read.csv("manip.csv", header=TRUE, na.strings="NA", as.is=TRUE)
man$spp<-with(man, paste(genus, species, sep=" "))

man_day<-subset(man, factormanipulated!="temperature")

int<-intxn
int<-int[order(int$studyid),]
int$intid<-1:nrow(int)

# take mean per shrub type - HMK 003- use spp_repl
stop(print('stopping code here, rest of code looks to be plotting ...'))

##### extra
#just slopes
sums<-unique(total2[,c("intid","studyid","neglatbi","poslatbi","neg_slope","pos_slope")])
sums$diff<-with(sums, neg_slope-pos_slope)
sums2<-subset(sums, neg_slope!="NA")

year1<-with(total2, aggregate(year, by=list(studyid), FUN=min))
names(year1)[1]<-"studyid"; names(year1)[2]<-"minyear"; 
#year_pheno<-merge(year1, obs[,c("studyid","predictvalue","phenovalue")], by=c("studyid","year")
year2<-with(total2, aggregate(year, by=list(studyid), FUN=max))
names(year2)[1]<-"studyid"; names(year2)[2]<-"maxyear"; 
year3<-merge(year1, year2, by="studyid")
year3$range<-with(year3, maxyear-minyear)

year1<-with(obs, aggregate(predictvalue, by=list(studyid), FUN=min))
names(year1)[1]<-"studyid"; names(year1)[2]<-"minyear"; 
year2<-with(obs, aggregate(predictvalue, by=list(studyid), FUN=max))
names(year2)[1]<-"studyid"; names(year2)[2]<-"maxyear"; 
year3<-merge(year1, year2, by="studyid")
year3$range<-with(year3, maxyear-minyear)

sums2<-merge(sums, year3, by="studyid")

obs2<-merge(roughintxn, obs, by="studyid")
obs3<-subset(obs3, predictunits=="year")
mat<-obs3[,c("studyid","predictvalue","phenovalue","role")]; mat<-na.omit(mat)
mat$predictvalue<-as.numeric(mat$predictvalue)
mat$role<-as.factor(mat$role)

model<-lme(phenovalue~predictvalue+role, random=~1|studyid, data=mat)

library(ggplot2)
sub<-subset(obs3, studyid=="HMK004" & pheno=="start")
ggplot(sub, aes(x=predictvalue, y=phenovalue, colour=taxa)) +
geom_point(size=3)+
geom_line(method="lm", size=1)+theme_bw()

# only studies that looked across SPACE not TIME e.g. HMK009
obs2<-subset(obs, studyid=="HMK006") #multiple sites/year
y2<-with(obs2, aggregate(phenovalue, by=list(spp,predictvalue), FUN=mean, na.rm=TRUE))
names(y2)[1]<-"spp"; names(y2)[2]<-"predictvalue_new"; names(y2)[3]<-"phenovalue_new"
obs3<-merge(y2,obs2[,c("spp","studyid","")], by="spp")
obs<-subset(obs, site!="NA")


###########
man<-read.csv("rough_man.csv", header=TRUE, na.strings="NA", as.is=TRUE)
sub<-subset(man, genus=="Lymantria")

sub<-subset(man, studyid=="HMK002" & fitnesscomp!="fecundity")
sub$fitnessvalue<-as.numeric(sub$fitnessvalue)
m1<-lm(fitnessvalue~treatchangeday, data=sub)
