rm(list=ls())
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis")
--------------------------------------------
				 ENV DATA
--------------------------------------------
## HAS TEMP CHANGED?
clim<-read.csv("climate2.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
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

std <- function(x) sd(x)/sqrt(length(x))
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


# ONLY FOR YEARS WITH SPECIES DATA
m1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype=="air" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & studyid!="HMK028"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-year); anova(m3,m1)
#water
m1<-lme(envvalue~year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype=="water" & envunits!="doy" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & species!="Keratella2 cochlearis" & species!="Leptodiaptomus ashlandi"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-year); anova(m3,m1)

m1<-lme(envvalue~envtype+year, random=~1|studyid/species, method="ML", data=subset(clim3, envfactor=="temperature" & envtype!="ground" & species!="Thermocyclops oithonoides" & species!="Parus3 major" & species!="Sitta europaea" & species!="Ficedula2 albicollis" & species!="Perca fluviatillis" & envunits!="doy"), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-envtype); anova(m3,m1)

################################################################
#### For interactions with same cue, does sensitivity differ?###
################################################################
so<-read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/spp_phenodata2.csv", header=TRUE, na.strings="NA", as.is=TRUE)
#so<-so[order(so$studyid),]
#bo<-merge(lets[,c("studyid","year","envvalue","z")], so, by=c("studyid","year"));
#bo<-bo[order(bo$studyid, bo$studyid),]

#temp only
HMK029- stratification better predictor than temp # NOW INCLUDED_  nov 2014
HMK038- diff temp cues but z scores are only r=0.49 correlated
HMK035- no correlation in paper and temp measured as doy NOT C
HMK031- temp measured as doy (just for diatom)
SET005- temp measured as score

bo<-merge(tot[,c("studyid","year","envvalue","z")], so, by=c("studyid","year")); #changed march 2015
#clim8<-subset(tot, envfactor=="temperature" | studyid=="HMK016" | studyid=="HMK034");
#bo<-merge(clim8[,c("studyid","year","envvalue","z")], so, by=c("studyid","year"));
#clim8<-clim8[order(clim8$studyid, clim8$species),]
#bo2<-unique(bo[,c("studyid","year","envvalue","z","species","phenovalue","intid","spp")])
#write.csv(clim8, "/users/kharouba/google drive/UBC/synchrony project/analysis/tempsens_nov2.csv")
write.csv(bo, "/users/kharouba/google drive/UBC/synchrony project/analysis/tempsens_nov3.csv")

#***********start here for stats
mom<-read.csv("tempsens_nov3.csv", header=TRUE, na.strings="NA", as.is=TRUE)
mom<-na.omit(mom)
mom$spp<-as.factor(mom$spp)
mom2<-subset(mom, studyid!="HMK033" & studyid!="HMK030" & studyid!="SET005" & studyid!="HMK035" & intid!="182" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188")
mom<-mom2
mom2<-unique(mom[,c("studyid","year","envvalue","z","species","phenovalue","slope","intid","int_type","spp")])
mom<-mom2

#load yo3 from stats_2015.R
yumm<-merge(mom, unique(yo3[,c("studyid","intid","year")]), by=c("studyid","intid","year"))
mom<-yumm

#mom<-read.csv("tempsens.csv", header=TRUE, na.strings="NA", as.is=TRUE)
#mom[mom$intid=="202","intid"]<-"207"
#mom[mom$intid=="203","intid"]<-"208"
#mom[mom$intid=="204","intid"]<-"209"

#study characteristics
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(mom, stud[,c("studyid","short_site","biome","terrestrial","big_interaction")], by="studyid")
mom<-yep

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

## does temp sensitivity differ for interacting species
m1<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m1)
m2<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK035" & intid!="182" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m2)
AIC(m1,m2)
m3<-update(m1, ~.-z:spp); anova(m1,m3)

m2<-lme(phenovalue~z*spp*terrestrial, random=~1|studyid/intid, method="ML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m2)
m3<-update(m2, ~.-z:spp:terrestrial); anova(m2,m3)

m1<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="aquatic"), na.action=na.omit); summary(m1)
m2<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="aquatic"), na.action=na.omit); summary(m2)
AIC(m1,m2)
m3<-update(m2, ~.-z:spp); anova(m2,m3)

m1<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="terrestrial"), na.action=na.omit); summary(m1)
m2<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="terrestrial"), na.action=na.omit); summary(m2)
AIC(m1,m2)
m3<-update(m2, ~.-z:spp); anova(m2,m3)

#extra
mom2<-subset(mom, int_type=="predation" | int_type=="herbivory")
m2<-lme(phenovalue~z*spp*int_type, random=~1|studyid/intid, method="ML", weights=varIdent(form=~1|spp), data=subset(mom2, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188"), na.action=na.omit); summary(m2)

m2<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", weights=varIdent(form=~1|spp), data=subset(mom2, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188" & int_type=="herbivory"), na.action=na.omit); summary(m2)
m3<-update(m2, ~.-z:spp); anova(m2,m3)

m2<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", weights=varIdent(form=~1|spp), data=subset(mom2, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188" & int_type=="predation"), na.action=na.omit); summary(m2)
m3<-update(m2, ~.-z:spp); anova(m2,m3)

####################################
## what is mean temp sensitivity?, ind't of interaction
only include unique spp! i.e. remove intid as random effect
since some species repeated e.g. Diatom 4 spp.
mom2<-subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="182" & intid!="184" & intid!="186" & intid!="187" & intid!="188")
mom2<-unique(mom2[c("studyid","species","year","z","phenovalue","biome","terrestrial","big_interaction")])
m2<-lme(phenovalue~z, random=~1|studyid/species, method="ML", mom2, na.action=na.omit); summary(m2)
m3<-update(m2, ~.-z); anova(m2,m3)
#remove HMK038
m2<-lme(phenovalue~z+spp, random=~1|studyid/intid, method="REML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK038"), na.action=na.omit); summary(m2)

m1<-lme(phenovalue~z*terrestrial, random=~1|studyid/species, method="ML", mom2, na.action=na.omit); summary(m1)
m2<-lme(phenovalue~z*terrestrial, random=~1|studyid/species, method="ML", weights=varIdent(form=~1|terrestrial), mom2, na.action=na.omit); summary(m2)
AIC(m1,m2)
m3<-update(m2, ~.-z:terrestrial); anova(m2,m3)


ggplot(mom2, aes(x=z, y=phenovalue, colour=factor(terrestrial))) +  geom_point(aes(shape=factor(terrestrial))) + geom_smooth(method="lm", se=FALSE, aes(fill = factor(terrestrial)))+theme_bw()

#extra
m2<-lme(phenovalue~spp+z, random=~1|studyid/species, method="ML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="terrestrial"), na.action=na.omit); summary(m2)
#study HMK018 outlier becasuse really late doy

m2<-lme(phenovalue~spp+z, random=~1|studyid/species, method="ML", weights=varIdent(form=~1|spp), data=subset(mom, studyid!="HMK035" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188" & terrestrial=="aquatic"), na.action=na.omit); summary(m2)

########################################
## can we predict shifts in synchrony?
#######################################

## do interacting species have different temperature sensitivities?
Bgroups<-unique(mom$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0, c(25, 3)))
names(new)[1] = "studyid"; names(new)[2]<-"intid"; names(new)[3]<-"p"; 
for(i in 1:length(b)) { 
yo<-mom[mom$intid==b[i],]
new[i,1]<-yo[1,c("studyid")]
new[i,2]<-yo[1,c("intid")]
m1<-with(yo, lm(phenovalue~z*spp))
new[i,3]<-summary(m1)$coefficients[16]
}
new<-subset(new, studyid!="0")
sig<-subset(new, p<=0.05)

new$sig<-c(1,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,1,1,0,0,0)
stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(new, stud[,c("studyid","terrestrial","big_interaction")], by="studyid")

m4<-glmer(sig~terrestrial+(1|studyid), family=binomial, yep, na.action=na.omit); summary(m4)
m2<-update(m4, ~.-terrestrial); anova(m4,m2)

#to get species' temp sensitivity
Bgroups<-unique(mom$species); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0, c(25, 5)))
names(new)[1] = "studyid"; names(new)[2]<-"intid"; names(new)[3]<-"species"; names(new)[4]<-"tempsens"; names(new)[5]<-"p" 
for(i in 1:length(b)) { 
yo<-mom[mom$species==b[i],]
new[i,1]<-yo[1,c("studyid")]
new[i,2]<-yo[1,c("intid")]
new[i,3]<-b[i]
m1<-with(yo, lm(phenovalue~z))
new[i,4]<-summary(m1)$coefficients[2]
new[i,5]<-summary(m1)$coefficients[8]
}

m1<-lme(tempsens~1, random=~1|studyid, method="ML", subset(new, species!="Diatom2 spp."), na.action=na.omit); summary(m1) #only species with temperature as main cue
m4<-update(m1, ~.-1); anova(m1,m4)

stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(new, stud[,c("studyid","terrestrial","big_interaction")], by="studyid")

m1<-lme(tempsens~terrestrial, random=~1|studyid, method="ML", subset(yep, species!="Diatom2 spp."), na.action=na.omit); summary(m1) #only species with temperature as main cue
m4<-update(m1, ~.-terrestrial); anova(m1,m4)

#Prep for overall model:
#Step 1 get temp sens per species
rowcount<-1
Bgroups<-unique(mom$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0, c(55, 8)))
names(new)[1] = "studyid"; names(new)[2]<-"species"; names(new)[3]<-"spp"; names(new)[4]<-"intid"; names(new)[5]<-"doybytemp"; names(new)[6]<-"tempsens_coef"; names(new)[7]<-"rsquared"; names(new)[8]<-"quad_p"
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
m3<-lm(phenovalue ~ poly(z,2), na.action=na.omit, data=spp) 
new[rowcount+length(c),5]<-summary(m1)$coefficients[2]
new[rowcount+length(c),6]<-summary(m1)$coefficients[4]
new[rowcount+length(c),7]<-summary(m1)$r.squared
new[rowcount+length(c),8]<-summary(m3)$coefficients[12]
rowcount<-rowcount+1
}
}
new<-subset(new, studyid!="0")
nl<-subset(new, quad_p<=0.06)


ggplot(data=subset(mom, intid=="1"), aes(z, phenovalue, colour=factor(species)))+geom_point()+geom_smooth(method="lm", formula=y ~ poly(x, 2), se=FALSE,size=0.5, aes(colour=factor(species))) + theme_bw()

Bgroups<-unique(new$intid); b<-Bgroups; b<-as.character(b)
no<-data.frame(array(0, c(length(b), 7)))
names(no)[1] = "studyid"; names(no)[2]<-"intid"; names(no)[3]<-"doybytemp";  names(no)[4]<-"rsquarediff"; names(no)[5]<-"proprsquare"; names(no)[6]<-"medrsquare"; names(no)[7]<-"minrsquare"
rowcount<-1
for(i in 1:length(b)) { 
yo<-new[new$intid==b[i],]
asdf<-(rowcount+nrow(yo))-1
no[rowcount:asdf,1]<-yo[1,c("studyid")]
no[rowcount:asdf,2]<-yo[1,c("intid")]
spp1<-subset(yo, spp=="1")
spp2<-subset(yo, spp=="2")
no[rowcount:asdf,3]<-spp1[1,5]-spp2[1,5]
no[rowcount:asdf,4]<-spp1[1,7]-spp2[1,7]
no[rowcount:asdf,6]<-median(c(spp1[1,7],spp2[1,7]))
no[rowcount:asdf,7]<-min(c(spp1[1,7],spp2[1,7]))
if(spp1[1,7]>spp2[1,7]){
no[rowcount:asdf,5]<-spp2[1,7]/spp1[1,7]	
}
if(spp2[1,7]>spp1[1,7]){
no[rowcount:asdf,5]<-spp1[1,7]/spp2[1,7]	
}
rowcount<-rowcount+nrow(yo)
}
no2<-unique(no[,c("studyid","intid","doybytemp","rsquarediff","proprsquare","medrsquare","minrsquare")])


NEED TO SQUARE ROOT
yum<-merge(no, yo3, by=c("studyid","intid"))
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
names(sup)[1] = "studyid"; names(sup)[2]<-"species"; names(sup)[3]<-"spp"; names(sup)[4]<-"intid"; names(sup)[5]<-"tempbyyr"; 
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
m3<-lm(z ~ poly(year,2), na.action=na.omit, data=spp) 
sup[rowcount+length(c),5]<-summary(m1)$coefficients[2]
#sup[rowcount+length(c),6]<-summary(m3)$coefficients[12]
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
m1<-with(sub, lm(phenodiff~year))
yes[i,3]<-summary(m1)$coefficients[2]
yes[i,4]<-summary(m1)$coefficients[8]
rowcount<-rowcount+1
}

lala<-merge(yes, unique(yo3[,c("studyid","intid","phenofreq","length")]), by=c("studyid","intid"))
lala2<-subset(lala, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")
sig<-subset(lala2, p<0.07)
pos<-subset(lala2, p<0.07 & phenodiffbyyr>0)

#get temp change per interaction
env<-merge(unique(yo3[,c("studyid","intid","year")]), unique(mom[,c("studyid","intid","year","envvalue","z")]), by=c("studyid","intid","year"))
Bgroups<-unique(env$intid); b<-Bgroups; b<-as.character(b)
soop<-data.frame(array(0, c(length(b), 5)))
names(soop)[1] = "studyid"; names(soop)[2]<-"intid"; names(soop)[3]<-"tempbyyr"; names(soop)[4]<-"tempp"; names(soop)[5]<-"quad_p";
for(i in 1:length(b)) { 
sub<-env[env$intid==b[i],]
soop[i,1]<-sub[1,c("studyid")]
soop[i,2]<-sub[1,c("intid")]
#m1<-with(sub, lm(envvalue~year))
m1<-with(sub, lm(z~year))
m3<-lm(z ~ poly(year,2), na.action=na.omit, data=sub) 
soop[i,3]<-summary(m1)$coefficients[2]
soop[i,4]<-summary(m1)$coefficients[8]
soop[i,5]<-summary(m3)$coefficients[12]
rowcount<-rowcount+1
}

nl<-subset(soop, quad_p<=0.06)
nl<-subset(soop, tempp<=0.06)

ggplot(data=subset(mom, intid=="191"), aes(year,envvalue))+geom_point()+geom_smooth(method="lm", formula=y ~ poly(x, 2), se=FALSE,size=0.5) + theme_bw()

all<-merge(no, yes, by=c("studyid","intid"))
all2<-merge(all, soop, by=c("studyid","intid"))

stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(all2, stud[,c("studyid","terrestrial","big_interaction")], by="studyid")
all2<-yep
all3<-unique(all2[,c("studyid","intid","doybytemp","rsquarediff","proprsquare","medrsquare","minrsquare","phenodiffbyyr","p","tempbyyr","tempp","terrestrial","big_interaction")])
all2<-all3

#to get temp change per spp
sup2<-merge(sup, all2[,c("studyid","intid","terrestrial")])
io<-subset(sup2, terrestrial=="aquatic")
io2<-unique(io[,c("studyid","species","tempbyyr")])
with(subset(io2, studyid!="EMW001" & studyid!="HMK035"), mean(tempbyyr))

# get z per year
#mom[mom$intid=="202","intid"]<-"207"
#mom[mom$intid=="203","intid"]<-"208"
#mom[mom$intid=="204","intid"]<-"209"
yo3$intid<-as.factor(yo3$intid)
#solo$intid<-as.factor(solo$intid)
solo<-merge(yo3, mom[,c("studyid","year","intid","envvalue","z")], by=c("studyid","year","intid"))
solo2<-unique(solo[,c("studyid","year","intid","interaction","phenodiff","phenodiff_base","length","phenofreq","envvalue","z")])
solo2$phenofreq<-as.factor(solo2$phenofreq)
m1<-lme(phenodiff~z, random=~1|studyid/intid, method="ML", data=subset(solo2, length>6 & studyid!="EMW001" & studyid!="HMK035" & studyid!="175")); summary(m1)
m3<-update(m1, ~.-z); anova(m1,m3)

##tempchange##
#based on envvalue
m1<-lme(tempbyyr~1, random=~1|studyid, method="ML", subset(soop, studyid!="HMK028"), na.action=na.omit); summary(m1) #excluding HMK028 because sum of temp and want to get C/year
m1<-lme(tempbyyr~1, random=~1|studyid, method="ML", subset(soop, studyid!="HMK028" & intid!="170" & intid!="171" & intid!="180" & intid!="191" & intid!="155"), na.action=na.omit); summary(m1) #without non-linear
m4<-update(m1, ~.-1); anova(m1,m4)
m1<-lme(tempbyyr~terrestrial, random=~1|studyid, method="ML", subset(all2, studyid!="HMK028"), na.action=na.omit); summary(m1)
m1<-lme(tempbyyr~terrestrial, random=~1|studyid, method="ML", subset(all2, studyid!="HMK028" & intid!="170" & intid!="171" & intid!="180" & intid!="191" & intid!="155"), na.action=na.omit); summary(m1) #without non-linear
m4<-update(m1, ~.-terrestrial); anova(m1,m4)

#based on z
m1<-lme(abs(phenodiffbyyr)~abs(tempbyyr), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1) #& studyid!="175"
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)

#without non-linearities
m1<-lme(abs(phenodiffbyyr)~abs(tempbyyr), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035" & intid!="155" & intid!="170" & intid!="171" & intid!="180" & intid!="191" & intid!="197"), method="ML"); summary(m1) #& studyid!="175"
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)
sig<-subset(all2, p<0.07)

#only sig
m1<-lme(abs(phenodiffbyyr)~abs(tempbyyr), random=~1|studyid, data=subset(all2, p<0.07 & studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1) #& studyid!="175"
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)
sig<-subset(all2, p<0.07)

##total difference in temp sensitivity (magnitude + direction)##
m1<-lme(abs(phenodiffbyyr)~abs(doybytemp), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp)); anova(m1,m3)

#without non-linear
m1<-lme(abs(phenodiffbyyr)~abs(doybytemp), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035" & intid!="165" & intid!="180" & intid!="194" & intid!="197" & intid!="155"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp)); anova(m1,m3)

#only sig
m1<-lme(abs(phenodiffbyyr)~abs(doybytemp), random=~1|studyid, data=subset(all2, p<0.07 & studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp)); anova(m1,m3)

m1<-lme(abs(phenodiffbyyr)~abs(rsquarediff), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(rsquarediff)); anova(m1,m3)

#without non-linear
m1<-lme(abs(phenodiffbyyr)~abs(rsquarediff), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035" & intid!="165" & intid!="180" & intid!="194" & intid!="197" & intid!="155"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(rsquarediff)); anova(m1,m3)

m1<-lme(sqrt(abs(phenodiffbyyr))~minrsquare, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-minrsquare); anova(m1,m3)

## aquatic vs. terrestrial
m1<-lme(abs(phenodiffbyyr)~abs(tempbyyr)*terrestrial, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1) #& studyid!="175"
m3<-update(m1, ~.-abs(tempbyyr):terrestrial); anova(m1,m3)

m1<-lme(abs(phenodiffbyyr)~abs(doybytemp)*terrestrial, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp):terrestrial); anova(m1,m3)

#can temp sens predict spp' pheno shift?
son2<-unique(son[c("studyid","species","doybyyr","doybytemp")])
m1<-lme(sqrt(abs(doybyyr))~abs(doybytemp), random=~1|studyid, data=subset(son2, species!="Keratella1 cochlearis" & studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp)); anova(m1,m3)

#can temp change predict spp' pheno shift?
son2<-unique(son[c("studyid","species","doybyyr","tempbyyr")])
m1<-lme(sqrt(abs(doybyyr))~abs(tempbyyr), random=~1|studyid, data=subset(son2, species!="Keratella1 cochlearis" & studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)

m1<-lme(rsquarediff~1, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-1); anova(m1,m3)

m1<-lme(proprsquare~1, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-1); anova(m1,m3)

m1<-lme(sqrt(abs(phenodiffbyyr))~rsquarediff, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-1); anova(m1,m3)

#just those with C
m1<-lme(sqrt(abs(phenodiff_base))~phenofreq+abs(envvalue), random=~1|studyid/intid, weights=varPower(form=~envvalue), data=subset(solo2, length>6 & studyid!="EMW001" & studyid!="HMK035" & studyid!="175" & studyid!="HMK028")); summary(m1)

#overall z change
solo2$z<-as.numeric(solo2$z)
m1<-lme(z~year, random=~1|studyid/intid,data=subset(solo2, length>6 & studyid!="EMW001" & studyid!="HMK035" & studyid!="175")); summary(m1)

ggplot(data=all2, aes(abs(tempbyyr), abs(phenodiffbyyr))) +  geom_point(size=3) +
geom_smooth(method="lm", se=FALSE)+theme_bw()

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