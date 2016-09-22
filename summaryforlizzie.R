################################################################
#### For interactions with same cue, does sensitivity differ?###
################################################################
mom<-read.csv("tempsens_nov2.csv", header=TRUE, na.strings="NA", as.is=TRUE)

## does temp sensitivity differ for interacting species (ignoring ecosystem), *spp has 2 levels- spp1 or spp2, z is the standardized temp value
m1<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", data=mom), na.action=na.omit); summary(m1)
m3<-update(m1, ~.-z:spp); anova(m1,m3)

## does the strength of temp sensitivity differ for interacting species DIFFERENTLy for terrestrial vs. aquatic
m2<-lme(phenovalue~z*spp*terrestrial, random=~1|studyid/intid, method="ML", mom, na.action=na.omit); summary(m2)
m3<-update(m2, ~.-z:spp:terrestrial); anova(m2,m3)

## how do the temp sensitivities of spp1 vs. spp2 compare in aquatic
m1<-lme(phenovalue~z*spp, random=~1|studyid/intid, method="ML", subset(mom,) terrestrial=="aquatic"), na.action=na.omit); summary(m1)
m3<-update(m2, ~.-z:spp); anova(m2,m3)

#####################################
## what is mean temp sensitivity of species?, ind't of interaction, *species is the ## actual scientific name
##############
m2<-lme(phenovalue~z, random=~1|studyid/species, method="ML", mom2, na.action=na.omit); summary(m2)
m3<-update(m2, ~.-z); anova(m2,m3)


########################################
## can we predict shifts in synchrony?
#######################################
# with this loop below which is the 'individual model' version of the group lme model above (the first model listed), there is no sig interaction for most interactions. There are 3/22 interactions with a sig interaction
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


#Prep for overall predictive models model:
#Step 1 get temp sens per species
rowcount<-1
Bgroups<-unique(mom$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0, c(55, 7)))
names(new)[1] = "studyid"; names(new)[2]<-"species"; names(new)[3]<-"spp"; names(new)[4]<-"intid"; names(new)[5]<-"doybytemp"; names(new)[6]<-"tempsens_coef"; names(new)[7]<-"rsquared"
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
new[rowcount+length(c),6]<-summary(m1)$coefficients[4]
new[rowcount+length(c),7]<-summary(m1)$r.squared
rowcount<-rowcount+1
}
}
new<-subset(new, studyid!="0")

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

yum<-merge(no, yo3, by=c("studyid","intid"))

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
m1<-with(sub, lm(phenodiff~year))
yes[i,3]<-summary(m1)$coefficients[2]
yes[i,4]<-summary(m1)$coefficients[8]
rowcount<-rowcount+1
}

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

stud<-read.csv("studies.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
yep<-merge(all2, stud[,c("studyid","terrestrial","big_interaction")], by="studyid")
all2<-yep
all3<-unique(all2[,c("studyid","intid","doybytemp","rsquarediff","proprsquare","medrsquare","minrsquare","phenodiffbyyr","p","tempbyyr","tempp","terrestrial","big_interaction")])
all2<-all3

#Step 3-predictive models
#tempchange
m1<-lme(abs(phenodiffbyyr)~abs(tempbyyr), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1) #& studyid!="175"
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)

#total difference in temp sensitivity (magnitude)
m1<-lme(abs(phenodiffbyyr)~abs(doybytemp), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp)); anova(m1,m3)

# 'ability' of temperature to predict phenology (i.e. R2)
m1<-lme(abs(phenodiffbyyr)~abs(rsquarediff), random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-minrsquare); anova(m1,m3)

m1<-lme(sqrt(abs(phenodiffbyyr))~minrsquare, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-minrsquare); anova(m1,m3)

## aquatic vs. terrestrial
m1<-lme(abs(phenodiffbyyr)~abs(tempbyyr)*terrestrial, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1) #& studyid!="175"
m3<-update(m1, ~.-abs(tempbyyr):terrestrial); anova(m1,m3)

m1<-lme(abs(phenodiffbyyr)~abs(doybytemp)*terrestrial, random=~1|studyid, data=subset(all2, studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp):terrestrial); anova(m1,m3)

#can temp SENSITIVITY predict spp' pheno shift?
son2<-unique(son[c("studyid","species","doybyyr","doybytemp")])
m1<-lme(sqrt(abs(doybyyr))~abs(doybytemp), random=~1|studyid, data=subset(son2, species!="Keratella1 cochlearis" & studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(doybytemp)); anova(m1,m3)

#can temp CHANGE predict spp' pheno shift?
son2<-unique(son[c("studyid","species","doybyyr","tempbyyr")])
m1<-lme(sqrt(abs(doybyyr))~abs(tempbyyr), random=~1|studyid, data=subset(son2, species!="Keratella1 cochlearis" & studyid!="EMW001" & studyid!="HMK035"), method="ML"); summary(m1)
m3<-update(m1, ~.-abs(tempbyyr)); anova(m1,m3)
