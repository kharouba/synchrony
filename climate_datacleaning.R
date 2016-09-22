rm(list=ls())
row.names=FALSE

SEE BELOW FOR CLEANING INDIV STUDIES !!!!

clim<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/climate.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
clim<-subset(clim, phenophase!="start" & extra!="10")
clim$envvalue<-as.numeric(clim$envvalue)
clim$envfactor[clim$envfactor=="temperaure"] <- "temperature"
clim$stud_block<-with(clim, paste(studyid,"_",seasonalblock))

#env change
new<-data.frame(array(0, c(nrow(clim), 7)))
names(new)[1] = "studyid"; names(new)[2] = "envfactor"; names(new)[3] = "envunits"; names(new)[4] = "species"; names(new)[5]<-"year"; names(new)[6]<-"envtype"; names(new)[7] = "envchange";
rowcount<-1
Bgroups<-unique(clim$studyid); b<-Bgroups; b<-as.character(b)
for(i in 1:length(b)) { 
yo<-clim[clim$studyid==b[i],]
if(yo$species=="all"){ #if there is only 1 phenophase/study
asdf<-rowcount+(nrow(yo)-1)
new[rowcount:asdf,1]<-yo[,1]
new[rowcount:asdf,2]<-yo[,c("envfactor")]
new[rowcount:asdf,3]<-yo[,c("envunits")]
new[rowcount:asdf,4]<-yo[,c("species")]
new[rowcount:asdf,5]<-yo[,c("year")]
new[rowcount:asdf,6]<-yo[,c("envtype")]
lm<-with(yo, lm(envvalue~year))
new[rowcount:asdf,7]<-summary(lm)$coefficients[2]
rowcount<-rowcount+nrow(yo)
}
if(yo$species!="all"){ #if there is only 1 phenophase/study
Cgroups<-unique(yo$species); c<-Cgroups; c<-as.character(c)
	for(j in 1:length(c)) {
		yo2<-yo[yo$species==c[j],] 
	asdf<-rowcount+(nrow(yo2)-1)
	new[rowcount:asdf,1]<-yo2[,1]
	new[rowcount:asdf,2]<-yo2[,c("envfactor")]
	new[rowcount:asdf,3]<-yo2[,c("envunits")]
	new[rowcount:asdf,4]<-yo2[,c("species")]
	new[rowcount:asdf,5]<-yo2[,c("year")]
	new[rowcount:asdf,6]<-yo2[,c("envtype")]
	lm<-with(yo2, lm(envvalue~year))
	new[rowcount:asdf,7]<-summary(lm)$coefficients[2]
	rowcount<-rowcount+nrow(yo2)
}
}
}
new2<-subset(new, studyid!="0")

write.csv(new2, "/Volumes/Music/UBC/synchrony project/analysis/envchange.csv")


## 

##standardizing IGNORING SITE
## taking avg env factor across site
sites<-subset(clim, studyid=="HMK018" | studyid=="HMK019" | studyid=="HMK023" & site!="tomakomai")
sites<-sites[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")]
nosites<-subset(clim, site=="tomakomai" | studyid!="HMK018" & studyid!="HMK019" & studyid!="HMK023" & stud_block!="HMK029 _ april" & stud_block!="HMK029 _ may" & stud_block!="HMK029 _ june")
nosites<-nosites[,c("studyid","envfactor","envunits","envtype","year","species","envvalue")]

new<-with(sites, aggregate(envvalue, by=list(studyid, year, species), FUN=mean, na.rm=T)) # across sites
names(new)[1]<-"studyid"; names(new)[2]<-"year"; names(new)[3]<-"species"; names(new)[4]<-"envvalue"
new<-new[order(new$studyid),]
new2<-merge(sites[,c("studyid","envfactor","envunits","envtype","year","species")], new, by=c("studyid","year","species"))
new2<-unique(new2[,c("studyid","year","species","envfactor","envunits","envtype","envvalue")])

clim2<-rbind(nosites, new2)
lol<-with(clim2, aggregate(envvalue, by=list(studyid, species,envfactor), FUN=mean, na.rm=T)) # across years
names(lol)[1]<-"studyid"; names(lol)[2]<-"species"; names(lol)[3]<-"envfactor"; names(lol)[4]<-"env_mean"
clim3<-merge(clim2, lol, by=c("studyid","species","envfactor"))

lol<-with(clim2, aggregate(envvalue, by=list(studyid, species, envfactor), FUN=sd, na.rm=T))
names(lol)[1]<-"studyid"; names(lol)[2]<-"species"; names(lol)[3]<-"envfactor"; names(lol)[4]<-"env_sd"
clim4<-merge(clim3, lol, by=c("studyid","species","envfactor"))

clim4$z<-with(clim4, (envvalue-env_mean)/env_sd)
clim5<-clim4[,c("studyid","species","year","envvalue","z", "envfactor")]
clim6<-merge(clim5, clim[,c("studyid","species","year")], by=c("studyid","species","year"))
clim7 <- unique(clim6[,c("studyid","species","year","envvalue","z","envfactor")])

#standardizing CONSIDERING SITE
lol<-with(clim, aggregate(envvalue, by=list(studyid,site,species), FUN=mean, na.rm=T))
names(lol)[1]<-"studyid"; names(lol)[2]<-"site"; names(lol)[3]<-"species"; names(lol)[4]<-"env_mean"
clim2<-merge(clim, lol, by=c("studyid","site","species"))

lol<-with(clim, aggregate(envvalue, by=list(studyid,site,species), FUN=sd, na.rm=T))
names(lol)[1]<-"studyid"; names(lol)[2]<-"site"; names(lol)[3]<-"species"; names(lol)[4]<-"env_sd"
clim3<-merge(clim2, lol, by=c("studyid","site","species"))

clim3$z<-with(clim3, (envvalue-env_mean)/env_sd)
clim4<-clim3[,c("studyid","site","species","year","z")]

## MERGE CLIMATE TO SPECIES PHENO DATA
all<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/spp_phenodata2.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)
int<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)
intxn<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/intxn.csv", header=TRUE, na.strings="<NA>", as.is=TRUE)

#For studies with SAME factors for diff species THEREFORE interaction gets ONE cue
lets<-subset(clim7, species=="all" | studyid=="HMK016" | studyid=="HMK030" | studyid=="HMK034" | studyid=="SET005")
lets2<-merge(lets[,c("studyid","year","envvalue","z")], int, by=c("studyid","year")); lets2<-na.omit(lets2)

#For studies with diff factors for diff species i.e. interactions have multiple cues
sun<-subset(clim7, species!="all")
clim_spp<-merge(sun, all, by=c("studyid","species","year"))

sun2<-merge(clim_spp, int, by=c("studyid","intid","year"))
sun2<-na.omit(sun2); sun2$phenodiff<-as.numeric(sun2$phenodiff)
write.csv(sun2, "/Volumes/Music/UBC/synchrony project/analysis/clim_spp_mult.csv")


mult<-read.csv("/Volumes/Music/UBC/synchrony project/analysis/clim_multcues.csv", header=TRUE, na.strings="NA", as.is=TRUE)
mix<-merge(mult, int[,c("studyid","intid","year","phenodiff")], by=c("studyid","intid","year"))
mix$phenodiff<-as.numeric(mix$phenodiff)
mix$z_spp1<-as.numeric(mix$z_spp1)
mix$z_spp2<-as.numeric(mix$z_spp2)

**** INTID WILL NO LONGER MATCH BECAUSE INTERACTION ADDED (MAY 2014)**
sun3<-subset(mix, studyid=="HMK011"); sun3<-na.omit(sun3); with(sun3, cor(z_spp1, z_spp2))
m1<-lm(abs(phenodiff)~z_spp1, data=sun3, na.action=na.omit); summary(m1)
m2<-lm(abs(phenodiff)~z_spp2, data=sun3, na.action=na.omit); summary(m2)
#spp1 best
#cor between z scores =0.81

sun3<-subset(mix, studyid=="HMK019" & intid=="170"); sun3<-na.omit(sun3); with(sun3, cor(z_spp1, z_spp2))
m1<-lm(abs(phenodiff)~z_spp1, data=sun3, na.action=na.omit); summary(m1)
m2<-lm(abs(phenodiff)~z_spp2, data=sun3, na.action=na.omit); summary(m2)
AIC(m1,m2)
#spp2 but NS
#cor between z scores =0.73

#intid=171

sun3<-subset(mix, studyid=="HMK031" & intid=="185");sun3<-na.omit(sun3); with(sun3, cor(z_spp1, z_spp2))
m1<-lm(abs(phenodiff)~z_spp1, data=sun3, na.action=na.omit); summary(m1)
m2<-lm(abs(phenodiff)~z_spp2, data=sun3, na.action=na.omit); summary(m2)
AIC(m1,m2)
#spp1 but NS
#cor between z scores =-0.76

sun3<-subset(mix, studyid=="HMK031" & intid=="188");sun3<-na.omit(sun3); with(sun3, cor(z_spp1, z_spp2))
m1<-lm(abs(phenodiff)~z_spp1, data=sun3, na.action=na.omit); summary(m1)
su<-subset(sun2, species=="Thermocyclops oithonoides" & year!="1988" & year!="1989" & year!="1990" & year!="2001")
m2<-lm(abs(phenodiff)~z, data=su, na.action=na.omit); summary(m2)
AIC(m1,m2)
#spp1
#cor between z scores =-0.79

sun3<-subset(mix, studyid=="HMK038" & intid=="197"); sun3<-na.omit(sun3); with(sun3, cor(z_spp1, z_spp2))
m1<-lm(abs(phenodiff)~z_spp1, data=sun3, na.action=na.omit); summary(m1)
m2<-lm(abs(phenodiff)~z_spp2, data=sun3, na.action=na.omit); summary(m2)
AIC(m1,m2)
#spp2 but NS
#cor between z scores =0.49

sun3<-subset(mix, studyid=="HMK038" & intid=="198"); sun3<-na.omit(sun3); with(sun3, cor(z_spp1, z_spp2))
m1<-lm(abs(phenodiff)~z_spp1, data=sun3, na.action=na.omit); summary(m1)
m2<-lm(abs(phenodiff)~z_spp2, data=sun3, na.action=na.omit); summary(m2)
AIC(m1,m2)
#spp2

sun3<-subset(mix, studyid=="HMK038" & intid=="199"); sun3<-na.omit(sun3); with(sun3, cor(z_spp1, z_spp2))
m1<-lm(abs(phenodiff)~z_spp1, data=sun3, na.action=na.omit); summary(m1)
m2<-lm(abs(phenodiff)~z_spp2, data=sun3, na.action=na.omit); summary(m2)
AIC(m1,m2)
spp2

sun3<-subset(mix, studyid=="HMK038" & intid=="200"); sun3<-na.omit(sun3);
m1<-lm(abs(phenodiff)~z_spp1, data=sun3, na.action=na.omit); summary(m1)
m2<-lm(abs(phenodiff)~z_spp2, data=sun3, na.action=na.omit); summary(m2)
AIC(m1,m2)
#spp2

lets3<-subset(clim7, species=="Epirrita autumnata" | species=="Daphnia3 spp." | species=="Diatom2 spp." | species=="Glis glis")
lets4<-merge(lets3[,c("studyid","year","envvalue","z")], int, by=c("studyid","year")); lets3<-na.omit(lets3)

tot<-rbind(lets, lets3); tot<-na.omit(tot) # use this one for temp sens
tot<-rbind(lets2, lets4); tot<-na.omit(tot)

## Standardizing phenodiff
lol<-with(tot, aggregate(phenodiff, by=list(studyid, intid), FUN=mean, na.rm=T)) # across years
names(lol)[1]<-"studyid"; names(lol)[2]<-"intid"; names(lol)[3]<-"phenodiff_mean"
tot2<-merge(tot, lol, by=c("studyid","intid"))

lol<-with(tot, aggregate(phenodiff, by=list(studyid, intid), FUN=sd, na.rm=T))
names(lol)[1]<-"studyid"; names(lol)[2]<-"intid"; names(lol)[3]<-"phenodiff_sd"
tot3<-merge(tot2, lol, by=c("studyid","intid"))

tot3$z_pheno<-with(tot3, (phenodiff-phenodiff_mean)/phenodiff_sd)
tot4<-tot3[,c("studyid","intid","year","z","spp1","spp2","interaction","phenodiff","z_pheno")]
write.csv(tot4, "/Volumes/Music/UBC/synchrony project/analysis/int_phenoenv.csv")

----------------------------------------------
INDIVIDUAL STUDIES
---------------------------------------------
#extract only march-may daily temp from HMK011
stud<-read.csv("/Volumes/Music/UBC/synchrony project/data/HMK011_dailytemp.csv", header=TRUE, na.strings="NA", as.is=TRUE)

stud2<-subset(stud, newmonth==3 & DAY>=27)
stud3<-subset(stud, newmonth==4)
stud4<-subset(stud, newmonth==5 & DAY<=6)
stud5<-rbind(stud2, stud3, stud4)

lol<-with(stud5, aggregate(mean.T, by=list(YEAR), FUN=mean, na.rm=T))
names(lol)[1]<-"year"; names(lol)[2]<-"meantemp"
lol$spp<-"willowtit"

stud2<-subset(stud, newmonth==3 & DAY>=22)
stud3<-subset(stud, newmonth==4)
stud4<-subset(stud, newmonth==5 & DAY<=27)
stud5<-rbind(stud2, stud3, stud4)
lol2<-with(stud5, aggregate(mean.T, by=list(YEAR), FUN=mean, na.rm=T))
names(lol2)[1]<-"year"; names(lol2)[2]<-"meantemp"
lol2$spp<-"initialcat"

stud2<-subset(stud, newmonth==3 & DAY>=13)
stud3<-subset(stud, newmonth==4)
stud4<-subset(stud, newmonth==5 & DAY<=25)
stud5<-rbind(stud2, stud3, stud4)
lol3<-with(stud5, aggregate(mean.T, by=list(YEAR), FUN=mean, na.rm=T))
names(lol3)[1]<-"year"; names(lol3)[2]<-"meantemp"
lol3$spp<-"peakcat"

lol4<-rbind(lol,lol2,lol3)

write.csv(lol4, "/Volumes/Music/UBC/synchrony project/data/HMK011_meantemp.csv")

## HMK043- take mean per year and then year-to-year variation
stud<-read.csv("/Volumes/Music/UBC/synchrony project/data/HMK043_temp.csv", header=TRUE, na.strings="NA", as.is=TRUE)
lol1<-with(stud, aggregate(TEMP, by=list(YEAR), FUN=mean, na.rm=T))
write.csv(lol1, "/Volumes/Music/UBC/synchrony project/data/HMK043_temp2.csv")

## HMK019- take mean per year and per site for each species
stud<-read.csv("/Volumes/Music/UBC/synchrony project/data/HMK019_temp.csv", header=TRUE, na.strings="NA", as.is=TRUE)
sub<-subset(stud, seasonal.block=="march" | seasonal.block=="april")
lol1<-with(sub, aggregate(envvalue, by=list(year, site), FUN=mean, na.rm=T))
names(lol1)[1]<-"year"; names(lol1)[2]<-"site"; names(lol1)[3]<-"envvalue"
lol1$species<-"phytoplankton"

sub<-subset(stud, seasonal.block=="april" | seasonal.block=="may")
lol2<-with(sub, aggregate(envvalue, by=list(year, site), FUN=mean, na.rm=T))
names(lol2)[1]<-"year"; names(lol2)[2]<-"site"; names(lol2)[3]<-"envvalue"
lol2$species<-"zooplankton_fish"
new<-rbind(lol1,lol2)
write.csv(new, "/Volumes/Music/UBC/synchrony project/data/HMK019_temp2.csv")

#HMK023- take mean per year
stud<-read.csv("/Volumes/Music/UBC/synchrony project/data/HMK023_tomaktemp.csv", header=TRUE, na.strings="NA", as.is=TRUE)
lol1<-with(stud, aggregate(temp, by=list(year), FUN=mean, na.rm=T))
names(lol1)[1]<-"year"; names(lol1)[2]<-"temp"
write.csv(lol1, "/Volumes/Music/UBC/synchrony project/data/HMK023_tomaktemp2.csv")

#HMK037- take mean for April and across depth, down to 10m per year
stud<-read.csv("/users/kharouba/google drive/UBC/synchrony project/data/HMK037_madison.csv", header=TRUE, na.strings="NA", as.is=TRUE)
names(stud)[7]<-"temp"
lol1<-with(stud, aggregate(temp, by=list(year), FUN=mean, na.rm=T))
names(lol1)[1]<-"year"; names(lol1)[2]<-"temp"
write.csv(lol1, "/users/kharouba/google drive/UBC/synchrony project/data/HMK037_madison2.csv")

