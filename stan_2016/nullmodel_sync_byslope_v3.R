## USE THIS ONE!
#Nov 2016- Look at workflow from Lizzie- 2016Nov18_nullmodel.pdf
## Null model used for synchrony change


## Some notes and updates by Lizzie on 27 Feb 2017 ##
## Note: This code seems to jump in with raw.long.tot already started ... #
# So I made some adjustments here and in syncmodels.R first! ##

# Set working directory: 
if(length(grep("Lizzie", getwd())>0)) {setwd("~/Documents/git/projects/trophsynch/synchrony/stan_2016") 
} else 
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")

rm(list=ls()) 
options(stringsAsFactors = FALSE)


library(ggplot2)
library(rstan)
library(shinystan)
library(grid)
library(nlme)
library(dplyr)
library(ggrepel)
library(reshape)
set_cppo("fast")  # for best running speed
source("/users/kharouba/google drive/UBC/multiplot.R")
#library(reshape)
# library(lme4)

#Step 1- create pre_climate change dataset
rawlong <- read.csv("input/rawlong2.csv", header=TRUE)
source("input/datacleaningmore.R")

#Step 1- create pre_climate change dataset
rawlong.tot2<-unique(rawlong.tot[,c("studyid","species","phenovalue","year","yr1981")]) 

rawlong.tot2$count<-1
pre<-subset(rawlong.tot2, year<=1981)

sss<- aggregate(pre["count"], pre[c("studyid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas

pre_cc<-merge(rawlong.tot2, sss2, by=c("studyid", "species"))
pre_cc<-subset(pre_cc, year<=1981)


#Step 2- Fit stan on each dataset
pre_cc  <- pre_cc[with(pre_cc, order(species, year)),]
N <- nrow(pre_cc)
y <- pre_cc$phenovalue
Nspp <- length(unique(pre_cc$species))
J <- length(unique(pre_cc$species))
species <- as.numeric(as.factor(pre_cc$species))
pre_cc$year2<-with(pre_cc, year-1956)
year <- pre_cc$year2
#16 interactions from pre_cc
#32 species from pre_cc (32 spp all together, 2 spp repeat across intxns)

null.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=6000, chains=4)
print(null.model, pars = c("mu_b", "sigma_y", "a", "b"))
goo <- extract(null.model)
# some very low n_eff, but then there is not much data (see below) *but* the intercept values seem really crazy ... in the thousands; may want to add a prior? 
dim(pre_cc)
dim(rawlong.tot2)
length(unique(pre_cc$studyid))
length(unique(rawlong.tot2$studyid))

#
fas<-extract(null.model)
it1000 <- matrix(0, ncol=3000, nrow=1)
for (i in 3000:6000){ # 3000 iterations?
        it1000[,(i-3000)] <- fas$mu_b[i]
}
mean(rowMeans(it1000, na.rm=TRUE))
sem<-sd(it1000)/sqrt(length(it1000)); sem
#95% confidence intervals of the mean
c(mean(it1000)-2*sem,mean(it1000)+2*sem)


#individual species
summ_studyspp <- unique(pre_cc[,c("studyid", "species")]);
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]
it1000 <- matrix(0, ncol=3000, nrow=Nspp) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    summ_studyspp$model <- goo$b[i,]
   it1000[,(i-3000)] <- summ_studyspp$model
}
#spp.model<-it1000
hist(rowMeans(it1000, na.rm=TRUE)); median(rowMeans(it1000, na.rm=TRUE))
sem<-sd(rowMeans(it1000, na.rm=TRUE))/sqrt(length(rowMeans(it1000, na.rm=TRUE))); sem
#95% confidence intervals of the mean
c(mean(rowMeans(it1000, na.rm=TRUE))-2*sem,mean(rowMeans(it1000, na.rm=TRUE))+2*sem)

# calculate observed sync change ONLY ONCE
specieschar.formodel <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "species")], FUN=length) 
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species"))
specieschar.formodel.sm  <- specieschar.formodel.sm[with(specieschar.formodel.sm, order(species)),]
intid_full <- read.csv("input/raw_oct.csv", header=TRUE)
sub<-subset(intid_full, intid=="194");
sub[,c("spp1")]<-"Diatom4a spp."
intid_full<-rbind(intid_full, sub)
sub<-subset(intid_full, intid=="195");
sub[,c("spp1")]<-"Diatom4b spp."
intid_full<-rbind(intid_full, sub)
sub<-subset(intid_full, intid=="196");
sub[,c("spp1")]<-"Diatom4c spp."
intid_full<-rbind(intid_full, sub)
intid_full<-subset(intid_full, spp1!="Diatom4 spp.")
intid2<-merge(intid_full, specieschar.formodel.sm, by.x=c("studyid", "spp1"), by.y=c("studyid", "species"),)
intid3<-merge(intid_full, specieschar.formodel.sm, by.x=c("studyid", "spp2"), by.y=c("studyid", "species"),)
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]

intid.nodups[11,2]<-"Diatom4a spp."
intid.nodups[12,2]<-"Diatom4b spp."
intid.nodups[13,2]<-"Diatom4c spp."

summ_studyspp <- unique(pre_cc[,c("studyid", "species")]);
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]

#For interactions AND synchrony change
it1000 <- matrix(0, ncol=3000, nrow=length(unique(intid.sm$intid))) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    summ_studyspp$model <- goo$b[i,]
    andtheanswer <- merge(intid.nodups, summ_studyspp, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, summ_studyspp, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
    it1000[,(i-3000)] <- andtheanswer$model.x-andtheanswer$model.y #model.x=spp1
}
synch.model<-it1000
meanchange <- rowMeans(it1000, na.rm=TRUE); mean(meanchange)

nulldata<-melt(it1000[,2001:3000])
names(nulldata)[1]<-"id"; names(nulldata)[2]<-"iteration"; names(nulldata)[3]<-"sync.change"; 

sem<-sd(rowMeans(it1000, na.rm=TRUE))/sqrt(length(rowMeans(it1000, na.rm=TRUE))); sem
#95% confidence intervals of the mean
c(mean(rowMeans(it1000, na.rm=TRUE))-2*sem,mean(rowMeans(it1000, na.rm=TRUE))+2*sem)

#to explore changes
andtheanswer$meanchange<-meanchange
pre_cc_int<-merge(rawlong.tot, sss2, by=c("studyid", "species"))

unis<-unique(pre_cc_int[,c("studyid","intid","species")])
pre_cc2<-merge(pre_cc, unis, by=c("studyid","species"))
ggplot(pre_cc2, aes(x=year, y=phenovalue, colour=factor(species)))+geom_point()+facet_wrap(~intid)+theme(legend.position="false")

******** REPEAT FROM HERE 5X **********
source("source/nullmodel_intxns.R")


# add with simple pair-wise, non-repeating spp, interactions
#total<-rbind(new, newer)

## STAN
total2<-unique(total[,c("species","y_null","year_null")])
total2<-total2[with(total2, order(species)),]
N <- nrow(total2)
y <- total2$y_null
Nspp <- length(unique(total2$species))
J <- length(unique(total2$species))
species <- as.numeric(as.factor(total2$species))
year <- total2$year_null

# new stan
new.model<-stan("stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=12000, chains=4)
print(new.model, pars = c("mu_b", "sigma_y", "a", "b"))

gooey <- extract(new.model)
summ_studyspp  <- unique(total[c("intid","species")])
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]
spponly<-data.frame(unique(summ_studyspp$species)); names(spponly)[1]<-"species" 

ints  <- summ_studyspp[with(summ_studyspp , order(intid, species)),]
ints$count<-1:nrow(ints)
odd <- seq_len(nrow(ints)) %% 2
first<-subset(ints, odd==1); names(first)[2]<-"spp1"
second<-subset(ints, odd==0); names(second)[2]<-"spp2"
#first<-ints[ints$species %in% odd_index,]; names(first)[2]<-"spp1"
#second<-ints[ints$species %in% even_index,]; names(second)[2]<-"spp2"
ints<-merge(first[,1:2], second[,1:2], by=c("intid"))


#For interactions AND synchrony change
it1000 <- matrix(0, ncol=3000, nrow=length(unique(ints$intid))) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    spponly$model <- gooey$b[i,]
    andtheanswer <- merge(ints, spponly, by.x=c("spp1"), by.y=c("species"), all.x=TRUE);
    names(andtheanswer)[4]<-"model.x"
    andtheanswer2 <- merge(ints, spponly, by.x=c("spp2"), by.y=c("species"), all.x=TRUE);
    names(andtheanswer2)[4]<-"model.y"
    it1000[,(i-3000)] <- andtheanswer$model.x-andtheanswer2$model.y #model.x=spp1
}
nots<-data.frame(array(0,c(54,2)))
nots$id<-1:54
nots$meanchange <- rowMeans(it1000, na.rm=TRUE) #mean 
write.csv(nots, "nullmodel_sync2.csv") 

### RESULTS:
*****only create data frame for the first simulation, 
distn<-data.frame(array(0,c(5,7)))
names(distn)[1]<-"id"; names(distn)[2]<-"mean"; names(distn)[3]<-"low_CI"; names(distn)[4]<-"high_CI"; names(distn)[5]<-"magnitude_mean"; names(distn)[6]<-"mag_low_CI"; names(distn)[7]<-"mag_high_CI";
distn$id<-1:5
*********
Now add to it manually:
mistake- replaced magnitude mean of 2nd sim with 3rd



distn[5,"mean"]<-mean(rowMeans(it1000, na.rm=TRUE))

sem<-sd(rowMeans(it1000, na.rm=TRUE))/sqrt(length(rowMeans(it1000, na.rm=TRUE))); sem
#95% confidence intervals of the mean
ci<-c(mean(rowMeans(it1000, na.rm=TRUE))-2*sem,mean(rowMeans(it1000, na.rm=TRUE))+2*sem)

distn[5,"low_CI"]<-ci[1]
distn[5,"high_CI"]<-ci[2]


#MAGNITUDE
distn[5,"magnitude_mean"]<-mean(abs(rowMeans(it1000, na.rm=TRUE)))
sem<-sd(abs(rowMeans(it1000, na.rm=TRUE)))/sqrt(length(abs(rowMeans(it1000, na.rm=TRUE)))); sem
#95% confidence intervals of the mean
ci<-c(mean(abs(rowMeans(it1000, na.rm=TRUE)))-2*sem,mean(abs(rowMeans(it1000, na.rm=TRUE)))+2*sem)

distn[5,"mag_low_CI"]<-ci[1]
distn[5,"mag_high_CI"]<-ci[2]

write.csv(distn, "nullmodel_sync.csv")

#FIGURE 
tog$null<-null[1:54]
text_high <- textGrob("Closer together", gp=gpar(fontsize=13, fontface="bold"))
text_low <- textGrob("Further apart", gp=gpar(fontsize=13, fontface="bold"))
ggplot(tog, aes(x=null))+geom_histogram(binwidth=.5, alpha=.5, position="identity", colour="black")+theme_bw()+xlim(-1.5, 1.4)+annotation_custom(text_high, xmin=1.0, xmax=1.0, ymin=-0.5, ymax=-0.5)+annotation_custom(text_low, xmin=-1.2, xmax=-1.2, ymin=-0.5, ymax=-0.5)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/year")+ annotation_custom(grob = textGrob(label = "a", hjust = 0, gp = gpar(cex = 1.5)), ymin = 15, ymax = 15, xmin = -2, xmax = -2)


nulldata<-melt(it1000[,2001:3000])
names(nulldata)[1]<-"id"; names(nulldata)[2]<-"iteration"; names(nulldata)[3]<-"sync.change"; 

c<-ggplot(nots, aes(x=meanchange))+geom_histogram(colour="black")+theme_bw()+xlim(range(tog$meanchange))+annotation_custom(text_high, xmin=1.0, xmax=1.0, ymin=-0.5, ymax=-0.5)+annotation_custom(text_low, xmin=-1.2, xmax=-1.2, ymin=-0.5, ymax=-0.5)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/year")+ annotation_custom(grob = textGrob(label = "a", hjust = 0, gp = gpar(cex = 1.5)), ymin = 15, ymax = 15, xmin = -2, xmax = -2)



#old

#For interactions AND synchrony change
it1000 <- matrix(0, ncol=3000, nrow=(length(unique(new$species))/2)) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    summ_studyspp$model <- gooey$b[i,]
    first<-summ_studyspp[summ_studyspp$species %in% spp1,]
    sec<-summ_studyspp[summ_studyspp$species %in% spp2,]
    it1000[,(i-3000)] <- first$model-sec$model #model.x=spp1
}
synch.model<-it1000
mean(rowMeans(it1000, na.rm=TRUE))
range(rowMeans(it1000, na.rm=TRUE))

