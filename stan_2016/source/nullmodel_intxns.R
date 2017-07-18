### SIMULATE DATA FOR SPECIES WITH SAME SPP INTERACTION STRUCTURE AS FULL DATASET ####
#### updated April 2017
# Start with 29 interactions (57 species) where BOTH partners are not found in other intxns (n=54species)

new<-data.frame(array(0,c(nrow(rawlong.tot2),4)))
spp1<-seq(1, 58, 2) #spp1<-seq(1, 53, 2)
spp2<-seq(2, 58, 2) #spp2<-seq(2, 54, 2) 
int<-seq(1, 29, 1) #need 27 interactions

rowcount<-1
for(i in 1:54){
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb
nice<-sample(1:31, size=1, replace=FALSE) #choose species from dataset for time series, n=31 species
rands<-subset(pre_cc, speciesid==nice)
year_null<-rands$year2

#for sppA
ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(rands), mean(ypred_a), sample(goo$sigma_y, size=1));

asdf<-rowcount+(nrow(rands)-1)
new[rowcount:asdf,1]<-spp1[i]
new[rowcount:asdf,2]<-int[i]
new[rowcount:asdf,3]<-y_a
new[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#for sppb
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

new[rowcount:asdf,1]<-spp2[i]
new[rowcount:asdf,2]<-int[i]
new[rowcount:asdf,3]<-y_b
new[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)
}
names(new)[1]<-"species"; names(new)[2]<-"intid";names(new)[3]<-"y_null"; names(new)[4]<-"year_null"
new<-na.omit(new)
new<-subset(new, year_null!=0)

# Next with 14 (not 20 because of HMK048) interactions where ONE partner is found in other intxns )
# use intid.nodups to look at interaction structure
#For two interactions, excluded double partner in this count e.g. Accipiter nisus is in 4 interactions but one is double so only 3 created here
#there are 5 interactions where both partner is repeating e.g. Parus2 caeruleus, Accipiter nisus, see next look below

newer<-data.frame(array(0,c(nrow(rawlong.tot2),4)))
repeatspp<-c(rep(59, 5), rep(60,4), rep(61,3), rep(62,2)) #spp b always repeats; spp a always unique
nonrepeatspp<-seq(63, 76, 1) #38
int<-seq(30, 43, 1)
listofsppa_b<-rep(999, max(repeatspp))

rowcount<-1
k<-1
l<-1
for(i in 59:64){ #count for sppa
	yas<-repeatspp==i
	newcount<-which(yas, TRUE)
#only sample spp A ONCE for each interaction
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb
listofsppa_b[i]<-sppa_b
nice<-sample(1:31, size=1, replace=FALSE) #choose species from dataset for time series, n=31
rands<-subset(pre_cc, speciesid==nice)
year_null<-rands$year2

#for sppA
ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(rands), mean(ypred_a), sample(goo$sigma_y, size=1));

asdf<-rowcount+(nrow(rands)-1)
newer[rowcount:asdf,1]<-i
newer[rowcount:asdf,2]<-int[k]
newer[rowcount:asdf,3]<-y_a
newer[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#for FIRST sppb
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newer[rowcount:asdf,1]<-nonrepeatspp[l]
newer[rowcount:asdf,2]<-int[k]
newer[rowcount:asdf,3]<-y_b
newer[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

k<-k+1
nonrepeatspp<-nonrepeatspp+1

#for OTHER sppb
for(j in 2:length(newcount)){
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

#first re-write sppa
asdf<-rowcount+(nrow(rands)-1)
newer[rowcount:asdf,1]<-i
newer[rowcount:asdf,2]<-int[k]
newer[rowcount:asdf,3]<-y_a
newer[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

newer[rowcount:asdf,1]<-nonrepeatspp[l]
newer[rowcount:asdf,2]<-int[k]
newer[rowcount:asdf,3]<-y_b
newer[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)
k<-k+1
nonrepeatspp<-nonrepeatspp+1
}
}

names(newer)[1]<-"species"; names(newer)[2]<-"intid";names(newer)[3]<-"y_null"; names(newer)[4]<-"year_null"
newer<-na.omit(newer)
newer<-subset(newer, year_null!=0)


### FIX 5 LAST INTERACTIONS- interaction by interaction
newest<-data.frame(array(0,c(nrow(rawlong.tot2),4)))
frame<-data.frame(array(0,c(6,2)))
repeatspp<-c(rep(102, 2), rep(103,2), rep(104,2)) #first 3 interactions are competing penguins
repeatspp<-seq(102,104,1)
#nonrepeatspp<-c(103,104,102)
int<-c(50,51,52,51,52,50)
frame[1:6,1]<-repeatspp
frame[1:6,2]<-int
names(frame)[1]<-"species"; names(frame)[2]<-"intid";

names(newest)[1]<-"species"; names(newest)[2]<-"intid";names(newest)[3]<-"y_null"; names(newest)[4]<-"year_null"

rowcount<-1
#First interaction
yas<-subset(frame, intid==50)
	
#only sample spp A ONCE for each interaction
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb
#listofsppa_b[i]<-sppa_b
nice<-sample(1:31, size=1, replace=FALSE) #choose species from dataset for time series, n=31
rands<-subset(pre_cc, speciesid==nice)
year_null<-rands$year2

#for sppA
ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(rands), mean(ypred_a), sample(goo$sigma_y, size=1));

asdf<-rowcount+(nrow(rands)-1)
newest[rowcount:asdf,1]<-yas[1,1]
newest[rowcount:asdf,2]<-yas$intid[1]
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#for FIRST sppb
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-yas[2,1]
newest[rowcount:asdf,2]<-yas$intid[1]
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#Second interaction
yas<-subset(frame, intid==51)
spp1<-subset(newest, species==102)

newest[rowcount:asdf,1]<-spp1$species
newest[rowcount:asdf,2]<-51
newest[rowcount:asdf,3]<-spp1$y_null
newest[rowcount:asdf,4]<-spp1$year_null
rowcount<-rowcount+nrow(spp1)
asdf<-rowcount+(nrow(spp1)-1)


ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(spp1), mean(ypred_a), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-yas[1,1]
newest[rowcount:asdf,2]<-51
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(spp1)
asdf<-rowcount+(nrow(spp1)-1)

#Third interaction
yas<-subset(frame, intid==52)
spp1<-subset(newest, species==103)

newest[rowcount:asdf,1]<-spp1$species
newest[rowcount:asdf,2]<-52
newest[rowcount:asdf,3]<-spp1$y_null
newest[rowcount:asdf,4]<-spp1$year_null
rowcount<-rowcount+nrow(spp1)
asdf<-rowcount+(nrow(spp1)-1)

spp2<-subset(newest, species==104)
newest[rowcount:asdf,1]<-spp2$species
newest[rowcount:asdf,2]<-52
newest[rowcount:asdf,3]<-spp2$y_null
newest[rowcount:asdf,4]<-spp2$year_null
rowcount<-rowcount+nrow(spp2)
asdf<-rowcount+(nrow(spp2)-1)


#For fourth and fifth interaction (53-54), 8 affected interactions because all need same years (110 is Parus2 caeruleus)
frame<-data.frame(array(0,c(16,2)))
repeatspp<-c(rep(105, 4), rep(106,4)) #
nonrepeatspp<-c(107,108,109,110,110,111,112,113)
int<-c(53:60)
frame[1:8,1]<-repeatspp
frame[9:16,1]<-nonrepeatspp
frame[1:8,2]<-int
frame[9:16,2]<-int
names(frame)[1]<-"species"; names(frame)[2]<-"intid";

#only sample spp A (spp105) ONCE for each interaction
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb
nice<-sample(1:31, size=1, replace=FALSE) #choose species from dataset for time series, n=31
rands<-subset(pre_cc, speciesid==nice)
year_null<-rands$year2

#for sppA- intid 105
ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(rands), mean(ypred_a), sample(goo$sigma_y, size=1));

asdf<-rowcount+(nrow(rands)-1)
newest[rowcount:asdf,1]<-105
newest[rowcount:asdf,2]<-53
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-107
newest[rowcount:asdf,2]<-53
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#intid 54
newest[rowcount:asdf,1]<-105
newest[rowcount:asdf,2]<-54
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-108
newest[rowcount:asdf,2]<-54
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#intid55
newest[rowcount:asdf,1]<-105
newest[rowcount:asdf,2]<-55
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-109
newest[rowcount:asdf,2]<-55
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#intid58, new spp (106) but same years as 105
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb

#for sppA-
ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(rands), mean(ypred_a), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-106
newest[rowcount:asdf,2]<-58
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-111
newest[rowcount:asdf,2]<-58
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#intid59
newest[rowcount:asdf,1]<-106
newest[rowcount:asdf,2]<-59
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-112
newest[rowcount:asdf,2]<-59
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#intid60
newest[rowcount:asdf,1]<-106
newest[rowcount:asdf,2]<-60
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-113
newest[rowcount:asdf,2]<-60
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#intid 57 # i.e. parus2 caeruleus
newest[rowcount:asdf,1]<-106
newest[rowcount:asdf,2]<-57
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

newest[rowcount:asdf,1]<-110
newest[rowcount:asdf,2]<-57
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

#intid 56 # i.e. parus2 caeruleus
newest[rowcount:asdf,1]<-105
newest[rowcount:asdf,2]<-56
newest[rowcount:asdf,3]<-y_a
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)

newest[rowcount:asdf,1]<-110
newest[rowcount:asdf,2]<-56
newest[rowcount:asdf,3]<-y_b
newest[rowcount:asdf,4]<-year_null
rowcount<-rowcount+nrow(rands)
asdf<-rowcount+(nrow(rands)-1)


newest<-subset(newest, year_null!=0)
total<-rbind(new, newer,newest)