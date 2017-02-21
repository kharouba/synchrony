#Nov 2016- Look at workflow from Lizzie- 2016Nov18_nullmodel.pdf
## Null model used for synchrony change

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#Step 1- create pre_climate change dataset
#load rawlong.tot
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
year <- pre_cc$year
#16 interactions from pre_cc
#32 species from pre_cc (32 spp all together, 2 spp repeat across intxns)

null.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=12000, chains=4)
print(null.model)
goo <- extract(null.model)

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

#individaul species
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

sem<-sd(rowMeans(it1000, na.rm=TRUE))/sqrt(length(rowMeans(it1000, na.rm=TRUE))); sem
#95% confidence intervals of the mean
c(mean(rowMeans(it1000, na.rm=TRUE))-2*sem,mean(rowMeans(it1000, na.rm=TRUE))+2*sem)

#to explore changes
andtheanswer$meanchange<-meanchange
pre_cc_int<-merge(rawlong.tot, sss2, by=c("studyid", "species"))

unis<-unique(pre_cc_int[,c("studyid","intid","species")])
pre_cc2<-merge(pre_cc, unis, by=c("studyid","species"))
ggplot(pre_cc2, aes(x=year, y=phenovalue, colour=factor(species)))+geom_point()+facet_wrap(~intid)+theme(legend.position="false")


### SIMULATE DATA FOR SPECIES WITH SAME SPP INTERACTION STRUCTURE AS FULL DATASET ####
#### NEW #### Dec 7 2016
# Start with 28 interactions where BOTH partners are not found in other intxns (n=54species)
new<-data.frame(array(0,c(nrow(rawlong.tot2),4)))
spp1<-seq(39, 94, 2) #spp1<-seq(1, 53, 2)
spp2<-seq(40, 95, 2) #spp2<-seq(2, 54, 2) 
int<-seq(1, 27, 1) #need 27 interactions

rowcount<-1
for(i in 1:54){
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb
nice<-sample(1:31, size=1, replace=FALSE) #choose species from dataset for time series, n=31 species
rands<-subset(pre_cc, speciesid==nice)
year_null<-rands$year

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

# Next with 22 interactions where at least ONE partner is found in other intxns (n=32 species (there are 21 species whose partner is repeating but not it e.g. Sitta europaea, there are 11 species that repeat e.g. Accipiter nisus))
#there are 5 interactions where both partner is repeating e.g. Parus2 caeruleus, Accipiter nisus
# for each repeat, same time series, FIX AFTER
### LOOP WILL HAVE ERROR BUT OK, PROCEED!!! (Error in (1 + (newnonrepeatspp[1] - nonrepeatspp[1])):21 : NA/NaN argument)
newer<-data.frame(array(0,c(nrow(rawlong.tot2),4)))
repeatspp<-c(rep(1, 5), rep(2,4), rep(3,4),rep(4,4),rep(5,3), rep(6,3), rep(7,2), rep(8,2), rep(9,2), rep(10,2), rep(11,2))
nonrepeatspp<-seq(12, 38, 1) #38
int<-seq(29, 61, 1)
listofsppa_b<-1:length(unique(repeatspp))

rowcount<-1
k<-1
l<-1
for(i in 1:length(unique(repeatspp))){ #count for sppa
	yas<-repeatspp==i
	newcount<-which(yas, TRUE)
#only sample spp A ONCE for each interaction
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb
listofsppa_b[i]<-sppa_b
nice<-sample(1:31, size=1, replace=FALSE) #choose species from dataset for time series, n=31
rands<-subset(pre_cc, speciesid==nice)
year_null<-rands$year

#for sppA
ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(rands), mean(ypred_a), sample(goo$sigma_y, size=1));

asdf<-rowcount+(nrow(rands)-1)
newer[rowcount:asdf,1]<-unique(repeatspp)[i]
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

int<-int+1
nonrepeatspp<-nonrepeatspp+1

#for OTHER sppb
for(j in 2:length(newcount)){
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

#first re-write sppa
asdf<-rowcount+(nrow(rands)-1)
newer[rowcount:asdf,1]<-unique(repeatspp)[i]
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
int<-int+1
nonrepeatspp<-nonrepeatspp+1
}
newint<-int
newnonrepeatspp<-nonrepeatspp

nonrepeatspp<-seq(12, 38, 1)
int<-seq(28, 61, 1)

int<-int[(1+(newint[1]-int[1])):33]
nonrepeatspp<-nonrepeatspp[(1+(newnonrepeatspp[1]-nonrepeatspp[1])):21]
}
names(newer)[1]<-"species"; names(newer)[2]<-"intid";names(newer)[3]<-"y_null"; names(newer)[4]<-"year_null"
newer<-na.omit(newer)
newer<-subset(newer, year_null!=0)


### FIX 5 LAST INTERACTIONS
#need to fix spp 7,8,9 create 9,10,11
newer<-subset(newer, intid!="34")
newer<-subset(newer, intid!="38")
newer<-subset(newer, intid!="39")
newer<-subset(newer, intid!="40")
newer<-subset(newer, intid!="41")
#newer<-subset(newer, intid!="51")
newer<-subset(newer, intid!="52")
newer<-subset(newer, intid!="53")
#newer<-subset(newer, intid!="55")
newer<-subset(newer, intid!="56")
newer<-subset(newer, intid!="57")
#newer<-subset(newer, intid!="58")

#newer[newer$species==17,c("species")]<-7

#for spp7 (as sppb first, for intxn with spp2, intid=="34")
rands<-subset(newer, species=="2" & intid=="35")
sppa_b<-listofsppa_b[2] #need to pull sppa_b from when first did spp2
year_null<-rands$year_null

spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

sub<-data.frame(array(0,c(nrow(rands),4)))
names(sub)[1]<-"species"; names(sub)[2]<-"intid";names(sub)[3]<-"y_null"; names(sub)[4]<-"year_null"
sub[,1]<-"7"
sub[,2]<-"34"
sub[,3]<-y_b
sub[,4]<-year_null
newer<-rbind(newer, sub)

sub2<-subset(newer, species=="2" & intid=="35")
sub2[,2]<-"34"
newer<-rbind(newer, sub2)


#for spp3, same years as spp2 because from same study
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb
rands<-subset(newer, species=="2" & intid=="35")
year_null<-rands$year
ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(rands), mean(ypred_a), sample(goo$sigma_y, size=1));

sub<-data.frame(array(0,c(nrow(rands),4)))
names(sub)[1]<-"species"; names(sub)[2]<-"intid";names(sub)[3]<-"y_null"; names(sub)[4]<-"year_null"
sub[,1]<-"3"
sub[,2]<-"38"
sub[,3]<-y_b
sub[,4]<-year_null
newer<-rbind(newer, sub)

#for spp7 (as sppb, for intxn with spp3, intid=="38"); to add spp7 again
sub2<-subset(newer, species=="7")
sub2[,2]<-"38"
newer<-rbind(newer, sub2)

# for spp22, (as sppb for intxn with spp3, intxn=39)
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

sub<-data.frame(array(0,c(nrow(rands),4)))
names(sub)[1]<-"species"; names(sub)[2]<-"intid";names(sub)[3]<-"y_null"; names(sub)[4]<-"year_null"
sub[,1]<-"22"
sub[,2]<-"39"
sub[,3]<-y_b
sub[,4]<-year_null
newer<-rbind(newer, sub)
#to add spp3 again so that spp22 has partner
sub2<-subset(newer, species=="3" & intid=="38")
sub2[,2]<-"39"
newer<-rbind(newer, sub2)

# for spp23, (as sppb for intxn with spp3, intxn=40)
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

sub<-data.frame(array(0,c(nrow(rands),4)))
names(sub)[1]<-"species"; names(sub)[2]<-"intid";names(sub)[3]<-"y_null"; names(sub)[4]<-"year_null"
sub[,1]<-"23"
sub[,2]<-"40"
sub[,3]<-y_b
sub[,4]<-year_null
newer<-rbind(newer, sub)

sub2<-subset(newer, species=="3" & intid=="38")
sub2[,2]<-"40"
newer<-rbind(newer, sub2)

# for spp24, (as sppb for intxn with spp3, intxn=41)
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

sub<-data.frame(array(0,c(nrow(rands),4)))
names(sub)[1]<-"species"; names(sub)[2]<-"intid";names(sub)[3]<-"y_null"; names(sub)[4]<-"year_null"
sub[,1]<-"24"
sub[,2]<-"41"
sub[,3]<-y_b
sub[,4]<-year_null
newer<-rbind(newer, sub)

sub2<-subset(newer, species=="3" & intid=="38")
sub2[,2]<-"41"
newer<-rbind(newer, sub2)

#for spp 9 (as spp a, for intxn with spp10, intxn 56)
sppa_b<-sample(goo$mu_b, size=1) # need value for sppb
nice<-sample(1:31, size=1, replace=FALSE) #choose species from dataset for time series, n=31 species
rands<-subset(pre_cc, speciesid==nice)
year_null<-rands$year
ypred_a<-(sample(goo$a, size=1))+(sppa_b*year_null)
y_a <- rnorm(nrow(rands), mean(ypred_a), sample(goo$sigma_y, size=1));

sub<-data.frame(array(0,c(nrow(rands),4)))
names(sub)[1]<-"species"; names(sub)[2]<-"intid";names(sub)[3]<-"y_null"; names(sub)[4]<-"year_null"
sub[,1]<-"9"
sub[,2]<-"56"
sub[,3]<-y_a
sub[,4]<-year_null
newer<-rbind(newer, sub)

#for spp 10 (as sppb, for intxn with spp9, intxn 56)
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

sub<-data.frame(array(0,c(nrow(rands),4)))
names(sub)[1]<-"species"; names(sub)[2]<-"intid";names(sub)[3]<-"y_null"; names(sub)[4]<-"year_null"
sub[,1]<-"10"
sub[,2]<-"56"
sub[,3]<-y_b
sub[,4]<-year_null
newer<-rbind(newer, sub)

#for spp 11 (as sppb, for intxn with spp9, intxn 57)
spp_b<-sppa_b-(sample(meanchange, size=1)) # calculate b for sppb from spp a and sync change
ypred_b<-(sample(goo$a, size=1))+(spp_b*year_null)
y_b <- rnorm(nrow(rands), mean(ypred_b), sample(goo$sigma_y, size=1));

sub<-data.frame(array(0,c(nrow(rands),4)))
names(sub)[1]<-"species"; names(sub)[2]<-"intid";names(sub)[3]<-"y_null"; names(sub)[4]<-"year_null"
sub[,1]<-"11"
sub[,2]<-"57"
sub[,3]<-y_b
sub[,4]<-year_null
newer<-rbind(newer, sub)

sub2<-subset(newer, species=="9" & intid=="56")
sub2[,2]<-"57"
newer<-rbind(newer, sub2)

#for intxn 58
sub<-subset(newer, species=="10")
sub2<-subset(newer, species=="11")
sub[,2]<-"58"
sub2[,2]<-"58"
newer<-rbind(newer, sub, sub2)

# add with simple pair-wise, non-repeating spp, interactions
total<-rbind(new, newer)

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
new.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=8000, chains=4)

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
mos<-mean(rowMeans(it1000, na.rm=TRUE)); null<-rowMeans(it1000, na.rm=TRUE)

sem<-sd(rowMeans(it1000, na.rm=TRUE))/sqrt(length(rowMeans(it1000, na.rm=TRUE))); sem
#95% confidence intervals of the mean
c(mean(rowMeans(it1000, na.rm=TRUE))-2*sem,mean(rowMeans(it1000, na.rm=TRUE))+2*sem)


#MAGNITUDE
hist(abs(rowMeans(it1000, na.rm=TRUE)))
median(abs(rowMeans(it1000, na.rm=TRUE)))

sem<-sd(abs(rowMeans(it1000, na.rm=TRUE)))/sqrt(length(abs(rowMeans(it1000, na.rm=TRUE)))); sem
#95% confidence intervals of the mean
c(mean(abs(rowMeans(it1000, na.rm=TRUE)))-2*sem,mean(abs(rowMeans(it1000, na.rm=TRUE)))+2*sem)

#FIGURE 
tog$null<-null[1:54]
text_high <- textGrob("Closer together", gp=gpar(fontsize=13, fontface="bold"))
text_low <- textGrob("Further apart", gp=gpar(fontsize=13, fontface="bold"))
ggplot(tog, aes(x=null))+geom_histogram(binwidth=.5, alpha=.5, position="identity", colour="black")+theme_bw()+xlim(-1.5, 1.4)+annotation_custom(text_high, xmin=1.0, xmax=1.0, ymin=-0.5, ymax=-0.5)+annotation_custom(text_low, xmin=-1.2, xmax=-1.2, ymin=-0.5, ymax=-0.5)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/year")+ annotation_custom(grob = textGrob(label = "a", hjust = 0, gp = gpar(cex = 1.5)), ymin = 15, ymax = 15, xmin = -2, xmax = -2)



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

