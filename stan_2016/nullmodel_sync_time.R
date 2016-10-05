#Aug 2016- Look at workflow from Lizzie- 2016Aug22_nullmodel.pdf
## Null model used for synchrony change

#Step 1- create pre_climate change dataset
#load rawlong.tot
interact <- read.csv("input/raw_oct.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue) #spp1-spp2


rawlong.tot2$count<-1
pre<-subset(rawlong.tot2, year<=1981)
pres<-merge(pre, unique(rawlong.tot[,c("studyid","intid","species")]), by=c("studyid","species"))
sss<- aggregate(pres["count"], pres[c("studyid", "intid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas

pre_cc<-merge(pres[,c(1:4,6)],  sss2, by=c("studyid", "intid", "species"))

#Step 2- Create distribution of means and sd (2)
means <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "intid")], FUN=mean)
names(means)[3]<-"mean_doy" #mean of each dataset
sds <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "intid")], FUN=sd)
names(sds)[3]<-"sd_doy" #sd of each dataset

# For each length of time series

timeseries<-c(5,10,15,20,25,30,35,40)
rowcount2<-1
asdf2<-rowcount2+(54-1) #number of interactions
output<-data.frame(array(0,c(8*54,4)))

for(le in 1:length(timeseries)){ #length(timeseries)


# FOR EACH INTERACTION
Bgroups<-unique(rawlong.tot$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0,c((timeseries[le]*2)*54,5)))
rowcount<-1
for(i in 1:length(b)){
int_long<-subset(rawlong.tot, intid==b[i])
asdf<-rowcount+((timeseries[le]*2)-1) #length of NEW time series*2 species

#Step 3- Create new distribution for spp1 from sample mean and sd
mu_spp1<-sample(means$mean_doy, size=1) #draw one mean from means
sigma_spp1<-sample(sds$sd_doy, size=1) #draw one sd from distribution of SDs
a<-rnorm(1000, mu_spp1, sigma_spp1); #create distribution of phenovalues based on moean and sd


#Step 3- Create new time-series for spp1 of X length based on normal distribution
spp1_data<-subset(int_long, spp=="spp1")
time<-spp1_data$yr1981
time2<-rnorm(1000, mean(time), sd(time)) # create distribution of years based on observed
spp1_datab<-data.frame(array(0,c(timeseries[le],3))) #make new dataframe based on new length
spp1_datab[,1]<-rep(unique(spp1_data$studyid), timeseries[le])
spp1_datab[,2]<-rep(unique(spp1_data$intid), timeseries[le])
spp1_datab[,3]<-rep(unique(spp1_data$species), timeseries[le])

spp1_datab$ynull<-sample(a, size=timeseries[le]); #start with a 5 year time series
spp1_datab$yr1981<-sample(time2, size=timeseries[le]); #sample X number of years based on length of new time series

#Step 4- Create new distribution for spp2
int_short<-subset(interact, intid==b[i])
mu_spp2<-mu_spp1-mean(int_short$newphenodiff) #adjust mean doy based on relative timing
sigma_spp2<-sample(sds$sd_doy, size=1)

bdist<-rnorm(1000, mu_spp2, sigma_spp2); #create distribution for spp2 based on moean and sd
spp2_data<-subset(int_long, spp=="spp2")
spp2_datab<-data.frame(array(0,c(timeseries[le],3))) #make new dataframe based on new length
spp2_datab[,1]<-rep(unique(spp2_data$studyid), timeseries[le])
spp2_datab[,2]<-rep(unique(spp2_data$intid), timeseries[le])
spp2_datab[,3]<-rep(unique(spp2_data$species), timeseries[le])
spp2_datab$ynull<-sample(bdist, size=timeseries[le]); #start with a 5 year time series
spp2_datab$yr1981<-sample(time2, size=timeseries[le]); #CHOOSE SAME YEARS AS SPP1


new[rowcount:asdf,]<-rbind(spp1_datab, spp2_datab)
rowcount<-rowcount+(timeseries[le]*2)
asdf<-rowcount+((timeseries[le]*2)-1)
}
names(new)[1]<-"studyid"; names(new)[2]<-"intid"; names(new)[3]<-"species"; names(new)[4]<-"ynull"; names(new)[5]<-"yr1981"

# Step 5- Run Stan
new <- new[with(new, order(species)),]
# clean to match syncmodels code, gets rid of non-unique data for repeating species within a study
Bgroups<-unique(new$species); b<-Bgroups; b<-as.character(b)
evennewer<-data.frame(array(0,c((timeseries[le]*2)*54,5)))
rowcount<-1
asdf<-rowcount+(timeseries[le]-1)
for(i in 1:length(b)){
spp<-subset(new, species==b[i])

if(length(unique(spp$intid))==1){
	spp<-subset(new, species==b[i])
	evennewer[rowcount:asdf,]<-spp
	rowcount<-rowcount+(timeseries[le])
	asdf<-rowcount+(timeseries[le]-1)
}
if(length(unique(spp$intid))>1){
	cgroups<-unique(spp$intid); c<-cgroups; c<-as.character(c)
	spp<-subset(new, species==b[i] & intid==c[sample(length(unique(spp$intid)), size=1)])
	evennewer[rowcount:asdf,]<-spp
rowcount<-rowcount+(timeseries[le])
asdf<-rowcount+(timeseries[le]-1)
}
}
names(evennewer)<-names(new)
evennewer2<-subset(evennewer, studyid!=0)
N <- nrow(evennewer2)
y <- evennewer2$ynull
Nspp <- length(unique(evennewer2$species))
J <- length(unique(evennewer2$species))
species <- as.numeric(as.factor(evennewer2$species))
year <- evennewer2$yr1981

null.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)



fh.sim <- extract(null.model) 
specieschar.formodel <- aggregate(rawlong.nodups["phenovalue"], rawlong.nodups[c("studyid", "species", "intid", "terrestrial","spp")], FUN=length) 
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species","intid"))
specieschar.formodel.sm  <- specieschar.formodel.sm [with(specieschar.formodel.sm , order(intid)),]
intid <- read.csv("input/raw_oct.csv", header=TRUE)
lal<-unique(rawlong.tot[,c("intid","terrestrial")])
intid2<-merge(intid, lal, by=c("intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction","terrestrial"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
#sync_int<-intid.nodups #synchrony change interactions


summ_studyspp <- subset(specieschar.formodel, select=c("studyid", "species")); summ_studyspp<-unique(summ_studyspp)
summ_studyspp  <- summ_studyspp[with(summ_studyspp , order(species)),]

#For interactions AND synchrony change
it1000 <- matrix(0, ncol=3000, nrow=length(unique(intid.sm$intid))) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 2000 iterations?
    summ_studyspp$model <- fh.sim$b[i,]
    andtheanswer <- merge(intid.nodups, summ_studyspp, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, summ_studyspp, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
    it1000[,(i-3000)] <- andtheanswer$model.x-andtheanswer$model.y #model.x=spp1
}
synch.model<-it1000


# FOR INTERPRETATION
#spp1-spp2=#resource-consumer
meanchange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH INTXN
andtheanswer$meanchange<-meanchange

#interact <- read.csv("input/raw_april.csv", header=TRUE)
#interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue) #spp1-spp2
#positive phenodiff= consumer emerges BEFORE resource
#negative phenodiff= resource emerges BEFORE consumer
season <- aggregate(interact["newphenodiff"], interact[c("studyid", "intid")], FUN=mean)

andtheanswer2<-merge(andtheanswer, season, by=c("studyid","intid"))

#increasing synchrony (i.e. species are getting closer together)
inc1<-subset(andtheanswer2, newphenodiff>0 & meanchange<0)
inc2<-subset(andtheanswer2, newphenodiff<0 & meanchange>0)
inc<-rbind(inc1, inc2)
inc$meanchange<-abs(inc$meanchange)
# decreasing synchrony (i.e. speceis are getting farther apart)
dec1<-subset(andtheanswer2, newphenodiff>0 & meanchange>0)
dec2<-subset(andtheanswer2, newphenodiff<0 & meanchange<0)
dec<-rbind(dec1, dec2)
dec$meanchange<-abs(dec$meanchange)
dec$meanchange<--dec$meanchange
tog<-rbind(dec, inc)


#to calculate mismatch
interact$count<-1
length<-aggregate(interact["count"], interact[c("studyid", "intid")], FUN=sum)
names(length)[3]<-"total"
neg<-subset(interact, newphenodiff<0)
length_neg<-aggregate(neg["count"], neg[c("studyid", "intid")], FUN=sum)
close<-merge(length, length_neg, by=c("studyid","intid"))
all<-merge(andtheanswer2, close, by=c("studyid","intid"))
all$count_diff<-with(all, count-total)

#fix mismatchs
sub<-subset(tog, intid=="1" | intid=="175" | intid=="235")
sub$meanchange<--sub$meanchange
sub2<-subset(tog, intid!="1" & intid!="175" & intid!="235")
tog<-rbind(sub, sub2)


# Magnitude only
mean(abs(tog$meanchange), na.rm=TRUE) #mean difference across interactions; they drift apart by half a day a decade

output[rowcount2:asdf2,1]<-timeseries[le]
output[rowcount2:asdf2,2:4]<-tog[,c("studyid","intid","meanchange")]

rowcount2<-rowcount2+nrow(tog)
asdf2<-rowcount2+(nrow(tog)-1)
}
names(output)[1]<-"length"; names(output)[2:4]<-c("studyid","intid","meanchange")

ggplot(subset(output, abs(meanchange)<0.4), (aes(x=length, y=abs(meanchange))))+geom_point()

sem<-sd(abs(tog$meanchange)/sqrt(length(abs(tog$meanchange))); sem
#95% confidence intervals of the mean
c(mean(abs(tog$meanchange))-2*sem,mean(abs(tog$meanchange))+2*sem)