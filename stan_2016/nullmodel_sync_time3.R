#Aug 2016- Look at workflow from Lizzie- 2016Aug22_nullmodel.pdf
## Null model used for synchrony change, to explore effect of length of time series on ability to detect sync change

#Step 1- create pre_climate change dataset
#load rawlong.tot
interact <- read.csv("input/raw_oct.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue) #spp1-spp2


rawlong.tot$count<-1
pre<-subset(rawlong.tot, year<=1981)
sums<-unique(rawlong.tot[,c("studyid","intid","species")])

#some species have non-unique data e.g. Glis glis, C

pre_uni<-unique(pre[,c("studyid","year","species","phenovalue","yr1981","count")])

#pres<-merge(pre, unique(rawlong.tot[,c("studyid","intid","species")]), by=c("studyid","species"))
sss<- aggregate(pre_uni["count"], pre_uni[c("studyid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas

pre_cc<-merge(pre_uni[,1:5], sss2, by=c("studyid","species"))
#pre_cc2<-merge(pre_cc, unique(rawlong.tot[,c("studyid","intid","species")]), by=c("studyid","species"))
#pre_cc2 <- pre_cc2[with(pre_cc2, order(species, year)),]
#pre_cc<-merge(pres[,c(1:4,6)],  sss2, by=c("studyid", "intid", "species"))

OR Use all data
rawlong.tot2<-unique(rawlong.tot[,c("studyid","species","phenovalue","yr1981")]) #CLEAN UP so only unique values across repeating species within studoes
pre_cc<-rawlong.tot2


!!! Manual up to here!!


#Step 2- Create distribution of means and sd (2)
means <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "species")], FUN=mean) #(not intid)
names(means)[3]<-"mean_doy" #mean of each dataset
sds <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "species")], FUN=sd)
names(sds)[3]<-"sd_doy" #sd of each dataset

# For each length of time series

timeseries<-c(5,10,15,20,25,30,35,40)
rowcount2<-1
asdf2<-rowcount2+(86-1) #number of interactions
output<-data.frame(array(0,c(8*33,5)))

for(le in 1:length(timeseries)){ #length(timeseries)


# FOR EACH INTERACTION
Bgroups<-unique(rawlong.tot$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0,c((timeseries[le]*2)*54,6)))
rowcount<-1
for(i in 1:length(b)){
int_long<-subset(rawlong.tot, intid==b[i])
asdf<-rowcount+((timeseries[le]*2)-1) #length of NEW time series*2 species

#Step 3- Create new distribution for spp1 from sample mean and sd
mu_spp1<-sample(means$mean_doy, size=1) #draw one mean from means
sigma_spp1<-sample(sds$sd_doy, size=1) #draw one sd from distribution of SDs
adist<-rnorm(1000, mu_spp1, sigma_spp1); #create distribution of phenovalues based on moean and sd


#Step 3- Create new time-series for spp1 of X length based on normal distribution
spp1_data<-subset(int_long, spp=="spp1")
time<-spp1_data$yr1981
time2<-rnorm(1000, mean(time), sd(time)) # create distribution of years based on observed
spp1_datab<-data.frame(array(0,c(timeseries[le],4))) #make new dataframe based on new length
spp1_datab[,1]<-rep(unique(spp1_data$studyid), timeseries[le])
spp1_datab[,2]<-rep(unique(spp1_data$intid), timeseries[le])
spp1_datab[,3]<-rep(unique(spp1_data$species), timeseries[le])
spp1_datab[,4]<-"spp1"
spp1_datab$ynull<-sample(adist, size=timeseries[le]); #start with a 5 year time series
spp1_datab$yr1981<-sample(time2, size=timeseries[le]); #sample X number of years based on length of new time series

#Step 4- Create new distribution for spp2 #NOT IN REFERENCE TO SPP1
mu_spp2<-sample(means$mean_doy, size=1) #draw one mean from means
sigma_spp2<-sample(sds$sd_doy, size=1) #draw one sd from distribution of SDs
bdist<-rnorm(1000, mu_spp2, sigma_spp2); #create distribution of phenovalues based on moean and sd


spp2_data<-subset(int_long, spp=="spp2")
time<-spp2_data$yr1981
time2<-rnorm(1000, mean(time), sd(time)) # create distribution of years based on observed

spp2_datab<-data.frame(array(0,c(timeseries[le],4))) #make new dataframe based on new length
spp2_datab[,1]<-rep(unique(spp2_data$studyid), timeseries[le])
spp2_datab[,2]<-rep(unique(spp2_data$intid), timeseries[le])
spp2_datab[,3]<-rep(unique(spp2_data$species), timeseries[le])
spp2_datab[,4]<-"spp2"
spp2_datab$ynull<-sample(bdist, size=timeseries[le]); #start with a 5 year time series
spp2_datab$yr1981<-spp1_datab$yr1981; #CHOOSE SAME YEARS AS SPP1


new[rowcount:asdf,]<-rbind(spp1_datab, spp2_datab)
rowcount<-rowcount+(timeseries[le]*2)
asdf<-rowcount+((timeseries[le]*2)-1)
}
names(new)[1]<-"studyid"; names(new)[2]<-"intid"; names(new)[3]<-"species"; names(new)[4]<-"spp"; names(new)[5]<-"ynull"; names(new)[6]<-"yr1981"

#Step 5- clean to match syncmodels code, gets rid of non-unique data for repeating species within a study (i.e. those species across multiple intxns, e.g. Accipiter nisus) so that models for unique species are built, not unique species-intxns ; Number of unique species should match syncmodels i.e. n=85
new <- new[with(new, order(species)),]
Bgroups<-unique(new$species); b<-Bgroups; b<-as.character(b)
evennewer<-data.frame(array(0,c((timeseries[le]*2)*54,6)))
rowcount<-1
asdf<-rowcount+(timeseries[le]-1)
for(i in 1:length(b)){
sun<-subset(new, species==b[i])

if(length(unique(sun$intid))==1){
	#sun<-subset(new, species==b[i])
	evennewer[rowcount:asdf,]<-sun
	rowcount<-rowcount+(timeseries[le])
	asdf<-rowcount+(timeseries[le]-1)
}
if(length(unique(sun$intid))>1){
	cgroups<-unique(sun$intid); c<-cgroups; c<-as.character(c)
	sun<-subset(new, species==b[i] & intid==c[sample(length(unique(sun$intid)), size=1)])
	evennewer[rowcount:asdf,]<-sun
rowcount<-rowcount+(timeseries[le])
asdf<-rowcount+(timeseries[le]-1)
}
}
names(evennewer)<-names(new)
evennewer2<-subset(evennewer, studyid!=0)

# Step 6- Run Stan

N <- nrow(evennewer2)
y <- evennewer2$ynull
Nspp <- length(unique(evennewer2$species))
J <- length(unique(evennewer2$species))
species <- as.numeric(as.factor(evennewer2$species))
year <- evennewer2$yr1981

null.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)


trial<-summary(null.model, pars="b")
#to get median coefficients from SUMMARY
median<-trial[[1]][1:86]; #new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) 
min<-trial[[1]][259:344]; e<-data.frame(min=unlist(min)) #-1.129 to -1.022
max<-trial[[1]][603:688]; f<-data.frame(max=unlist(max)) #-0.3077 to 0.14
d$min<-e$min; d$max<-f$max;

#uni<-unique(evennewer2[,c("studyid","intid","species","spp")])
#tot<-cbind(uni,d)
#tot <- tot[with(tot, order(intid, spp)),]
#spp1s<-subset(tot, spp=="spp1"); names(spp1s)[5]<-"spp1_pheno"
#spp2s<-subset(tot, spp=="spp2"); names(spp2s)[5]<-"spp2_pheno"
#tot<-merge(spp1s, spp2s[,c("intid","spp2_pheno")], by=c("intid"))

output[rowcount2:asdf2,1]<-timeseries[le]
output[rowcount2:asdf2,2:5]<-d

rowcount2<-rowcount2+nrow(d)
asdf2<-rowcount2+(nrow(d)-1)

}

names(output)[1]<-"length"; names(output)[2:5]<-c("phenochange","species","min","max")
output$id<-with(output, paste(length, species, sep="."))
output$id<-as.numeric(output$id)

ggplot(output, (aes(x=factor(id), y=phenochange)))+geom_errorbar(aes(ymin=min, ymax=max, width=.0025, colour="black"))+geom_point()+theme(legend.position="false")+geom_hline(yintercept=0, linetype="dashed")

lol<- aggregate(output["phenochange"], output[c("length")], FUN=mean); names(lol)[2]<-"mean"; 
st.err <- function(x) {sd(x)/sqrt(length(x))}
lol2<- aggregate(output["phenochange"], output[c("length")], FUN=st.err);
names(lol2)[2]<-"se"
lols<-merge(lol, lol2, by="length")

ggplot(lols, (aes(x=factor(length), y=mean)))+geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=.0025, colour="black"))+geom_point()+theme(legend.position="false")+geom_hline(yintercept=0, linetype="dashed")

# attempt to link with synchrony change
spp1<-subset(output, species<44)
names(spp1)[2]<-"spp1_phenochange"; names(spp1)[3]<-"spp1_species";names(spp1)[4]<-"spp1_min"; names(spp1)[5]<-"spp1_max"; names(spp1)[6]<-"spp1_id"
spp2<-subset(output, species>43)
names(spp2)[2]<-"spp2_phenochange"; names(spp2)[3]<-"spp2_species"; names(spp2)[4]<-"spp2_min"; names(spp2)[5]<-"spp2_max"; names(spp2)[6]<-"spp2_id"

tot<-cbind(spp1, spp2)
tot$sync<-with(tot, spp1_phenochange-spp2_phenochange)

ggplot(tot, (aes(x=factor(spp1_id), y=sync, colour=factor(length))))+geom_point()+theme(legend.position="false")+geom_hline(yintercept=0, linetype="dashed")

lol<- aggregate(tot["sync"], tot[c("length")], FUN=mean); names(lol)[2]<-"mean"; 
st.err <- function(x) {sd(x)/sqrt(length(x))}
lol2<- aggregate(tot["sync"], tot[c("length")], FUN=st.err);
names(lol2)[2]<-"se"
lols<-merge(lol, lol2, by="length")


ggplot(lols, (aes(x=factor(length), y=mean)))+geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=.0025, colour="black"))+geom_point()+theme(legend.position="false")+geom_hline(yintercept=0, linetype="dashed")+ggtitle("Mean synchrony change +/-SE")
