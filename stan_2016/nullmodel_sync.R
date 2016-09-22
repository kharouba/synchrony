#Aug 2016- Look at workflow from Lizzie- 2016Aug22_nullmodel.pdf
## Null model used for synchrony change


#Step 1- create pre_climate change dataset
#get interactions ready
interact <- read.csv("input/raw_april.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue)

#load rawlong.tot
rawlong.tot$count<-1
pre<-subset(rawlong.tot, year<=1981)
sss<- aggregate(pre["count"], pre[c("studyid", "intid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas

pre_cc<-merge(rawlong.tot, sss2, by=c("studyid", "intid", "species"))

#Step 2- Create distribution of means and sd (2)
means <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "intid", "species")], FUN=mean)
names(means)[4]<-"mean_doy" #mean of each dataset
sds <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "intid", "species")], FUN=sd)
names(sds)[4]<-"sd_doy" #sd of each dataset

# FOR EACH INTERACTION
Bgroups<-unique(rawlong.tot$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0,c(nrow(rawlong.tot),19)))
rowcount<-1
for(i in 1:length(b)){
int_long<-subset(rawlong.tot, intid==b[i])
asdf<-rowcount+(nrow(int_long)-1)
#Step 3- Create new distribution for spp1 from sample mean and sd
mu_spp1<-sample(means$mean_doy, size=1) #draw one mean from means
sigma_spp1<-sample(sds$sd_doy, size=1) #draw one sd from distribution of SDs
a<-rnorm(1000, mu_spp1, sigma_spp1); #create distribution based on moean and sd

#Step 3- Create new time-series for spp1 of X length based on normal distribution
spp1_data<-subset(int_long, spp=="spp1")
spp1_data$ynull<-sample(a, size=nrow(spp1_data)); 

#Step 4- Create new distribution for spp2
int_short<-subset(interact, intid==b[i])
mu_spp2<-mu_spp1-mean(int_short$newphenodiff) #adjust mean doy based on relative timing
sigma_spp2<-sample(sds$sd_doy, size=1)

bdist<-rnorm(1000, mu_spp2, sigma_spp2); #create distribution for spp2 based on moean and sd
spp2_data<-subset(int_long, spp=="spp2")
spp2_data$ynull<-sample(bdist, size=nrow(spp2_data))
new[rowcount:asdf,]<-rbind(spp1_data, spp2_data)
rowcount<-rowcount+nrow(int_long)
asdf<-rowcount+(nrow(int_long)-1)
}
names(new)[19]<-"ynull"
names(new)[1:18]<-names(rawlong.tot)[1:18]

# Step 5- Run Stan
new <- new[with(new, order(species)),]
N <- nrow(new)
y <- new$ynull
specieschar.hin<- aggregate(new["ynull"], new[c("studyid", "species", "int_type", "terrestrial", "spp")], FUN=length) #number of years per species
specieschar.hin <- specieschar.hin[with(specieschar.hin, order(species)),]
Nspp <- nrow(specieschar.hin)
J <- nrow(specieschar.hin)
species <- as.numeric(as.factor(new$species))
new$yr1981 <- new$newyear-1981
year <- new$yr1981

null.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

goo <- extract(null.model) 
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species"))
specieschar.formodel.sm <- specieschar.formodel.sm[with(specieschar.formodel.sm, order(species)),]
intid <- read.csv("input/raw_april.csv", header=TRUE)
lal<-unique(rawlong.tot[,c("intid","terrestrial")])
intid2<-merge(intid, lal, by=c("intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction","terrestrial"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
intid.nodups <- intid.nodups[with(intid.nodups, order(species)),]
intid.nodups<-na.omit(intid.nodups)
sync_int<-intid.nodups #synchrony change interactions

summ_studyspp <- subset(specieschar.formodel, select=c("studyid", "species")); summ_studyspp<-unique(summ_studyspp)

## To estimate pheno change per species
#it1000 <- matrix(0, ncol=3000, nrow=Nspp)
#for (i in 3000:6000){ # 3000 iterations?
#    summ_studyspp$model <- goo$b[i,]
#    it1000[,(i-3000)] <- goo$b[i,]
#}
#summ_studyspp$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP
#mean(summ_studyspp$stanfit) 

## TO measure synchrony change
it1000 <- matrix(0, ncol=3000, nrow=length(unique(intid.sm$intid))) #2000 iterations for 53 interactions;
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- goo$b[i,]
    andtheanswer <- merge(intid.nodups, summ_studyspp, by.x=c("studyid", "spp1"),
        by.y=c("studyid", "species"), all.x=TRUE)
    andtheanswer <- merge(andtheanswer, summ_studyspp, by.x=c("studyid", "spp2"),
        by.y=c("studyid", "species"), all.x=TRUE)
    it1000[,(i-3000)] <- andtheanswer$model.x-andtheanswer$model.y #model.x=spp1
}

synch.model<-it1000

# FOR INTERPRETATION
meanchange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH INTXN
andtheanswer$meanchange<-meanchange
#Negative difference= decreasing synchrony ()
#Positive difference= increasing synchrony ()
# spp1 = negative species in trophic interaction and spp2= positive species

interact <- read.csv("input/raw_april.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue) #spp1-spp2
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

mean(tog$meanchange, na.rm=TRUE) #mean difference across interactions; they drift apart by half a day a decade
#computation of the standard error of the mean
sem<-sd(tog$meanchange)/sqrt(length(tog$meanchange)); sem
#95% confidence intervals of the mean
c(mean(tog$meanchange)-2*sem,mean(tog$meanchange)+2*sem)

### Extra
sub<-subset(rawlong, intid=="1")
ggplot(sub, aes(x=year, y=phenovalue))+geom_point()+facet_wrap(~spp)

1. Randomly choose pre-climate change dataset(<1981, n=11)
2. Run stan on this dataset
3. Simulate new data with this model for 'new' datset i.e. all years
4. Run stan again but now across all years of dataset
5. Repeat steps above 1000x
6. Match mean slope with interacting species to calcuate sync change

#predict
sub<-singlespeciesdata
sub2<-subset(singlespeciesdata, year<=1981)
m1<-lm(phenovalue~year, data=sub2)
newdata<-data.frame(phenovalue=sub$phenovalue, year=sub$newyear)
sub$ypred<-predict(m1, newdata)

#stan predict (i.e. posterior simulations)
In Stan, posterior simulations can be generated in two ways. The first approach is to treat the predicted variables as parameters and then define their distributions in the model block. The second approach, which also works for discrete variables, is to generate replicated data using random-number generators in the generated quantities block.

add:
data{
	real new_x
}

generated quantities {
  real y_sim;
  y_sim <- normal(alpha+new_x * beta, sigma);
}

READ: https://groups.google.com/forum/#!topic/stan-users/5GJRlkTQ3Zk


Overall approach:
Randomize phenovalues for all species. Keep same data structure as full dataset (i.e same number of years per species (all years included).

Work flow:

For each of 92 species, randomly select one of the pre datasets. Run bootstrap (n=1000) on phenovalues of that selected dataset. Use mean and sd from bootstrap to generate distribution (rnorm) to generate new phenovalues of the same length as full dataset of species x.


load rawlong.tot

rawlong.tot$count<-1
pre<-subset(rawlong.tot, year<=1981)
sss<- aggregate(pre["count"], pre[c("studyid", "intid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datasets

b<-unique(rawlong.tot$species)

for(i 1:length(unique(rawlong.tot$species)))
sub<-subset(rawlong.tot, species=="Acartia hudsonica")
newsub<-subset(sss2, speciesid==sample(32, 1)) #sample 1 species/dataset from 32 datasets
newpre<-subset(pre, species==newsub$species) #isolate pre dataset for this species
bootmean<-boot(newpre$phenovalue, function(x, index) mean(x[index]), R=1000) #bootstrap to generate mean and sd
newphenovalue<-rnorm(nrow(sub), mean(bootmean$t), sd(bootmean$t)) #sample from normal dsitribution based on this mean and sd