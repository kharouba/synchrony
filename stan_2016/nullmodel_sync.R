#Aug 2016- Look at workflow from Lizzie- 2016Aug22_nullmodel.pdf
## Null model used for synchrony change

#Step 1- create pre_climate change dataset
#get interactions ready
interact <- read.csv("input/raw_oct.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue)

#load rawlong.tot
rawlong.tot2<-unique(rawlong.tot[,c("studyid","species","phenovalue","year","yr1981")]) 
rawlong.tot2$count<-1
pre<-subset(rawlong.tot2, year<=1981)
sss<- aggregate(pre["count"], pre[c("studyid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas

pre_cc<-merge(rawlong.tot2, sss2, by=c("studyid", "species"))
pre_cc<-subset(pre_cc, year<=1981)


#Step 2- Create distribution of means and sd (2)
means <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "species")], FUN=mean)
names(means)[3]<-"mean_doy" #mean of each dataset
sds <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "species")], FUN=sd)
names(sds)[3]<-"sd_doy" #sd of each dataset

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

# Step 5- clean to match syncmodels code, gets rid of non-unique data for repeating species within a study (i.e. those species across multiple intxns, e.g. Accipiter nisus) so that models for unique species are built, not unique species-intxns ; Number of unique species should match syncmodels i.e. n=86
new <- new[with(new, order(species)),]
Bgroups<-unique(new$species); b<-Bgroups; b<-as.character(b)
evennewer<-data.frame(array(0,c(nrow(new)*54,19)))
rowcount<-1

for(i in 1:length(b)){
spp<-subset(new, species==b[i])
asdf<-rowcount+(nrow(spp)-1)

if(length(unique(spp$intid))==1){
	evennewer[rowcount:asdf,]<-spp
	rowcount<-rowcount+nrow(spp)
	#asdf<-rowcount+(nrow(spp)-1)
}
if(length(unique(spp$intid))>1){
	cgroups<-unique(spp$intid); c<-cgroups; c<-as.character(c)
	j<-sample(length(unique(spp$intid)), size=1)
	spp2<-subset(spp, intid==c[j])
	asdf<-rowcount+(nrow(spp2)-1)
	evennewer[rowcount:asdf,]<-spp2
	rowcount<-rowcount+nrow(spp2)
}
}
names(evennewer)<-names(new)
evennewer2<-subset(evennewer, studyid!=0)


# step 6- Run Stan
N <- nrow(evennewer2)
y <- evennewer2$ynull
Nspp <- length(unique(evennewer2$species))
J <- length(unique(evennewer2$species))
species <- as.numeric(as.factor(evennewer2$species))
year <- evennewer2$yr1981


null.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

goo <- extract(null.model) 
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

interact <- read.csv("input/raw_oct.csv", header=TRUE)
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
tog<-rbind(dec, inc); nrow(tog)

!! IGNORING MISMATCHES BECAUSE SHOULDNT BE ANY 
!! double check sample size !!

mean(tog$meanchange, na.rm=TRUE) #mean difference across interactions; they drift apart by half a day a decade
#computation of the standard error of the mean
sem<-sd(tog$meanchange)/sqrt(length(tog$meanchange)); sem
#95% confidence intervals of the mean
c(mean(tog$meanchange)-2*sem,mean(tog$meanchange)+2*sem)

## Magnitude
mean(abs(tog$meanchange), na.rm=TRUE) #
sem<-sd(abs(tog$meanchange))/sqrt(length(abs(tog$meanchange))); sem
#95% confidence intervals of the mean
c(mean(abs(tog$meanchange))-2*sem,mean(abs(tog$meanchange))+2*sem)

Histogram figures
text_high <- textGrob("Closer together", gp=gpar(fontsize=13, fontface="bold"))
text_low <- textGrob("Further apart", gp=gpar(fontsize=13, fontface="bold"))
ggplot(tog, aes(x=meanchange))+geom_histogram(binwidth=.5, alpha=.5, position="identity", colour="black")+theme_bw()+geom_vline(xintercept=-0.054, linetype="solid",size=1)+geom_vline(xintercept=0.0022, linetype=2,size=1)+geom_xlim(-1.5, 1.5)+annotation_custom(text_high, xmin=1.5, xmax=1.5, ymin=-0.5, ymax=-0.5)+annotation_custom(text_low, xmin=-1.5, xmax=-1.5, ymin=-0.5, ymax=-0.5)+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/year")


ggplot(tog, aes(x=meanchange))+geom_histogram(position="identity", colour="black")+theme_bw()+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of interactions")+xlab("Change in number of days/year")+xlim(-1.5, 1.5)+geom_vline(xintercept=0, linetype="solid",size=1)


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