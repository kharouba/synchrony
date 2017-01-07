#Dec 2016- Look at workflow from Lizzie- 2016Nov18_nullmodel.pdf
## Null model used for synchrony change, to explore effect of length of time series on ability to detect sync change
# Calculates slopes first

st.err <- function(x) {sd(x)/sqrt(length(x))}

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
#31 species from pre_cc (31 spp all together, 2 spp repeat across intxns)

null.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)
faa <- extract(null.model)
solo<-unique(pre_cc[,c("studyid","species")])

it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    solo$model <- faa$b[i,]
    it1000[,(i-3000)] <- faa$b[i,]
}
solo$phenochange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EAC

#mean phenochange=-0.26

# Step 3- Simulate data for all speices and re-run stan

timeseries<-c(5,10,15,20,25,30,35,40)
rowcount<-1
rowcount2<-1
asdf2<-rowcount2+(88-1)
evens<-data.frame(array(0,c(88*8, 5)))
names(evens)[1]<-"species"; names(evens)[2]<-"time_length";names(evens)[3]<-"phenochange"; names(evens)[4]<-"min"; names(evens)[5]<-"max"

for(le in 1:length(timeseries)){ #length(timeseries)
asdf<-rowcount+(timeseries[le]-1) #number of interactions
new<-data.frame(array(0,c(20000,4)))
names(new)[1]<-"species"; names(new)[2]<-"time_length";names(new)[3]<-"y_null"; names(new)[4]<-"year_null"

for(i in 1:88){ # for each species
year_null<-sample(rawlong.tot2$year, size=timeseries[le], replace=FALSE) #choose years from observed distribution of years from full dataset

sppa_b<-sample(faa$mu_b, size=1)
ypred_a<-(sample(faa$a, size=timeseries[le]))+(sppa_b*year_null)
y_a <- rnorm(timeseries[le], mean(ypred_a), sample(faa$sigma_y, size=1));

new[rowcount:asdf,1]<-i
new[rowcount:asdf,2]<-timeseries[le]
new[rowcount:asdf,3]<-y_a
new[rowcount:asdf,4]<-year_null

rowcount<-rowcount+timeseries[le]
asdf<-rowcount+(timeseries[le]-1)

}
new<-subset(new, species!=0)

N <- nrow(new)
y <- new$y_null
Nspp <- length(unique(new$species))
species <- as.numeric(as.factor(new$species))
year <- new$year_null

new.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)


asdf<-summary(new.model, pars="b")
#to get median coefficients from SUMMARY
median<-asdf[[1]][1:88]; #new<-as.data.frame(y); #number of species =91
d<-data.frame(y=unlist(median), grp=1:length(median)) 
min<-asdf[[1]][265:352]; e<-data.frame(min=unlist(min)) #-0.026 -0.044
max<-asdf[[1]][617:704]; f<-data.frame(max=unlist(max)) #0.018 0.0084
d$min<-e$min; d$max<-f$max;


evens[rowcount2:asdf2,1]<-unique(species)
evens[rowcount2:asdf2,2]<-timeseries[le]
evens[rowcount2:asdf2,3]<-d$y
evens[rowcount2:asdf2,4]<-d$min
evens[rowcount2:asdf2,5]<-d$max

rowcount2<-rowcount2+88
asdf2<-rowcount2+(88-1)

}


evens2<-subset(evens, species!="0")
evens2$names<-with(evens2, paste(time_length, species,sep="."))
evens2$names<-as.numeric(evens2$names)
evens2  <- evens2[with(evens2, order(names)),]

write.csv(evens2, "/output/nullmodel_pheno_time_byslope.csv")

evens2$names2<-evens2$names
sub<-subset(evens2, species=="1");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="1");
evens2<-rbind(evens2, sub)

sub<-subset(evens2, species=="2");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="2");
evens2<-rbind(evens2, sub)

sub<-subset(evens2, species=="3");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="3");
evens2<-rbind(evens2, sub)

sub<-subset(evens2, species=="4");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="4");
evens2<-rbind(evens2, sub)

sub<-subset(evens2, species=="5");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="5");
evens2<-rbind(evens2, sub)

sub<-subset(evens2, species=="6");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="6");
evens2<-rbind(evens2, sub)

sub<-subset(evens2, species=="7");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="7");
evens2<-rbind(evens2, sub)

sub<-subset(evens2, species=="8");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="8");
evens2<-rbind(evens2, sub)

sub<-subset(evens2, species=="9");
sub$names2<-with(sub, names+0.005)
evens2<-subset(evens2, species!="9");
evens2<-rbind(evens2, sub)

evens2  <- evens2[with(evens2, order(names2)),]

#ggtitle("Mean phenological change with 95% CI")
ggplot(evens2, (aes(x=factor(names2), y=phenochange, colour=factor(time_length))))+geom_errorbar(aes(ymin=min, ymax=max, width=.0025, colour=factor(time_length)))+geom_point(colour="black")+theme(legend.position="false")+geom_hline(yintercept=0, linetype="dashed", size=1)+geom_vline(xintercept=264, linetype="dotted", size=0.75)+geom_vline(xintercept=352, linetype="dotted", size=0.75)+theme(axis.title.x = element_text(size=15), axis.text.x=element_blank(), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Phenological change (days/yr)")+xlab("Time series length")


##
baas <- extract(new.model)
boo <-as.data.frame(unique(new[,c("species")])); 

it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    boo$model <- baas$b[i,]
    it1000[,(i-3000)] <- baas$b[i,]
}
boo$phenochange <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EAC
names(boo)[1]<-"species"

