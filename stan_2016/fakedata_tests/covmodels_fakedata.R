For 3-level random effects model- matches threelevelrandomeffects.stan

#lm for predictive check- for interactions
uniquespp<-unique(newframe2$intid)
slopefits<-rep(NA< length(uniquespp))
varfits<-rep(NA< length(uniquespp))
intfits<-rep(NA< length(uniquespp))
for(eachsp in 1:length(uniquespp)){
	lmhere<-lm(sync_change~temp_change, data=subset(newframe2, species==uniquespp[eachsp]))
	slopefits[eachsp]<-coef(lmhere)[2]
	varfits[eachsp]<-(summary(lmhere)$sigma)**2
	intfits[eachsp]<-coef(lmhere)[1]
}
dr<-data.frame(cbind(uniquespp, slopefits, varfits, intfits))
mean(dr$slopefits); sd(dr$slopefits)


# True parameter values
mu_a<-0.1719781 # overall intercept: mean response of raw data, here mean sync_change across all iterations (e.g. doy) 
mu_b<--4.914741 # overall slope across all data (i.e interactions)- sync_change~temp_change (like mu_b)
sigma_y<-0.9932294 #sd associated with response

sigma_a<-0.05
a_spp<-rnorm(Nj, mu_a, sigma_a) # create intercepts for each interaction


sigma_b<-1.42784 #sd of mean slope across interactions
b_spp<-rnorm(Nj, mu_b, sigma_b) # create slopes for each interaction


#Simulate/create the data
Ni<-2200
Nj<-22
Nk<-13
n_data_per_int <- round(runif(Nj, 100, 100)) # how many iterations per interaction
intid<-rep(1:Nj, n_data_per_int) #adds intid
year<-rep(NA, Ni) #sub for temp_change

for(j in 1:Nj){
	year[intid==j]<-rnorm(n_data_per_int[j], mean(newframe2$temp_change), sd(newframe2$temp_change)) #for each iteration for each interaction, pull temp_change value from distribution
}
desMat <- model.matrix(object = ~ 1 + year)

p<-ncol(desMat)

ypred<-length(Ni)
for (i in 1:Ni){ 
s <- intid[i] #sppid for each row
ypred[i]<-beta_0jk[intid[s]]+beta_1jk[intid[s]]*desMat[i,2];
}
y<-rnorm(Ni,ypred,sigma_y)

#intid <- as.integer(as.factor(newframe2$intid))
studyid <- as.integer(as.factor(newframe2$studyid))
studyLookup <- unique(newframe2[c("intid","studyid")])[,"studyid"] # C
#nVars<-1
#Imat <- diag(1, nVars)

fit_simple<-stan("threelevelrandomeffects.stan", data=c("Ni","Nj","Nk","p","intid","studyLookup","year","y"), iter=2000, chains=4)

print(fit_simple, pars=c("mu_a", "mu_b", "sigma_y", "sigma_a_spp","sigma_b_spp"))

----
# Create the interaction-level parameters
beta_0jk<-rnorm(J, 0, 1) #mean intercept across interactions
sigma_0jk<-
lmfits <- read.csv("lmfits.csv", header=TRUE); names(lmfits)[2]<-"species"
spp1<-merge(lmfits[,2:4], intid.nodups[,c("studyid","spp1","intid")], by.x=c("species"), by.y="spp1")
names(spp1)[1]<-"spp1"; names(spp1)[2]<-"spp1_lmfits"
spp2<-merge(lmfits[,2:4], intid.nodups[,c("studyid","spp2","intid")], by.x=c("species"), by.y="spp2")
names(spp2)[1]<-"spp2"; names(spp2)[2]<-"spp2_lmfits"
spp<-merge(spp1, spp2, by=c("studyid","intid"))
#join species to form interactions
beta_1k<-rnorm(J, 0, 1) #mean slope across interactions

#Create study-level parameters
Nk<- 13 #studies
aggregate(lmfits2["lmfits"],)

beta_0k<-rnorm(J, 0, 1) #mean intercept across interactions
beta_1k<-rnorm(J, 0, 1) #mean slope across interactions
