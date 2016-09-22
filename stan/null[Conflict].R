# Simulate fake data

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")


# RANDOM SLOPE AND INTERCEPT MODEL #
# Create the species-level parameters
J <- 71 # number of species
#N<- 1950 # number of data pots
species_level <- rnorm(J, 0, 1) #species error
species_trend <- rnorm(J, 0, 1)

year_0 <- 1976
sigma_y <- 5 #measurement error,  noise

#create hyperparameters
mu_a <- c(125, 125) # mean intercept across species; based on intercepts from latest stan model, mean is 141
sigma_a <- c(25, 25) # variance around mean intercept; based on latest stan model, mean SE is 0.11, mean SD is 5.49
mu_b <- c(-.1, -.2) # mean slope?; based on latest stan model, range is -1.74 to 1.21, mean is -0.31
sigma_b <- c(.2, .2)  # variance around slope?  0.0020 to 0.0103



# Create the data
n_data_per_species <- round(runif(J, 5, 40)) #median is 24, max is 37
species <- rep(1:J, n_data_per_species)
N <- length(species)
year <- rep(NA, N)
for (j in 1:J){
  year[species==j] <- rev(2009 - 1:(n_data_per_species[j])) - year_0
}


y<-rnorm(N, mean=134, sd=48.28) #mean phenovalue and sd
year<-rnorm(N, mean=16, sd=11.02) #median year is 1992, 16 years from 1976
nVars<-1
Imat <- diag(1, nVars)
fit.notype <- stan("synchrony1_notype_hk.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=3)

[IDEA IS TO RANDOMLY ASSIGN PAIR-WISE INTERACTIONS TO 54 SPECIES, TAKE DIFFERENCE BETWEEN SPECIES, TAKE MEDIAN DIFFERENCE ACROSS INTERACTIONS, REPEAT 999X; problem- 54 unique species but 82 species]

raw <- read.csv("raw.csv", header=TRUE)
rawlong <- read.csv("indiv_sppdata2.csv", header=TRUE)
rawlong <- rawlong[with(rawlong, order(intid)),]
rawlong$int_type[which(rawlong$int_type=="poll")] <- "pollination"
rawlong.X <- NULL

# DON'T CHANGE NAMES AS BELOW, INSTEAD USE RAWLONG AFTER EDITS IN SYNCMODELS_LIZZIE3 #
rawlong$species[rawlong$species=="Ficedula1 albicollis"] <- "Ficedula albicollis"
rawlong$species[rawlong$species=="Ficedula2 albicollis"] <- "Ficedula albicollis"
rawlong$species[rawlong$species=="Keratella2 cochlearis"]<-"Keratella cochlearis"
rawlong$species[rawlong$species=="Parus1 major"]<-"Parus major"
rawlong$species[rawlong$species=="Parus2 major"]<-"Parus major"
rawlong$species[rawlong$species=="Parus3 major"]<-"Parus major"
rawlong$species[rawlong$species=="Parus4 major"]<-"Parus major"
rawlong$species[rawlong$species=="Quercus1 robur"]<-"Quercus robur"
rawlong$species[rawlong$species=="Quercus2 robur"]<-"Quercus robur"
rawlong$species[rawlong$species=="Quercus3 robur"]<-"Quercus robur"
rawlong$species[rawlong$species=="Mnemiopsis2 leidyi"]<-"Mnemiopsis leidyi"
rawlong$species[rawlong$species=="Mnemiopsis1 leidyi"]<-"Mnemiopsis leidyi"
rawlong$species[rawlong$species=="Rissa1 tridactyla"]<-"Rissa tridactyla"
rawlong$species[rawlong$species=="Rissa2 tridactyla"]<-"Rissa tridactyla"
rawlong$species[rawlong$species=="Parus2 caeruleus"]<-"Parus caeruleus"
rawlong$species[rawlong$species=="Caterpillar2 spp."]<-"Caterpillar spp."
rawlong$species[rawlong$species=="Copepod1 spp."]<-"Copepod spp."
rawlong$species[rawlong$species=="Copepod2 spp."]<-"Copepod spp."
rawlong$species[rawlong$species=="Daphnia1 spp."]<-"Daphnia spp."
rawlong$species[rawlong$species=="Daphnia2 spp."]<-"Daphnia spp."
rawlong$species[rawlong$species=="Daphnia3 spp."]<-"Daphnia spp."
rawlong$species[rawlong$species=="Diatom1 spp."]<-"Diatom spp."
rawlong$species[rawlong$species=="Diatom2 spp."]<-"Diatom spp."
rawlong$species[rawlong$species=="Diatom3 spp."]<-"Diatom spp."
rawlong$species[rawlong$species=="Diatom4 spp."]<-"Diatom spp."
rawlong$species[rawlong$species=="Phytoplankton1 spp."]<-"Phytoplankton spp."
rawlong$species[rawlong$species=="Phytoplankton2 spp."]<-"Phytoplankton spp."
rawlong$species[rawlong$species=="Plant2 spp."]<-"Plant spp."
rawlong$species[rawlong$species=="Plant3 spp."]<-"Plant spp."

#randomly order dataframe to assign interactions
df2 <- rawlong[sample(nrow(rawlong)),]; df2$count<-1 # use rawlong  not rawlong.nodups because not all interactions have 2 spp
summ<-unique(rawlong[,c("intid","species")]); summ$count<-1 # use rawlong not rawlong.nodups because not all interactions have 2 spp
summ <- summ[with(summ, order(species)),]

summ2 <- aggregate(summ["count"], summ[c("species")], FUN=length)
#tab<-with(summ2, table(species,count))
#vec<-rep()
#vec <- rep(tab$value, tab$freq)

sam<-data.frame(array(0, c(sum(summ2$count), 2))); #of unique spp x #unique intxns (nrow summ)
#vec<-c(rep(1:55,1), rep(56:64,2), rep(65:67,3), rep(68:69,4), rep(70:71,5)); #9 spp repeated twice- because in two intxns, 3 species repeated 3x, 2 species repeated 4x, 2 spp repeated 5x
#vec<-c(rep(1:32,1), rep(33:45,2), rep(46:47,3), rep(48:50,4), rep(51:52,5), rep(53, 6), rep(54, 8)); #13 spp repeated twice- because in two intxns, 2 species repeated 3x, 3 species repeated 4x, 2 spp repeated 5x, 1 spp repeated 6x, 1 spp repeated 8x
OR
#IF each 'unique' spp included
vec<-c(rep(1:73,1), rep(74:77,2), rep(78:79,3), rep(80:81,4), rep(82,5));  #4 spp repeated 2x, 2 spp repeated 3x; 2 spp repeated 4x, 1 spp repeated 5x
sam$sppid<-vec
sam <- sam[with(sam, order(sppid)),]

vec<-1:999
new<-data.frame(array(0, c(500, 2))); names(new)[1]<-"iteration"; names(new)[2]<-"null_diff"
for(j in 1:999){
null <- sam[sample(nrow(sam)),]; #randomly order dataframe to assign interactions
null$intid2<-rep(1:50,2)

sam2<-merge(null[,c(3:4)], fits3, by=c("sppid")) #add fits from stan model
sam3<-unique(sam2[,c("sppid","intid2","spptrends","sppSE")])

best<-data.frame(array(0, c(length(unique(sam3$intid2)), 2)))
names(best)[1]<-"intid"; names(best)[2]<-"diff"
Bgroups<-unique(sam3$intid2); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-sam3[sam3$intid2==b[i],]
best[i,1]<-b[i]
best[i,2]<-spp[1,c("spptrends")]-spp[2,c("spptrends")]
}
new[j,1]<-j
new[j,2]<-median(best$diff, na.rm=TRUE)
}


//level=intercept=a
//trend=slope=b
//sigma=variance