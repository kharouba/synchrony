 ### 29 April 2013 ###
### By Lizzie ###

## It's time to start the code for Heather's synchrony project ##
## Let's go! ##

options(stringsAsFactors=FALSE)
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis")

# get data
taxer <- read.csv("taxa.csv", header=TRUE)
#obsraw <- read.csv("obs.csv", header=TRUE)

# data clean up (more needed here!)
taxer$latbi <- paste(taxer$genus, taxer$species)
taxer$role <- tolower(taxer$role)
taxer$role[taxer$role=="pos"] <- "positive"
taxer$role[taxer$role=="neg"] <- "negative"
#taxer$interaction[taxer$interactio=="poll"] <- "pollination"

# the f(x)s I built require you to divide out trophic intxns
# from other interactions so that's what I do next

taxer.trophic <- subset(taxer, interaction != "mutualist" & interaction != "pollination" & interaction != "competition" & studyid!="HMK027" & studyid!="HMK043" & studyid!="HMK019" & studyid!="SET006" & studyid!="HMK048")

# interaction != "mutualist" & interaction != "pollination" : trophic! BUT keep statements in because otherwise crash- different coding needed for mutualism

taxer.mutual <- subset(taxer, interaction=="mutualist" | interaction == "pollination") #OR
taxer.comp <- subset(taxer, interaction=="competition")

######
## f(x) to build all possible POSITIVE same-same interactions
######
# note:
# (1) make sure you have a latbi column
# (2) make sure theses really are non-trophic interactions!
buildintxns.nontrophic <- function(dater, idcol) {
intxndf.nt <- lapply(unique(dater[[idcol]]), function(uniquesite){
	subby <- subset(dater, dater[[idcol]]==uniquesite)
	combos <- combn(subby$latbi, 2, simplify=TRUE)
	intxnframe <- data.frame(
	studyid = uniquesite,                     
	poslatbi1 = combos[1,],
	poslatbi2 = combos[2,]
	)
})
do.call("rbind", intxndf.nt)
}

######
## f(x) to build all possible NEGATIVE same-same interactions
######
# note:
# (1) make sure you have a latbi column
# (2) make sure theses really are non-trophic interactions!
buildintxns.comp <- function(dater, idcol) {
intxndf.nt <- lapply(unique(dater[[idcol]]), function(uniquesite){
	subby <- subset(dater, dater[[idcol]]==uniquesite)
	combos <- combn(subby$latbi, 2, simplify=TRUE)
	intxnframe <- data.frame(
	studyid = uniquesite,                     
	neglatbi1 = combos[1,],
	neglatbi2 = combos[2,]
	)
})
do.call("rbind", intxndf.nt)
}

######
## f(x) to build all possible positive-negative interactions
######
# note:
# (1) make sure you have a 'role' column with positive or negative as values
# (2) make sure you have a latbi column

buildintxns <- function(dater, idcol) {
intxndf <- lapply(unique(dater[[idcol]]), function(uniquesite){
	subby <- subset(dater, dater[[idcol]]==uniquesite)
	allgood <- subby[subby$role == "positive",]
	ohdear <- subby[subby$role == "negative",]
	combinate <- expand.grid(ohdear$latbi, allgood$latbi)
	badgood <- data.frame(
	studyid = uniquesite,
	neglatbi = as.character(combinate$Var1),
	poslatbi = as.character(combinate$Var2)
)
})
do.call("rbind", intxndf)
}

## run the f(x)s and say 'ahh':
#buildintxns.nontrophic(taxer.mutual, "studyid")
#buildintxns(taxer.trophic, "studyid")

# in reality though you will more often run them and assign them to an object, like this:
intxn_trophic <- buildintxns(taxer.trophic, "studyid")
names(intxn_trophic)[2]<-"latbi"
intxn_trophic2<-merge(intxn_trophic, unique(taxer.trophic[,c("latbi","interaction","taxatype")]),by="latbi")
names(intxn_trophic2)[1]<-"spp1"; names(intxn_trophic2)[3]<-"spp2"

intxn_mut<- buildintxns.nontrophic(taxer.mutual, "studyid")
names(intxn_mut)[2]<-"latbi"
intxn_mut2<-merge(intxn_mut, unique(taxer.mutual[,c("latbi","interaction","taxatype")]),by="latbi")
names(intxn_mut2)[1]<-"spp1"; names(intxn_mut2)[3]<-"spp2"

intxn_comp<- buildintxns.comp(taxer.comp, "studyid")
names(intxn_comp)[2]<-"latbi"
intxn_comp2<-merge(intxn_comp, unique(taxer.comp[,c("latbi","interaction","taxatype")]),by="latbi")
names(intxn_comp2)[1]<-"spp1"; names(intxn_comp2)[3]<-"spp2"

intxn<-rbind(intxn_trophic2,intxn_comp2, intxn_mut2)
intxn<-intxn[,c("studyid","spp1","spp2","interaction","taxatype")]

#add HMK048 manually because 4-trophic web (9 interactions)
new<-data.frame(array(0, c(9, 5)))
names(new)[1] = "studyid"; names(new)[2] = "spp1"; names(new)[3] = "spp2"; names(new)[4] = "interaction"; names(new)[5]="taxatype"
new[,1]<-"HMK048"
new[,2]<-c("Ficedula hypoleuca","Parus ater", "Parus2 caeruleus","Parus4 major", "Caterpillar2 spp.", "Caterpillar2 spp.", "Caterpillar2 spp.", "Caterpillar2 spp.", "Quercus3 robur")
new[,3]<-c("Accipiter nisus","Accipiter nisus", "Accipiter nisus","Accipiter nisus", "Ficedula hypoleuca","Parus ater", "Parus2 caeruleus","Parus4 major", "Caterpillar2 spp.")
new[,4]<-c("predation","predation","predation","predation", "predation","predation","predation","predation", "herbivory")
new[,5]<-c("spp","spp","spp","spp", "genus","genus","genus","genus", "genus")

intxn2<-rbind(intxn, new)
intxn<-intxn2


#add HMK027 manually because 3-trophic web (4 interactions)
new<-data.frame(array(0, c(3, 5)))
names(new)[1] = "studyid"; names(new)[2] = "spp1"; names(new)[3] = "spp2"; names(new)[4] = "interaction"; names(new)[5] = "taxatype";
new[,1]<-"HMK027"
new[,2]<-c("Caterpillar spp.", "Caterpillar spp.", "Quercus2 robur")
new[,3]<-c("Parus1 major", "Ficedula1 albicollis", "Caterpillar spp.")
new[,4]<-c("predation","predation","herbivory")
new[,5]<-c("genus","genus","genus")

intxn2<-rbind(intxn, new)
intxn<-intxn2

#add HMK043 manually because 2 diff interactions
new<-data.frame(array(0, c(2, 5)))
names(new)[1] = "studyid"; names(new)[2] = "spp1"; names(new)[3] = "spp2"; names(new)[4] = "interaction";names(new)[5] = "taxatype";
new[,1]<-"HMK043"
new[,2]<-c("Copepod2 spp.", "Pleurobrachia pileus")
new[,3]<-c("Pleurobrachia pileus","Beroe gracilis")
new[,4]<-c("predation","predation")
new[,5]<-c("genus","spp")

intxn2<-rbind(intxn, new)
intxn<-intxn2

#add HMK019 manually because 2 diff interactions
new<-data.frame(array(0, c(2, 5)))
names(new)[1] = "studyid"; names(new)[2] = "spp1"; names(new)[3] = "spp2"; names(new)[4] = "interaction";names(new)[5] = "taxatype";
new[,1]<-"HMK019"
new[,2]<-c("Phytoplankton1 spp.", "Daphnia3 spp.")
new[,3]<-c("Daphnia3 spp.","Perca fluviatillis")
new[,4]<-c("herbivory","predation")
new[,5]<-c("genus","genus")

intxn2<-rbind(intxn, new)
intxn<-intxn2

#add SET006 manually because 2 diff interactions
new<-data.frame(array(0, c(2, 5)))
names(new)[1] = "studyid"; names(new)[2] = "spp1"; names(new)[3] = "spp2"; names(new)[4] = "interaction";names(new)[5] = "taxatype";
new[,1]<-"SET006"
new[,2]<-c("Quercus robur", "Operophtera brumata")
new[,3]<-c("Operophtera brumata","Parus major")
new[,4]<-c("herbivory","predation")
new[,5]<-c("spp","spp")

intxn2<-rbind(intxn, new)
intxn<-intxn2

#Parus major <> caterpillar
#Ficedula albicollis <> caterpillar
#caterpillar <> Quercus robur






