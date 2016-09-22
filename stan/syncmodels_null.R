## Started 21 August 2015 ##
## By Lizzie ##

## Trying to do a quick null model from Stan output ##
## Based off current code from syncmodels_PPchecks.R ##

setwd("~/Documents/git/projects/trophsynch/synchrony")
source("syncmodels.R")

# look at output 
goo <- extract(fit.notype)
hist(goo$mu_b) # example!

# extract means for now #
sigma_y <- mean(goo$sigma_y) # from stan output (which I named 'goo' for now)
sigma_a <- mean(goo$sigma_a) # from stan output 
sigma_b <- mean(goo$sigma_b) # from stan output 
mu_b <- mean(goo$mu_b) # from stan output 
mu_a <- mean(goo$mu_a) # from stan output 

a <- rnorm(J, mean=mu_a, sd=sigma_a) # alert! perhaps should not set sd to sigma exactly?
b <- rnorm(J, mean=mu_b, sd=sigma_b) # alert! perhaps should not set sd to sigma exactly?

# Create the data using new a and b for each of 71 species
m <- matrix(0, ncol = 20, nrow = nrow(rawlong.nodups))
ypreds <- data.frame(m)
ypred <- length(N)
slope<-length(J)
for (x in c(1:ncol(ypreds))){
    for (n in 1:N){
        s <- species[n]
        ypred[n] <- a[s] + b[s]*year[n]
        slope[s]<-b[s]
    }
    y <- rnorm(N, ypred, sigma_y)
    ypreds[,x] <- y
}
quicknull <- cbind(rawlong.nodups, ypreds)

best<-data.frame(array(0, c(1740, 5)))
names(best)[1]<-"iteration"; names(best)[2]<-"sppid"; names(best)[3]<-"intercepts";
names(best)[4]<-"spptrends"; names(best)[5]<-"sppSE"
mat<-quicknull[,12:15] #31
rowcount<-1
asdf<-rowcount+(length(sppid)-1)
for (x in c(1:ncol(mat))){
y <- mat[,x]
fit <- stan("synchrony1_notype_hk.stan", data=c("N","J","y","species","year", "nVars","Imat"), iter=2000, chains=3)
zz<-summary(fit)
intercepts <- as.vector(zz[[1]][,1])[1:87]
spptrends<- as.vector(zz[[1]][,1])[88:174]
sppSE<- as.vector(zz[[1]][,2])[88:174]
sppid<-1:87
best[rowcount:asdf,1]<-x
best[rowcount:asdf, 2]<-sppid
best[rowcount:asdf, 3]<-intercepts
best[rowcount:asdf, 4]<-spptrends
best[rowcount:asdf, 5]<-sppSE
rowcount<-rowcount+length(sppid)
asdf<-asdf+length(sppid)
}
