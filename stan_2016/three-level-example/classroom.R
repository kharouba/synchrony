library(magrittr)
library(rstan)
library(lme4)
library(dplyr)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/three-level-example/")
classroom <- read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/three-level-example/classroom.csv")

## Sort by class ID for easier handling in Stan
classroom <- arrange(classroom, classid, schoolid)

## Create a vector of school IDs where j-th element gives school ID for class ID j
schoolLookupVec <- unique(classroom[c("classid","schoolid")])[,"schoolid"]


#RANDOM INTERCEPT ONLY MODEL
## Combine as a stan dataset
dat <- with(classroom,
            list(Ni           = length(unique(childid)),
                 Nj           = length(unique(classid)),
                 Nk           = length(unique(schoolid)),
                 classid      = classid,
                 schoolid     = schoolid,
                 schoolLookup = schoolLookupVec,
                 mathgain     = mathgain))
                 
 
## Load Stan file
fileName <- "./multilevel3.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
cat(stan_code)

resStan <- stan(model_code = stan_code, data = dat,
                chains = 4, iter = 10000, warmup = 1000, thin = 10)



                
#RANDOM INTERCEPTS WITH FIXED EFFECT COVARIATES
desMat <- model.matrix(object = ~ 1 + mathkind + sex + minority + ses + housepov, data = classroom)
dat2 <- with(classroom,
            list(Ni           = length(unique(childid)),
                 Nj           = length(unique(classid)),
                 Nk           = length(unique(schoolid)),
                 p            = ncol(desMat),
                 desMat       = desMat,
                 classid      = classid,
                 schoolid     = schoolid,
                 schoolLookup = schoolLookupVec,
                 mathgain     = mathgain))
                 
fileName <- "./multilevel3.2.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
cat(stan_code)
resStan <- stan(model_code = stan_code, data = dat2, chains = 4, iter = 2000, warmup = 1000, thin = 10)

N           = length(unique(classroom$childid))
Nspp           = length(unique(classroom$classid))
Nstudy           = length(unique(classroom$schoolid))
species      = classroom$classid
studyid = schoolLookupVec
year<- classroom$mathgain

fit_simple<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/threelevelrandomintercept.stan", data=c("N","Nspp","Nstudy","species", "studyid","year"), iter=2000, chains=4)

fit_simple<-stan("multilevel3.2.stan", data=c("Ni","Nj","Nk","species", "studyid","year"), iter=2000, chains=4)
