data {
//int<lower=0> Ndata;
 int Jspp; // number of data sets (spp)
  int Ndata[Jspp]; //number of data points per dataset   
real y[Jspp,Ndata];		
   real year[Jspp, Ndata]; //x vector, year of data point
   }
parameters {
  real alpha[Jspp];                          // the intercept for each species
  real beta[Jspp];                          // the slope for each species
  real sigma[Jspp];                // measurement error, noise, etc.
//one model per species
}
model {
  for (j in 1:Jspp){
  y[j]~normal(alpha[j]+beta[j]*year[j], sigma[j]);
	}
}
  