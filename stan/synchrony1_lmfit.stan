data {
  int<lower=0> N;						
  vector [N] y; // y vector
  vector [N] year; //x vector
  }
parameters {
  real alpha; //intercept
  real beta; //slope
  real<lower=0> sigma; //regression variance, dispersion
}
//one big model: complete pooling
model {
//VECTORIZED model structure:
	y~normal(alpha+beta*year, sigma);
  }

//In this case, it is assumed that each data point was drawn from a normal distribution with mean equal to the regression model and standard deviation equal to sigma

  
//priors
//alpha ~ normal(0,100);
//beta ~ normal(0,100);
//sigma ~ uniform(0,1000);


  //UN-VECTORIZED model structure:
  // for (n in 1:N){
  //y[n]~normal(alpha+beta*year[n], sigma);
  // }
