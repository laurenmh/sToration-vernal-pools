
library(rstan)

BH_model_block <- "
data{
  int y; 
  int a; 
}

parameters{
  real <lower = 0, upper = 1> c;
  real <lower = 0, upper = 1> b;
}

model{
  y ~ binomial(a,b)+ poisson(c);
}
"
BH_model <- stan_model(model_code = BH_model_block)
BH_fit <- sampling(BH_model,
                   data = list(y = 5,
                               a = 2), 
                   iter= 1000)
