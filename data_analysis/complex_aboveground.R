# Model LACO without storage effect but with environmental fluctuations 

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

# simulated data
sim_n_pools <- 256 #number of pools
sim_n_years <- 18 #years of data

set.seed(124) #this helps create simulated values that are reproducible
sim_obs_EG <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 100), ncol=sim_n_years) #simulate poisson distributed observed EG density with mean of 150
set.seed(123)
sim_obs_ERVA <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 30), ncol=sim_n_years)
set.seed(122)
sim_obs_NF <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 50), ncol=sim_n_years)

sim_aii <- as.vector(rbeta(sim_n_years, 3, 1))
sim_a1 <- as.vector(rbeta(sim_n_years, 1, 2))
sim_a2 <- as.vector(rbeta(sim_n_years, 2, 5))
sim_a3 <- as.vector(rbeta(sim_n_years, 1, 5))
sim_lambda <- as.vector(rpois(sim_n_years, 40))

sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years)
sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years)

bh.sim <- function(n_pools, init, EG, ERVA, NF, aii, a1, a2, a3, lambda){
  sim_obs_LACO[,1]<- rbinom(n_pools,100,0.8)
  sim_mu[,1]<- init
  for(i in 1:nrow(sim_obs_LACO)){
    for(j in 2:ncol(sim_obs_LACO)){
      sim_mu[i,j] <- lambda[j-1]*sim_obs_LACO[i,j-1]/(1+sim_obs_LACO[i,j-1]*aii[j-1]+EG[i,j-1]*a1[j-1]+ERVA[i,j-1]*a2[j-1]+NF[i,j-1]*a3[j-1])
      sim_obs_LACO[i,j] <- rpois(1, lambda=sim_mu[i,j])
    }
  }
  return(sim_obs_LACO)
}

#list "true" lambda and alpha parameter values here. 
#check that the model estimates parameters close to these values.
sim_obs_LACO <- bh.sim(n_pools = sim_n_pools,
                       init = 100,
                       EG = sim_obs_EG,
                       ERVA = sim_obs_ERVA,
                       NF = sim_obs_NF,
                       aii = sim_aii,
                       a1 = sim_a1,
                       a2 = sim_a2, 
                       a3 = sim_a3,
                       lambda = sim_lambda)

hist(sim_obs_LACO)

# Stan model
BH_model_block <- "
data{
    int n_pools; // number of pools
    int n_years; // number of years
    int obs_LACO [n_pools, n_years] ; // observed LACO density // this is an array of integers because the model is a discrete poisson model (see below)
    matrix[n_pools, n_years] obs_EG; // exotic grass density
    matrix[n_pools, n_years] obs_ERVA; // ERVA density
    matrix[n_pools, n_years] obs_NF; // non-native forb density
}
parameters{
    vector<lower = 0>[n_years] lambda; // max growth rate of LACO in absence of competition
    vector<lower = 0, upper = 1>[n_years] alpha_LACO; // competition term for LACO-LACO
    vector<lower = 0, upper = 1>[n_years] alpha_EG; // competition term for LACO-exotic grass
    vector<lower = 0, upper = 1>[n_years] alpha_ERVA;// competition term for LACO-ERVA
    vector<lower = 0, upper = 1>[n_years] alpha_NF; // competition term for LACO-non-native forb
    real <lower = 0, upper = 1> germ_LACO; // germination rate of LACO on t=1
    real <lower = 0, upper = 0.1> sigma; // error term for expected value of LACO
}
transformed parameters{
    matrix [n_pools, n_years] mu_LACO;// expected value of LACO at time t
    for(i in 1:n_pools){  
          for(j in 1:1){
                mu_LACO[i,j] = 100;
          }
          for(j in 2:n_years){
                mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
                obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]); // Beverton Holt model
          }
    }
}
model{
    for(i in 1:n_pools){
          for(j in 1:1){
                obs_LACO[i,j] ~ binomial(100, germ_LACO); //the first year's obs_LACO is the initial germination of 100 seeds
          }
          for(j in 2:n_years){
                obs_LACO[i,j] ~ poisson(mu_LACO[i,j] + sigma); //the rest of the year's obs_LACO is from a poisson distribution of mu_LACO. 
          }
    }
    lambda ~ normal(40,10); //get partially-informed priors from lit
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
    sigma ~ normal(0,0.01);
}"

#run the model with simulated data
BH_model <- stan_model(model_code = BH_model_block)
BH_fit <- sampling(BH_model,
                   data = list(n_pools = sim_n_pools,
                               n_years = sim_n_years,
                               obs_LACO = sim_obs_LACO,
                               obs_EG = sim_obs_EG,
                               obs_ERVA = sim_obs_ERVA,
                               obs_NF = sim_obs_NF), 
                   iter= 1000)

#check if there is enough iteration
stan_trace(BH_fit, pars = c("lambda"))

#mean posterior estimates of parameters
get_posterior_mean(BH_fit, pars = c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF"))

#zoom into posterior distribution of parameters
plot(BH_fit, pars = c("alpha_NF"))


