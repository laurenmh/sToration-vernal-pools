# Model LACO without storage effect or environmental fluctuations 

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

# simulated data
sim_n_pools <- 250 #number of pools
sim_n_years <- 3 #years of data
set.seed(125) #this helps create simulated values that are reproducible
sim_obs_EG <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 100), ncol=sim_n_years) #simulate poisson distributed observed EG density with mean of 150
set.seed(126)
sim_obs_ERVA <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 30), ncol=sim_n_years)
set.seed(127)
sim_obs_NF <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 50), ncol=sim_n_years)

sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years)
sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years)

bh.sim <- function(n_pools, init, EG, ERVA, NF, aii, a1, a2, a3, lambda){
  sim_obs_LACO[,1]<- rbinom(n_pools,100,0.7)
  sim_mu[,1]<- init
  for(i in 1:nrow(sim_obs_LACO)){
    for(j in 2:ncol(sim_obs_LACO)){
      sim_mu[i,j] <- lambda*sim_obs_LACO[i,j-1]/(1+sim_obs_LACO[i,j-1]*aii+EG[i,j-1]*a1+ERVA[i,j-1]*a2+NF[i,j-1]*a3)
      sim_obs_LACO[i,j] <- rpois(1, lambda=sim_mu[i,j])
    }
  }
 return(sim_obs_LACO)
}

sim_obs_LACO <- bh.sim(n_pools = sim_n_pools,
                       init = 100,
                       EG = sim_obs_EG,
                       ERVA = sim_obs_ERVA,
                       NF = sim_obs_NF,
                       aii = 0.9,
                       a1 = 0.3,
                       a2 = 0.2, 
                       a3 = 0.1,
                       lambda = 40)
#list "true" lambda and alpha parameter values here. start with constant parameters. 
#check that the model estimates parameters close to these values.

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
    real <lower = 0> lambda; // max growth rate of LACO in absence of competition
    real <lower = 0, upper = 1> alpha_LACO; // competition term for LACO-LACO
    real <lower = 0, upper = 1> alpha_EG; // competition term for LACO-exotic grass
    real <lower = 0, upper = 1> alpha_ERVA;// competition term for LACO-ERVA
    real <lower = 0, upper = 1> alpha_NF; // competition term for LACO-non-native forb
    real <lower = 0, upper = 1> germ_LACO; // germination rate of LACO at t=1
    real <lower = 0, upper = 0.1> sigma; // error term for expected value of LACO
}
transformed parameters{
    matrix [n_pools, n_years] mu_LACO;// expected value of LACO at time t
    for(i in 1:n_pools){  
        for(j in 1:1){
            mu_LACO[i,j] = 100; // this column does not get used in the model
        }
        for(j in 2:n_years){
            mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda)./(1 + obs_LACO[i,j-1] * alpha_LACO + 
            obs_EG[i,j-1] * alpha_EG + obs_ERVA[i,j-1] * alpha_ERVA + obs_NF[i,j-1] * alpha_NF); // Beverton Holt model
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
    lambda ~ normal(40,1); //get partially-informed priors from lit
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
get_posterior_mean(BH_fit, pars = c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF", "germ_LACO"))

#zoom into posterior distribution of parameters
plot(BH_fit, pars = c("alpha_NF", "germ_LACO"))


