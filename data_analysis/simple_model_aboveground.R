# Model LACO without storage effect or environmental fluctuations 

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

# simulated data
sim_n_pools <- 3 #number of pools
sim_n_years <- 18 #years of data
set.seed(124) #this helps create simulated values that are reproducible
sim_obs_EG <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 150), ncol=sim_n_years) #simulate poisson distributed observed EG density with mean of 150
set.seed(123)
sim_obs_ERVA <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 10), ncol=sim_n_years)
set.seed(122)
sim_obs_NF <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 50), ncol=sim_n_years)

sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years)

bh.sim <- function(init, EG, ERVA, NF, aii, a1, a2, a3, lambda){
  sim_obs_LACO[,1]<- init
  for(i in 1:nrow(sim_obs_LACO)){
    for(j in 2:ncol(sim_obs_LACO)){
      mu <- lambda*sim_obs_LACO[i,j-1]/(1+sim_obs_LACO[i,j-1]*aii+EG[i,j-1]*a1+ERVA[i,j-1]*a2+NF[i,j-1]*a3)
      sim_obs_LACO[i,j] <- rpois(1, lambda=mu)
    }
  }
 return(sim_obs_LACO)
}

sim_obs_LACO <- bh.sim(init = 100,
                       EG = sim_obs_EG,
                       ERVA = sim_obs_ERVA,
                       NF = sim_obs_NF,
                       aii = 1,
                       a1 = 0.1,
                       a2 = 0.2, 
                       a3 = 0.5,
                       lambda = 40)
#list "true" lambda and alpha parameter values here. start with constant parameters. 
#check that the model estimates parameters close to these values.

hist(sim_obs_LACO)

# Stan model
BH_model <- "
data{
    int n_pools; // number of pools
    int n_years; // number of years
    matrix[n_pools, n_years] obs_LACO; // observed LACO density
    matrix[n_pools, n_years] obs_EG; // exotic grass density
    matrix[n_pools, n_years] obs_ERVA; // ERVA density
    matrix[n_pools, n_years] obs_NF; // non-native forb density
}
parameters{
    real <lower = 0, upper = 1> lambda; // max growth rate of LACO in absence of competition
    real <lower = 0, upper = 1> alpha_LACO; // competition term for LACO-LACO
    real <lower = 0, upper = 1> alpha_EG; // competition term for LACO-exotic grass
    real <lower = 0, upper = 1> alpha_ERVA;// competition term for LACO-ERVA
    real <lower = 0, upper = 1> alpha_NF; // competition term for LACO-non-native forb
    matrix [n_pools, n_years] N_LACO ;// total population (observed + unobserved) of LACO at time t
}
transformed parameters{
    matrix[n_pools, n_years] mu_LACO;// total population (observed + unobserved) of LACO at time t
    for(i in 1:n_pools){
      for(j in 2:n_years){
          mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda)./(1 + obs_LACO[i,j-1] * alpha_LACO + 
          obs_EG[i,j-1] * alpha_EG + obs_ERVA[i,j-1] * alpha_ERVA + obs_NF[i,j-1] * alpha_NF);
      }
    }
}
model{
    for(i in 1:n_pools){
      for(j in 1:1){
          N_LACO[i,j] = obs_LACO[i,1]; #SYNTAX ERROR here - look up how to set up the first column
      }
      for(j in 2:n_years){
          N_LACO[i,j] ~ poisson(mu_LACO[i,j]); #ERROR here too - poisson only for int but N_LACO is matrix of real numbers.
      }
    }
    lambda ~ normal(122.2561,83.95284); //get partially-informed priors from lit
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
}"

#run the model with simulated data
stan_fit <- stan(model_code = BH_model,
                 data = list(n_pools = sim_n_pools,
                             n_years = sim_n_years,
                             obs_LACO = sim_obs_LACO,
                             obs_EG = sim_obs_EG,
                             obs_ERVA = sim_obs_ERVA,
                             obs_NF = sim_obs_NF), 
                 iter= 1000)

#check if there is enough iteration
stan_trace(stan_fit, pars = c("lambda"))

#mean posterior estimates of parameters
get_posterior_mean(stan_fit, pars = c("lambda"))

#zoom into posterior distribution of parameters
plot(stan_fit, pars = c("lambda", "alpha_LACO", "s_LACO"))


