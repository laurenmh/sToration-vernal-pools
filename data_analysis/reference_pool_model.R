# POPULATION MODEL FOR LACO IN REFERENCE POOLS
# Model LACO stem counts with storage effect and environmental fluctuations 
# Recursive modeling of seedbank to account for years where aboveground LACO is zero but seedbank persists to next year
# Germination rate of LACO depends on the previous year's exotic grass cover

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

### CREATE SIMULATED DATA ###

sim_n_pools <- 152 #number of pools
sim_n_years <- 7 #years of data

set.seed(124) #this helps create simulated values that are reproducible
sim_ref_EG <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 35), ncol=sim_n_years) #simulate exotic grass(EG) cover 
set.seed(123)
sim_ref_ERVA <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 10), ncol=sim_n_years) #simulate ERVA cover 
set.seed(122)
sim_ref_NF <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 25), ncol=sim_n_years) #simulate native forb (NF) cover 

sim_ref_aii <- as.vector(rbeta(sim_n_years, 3, 1)) #simulate alpha_LACO from beta distribution
sim_ref_a1 <- as.vector(rbeta(sim_n_years, 1, 2)) #simulate alpha_EG from beta distribution
sim_ref_a2 <- as.vector(rbeta(sim_n_years, 2, 5)) #simulate alpha_ERVA from beta distribution
sim_ref_a3 <- as.vector(rbeta(sim_n_years, 1, 5)) #simulate alpha_NF from beta distribution
sim_ref_lambda <- as.vector(rpois(sim_n_years, 40)) #simulate lambda from poisson distribution

sim_ref_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of LACO stem counts
sim_ref_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of mean LACO stem counts

bh.formula <- function(sim_ref_LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g){
  sim_ref_LACO*lambda/(1+sim_ref_LACO*aii+EG*a1+ERVA*a2+NF*a3)+s*(1-g)*sim_ref_LACO/g
} #this is the modified Beverton-Holt model we'll use for LACO stem counts

#now simulate LACO stem counts
bh.sim <- function(n_pools, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
  for(i in 1:nrow(sim_ref_mu)){
    for(j in 1:1){
      sim_ref_mu[i,j] <- 100
      sim_ref_LACO[i,j] <- rpois(1, lambda = (sim_ref_mu[i,j]))
    }
    for(j in 2:3){
      if (EG[i,j-1]> 100){
        g = glow
      }
      else{g = g}
      sim_ref_mu[i,j] <- bh.formula(sim_ref_LACO = sim_ref_LACO[i,j-1],
                                EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                lambda = lambda[j-1], s = s, g = g)
      sim_ref_LACO[i,j] <- rpois(1, lambda = (sim_ref_mu[i,j]))
    }
    for(j in 4:ncol(sim_ref_mu)){
      if (EG[i,j-1]> 100){
        g = glow
      }
      else{g = g}
      if (sim_ref_LACO[i,j-1] > 0){
        sim_ref_mu[i,j] <- bh.formula(sim_ref_LACO = sim_ref_LACO[i,j-1],
                                  EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                  aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                  lambda = lambda[j-1], s = s, g = g)
      }
      else {
        sim_ref_mu[i,j] <- bh.formula(sim_ref_LACO = sim_ref_LACO[i,j-2]*lambda[j-2]/(1+sim_ref_LACO[i,j-2]*aii[j-2]+EG[i,j-2]*a1[j-2]+ERVA[i,j-2]*a2[j-2]+NF[i,j-2]*a3[j-2])+s*(1-g)*sim_obs_LACO[i,j-2]/g,
                                  EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                  aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                  lambda = lambda[j-1], s = s, g = g)
      }
      sim_ref_LACO[i,j] <- rpois(1, lambda = sim_ref_mu[i,j])
    }
  }
  return(sim_ref_LACO)
}

# List "true" lambda and alpha parameter values here. Start with constant parameters. 
# After running the model, check that the model outputs are close to these values.
sim_ref_LACO <- bh.sim(n_pools = sim_n_pools,
                       EG = sim_ref_EG,
                       ERVA = sim_ref_ERVA,
                       NF = sim_ref_NF,
                       aii = sim_ref_aii,
                       a1 = sim_ref_a1,
                       a2 = sim_ref_a2, 
                       a3 = sim_ref_a3,
                       lambda = sim_ref_lambda,
                       s = 0.2,
                       g = 0.7,
                       glow = 0.2)

hist(sim_ref_LACO) # Check distribution of simulated LACO

### CREATE A STAN MODEL ###

BH_ref_model_block <- "
data{
      int n_pools; // number of pools
      int n_years; // number of years
      int obs_LACO [n_pools, n_years] ; // observed LACO density // this is an array of integers because the model is a discrete poisson model (see below)
      matrix[n_pools, n_years] obs_EG; // exotic grass density
      matrix[n_pools, n_years] obs_ERVA; // ERVA density
      matrix[n_pools, n_years] obs_NF; // non-native forb density
      real low_germ_LACO; // germination rate of LACO in normal years
      real high_germ_LACO; // germination rate of LACO in high litter years
}
parameters{
      vector<lower = 0>[n_years-1] lambda; // max growth rate of LACO in absence of competition
      vector<lower = 0, upper = 1>[n_years-1] alpha_LACO; // competition term for LACO-LACO
      vector<lower = 0, upper = 1>[n_years-1] alpha_EG; // competition term for LACO-exotic grass
      vector<lower = 0, upper = 1>[n_years-1] alpha_ERVA;// competition term for LACO-ERVA
      vector<lower = 0, upper = 1>[n_years-1] alpha_NF; // competition term for LACO-non-native forb
      real <lower = 0, upper = 0.1> sigma; // error term for expected value of LACO
      real <lower = 0, upper = 1> survival_LACO; // survival rate of LACO seeds in the seedbank
}
transformed parameters{
      matrix [n_pools, n_years-1] mu_LACO;// mean expected value of LACO at time t from a discrete BH model
      for(i in 1:n_pools){  
          for(j in 1:1){
              mu_LACO[i,j] = 100;
          }
          for(j in 2:2){
              real germ_LACO;
              if (obs_EG[i,j-1] > 100)
                      germ_LACO = low_germ_LACO;
              else
                      germ_LACO = high_germ_LACO;
              mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
                              obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
                              survival_LACO * (1-germ_LACO) * obs_LACO[i,j-1] ./ germ_LACO; // modified Beverton-Holt model
          }
          for(j in 3:(n_years-1)){
              matrix [n_pools, n_years-3] int_LACO;// intermediate matrix of LACO at time t-1 estimated from values at t-2
              real germ_LACO;
              if (obs_EG[i,j-1] > 100)
                      germ_LACO = low_germ_LACO;
              else
                      germ_LACO = high_germ_LACO;
              if (obs_EG[i,j-1] > 0){
                      int_LACO[i,j-2] = 100;
                      mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
                                      obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
                                      survival_LACO * (1-germ_LACO) * obs_LACO[i,j-1] ./ germ_LACO; // modified Beverton-Holt model
              }
              else{
                      int_LACO[i,j-2] = (obs_LACO[i,j-2] * lambda[j-2])./(1 + obs_LACO[i,j-2] * alpha_LACO[j-2] + 
                                        obs_EG[i,j-2] * alpha_EG[j-2] + obs_ERVA[i,j-2] * alpha_ERVA[j-2] + obs_NF[i,j-2] * alpha_NF[j-2]) +
                                        survival_LACO * (1-germ_LACO) * obs_LACO[i,j-2] ./ germ_LACO; // modified Beverton-Holt model
                      mu_LACO[i,j] = (int_LACO[i,j-2] * lambda[j-1])./(1 + int_LACO[i,j-2] * alpha_LACO[j-1] + 
                                      obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
                                      survival_LACO * (1-germ_LACO) * int_LACO[i,j-2] ./ germ_LACO; // modified Beverton-Holt model
              }
          }
      }
}
model{
    for(a in 1:n_pools){
        for(b in 1:n_years-1){
             obs_LACO[a,b] ~ poisson(mu_LACO[a,b] + sigma); // obs_LACO is from a poisson distribution of mu_LACO.
        }
    }
    lambda ~ normal(40,10); //get partially-informed priors from lit
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
    sigma ~ normal(0,0.01);
    survival_LACO ~ normal(0,1);
}"

BH_ref_model <- stan_model(model_code = BH_ref_model_block)

### RUN THE MODEL ###

## Option 1: run with simulated data
BH_ref_fit <- sampling(BH_ref_model,
                   data = list(n_pools = sim_n_pools,
                               n_years = sim_n_years,
                               obs_LACO = sim_ref_LACO,
                               obs_EG = sim_ref_EG,
                               obs_ERVA = sim_ref_ERVA,
                               obs_NF = sim_ref_NF,
                               low_germ_LACO = 0.2,
                               high_germ_LACO = 0.7), 
                   iter= 1000)

## Option 2: run with real data 
## See data prep file before running this
BH_fit <- sampling(BH_ref_model,
                   data = list(n_pools = n_pools,
                               n_years = n_years,
                               obs_LACO = LACOdens,
                               obs_EG = sumEGcover,
                               obs_ERVA = ERVAdens,
                               obs_NF = sumNFcover,
                               low_germ_LACO = 0.2,
                               high_germ_LACO = 0.8), 
                   iter= 1000)