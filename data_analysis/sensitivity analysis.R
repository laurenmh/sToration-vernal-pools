# MODEL SENSITIVITY ANALYSIS
# Bootstrapping method to check for model sensitivity
# Remove one pool at a time

# Load packages
library(dplyr)
library(rstan)
library(StanHeaders)
library(tidyverse)

# Get data
source("data_wrangling/prep data before modeling.R") 
------------------------------------------------------

### SET UP REFERENCE MODEL ###

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
      real <lower = 0, upper = 1> survival_LACO; // survival rate of LACO seeds in the seedbank
}
transformed parameters{
      matrix [n_pools, n_years-1] mu_LACO;// mean expected value of LACO at time t from a discrete BH model
      matrix [n_pools, n_years-3] int_LACO;// intermediate matrix of LACO at time t-1 estimated from values at t-2
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
              real germ_LACO;
              if (obs_EG[i,j-1] > 100)
                      germ_LACO = low_germ_LACO;
              else
                      germ_LACO = high_germ_LACO;
              if (obs_LACO[i,j-1] > 0){
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
    real germ_LACO;
    matrix[n_pools, n_years] stems_LACO;
    for(a in 1:n_pools){
        for(b in 1:1){
             germ_LACO = high_germ_LACO;
             stems_LACO[a,b] = mu_LACO[a,b] * germ_LACO;
             if(stems_LACO[a,b] > 0)
             obs_LACO[a,b] ~ poisson(stems_LACO[a,b]); 
        }
        for(b in 2:(n_years-1)){
             if (obs_EG[a,b-1] > 100)
                      germ_LACO = low_germ_LACO;
              else
                      germ_LACO = high_germ_LACO;
             stems_LACO[a,b] = mu_LACO[a,b] * germ_LACO;
             if(stems_LACO[a,b] > 0)
             obs_LACO[a,b] ~ poisson(stems_LACO[a,b]); // obs_LACO is from a poisson distribution of mu_LACO.
        }
    }
    lambda ~ normal(60,20); //get partially-informed priors from lit
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
    survival_LACO ~ beta(0.5,0.5);
}"

BH_ref_model <- stan_model(model_code = BH_ref_model_block)

### Data WRANGLING ###

ref_LACOcover_test <- rbind(ref_LACOcover, ref_LACOcover)
ref_sumEGcover_test <- rbind(ref_sumEGcover, ref_sumEGcover)
ref_ERVAcover_test<- rbind(ref_ERVAcover, ref_ERVAcover)
ref_sumNFcover_test <- rbind(ref_sumNFcover, ref_sumNFcover)
ref_n_pools_test <- ref_n_pools-1

### RUN THE MODEL ###

for(i in 1:9){
BH_ref_fit_ <- sampling(BH_ref_model,
                       data = list(n_pools = ref_n_pools_test,
                                   n_years = ref_n_years,
                                   obs_LACO = ref_LACOcover_test[i:(i+7),],
                                   obs_EG = ref_sumEGcover_test[i:(i+7),],
                                   obs_ERVA = ref_ERVAcover_test[i:(i+7),],
                                   obs_NF = ref_sumNFcover_test[i:(i+7),],
                                   low_germ_LACO = 0.2,
                                   high_germ_LACO = 0.7), 
                       iter= 1000)
assign(paste0("BH_ref_fit_", i), BH_ref_fit_)
}

### CHECK ESTIMATES ###

get_posterior_mean(BH_ref_fit_1, pars = c("alpha_LACO"))
get_posterior_mean(BH_ref_fit_2, pars = c("alpha_LACO"))
get_posterior_mean(BH_ref_fit_9, pars = c("alpha_LACO"))
get_posterior_mean(BH_ref_fit_1, pars = c("lambda"))
get_posterior_mean(BH_ref_fit_2, pars = c("lambda"))
get_posterior_mean(BH_ref_fit_9, pars = c("lambda"))

### EXTRACT ESTIMATES ###


### SET UP CONSTRUCTED MODEL ###

BH_model_block <- "
data{
    int n_pools; // number of pools
    int n_years; // number of years
    int obs_LACO [n_pools, n_years] ; // observed LACO density // this is an array of integers because the model is a discrete poisson model (see below)
    matrix[n_pools, n_years] obs_EG; // exotic grass density
    matrix[n_pools, n_years] obs_ERVA; // ERVA density
    matrix[n_pools, n_years] obs_NF; // non-native forb density
    int seeds_added [n_pools, 3]; // number of seeds added in 1999-2001
    real low_germ_LACO; // germination rate of LACO in normal years
    real high_germ_LACO; // germination rate of LACO in high litter years
}
parameters{
    vector<lower = 0>[n_years-1] lambda; // max growth rate of LACO in absence of competition
    vector<lower = 0, upper = 1>[n_years-1] alpha_LACO; // competition term for LACO-LACO
    vector<lower = 0, upper = 1>[n_years-1] alpha_EG; // competition term for LACO-exotic grass
    vector<lower = 0, upper = 1>[n_years-1] alpha_ERVA;// competition term for LACO-ERVA
    vector<lower = 0, upper = 1>[n_years-1] alpha_NF; // competition term for LACO-non-native forb
    real <lower = 0, upper = 1> survival_LACO; // survival rate of LACO seeds in the seedbank
}
transformed parameters{
    matrix [n_pools, n_years-1] mu_LACO;// mean expected value of seed LACO at time t from a discrete BH model
    matrix [n_pools, n_years-3] int_LACO;// intermediate matrix of seed LACO at time t-1 estimated from values at t-2
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
        for(j in 3:n_years-1){
            real germ_LACO;
            if (obs_EG[i,j-1] > 100)
                germ_LACO = low_germ_LACO;
            else
                germ_LACO = high_germ_LACO;
            if (obs_LACO[i,j-1] > 0){
                int_LACO[i,j-2] = (obs_LACO[i,j-2] * lambda[j-2])./(1 + obs_LACO[i,j-2] * alpha_LACO[j-2] + 
                                  obs_EG[i,j-2] * alpha_EG[j-2] + obs_ERVA[i,j-2] * alpha_ERVA[j-2] + obs_NF[i,j-2] * alpha_NF[j-2]) +
                                  survival_LACO * (1-germ_LACO) * obs_LACO[i,j-2] ./ germ_LACO;
                mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
                                obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
                                survival_LACO * (1-germ_LACO) * obs_LACO[i,j-1] ./ germ_LACO; 
            }
            else{
                int_LACO[i,j-2] = (obs_LACO[i,j-2] * lambda[j-2])./(1 + obs_LACO[i,j-2] * alpha_LACO[j-2] + 
                                  obs_EG[i,j-2] * alpha_EG[j-2] + obs_ERVA[i,j-2] * alpha_ERVA[j-2] + obs_NF[i,j-2] * alpha_NF[j-2]) +
                                  survival_LACO * (1-germ_LACO) * obs_LACO[i,j-2] ./ germ_LACO; // LACO at t-2
                mu_LACO[i,j] = (int_LACO[i,j-2] * lambda[j-1])./(1 + int_LACO[i,j-2] * alpha_LACO[j-1] + 
                                obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
                                survival_LACO * (1-germ_LACO) * int_LACO[i,j-2] ./ germ_LACO; // plugging in LACO t-2 for obs LACO t-1
            }
        }
    }
}
model{
    real germ_LACO;
    matrix[n_pools, n_years] stems_LACO;
    for(a in 1:n_pools){
        for(b in 1:1){
            germ_LACO = high_germ_LACO;
            stems_LACO[a,b] = mu_LACO[a,b] * germ_LACO;
            obs_LACO[a,b] ~ binomial(seeds_added[a,b], germ_LACO); //the first year's obs_LACO is the initial germination of seeds added in 1999
        }
        for(b in 2:3){
            if (obs_EG[a,b-1] > 100)
                germ_LACO = low_germ_LACO;
            else
                germ_LACO = high_germ_LACO;
            stems_LACO[a,b] = mu_LACO[a,b] * germ_LACO;
            if(stems_LACO[a,b] >0)
            obs_LACO[a,b] ~ poisson(stems_LACO[a,b] + seeds_added[a,b] * germ_LACO); //the second and third year's obs_LACO is the sum of germination of seeds added in 2000 and 2001 and previous year's population. 
        }
        for(b in 4:(n_years-1)){
            if (obs_EG[a,b-1] > 100)
                germ_LACO = low_germ_LACO;
            else
                germ_LACO = high_germ_LACO;
            stems_LACO[a,b] = mu_LACO[a,b] * germ_LACO;
            if(stems_LACO[a,b] > 0)
            obs_LACO[a,b] ~ poisson(stems_LACO[a,b]); //the rest of the year's obs_LACO is from a poisson distribution of mu_LACO. 
        }
    }
    lambda ~ normal(60,20); //partially informed prior from literature 
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
    survival_LACO ~ beta(0.5,0.5);//Jeffery's prior
}"

BH_model <- stan_model(model_code = BH_model_block)

### DATA WRANGLING ###

LACOcover_test <- rbind(LACOdens, LACOdens)
sumEGcover_test <- rbind(sumEGcover, sumEGcover)
ERVAcover_test<- rbind(ERVAdens, ERVAdens)
sumNFcover_test <- rbind(sumNFcover, sumNFcover)
n_pools_test <- n_pools-1
seedtrt_test <- rbind(seedtrt, seedtrt)

### RUN THE MODEL ###

for(i in 3:20){
BH_fit_ <- sampling(BH_model,
                   data = list(n_pools = n_pools_test,
                               n_years = 16,
                               obs_LACO = LACOcover_test[i:(i+140),1:16],
                               obs_EG = sumEGcover_test[i:(i+140),1:16],
                               obs_ERVA = ERVAcover_test[i:(i+140),1:16],
                               obs_NF = sumNFcover_test[i:(i+140),1:16],
                               seeds_added = seedtrt_test[i:(i+140),4:6],
                               low_germ_LACO = 0.2,
                               high_germ_LACO = 0.7), 
                   iter= 1000)
assign(paste0("BH_fit_", i), BH_fit_)
}

### CHECK ESTIMATES ###
get_posterior_mean(BH_fit_1, pars = c("alpha_LACO"))
get_posterior_mean(BH_fit_2, pars = c("alpha_LACO"))
get_posterior_mean(BH_fit_1, pars = c("lambda"))
get_posterior_mean(BH_fit_2, pars = c("lambda"))

