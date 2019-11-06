# Lina's first attempt at modeling LACO population dynamics using the Beverton-Holt model
# Feel free to make changes on my code if you see fit 

# load packages
library(dplyr)
library(rstan)

# source the data
source("data_compiling/compile_composition.R")
LACOseeds <- read.csv(paste(datpath, "/Seed count data/Vernal Pool LACO seed count data.csv", sep="")) 

# get distribution the seed production of LACO to inform prior for lambda
restoredLACO <- LACOseeds %>%
  filter(Pool.Type == "Restored")
hist(restoredLACO$X..Potentially.Viable)
mean(restoredLACO$X..Potentially.Viable)
sd(restoredLACO$X..Potentially.Viable)

# Stan model
BH_model <- stan_model(model_code="
data{
    int n; // number of observations
    int n_pools; // number of pools
    vector[n_pools] pool_id; // pool id
    int n_years; // number of years
    vector[n_years] year; // year
    int obs_LACO [n_pools, n_years]; // LACOdens
    int obs_EG [n_pools, n_years]; // exotic grass density
    int obs_ERVA [n_pools, n_years]; // ERVAdens
    int obs_NF [n_pools, n_years]; // non-native forb density
    vector[n_years] g_LACO; // growth rate of LACO
    real s_LACO; // survival rate of LACO
}
parameters{
    matrix<lower = 0>[n_pools, n_years] lambda; // max growth rate of LACO in absence of competition
    matrix<lower = 0>[n_pools, n_years] alpha_LACO; // competition term for LACO-LACO
    matrix<lower = 0>[n_pools, n_years] alpha_EG; // competition term for LACO-exotic grass
    matrix<lower = 0>[n_pools, n_years] alpha_ERVA;// competition term for LACO-ERVA
    matrix<lower = 0>[n_pools, n_years] alpha_NF; // competition term for LACO-non-native forb
}
model{
    for(j in 1:n_pools){
        obs_LACO[j,] ~ poisson(100, g_LACO); // this g_LACO[1] is just first year's growth rate
    }
    for(i in 1:n.years){
        obs_LACO[,i] ~ (g_LACO[i-1] * obs_LACO[,i-1] * lambda[,i-1])/(1 + obs_LACO[,i-1] * alpha_LACO[,i-1] + 
        obs_EG[,i-1] * alpha_EG[,i-1] + obs_ERVA[,i-1] * alpha_ERVA[,i-1] + obs_NF[,i-1] * alpha_NF[,i-1]) + 
        s_LACO * (1 - g_LACO[i-1]) * obs_LACO[,i-1]);
    }
    lambda ~ normal(122.2561,83.95284);
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
}")

