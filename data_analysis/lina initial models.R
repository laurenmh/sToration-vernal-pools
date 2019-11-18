# Lina's first attempt at modeling LACO population dynamics using the Beverton-Holt model
# Feel free to make changes on my code if you see fit 

# load packages
library(dplyr)
library(rstan)

# source the data
source("data_compiling/compile_composition.R")
LACOseeds <- read.csv(paste(datpath, "/Seed count data/Vernal Pool LACO seed count data.csv", sep="")) 

# this is for later: seed production of LACO to check the model
restoredLACO <- LACOseeds %>%
  filter(Pool.Type == "Restored")
hist(restoredLACO$X..Potentially.Viable)
mean(restoredLACO$X..Potentially.Viable)
sd(restoredLACO$X..Potentially.Viable)

# set up data



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
    vector[n_years] g_LACO; // per seed germination rate of LACO
    real s_LACO; // seed survival of LACO seedbank
}
parameters{
    vector<lower = 0>[n_years] lambda; // max growth rate of LACO in absence of competition
    vector<lower = 0>[n_years] alpha_LACO; // competition term for LACO-LACO
    vector<lower = 0>[n_years] alpha_EG; // competition term for LACO-exotic grass
    vector<lower = 0>[n_years] alpha_ERVA;// competition term for LACO-ERVA
    vector<lower = 0>[n_years] alpha_NF; // competition term for LACO-non-native forb
}
model{
    for(j in 1:n_pools){
        obs_LACO[j,] ~ binomial(100, g_LACO); //
    }
    for(i in 1:n_years){
        obs_LACO[,i] ~ (g_LACO[i-1] * obs_LACO[,i-1] * lambda[i-1])/(1 + obs_LACO[,i-1] * alpha_LACO[i-1] + 
        obs_EG[,i-1] * alpha_EG[i-1] + obs_ERVA[,i-1] * alpha_ERVA[i-1] + obs_NF[,i-1] * alpha_NF[i-1]) + 
        s_LACO * (1 - g_LACO[i-1]) * obs_LACO[,i-1]);
    }
    lambda ~ normal(122.2561,83.95284); //get partially-informed priors from lit
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
}")

#run the model 
stan_fit <- rstan::optimizing(object = BH_model, 
                              data = list(n = ,
                                          n_pools = ,
                                          pool_id = ,
                                          n_years = ,
                                          year = ,
                                          obs_LACO = ,
                                          obs_EG = ,
                                          obs_ERVA = ,
                                          obs_NF = ,
                                          g_LACO = ,
                                          s_LACO = ), iter= 1000)
