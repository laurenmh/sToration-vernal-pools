# Lina's first attempt at modeling LACO population dynamics using the Beverton-Holt model
# Feel free to make changes on my code if you see fit 

# load packages
library(dplyr)
library(rstan)

# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# Stan model
BH_model <- stan_model(model_code="
data{
    int n; // number of observations
    int n.pools; // number of pools
    vector[n.pools] pool.id; // pool id
    int n.years; // number of years
    vector[n.years] year; // year
    int [n.pools, year] obs.LACO; // LACOdens
    int [n.pools, year] obs.EG; // exotic grass density
    int [n.pools, year] obs.ERVA; // ERVAdens
    int [n.pools, year] obs.NF; // non-native forb density
    vector[year] g.LACO; // growth rate of LACO
    real s.LACO; // survival rate of LACO
}
parameters{
    matrix<lower = 0>[n.pools, year] lambda; // max growth rate of LACO in absence of competition
    matrix<lower = 0>[n.pools, year] alpha.LACO; // competition term for LACO-LACO
    matrix<lower = 0>[n.pools, year] alpha.EG; // competition term for LACO-exotic grass
    matrix<lower = 0>[n.pools, year] alpha.ERVA;// competition term for LACO-ERVA
    matrix<lower = 0>[n.pools, year] alpha.NF; // competition term for LACO-non-native forb
}
model{
    for(j in 1:n.pools){
        obs.LACO[j,] ~ dpois(100, g.LACO[1]); // this g.LACO[1] is just first year's growth rate
    }
    for(i in 1:n.years){

    }
}")

