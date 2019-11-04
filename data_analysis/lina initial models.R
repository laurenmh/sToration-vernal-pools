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
    vector[n] obs.LACO; // LACOdens
    vector[n] g.LACO; // growth rate of LACO
    vector[n] s.LACO; // survival rate of LACO
    vector[n] obs.EG; // exotic grass density
    vector[n] obs.ERVA; // ERVAdens
    vector[n] obs.NF; // non-native forb density
}

parameters{
    vector<lower = 0>[pool.id, year] lambda; // max finite rate of LACO increase in absence of competition
    alpha.LACO;
    alpha.EG;
    alpha.ERVA;
    alpha.NF;
}
model{
  
}")

