# Lina's first attempt at modeling LACO population dynamics using the Beverton-Holt model
# Feel free to make changes on my code as you see fit 

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

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


# simulated data
sim_n_pools <- 3 #number of pools
sim_n_years <- 18 #years of data
set.seed(125) #this helps create simulated values that are reproducible
sim_obs_LACO <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 4), ncol=sim_n_years) #simulate poisson distributed LACO population
set.seed(124)
sim_obs_EG <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 4), ncol=sim_n_years)
set.seed(123)
sim_obs_ERVA <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 4), ncol=sim_n_years)
set.seed(122)
sim_obs_NF <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 4), ncol=sim_n_years)
set.seed(121)
sim_g_LACO <- as.vector(rbeta(sim_n_years, 2, 4)) #simulate germination rate of seedbank from a beta distribution

# Stan model
BH_model <- "
data{
    int n_pools; // number of pools
    int n_years; // number of years
    matrix[n_pools, n_years] obs_LACO; // observed LACO density
    matrix[n_pools, n_years] obs_EG; // exotic grass density
    matrix[n_pools, n_years] obs_ERVA; // ERVA density
    matrix[n_pools, n_years] obs_NF; // non-native forb density
    vector[n_years] g_LACO; // per seed germination rate of LACO
}
parameters{
    vector<lower = 0>[n_years] lambda; // max growth rate of LACO in absence of competition
    vector<lower = 0>[n_years] alpha_LACO; // competition term for LACO-LACO
    vector<lower = 0>[n_years] alpha_EG; // competition term for LACO-exotic grass
    vector<lower = 0>[n_years] alpha_ERVA;// competition term for LACO-ERVA
    vector<lower = 0>[n_years] alpha_NF; // competition term for LACO-non-native forb
    real<lower = 0, upper = 1> s_LACO; // seed survival of LACO seedbank
}
transformed parameters{
matrix[n_pools, n_years] N_LACO;// total population (observed + unobserved) of LACO at time t
for(i in 2:n_years){
        N_LACO[,i] = (g_LACO[i-1] * obs_LACO[,i-1] * lambda[i-1])./(1 + obs_LACO[,i-1] * alpha_LACO[i-1] * g_LACO[i-1]+ 
        obs_EG[,i-1] * alpha_EG[i-1] + obs_ERVA[,i-1] * alpha_ERVA[i-1] + obs_NF[,i-1] * alpha_NF[i-1]) + 
        s_LACO * (1 - g_LACO[i-1]) * obs_LACO[,i-1]  / g_LACO[i-1];
    }
}
model{
    lambda ~ normal(122.2561,83.95284); //get partially-informed priors from lit
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
    s_LACO ~ normal(0,1);
}"

#run the model with simulated data
stan_fit <- stan(model_code = BH_model,
                  data = list(n_pools = sim_n_pools,
                              n_years = sim_n_years,
                              obs_LACO = sim_obs_LACO,
                              obs_EG = sim_obs_EG,
                              obs_ERVA = sim_obs_ERVA,
                              obs_NF = sim_obs_NF,
                              g_LACO = sim_g_LACO), 
                  iter= 1000)

#check if there is enough iteration
stan_trace(stan_fit, pars = c("lambda"))

#mean posterior estimates of parameters
get_posterior_mean(stan_fit, pars = c("lambda"))

#zoom into posterior distribution of parameters
plot(stan_fit, pars = c("lambda", "alpha_LACO", "s_LACO"))


