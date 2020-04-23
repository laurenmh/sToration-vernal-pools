# Model LACO stem counts with storage effect and environmental fluctuations 
# Multiple year seeding addition included

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

### CREATE SIMULATED DATA ###

sim_n_pools <- 256 #number of pools
sim_n_years <- 10 #years of data

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

sim_seedtrt <- matrix(nrow = sim_n_pools, ncol = 3)
sim_seedtrt[,1] <- sample(c(100,100,100,100,100), sim_n_pools, replace = TRUE)
sim_seedtrt[,2] <- sample(c(0,0,100,100,100), sim_n_pools, replace = TRUE)
sim_seedtrt[,3] <- sample(c(0,0,0,0,100), sim_n_pools, replace = TRUE)

sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years)
sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years)

bh.sim <- function(n_pools, init, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g){
  sim_obs_LACO[,1]<- rbinom(n_pools,100,0.8)
  sim_mu[,1]<- init
  for(i in 1:nrow(sim_obs_LACO)){
    for(j in 2:ncol(sim_obs_LACO)){
      sim_mu[i,j] <- lambda[j-1]*sim_obs_LACO[i,j-1]/(1+sim_obs_LACO[i,j-1]*aii[j-1]+EG[i,j-1]*a1[j-1]+ERVA[i,j-1]*a2[j-1]+NF[i,j-1]*a3[j-1])+s*(1-g)*sim_obs_LACO[i,j-1]/g
      sim_obs_LACO[i,j] <- rpois(1, lambda=sim_mu[i,j])
    }
  }
  return(sim_obs_LACO)
}

# List "true" lambda and alpha parameter values here. Start with constant parameters. 
# After running the model, check that the model outputs are close to these values.
sim_obs_LACO <- bh.sim(n_pools = sim_n_pools,
                       init = 100,
                       EG = sim_obs_EG,
                       ERVA = sim_obs_ERVA,
                       NF = sim_obs_NF,
                       aii = sim_aii,
                       a1 = sim_a1,
                       a2 = sim_a2, 
                       a3 = sim_a3,
                       lambda = sim_lambda,
                       s = 0.5,
                       g = 0.7)

hist(sim_obs_LACO) # Check distribution of simulated LACO

### CREATE A STAN MODEL ###

BH_model_block <- "
data{
    int n_pools; // number of pools
    int n_years; // number of years
    int obs_LACO [n_pools, n_years] ; // observed LACO density // this is an array of integers because the model is a discrete poisson model (see below)
    matrix[n_pools, n_years] obs_EG; // exotic grass density
    matrix[n_pools, n_years] obs_ERVA; // ERVA density
    matrix[n_pools, n_years] obs_NF; // non-native forb density
    int seeds_added [n_pools, 3]; // number of seeds added in 1999-2001
    real germ_LACO; // germination rate of LACO
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
    matrix [n_pools, n_years-1] mu_LACO;// expected value of LACO at time t
    for(i in 1:n_pools){  
        for(j in 1:1){
            mu_LACO[i,j] = 100;
        }
        for(j in 2:(n_years-1)){
            mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
            obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
            survival_LACO * (1-germ_LACO) * obs_LACO[i,j-1] ./ germ_LACO; // modified Beverton-Holt model
        }
    }
}
model{
    for(a in 1:n_pools){
        for(b in 1:1){
            obs_LACO[a,b] ~ binomial(seeds_added[a,b], germ_LACO); //the first year's obs_LACO is the initial germination of seeds added in 1999
        }
        for(b in 2:3){
            obs_LACO[a,b] ~ poisson(mu_LACO[a,b] + sigma + seeds_added[a,b] * germ_LACO); //the second and third year's obs_LACO is the sum of germination of seeds added in 2000 and 2001 and previous year's population. 
        }
        for(b in 4:(n_years-1)){
            obs_LACO[a,b] ~ poisson(mu_LACO[a,b] + sigma); //the rest of the year's obs_LACO is from a poisson distribution of mu_LACO. 
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

BH_model <- stan_model(model_code = BH_model_block)

### RUN THE MODEL ###

## Option 1: run with simulated data
BH_fit <- sampling(BH_model,
                   data = list(n_pools = sim_n_pools,
                               n_years = sim_n_years,
                               obs_LACO = sim_obs_LACO,
                               obs_EG = sim_obs_EG,
                               obs_ERVA = sim_obs_ERVA,
                               obs_NF = sim_obs_NF,
                               seeds_added = sim_seedtrt,
                               germ_LACO = 0.7), 
                   iter= 1000)

## Option 2: run with real data 
## See data prep file before running this
BH_fit <- sampling(BH_model,
                   data = list(n_pools = n_pools,
                               n_years = n_years,
                               obs_LACO = LACOdens,
                               obs_EG = sumEGdens,
                               obs_ERVA = ERVAdens,
                               obs_NF = sumNFdens,
                               seeds_added = seedtrt[,4:6],
                               germ_LACO = 0.7), 
                   iter= 1000)

### EXTRACT MODEL OUTPUT ###

# Check there is enough iteration
stan_trace(BH_fit, pars = c("lambda"))

# Extract R hat (value greater than 1.1 means inadequate MCMC convergence)
summary(BH_fit)$summary[,"Rhat"]

# Mean posterior estimates of parameters
get_posterior_mean(BH_fit, pars = c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF", "survival_LACO"))

# Zoom into posterior distribution of parameters
plot(BH_fit, pars = c("lambda"))

# Are the parameters correlated?
list_of_draws <- extract(BH_fit) #extract the list of draws
#if tidyr is interfering with extract(), run this line of code: 
#.rs.unloadPackage("tidyr")
print(names(list_of_draws)) #see the names of parameters
head(list_of_draws$lambda) #see the first 6 draws of lambda
lambda_draws <- list_of_draws$lambda
alpha_LACO_draws <- list_of_draws$alpha_LACO

plot(lambda_draws[,1], alpha_LACO_draws[,1]) #lambda vs alpha_LACO at time 1
plot(lambda_draws[,2], alpha_LACO_draws[,2]) #lambda vs alpha_LACO at time 2
plot(lambda_draws[,3], alpha_LACO_draws[,3]) #lambda vs alpha_LACO at time 3
plot(lambda_draws[,4], alpha_LACO_draws[,4]) #lambda vs alpha_LACO at time 3

mean_lambda <- summary(BH_fit, pars = c("lambda"))$summary[,"mean"]
mean_alpha_LACO <- summary(BH_fit, pars = c("alpha_LACO"))$summary[,"mean"]
plot(mean_lambda, mean_alpha_LACO)


