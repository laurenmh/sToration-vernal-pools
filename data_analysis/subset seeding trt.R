# Model LACO stem counts with storage effect and environmental fluctuations 
# Multiple year seeding addition
# Recursive modeling of seedbank to account for years where aboveground LACO is zero but seedbank persists to next year
# Germination rate of LACO depends on the previous year's exotic grass cover
# Germination rate is now a constant input, not an output.

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of LACO stem counts
sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of mean LACO stem counts

bh.formula <- function(sim_obs_LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g){
  sim_obs_LACO*lambda/(1+sim_obs_LACO*aii+EG*a1+ERVA*a2+NF*a3)+s*(1-g)*sim_obs_LACO/g
} #this is the modified Beverton-Holt model we'll use for LACO stem counts

#now simulate LACO stem counts
bh.sim <- function(n_pools, seedtrt, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
  for(i in 1:nrow(sim_mu)){
    for(j in 1:1){
      sim_mu[i,j] <- 100
      sim_obs_LACO[i,j] <- rbinom(1,100,g)
    }
    for(j in 2:3){
      if (EG[i,j-1]> 100){
        g = glow
      }
      else{g = g}
      sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-1],
                                EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                lambda = lambda[j-1], s = s, g = g)
      sim_obs_LACO[i,j] <- rpois(1, lambda = (sim_mu[i,j] + seedtrt[i,j] * g))
    }
    for(j in 4:ncol(sim_mu)){
      if (EG[i,j-1]> 100){
        g = glow
      }
      else{g = g}
      if (sim_obs_LACO[i,j-1] > 0){
        sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-1],
                                  EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                  aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                  lambda = lambda[j-1], s = s, g = g)
      }
      else {
        sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-2]*lambda[j-2]/(1+sim_obs_LACO[i,j-2]*aii[j-2]+EG[i,j-2]*a1[j-2]+ERVA[i,j-2]*a2[j-2]+NF[i,j-2]*a3[j-2])+s*(1-g)*sim_obs_LACO[i,j-2]/g,
                                  EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                  aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                  lambda = lambda[j-1], s = s, g = g)
      }
      sim_obs_LACO[i,j] <- rpois(1, lambda = sim_mu[i,j])
    }
  }
  return(sim_obs_LACO)
}

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
  for(a in 1:n_pools){
      for(b in 1:1){
          real germ_LACO;
          germ_LACO = high_germ_LACO;
          obs_LACO[a,b] ~ binomial(seeds_added[a,b], germ_LACO); //the first year's obs_LACO is the initial germination of seeds added in 1999
      }
      for(b in 2:3){
          real germ_LACO;
          if (obs_EG[a,b-1] > 100)
              germ_LACO = low_germ_LACO;
          else
              germ_LACO = high_germ_LACO;
          obs_LACO[a,b] ~ poisson(mu_LACO[a,b] + sigma + seeds_added[a,b] * germ_LACO); //the second and third year's obs_LACO is the sum of germination of seeds added in 2000 and 2001 and previous year's population. 
      }
      for(b in 4:(n_years-1)){
          obs_LACO[a,b] ~ poisson(mu_LACO[a,b] + sigma); //the rest of the year's obs_LACO is from a poisson distribution of mu_LACO. 
      }
  }
  lambda ~ normal(60,20); //get partially-informed priors from lit
  alpha_LACO ~ normal(0,1);
  alpha_EG ~ normal(0,1);
  alpha_ERVA ~ normal(0,1);
  alpha_NF ~ normal(0,1);
  sigma ~ normal(0,0.01);
  survival_LACO ~ normal(0,1);
}"

BH_model <- stan_model(model_code = BH_model_block)

### RUN THE MODEL ###

## Run with subset of real data 
## See data prep file before running this
BH_fit_GrAB <- sampling(BH_model,
                   data = list(n_pools = n_pools_AB,
                               n_years = n_years_AB,
                               obs_LACO = LACOdens_AB,
                               obs_EG = sumEGcover_AB,
                               obs_ERVA = ERVAdens_AB,
                               obs_NF = sumNFcover_AB,
                               seeds_added = seedtrt_AB[,4:6],
                               low_germ_LACO = 0.2,
                               high_germ_LACO = 0.8), 
                   iter= 1000)

BH_fit_GrBA <- sampling(BH_model,
                        data = list(n_pools = n_pools_BA,
                                    n_years = n_years_BA,
                                    obs_LACO = LACOdens_BA,
                                    obs_EG = sumEGcover_BA,
                                    obs_ERVA = ERVAdens_BA,
                                    obs_NF = sumNFcover_BA,
                                    seeds_added = seedtrt_BA[,4:6],
                                    low_germ_LACO = 0.2,
                                    high_germ_LACO = 0.8), 
                        iter= 1000)

BH_fit_LALA <- sampling(BH_model,
                        data = list(n_pools = n_pools_LALA,
                                    n_years = n_years_LALA,
                                    obs_LACO = LACOdens_LALA,
                                    obs_EG = sumEGcover_LALA,
                                    obs_ERVA = ERVAdens_LALA,
                                    obs_NF = sumNFcover_LALA,
                                    seeds_added = seedtrt_LALA[,4:6],
                                    low_germ_LACO = 0.2,
                                    high_germ_LACO = 0.8), 
                        iter= 1000)

BH_fit_LANo <- sampling(BH_model,
                        data = list(n_pools = n_pools_LANo,
                                    n_years = n_years_LANo,
                                    obs_LACO = LACOdens_LANo,
                                    obs_EG = sumEGcover_LANo,
                                    obs_ERVA = ERVAdens_LANo,
                                    obs_NF = sumNFcover_LANo,
                                    seeds_added = seedtrt_LANo[,4:6],
                                    low_germ_LACO = 0.2,
                                    high_germ_LACO = 0.8), 
                        iter= 1000)

### EXTRACT MODEL OUTPUT ###

# Check there is enough iteration
stan_trace(BH_fit_GrAB, pars = c("lambda"))
stan_trace(BH_fit_GrBA, pars = c("lambda"))
stan_trace(BH_fit_LALA, pars = c("lambda"))
stan_trace(BH_fit_LANo, pars = c("lambda"))

# Extract R hat (value greater than 1.1 means inadequate MCMC convergence)
summary(BH_fit_GrAB)$summary[,"Rhat"]
summary(BH_fit_GrBA)$summary[,"Rhat"]
summary(BH_fit_LALA)$summary[,"Rhat"]
summary(BH_fit_LANo)$summary[,"Rhat"]

# Mean posterior estimates of parameters
get_posterior_mean(BH_fit_GrAB, pars = c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF", "survival_LACO"))
get_posterior_mean(BH_fit_GrBA, pars = c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF", "survival_LACO"))
get_posterior_mean(BH_fit_LALA, pars = c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF", "survival_LACO"))
get_posterior_mean(BH_fit_LANo, pars = c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF", "survival_LACO"))

# Zoom into posterior distribution of parameters
plot(BH_fit_GrAB, pars = c("lambda"))
plot(BH_fit_GrBA, pars = c("lambda"))
plot(BH_fit_LALA, pars = c("alpha_LACO"))
plot(BH_fit_LANo, pars = c("alpha_LACO"))

# extract mean estimates 
alpha_LACO_mean_GrAB <- as.data.frame(get_posterior_mean(BH_fit_GrAB, pars = c("alpha_LACO")))
alpha_EG_mean_GrAB <- as.data.frame(get_posterior_mean(BH_fit_GrAB, pars = c("alpha_EG")))
alpha_ERVA_mean_GrAB <- as.data.frame(get_posterior_mean(BH_fit_GrAB, pars = c("alpha_ERVA")))
alpha_NF_mean_GrAB <- as.data.frame(get_posterior_mean(BH_fit_GrAB, pars = c("alpha_NF")))
lambda_mean_GrAB <- as.data.frame(get_posterior_mean(BH_fit_GrAB, pars = c("lambda")))
s_mean_GrAB <- as.data.frame(get_posterior_mean(BH_fit_GrAB, pars = c("survival_LACO")))

alpha_LACO_mean_GrBA <- as.data.frame(get_posterior_mean(BH_fit_GrBA, pars = c("alpha_LACO")))
alpha_EG_mean_GrBA <- as.data.frame(get_posterior_mean(BH_fit_GrBA, pars = c("alpha_EG")))
alpha_ERVA_mean_GrBA <- as.data.frame(get_posterior_mean(BH_fit_GrBA, pars = c("alpha_ERVA")))
alpha_NF_mean_GrBA <- as.data.frame(get_posterior_mean(BH_fit_GrBA, pars = c("alpha_NF")))
lambda_mean_GrBA <- as.data.frame(get_posterior_mean(BH_fit_GrBA, pars = c("lambda")))
s_mean_GrBA <- as.data.frame(get_posterior_mean(BH_fit_GrBA, pars = c("survival_LACO")))

alpha_LACO_mean_LALA <- as.data.frame(get_posterior_mean(BH_fit_LALA, pars = c("alpha_LACO")))
alpha_EG_mean_LALA <- as.data.frame(get_posterior_mean(BH_fit_LALA, pars = c("alpha_EG")))
alpha_ERVA_mean_LALA <- as.data.frame(get_posterior_mean(BH_fit_LALA, pars = c("alpha_ERVA")))
alpha_NF_mean_LALA <- as.data.frame(get_posterior_mean(BH_fit_LALA, pars = c("alpha_NF")))
lambda_mean_LALA <- as.data.frame(get_posterior_mean(BH_fit_LALA, pars = c("lambda")))
s_mean_LALA <- as.data.frame(get_posterior_mean(BH_fit_LALA, pars = c("survival_LACO")))

alpha_LACO_mean_LANo <- as.data.frame(get_posterior_mean(BH_fit_LANo, pars = c("alpha_LACO")))
alpha_EG_mean_LANo <- as.data.frame(get_posterior_mean(BH_fit_LANo, pars = c("alpha_EG")))
alpha_ERVA_mean_LANo <- as.data.frame(get_posterior_mean(BH_fit_LANo, pars = c("alpha_ERVA")))
alpha_NF_mean_LANo <- as.data.frame(get_posterior_mean(BH_fit_LANo, pars = c("alpha_NF")))
lambda_mean_LANo <- as.data.frame(get_posterior_mean(BH_fit_LANo, pars = c("lambda")))
s_mean_LANo <- as.data.frame(get_posterior_mean(BH_fit_LANo, pars = c("survival_LACO")))

#########################################
# Plot lambda constructed vs. reference #
#########################################
# combine lambda estimates from constructed and reference models
const_lambda_trim_GrAB <- as.data.frame(lambda_mean_GrAB[-c(1,2,16,17),5])#trim 2001, 2002, 2016, and 2017 estimates
const_lambda_trim_GrBA <- as.data.frame(lambda_mean_GrBA[-c(1,2,16,17),5])
const_lambda_trim_LALA <- as.data.frame(lambda_mean_LALA[-c(1,2,16,17),5])
const_lambda_trim_LANo <- as.data.frame(lambda_mean_LANo[-c(1,2,16,17),5])
ref_lambda_trim <- as.data.frame(reflambda_mean[,5]) 
join_lambda <- cbind(ref_lambda_trim, const_lambda_trim_GrAB, const_lambda_trim_GrBA, const_lambda_trim_LALA, const_lambda_trim_LANo)
row.names(join_lambda) <- c(2003:2015)
colnames(join_lambda) <- c("reference", "GrAB", "GrBA", "LALA", "LANo")

ggplot(join_lambda, aes(y = GrAB, x = reference)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO lambda", x = "Reference LACO lambda") +
  annotate("text", label = "R^2 = 0.1045", x = 20, y = 50) +
  theme_bw()

