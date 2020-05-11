# Model LACO stem counts with storage effect and environmental fluctuations 
# Multiple year seeding addition
# Recursive modeling of seedbank to account for years where aboveground LACO is zero but seedbank persists to next year.

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

### CREATE SIMULATED DATA ###

sim_n_pools <- 152 #number of pools
sim_n_years <- 7 #years of data

set.seed(124) #this helps create simulated values that are reproducible
sim_obs_EG <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 100), ncol=sim_n_years) #simulate poisson distributed observed EG density with mean of 100
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
sim_seedtrt[,1] <- sample(c(100,100,100,100), sim_n_pools, replace = TRUE)
sim_seedtrt[,2] <- sample(c(0,100,100,100), sim_n_pools, replace = TRUE)
sim_seedtrt[,3] <- sample(c(0,0,0,100), sim_n_pools, replace = TRUE)

sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years)
sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years)

#fix this simulation to include seeding treatments and seedbank
bh.sim <- function(n_pools, seedtrt, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g){
  for(i in 1:nrow(sim_mu)){
    for(j in 1:1){
      sim_mu[i,j] <- 100
      sim_obs_LACO[i,j] <- rbinom(1,100,g)
    }
    for(j in 2:3){
      sim_mu[i,j] <- sim_obs_LACO[i,j-1]*lambda[j-1]/(1+sim_obs_LACO[i,j-1]*aii[j-1]+EG[i,j-1]*a1[j-1]+ERVA[i,j-1]*a2[j-1]+NF[i,j-1]*a3[j-1])+s*(1-g)*sim_obs_LACO[i,j-1]/g
      sim_obs_LACO[i,j] <- rpois(1, lambda = (sim_mu[i,j] + seedtrt[i,j] * g))
    }
    for(j in 4:ncol(sim_mu)){
      if (sim_obs_LACO[i,j-1] > 0){
        sim_mu[i,j] <- sim_obs_LACO[i,j-1]*lambda[j-1]/(1+sim_obs_LACO[i,j-1]*aii[j-1]+EG[i,j-1]*a1[j-1]+ERVA[i,j-1]*a2[j-1]+NF[i,j-1]*a3[j-1])+s*(1-g)*sim_obs_LACO[i,j-1]/g
      }
      else {
        sim_mu[i,j] <- (sim_obs_LACO[i,j-2]*lambda[j-2]*lambda[j-1]/(1+aii[i,j-2]*sim_obs_LACO[i,j-2]+EG[i,j-2]*a1[j-2]+ERVA[i,j-2]*a2[j-2]+NF[i,j-2]*a3[j-2])+
                          (s*(1-g)*sim_obs_LACO[i,j-2]*lambda[j-1])/g)/
          (1+(aii[j-1]*sim_obs_LACO[i,j-2]*lambda[j-2])/(1+aii[i,j-2]*sim_obs_LACO[i,j-2]+EG[i,j-2]*a1[j-2]+ERVA[i,j-2]*a2[j-2]+NF[i,j-2]*a3[j-2])+
             s*(1-g)*sim_obs_LACO[i,j-2]*aii[j-2]/g + EG[i,j-1]*a1[j-1]+ERVA[i,j-1]*a2[j-1]+NF[i,j-1]*a3[j-1]) +
          ((s*(1-g)*sim_obs_LACO[i,j-2]*lambda[j-2])+(s^2*(1-g)^2*sim_obs_LACO[i,j-2]/g)/g)
      }
      sim_obs_LACO[i,j] <- rpois(1, lambda = sim_mu[i,j])
    }
  }
  return(sim_obs_LACO)
}

# List "true" lambda and alpha parameter values here. Start with constant parameters. 
# After running the model, check that the model outputs are close to these values.
sim_obs_LACO <- bh.sim(n_pools = sim_n_pools,
                       seedtrt = sim_seedtrt,
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
}
parameters{
    vector<lower = 0>[n_years-1] lambda; // max growth rate of LACO in absence of competition
    vector<lower = 0, upper = 1>[n_years-1] alpha_LACO; // competition term for LACO-LACO
    vector<lower = 0, upper = 1>[n_years-1] alpha_EG; // competition term for LACO-exotic grass
    vector<lower = 0, upper = 1>[n_years-1] alpha_ERVA;// competition term for LACO-ERVA
    vector<lower = 0, upper = 1>[n_years-1] alpha_NF; // competition term for LACO-non-native forb
    real <lower = 0, upper = 0.1> sigma; // error term for expected value of LACO
    real <lower = 0, upper = 1> survival_LACO; // survival rate of LACO seeds in the seedbank
    real <lower = 0, upper = 1> germ_LACO; // germination rate of LACO
}
transformed parameters{
    matrix [n_pools, n_years-1] mu_LACO;// expected value of LACO at time t
    for(i in 1:n_pools){  
        for(j in 1:1){
            mu_LACO[i,j] = 100;
        }
        for(j in 2:2){
            mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
            obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
            survival_LACO * (1-germ_LACO) * obs_LACO[i,j-1] ./ germ_LACO; // modified Beverton-Holt model
        }
        for(j in 3:(n_years-1)){
        if (obs_LACO[i,j-1] > 0)
            mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
            obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
            survival_LACO * (1-germ_LACO) * obs_LACO[i,j-1] ./ germ_LACO; // modified Beverton-Holt model
        else
            mu_LACO[i,j] = ((obs_LACO[i,j-2] * lambda[j-2] * lambda[j-1]) ./ 
            (1 + alpha_LACO[j-2] * obs_LACO[i,j-2] + obs_EG[i,j-2] * alpha_EG[j-2] + obs_ERVA[i,j-2] * alpha_ERVA[j-2] + obs_NF[i,j-2] * alpha_NF[j-2]) + 
            (survival_LACO * (1-germ_LACO) * obs_LACO[i,j-2] * lambda[j-1]) ./ germ_LACO) ./
            (1 + (alpha_LACO[j-1] * obs_LACO[i,j-2] * lambda[j-2]) ./ 
            (1 + alpha_LACO[j-2] * obs_LACO[i,j-2] + obs_EG[i,j-2] * alpha_EG[j-2] + obs_ERVA[i,j-2] * alpha_ERVA[j-2] + obs_NF[i,j-2] * alpha_NF[j-2]) +
            (survival_LACO * (1-germ_LACO) * obs_LACO[i, j-2] * alpha_LACO[j-2]) ./ germ_LACO + 
            obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
            ((survival_LACO * (1-germ_LACO) * obs_LACO[i,j-2] * lambda[j-2]) + 
            (survival_LACO ^ 2 * (1-germ_LACO) ^ 2 * obs_LACO[i,j-2] ./ germ_LACO))./germ_LACO; // use t-2 pop data to model t pop
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
                               seeds_added = sim_seedtrt), 
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
                               seeds_added = seedtrt[,4:6]), 
                   iter= 1000)

### EXTRACT MODEL OUTPUT ###

# Check there is enough iteration
stan_trace(BH_fit, pars = c("lambda"))

# extract mean estimates 
alpha_LACO_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_LACO")))
alpha_EG_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_EG")))
alpha_ERVA_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_ERVA")))
alpha_NF_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_NF")))
lambda_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("lambda")))
s_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("survival_LACO")))
g_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("germ_LACO")))

### COMPARE OBSERVED AND PREDICTED LACO ###
library(tidyr)
library(ggplot2)

#Option 1: use simulated data
#make a table of predicted LACO from estimated parameters
predicted_LACO_sim <- bh.sim(n_pools = sim_n_pools,
                         seedtrt = sim_seedtrt,
                         EG = sim_obs_EG,
                         ERVA = sim_obs_ERVA,
                         NF = sim_obs_NF,
                         aii = alpha_LACO_mean[,5],
                         a1 = alpha_EG_mean[,5],
                         a2 = alpha_ERVA_mean[,5], 
                         a3 = alpha_NF_mean[,5],
                         lambda = lambda_mean[,5],
                         s = s_mean[,5],
                         g = g_mean[,5])

#plot simulated LACOdens vs predicted_LACO_sim to check model fit
colnames(predicted_LACO_sim) <- c(1:7)
predicted_LACO_sim <- as.data.frame(predicted_LACO_sim) %>% 
  mutate(Pool = row_number()) %>%
  gather(`1`,`2`,`3`,`4`,`5`,`6`,`7`, key = time, value = predicted_LACO)
colnames(sim_obs_LACO) <- c(1:7)
sim_obs_LACO <- as.data.frame(sim_obs_LACO) %>% 
  mutate(Pool = row_number()) %>%
  gather(`1`,`2`,`3`,`4`,`5`,`6`,`7`, key = time, value = sim_LACO)
join_sim_LACO <- left_join(predicted_LACO_sim, sim_obs_LACO, by = c("Pool", "time"))

summary(lm(predicted_LACO ~ sim_LACO, data = join_sim_LACO)) #R2 = 0.875
ggplot(join_sim_LACO, aes(x = sim_LACO, y = predicted_LACO)) +
  geom_point() +
  annotate("text", label = "R^2 = 0.875", x = 50, y = 150)#looks like a good fit

#plot timeseries of simulated LACOdens and predicted_LACO_sim
long_join_sim <- join_sim_LACO %>% gather(`predicted_LACO`, `sim_LACO`, key = type, value = LACO)
ggplot(long_join_sim, aes(x = time, y = LACO, color = type)) +
  geom_jitter() +
  labs(y = "LACO count") +
  scale_color_discrete(breaks = c("predicted_LACO", "sim_LACO"),
                       labels = c("predicted", "simulated"))

#Option 2: use real data
predicted_LACO <- bh.sim(n_pools = n_pools,
                      seedtrt = as.matrix(seedtrt[,4:6]),
                      EG = as.matrix(sumEGdens),
                      ERVA = as.matrix(ERVAdens),
                      NF = as.matrix(sumNFdens),
                      aii = alpha_LACO_mean[,5],
                      a1 = alpha_EG_mean[,5],
                      a2 = alpha_ERVA_mean[,5], 
                      a3 = alpha_NF_mean[,5],
                      lambda = lambda_mean[,5],
                      s = s_mean[,5],
                      g = g_mean[,5])

#plot LACOdens vs predicted_LACO for modelfit
colnames(predicted_LACO) <- c("2000", "2001", "2002", "2003", "2004", "2005", "2006")
predicted_LACO <- as.data.frame(predicted_LACO) %>% 
  mutate(Pool = row_number()) %>%
  gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, key = time, value = predicted_LACO)
obs_LACO <- LACOdens %>% 
  mutate(Pool = row_number()) %>%
  gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, key = time, value = observed_LACO)
join_real_LACO <- left_join(predicted_LACO, obs_LACO, by = c("Pool", "time"))

summary(lm(predicted_LACO ~ observed_LACO, data = join_real_LACO)) #R2 = 0.111
ggplot(join_real_LACO, aes(x = observed_LACO, y = predicted_LACO)) +
  geom_point()+
  annotate("text", label = "R^2 = 0.111", x = 3000, y = 1500)

#plot timeseries of LACOdens and predicted_LACO 
long_join_real <- join_real_LACO %>% gather(`predicted_LACO`, `observed_LACO`, key = type, value = LACO)
ggplot(long_join_real, aes(x = time, y = LACO, color = type)) +
  geom_jitter() +
  labs(y = "LACO count") +
  scale_color_discrete(breaks = c("predicted_LACO", "observed_LACO"),
                       labels = c("predicted", "observed"))



