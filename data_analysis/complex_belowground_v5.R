# Model LACO stem counts with storage effect and environmental fluctuations 
# Multiple year seeding addition
# Recursive modeling of seedbank to account for years where aboveground LACO is zero but seedbank persists to next year
# Germination rate of LACO depends on the previous year's exotic grass cover
# Germination rate is now a constant input, not an output.

# load packages
library(dplyr)
library(rstan)
library(StanHeaders)

# ### CREATE SIMULATED DATA ###
# 
# sim_n_pools <- 142 #number of pools
# sim_n_years <- 18 #years of data
# 
# set.seed(124) #this helps create simulated values that are reproducible
# sim_obs_EG <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 35), ncol=sim_n_years) #simulate exotic grass(EG) cover 
# set.seed(123)
# sim_obs_ERVA <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 10), ncol=sim_n_years) #simulate ERVA cover 
# set.seed(122)
# sim_obs_NF <- matrix(rpois(sim_n_pools*sim_n_years, lambda = 25), ncol=sim_n_years) #simulate native forb (NF) cover 
# 
# sim_aii <- as.vector(rbeta(sim_n_years, 3, 1)) #simulate alpha_LACO from beta distribution
# sim_a1 <- as.vector(rbeta(sim_n_years, 1, 2)) #simulate alpha_EG from beta distribution
# sim_a2 <- as.vector(rbeta(sim_n_years, 2, 5)) #simulate alpha_ERVA from beta distribution
# sim_a3 <- as.vector(rbeta(sim_n_years, 1, 5)) #simulate alpha_NF from beta distribution
# sim_lambda <- as.vector(rpois(sim_n_years, 40)) #simulate lambda from poisson distribution
# 
# sim_seedtrt <- matrix(nrow = sim_n_pools, ncol = 3) #dummy seeding treatment
# sim_seedtrt[,1] <- sample(c(100,100,100,100), sim_n_pools, replace = TRUE)
# sim_seedtrt[,2] <- sample(c(0,100,100,100), sim_n_pools, replace = TRUE)
# sim_seedtrt[,3] <- sample(c(0,0,0,100), sim_n_pools, replace = TRUE)
# 
# sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of LACO seed counts
# sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of mean LACO seed counts
# sim_stem_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) # empty metrix of LACO stem counts
# 
# bh.formula <- function(sim_obs_LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g){
#   sim_obs_LACO*lambda/(1+sim_obs_LACO*aii+EG*a1+ERVA*a2+NF*a3)+s*(1-g)*sim_obs_LACO/g
# } #this is the modified Beverton-Holt model we'll use for LACO stem counts
# 
# #now simulate LACO stem counts
# bh.sim <- function(n_pools, seedtrt, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
#   for(i in 1:nrow(sim_mu)){
#     for(j in 1:1){
#       sim_mu[i,j] <- 100
#       sim_stem_LACO[i,j] <- sim_mu[i,j] * g
#       sim_obs_LACO[i,j] <- rbinom(1,100,g)
#     }
#     for(j in 2:3){
#       if (EG[i,j-1]> 100){
#         g = glow
#       }
#       else{g = g}
#       sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-1],
#                                 EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
#                                 aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
#                                 lambda = lambda[j-1], s = s, g = g)
#       sim_stem_LACO[i,j] <- sim_mu[i,j] *g
#       if(sim_stem_LACO[i,j] > 0){
#         sim_obs_LACO[i,j] <- rpois(1, lambda = (seedtrt[i,j] * g + sim_stem_LACO[i,j]))
#       }
#       else{
#         sim_obs_LACO[i,j] <- rpois(1, lambda = seedtrt[i,j] * g)
#       }
#     }
#     for(j in 4:ncol(sim_mu)){
#       if (EG[i,j-1]> 100){
#         g = glow
#       }
#       else{g = g}
#       if (sim_obs_LACO[i,j-1] > 0){
#         sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-1],
#                                    EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
#                                    aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
#                                    lambda = lambda[j-1], s = s, g = g)
#       }
#       else {
#         sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-2]*lambda[j-2]/(1+sim_obs_LACO[i,j-2]*aii[j-2]+EG[i,j-2]*a1[j-2]+ERVA[i,j-2]*a2[j-2]+NF[i,j-2]*a3[j-2])+s*(1-g)*sim_obs_LACO[i,j-2]/g,
#                                   EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
#                                   aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
#                                   lambda = lambda[j-1], s = s, g = g)
#       }
#       sim_stem_LACO[i,j] <- sim_mu[i,j] *g
#       if(sim_stem_LACO[i,j] > 0){
#         sim_obs_LACO[i,j] <- rpois(1, lambda = sim_mu[i,j]*g)
#       }
#       else{
#         sim_obs_LACO[i,j] = 0
#       }
#     }
#   }
#   return(sim_obs_LACO)
# }
# 
# # List "true" lambda and alpha parameter values here. Start with constant parameters. 
# # After running the model, check that the model outputs are close to these values.
# sim_obs_LACO <- bh.sim(n_pools = sim_n_pools,
#                        seedtrt = sim_seedtrt,
#                        EG = sim_obs_EG,
#                        ERVA = sim_obs_ERVA,
#                        NF = sim_obs_NF,
#                        aii = sim_aii,
#                        a1 = sim_a1,
#                        a2 = sim_a2, 
#                        a3 = sim_a3,
#                        lambda = sim_lambda,
#                        s = 0.2,
#                        g = 0.7,
#                        glow = 0.2)
# 
# hist(sim_obs_LACO) # Check distribution of simulated LACO

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
### RUN THE MODEL ###

## Option 1: run with simulated data
# BH_fit <- sampling(BH_model,
#                    data = list(n_pools = sim_n_pools,
#                                n_years = sim_n_years,
#                                obs_LACO = sim_obs_LACO,
#                                obs_EG = sim_obs_EG,
#                                obs_ERVA = sim_obs_ERVA,
#                                obs_NF = sim_obs_NF,
#                                seeds_added = sim_seedtrt,
#                                low_germ_LACO = 0.2,
#                                high_germ_LACO = 0.7), 
#                    iter= 1000)

## Option 2: run with real data 
## See data prep file before running this
BH_fit <- sampling(BH_model,
                   data = list(n_pools = n_pools,
                               n_years = n_years,
                               obs_LACO = LACOdens,
                               obs_EG = sumEGcover,
                               obs_ERVA = ERVAdens,
                               obs_NF = sumNFcover,
                               seeds_added = seedtrt[,4:6],
                               low_germ_LACO = 0.2,
                               high_germ_LACO = 0.7), 
                   iter= 1000)

### EXTRACT MODEL OUTPUT ###

# Check there is enough iteration
stan_trace(BH_fit, pars = c("lambda"))

# Extract R hat (value greater than 1.1 means inadequate MCMC convergence)
summary(BH_fit)$summary[,"Rhat"]

# Mean posterior estimates of parameters
get_posterior_mean(BH_fit, pars = c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF", "survival_LACO"))

# Zoom into posterior distribution of parameters
plot(BH_fit, pars = c("alpha_LACO"))

# extract mean estimates 
alpha_LACO_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_LACO")))
alpha_EG_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_EG")))
alpha_ERVA_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_ERVA")))
alpha_NF_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_NF")))
lambda_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("lambda")))
s_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("survival_LACO")))

# ### COMPARE OBSERVED AND PREDICTED LACO ###
# library(tidyr)
# library(ggplot2)

# #Option 1: use simulated data
# #make a table of predicted LACO from estimated parameters
# predicted_LACO_sim <- bh.sim(n_pools = sim_n_pools,
#                          seedtrt = sim_seedtrt,
#                          EG = sim_obs_EG,
#                          ERVA = sim_obs_ERVA,
#                          NF = sim_obs_NF,
#                          aii = alpha_LACO_mean[,5],
#                          a1 = alpha_EG_mean[,5],
#                          a2 = alpha_ERVA_mean[,5], 
#                          a3 = alpha_NF_mean[,5],
#                          lambda = lambda_mean[,5],
#                          s = s_mean[,5],
#                          g = 0.7,
#                          glow = 0.2)
# 
# #plot simulated LACOdens vs predicted_LACO_sim to check model fit
# colnames(predicted_LACO_sim) <- c(1:18)
# predicted_LACO_sim <- as.data.frame(predicted_LACO_sim) %>% 
#   mutate(Pool = row_number()) %>%
#   gather(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`,`13`,`14`,`15`,`16`,`17`,`18`, key = time, value = predicted_LACO)
# colnames(sim_obs_LACO) <- c(1:18)
# sim_obs_LACO <- as.data.frame(sim_obs_LACO) %>% 
#   mutate(Pool = row_number()) %>%
#   gather(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`,`13`,`14`,`15`,`16`,`17`,`18`, key = time, value = sim_LACO)
# join_sim_LACO <- left_join(predicted_LACO_sim, sim_obs_LACO, by = c("Pool", "time"))
# 
# summary(lm(predicted_LACO ~ sim_LACO, data = join_sim_LACO)) #R2 = 0.8644
# ggplot(join_sim_LACO, aes(x = sim_LACO, y = predicted_LACO)) +
#   geom_point() +
#   annotate("text", label = "R^2 = 0.8644", x = 50, y = 150) + #looks like a good fit 
#   geom_smooth(method = "lm") +
#   labs(x = "simulated LACO counts", y = "predicted LACO counts")
# 
# #plot timeseries of simulated LACOdens and predicted_LACO_sim
# long_join_sim <- join_sim_LACO %>% gather(`predicted_LACO`, `sim_LACO`, key = type, value = LACO)
# ggplot(long_join_sim, aes(x = time, y = LACO, color = type)) +
#   geom_jitter() +
#   labs(y = "LACO count") +
#   scale_color_discrete(breaks = c("predicted_LACO", "sim_LACO"),
#                        labels = c("predicted", "simulated"))

# #Option 2: use real data
# predicted_LACO <- bh.sim(n_pools = n_pools,
#                       seedtrt = as.matrix(seedtrt[,4:6]),
#                       EG = as.matrix(sumEGcover),
#                       ERVA = as.matrix(ERVAdens),
#                       NF = as.matrix(sumNFcover),
#                       aii = alpha_LACO_mean[,5],
#                       a1 = alpha_EG_mean[,5],
#                       a2 = alpha_ERVA_mean[,5], 
#                       a3 = alpha_NF_mean[,5],
#                       lambda = lambda_mean[,5],
#                       s = s_mean[,5],
#                       g = 0.7,
#                       glow = 0.2)

# #plot LACOdens vs predicted_LACO for modelfit
# colnames(predicted_LACO) <- c("2000", "2001", "2002", "2003", "2004", "2005", "2006",
#                               "2007", "2008", "2009", "2010", "2011", "2012", "2013",
#                               "2014", "2015", "2016", "2017")
# predicted_LACO <- as.data.frame(predicted_LACO) %>% 
#   mutate(Pool = row_number()) %>%
#   gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, `2007`, `2008`, `2009`, `2010`,
#          `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = predicted_LACO)
# obs_LACO <- LACOdens %>% 
#   mutate(Pool = row_number()) %>%
#   gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`,`2007`,`2008`,`2009`,`2010`,
#          `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = observed_LACO)
# 
# join_real_LACO <- left_join(predicted_LACO, obs_LACO, by = c("Pool", "time")) %>%
#   mutate(log_predicted_LACO = log(predicted_LACO)) %>%
#   mutate(log_observed_LACO = log(observed_LACO)) %>%
#   mutate_if(is.numeric, ~replace(., is.infinite(.), 0))
# 
# summary(lm(predicted_LACO ~ observed_LACO, data = join_real_LACO)) #R2 = 0.1806
# ggplot(join_real_LACO, aes(x = observed_LACO, y = predicted_LACO)) +
#   geom_point()+
#   annotate("text", label = "R^2 = 0.1208", x = 3000, y = 2500) +
#   ylim(0, 4000) #most values are small, so try log scale
# 
# summary(lm(log_predicted_LACO ~ log_observed_LACO, data = join_real_LACO))
# ggplot(join_real_LACO, aes(x = log_predicted_LACO, y = log_observed_LACO)) +
#   geom_point()+
#   annotate("text", label = "R^2 = 0.3009", x = 1.5, y = 8) +
#   geom_abline(intercept = 0, slope =1) #there are a bunch of zeros that are predicted non-zero values.
# 
# #how many of the observed zeros have non-zero predicted values?
# nonzero <- join_real_LACO %>%
#   filter(observed_LACO == 0) %>%
#   filter(predicted_LACO > 0) #851 out of 2556 observations
# 
# #plot timeseries of LACOdens and predicted_LACO 
# long_join_real <- join_real_LACO %>% gather(`predicted_LACO`, `observed_LACO`, key = type, value = LACO)
# ggplot(long_join_real, aes(x = time, y = log(LACO), color = type)) +
#   geom_jitter() +
#   labs(y = "log LACO count") +
#   scale_color_discrete(breaks = c("predicted_LACO", "observed_LACO"),
#                        labels = c("predicted", "observed"))

#COMPARE SEEDING TREATMENTS
BH_fit_AB <- sampling(BH_model,
                   data = list(n_pools = n_pools_AB,
                               n_years = n_years_AB,
                               obs_LACO = LACOdens_AB,
                               obs_EG = sumEGcover_AB,
                               obs_ERVA = ERVAdens_AB,
                               obs_NF = sumNFcover_AB,
                               seeds_added = seedtrt_AB[,4:6],
                               low_germ_LACO = 0.2,
                               high_germ_LACO = 0.7), 
                   iter= 1000)

BH_fit_BA <- sampling(BH_model,
                      data = list(n_pools = n_pools_BA,
                                  n_years = n_years_BA,
                                  obs_LACO = LACOdens_BA,
                                  obs_EG = sumEGcover_BA,
                                  obs_ERVA = ERVAdens_BA,
                                  obs_NF = sumNFcover_BA,
                                  seeds_added = seedtrt_BA[,4:6],
                                  low_germ_LACO = 0.2,
                                  high_germ_LACO = 0.7), 
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
                                  high_germ_LACO = 0.7), 
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
                                    high_germ_LACO = 0.7), 
                        iter= 1000)
