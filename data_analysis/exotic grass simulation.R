# Simulate LACO populations with and without exotic grass removal

library(tidyverse)

# Simulation model:
sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of LACO stem counts
sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of mean LACO stem counts

bh.formula <- function(sim_obs_LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g){
  sim_obs_LACO*lambda/(1+sim_obs_LACO*aii+EG*a1+ERVA*a2+NF*a3)+s*(1-g)*sim_obs_LACO/g
} #this is the modified Beverton-Holt model we'll use for LACO stem counts

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

#Simulate without exotic grass removal
predicted_LACO <- bh.sim(n_pools = n_pools,
                         seedtrt = as.matrix(seedtrt[,4:6]),
                         EG = as.matrix(sumEGcover),
                         ERVA = as.matrix(ERVAdens),
                         NF = as.matrix(sumNFcover),
                         aii = alpha_LACO_mean[,5],
                         a1 = alpha_EG_mean[,5],
                         a2 = alpha_ERVA_mean[,5], 
                         a3 = alpha_NF_mean[,5],
                         lambda = lambda_mean[,5],
                         s = s_mean[,5],
                         g = 0.7,
                         glow = 0.2)

#Remove 50% of EG every year from 2001-2017
mult.5 <- function(x)(x*0.5)
reduced50EGcover <- sumEGcover %>%
  mutate_at(c("2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016", "2017"), mult.5)

#Remove 25% of EG every year from 2001-2017
mult.25 <- function(x)(x*0.75)
reduced25EGcover <- sumEGcover %>%
  mutate_at(c("2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016", "2017"), mult.25)

#Remove 75% of EG every year from 2001-2017
mult.75 <- function(x)(x*0.25)
reduced75EGcover <- sumEGcover %>%
  mutate_at(c("2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016", "2017"), mult.75)

#Simulate with exotic grass removal
reduced50EG_LACO <- bh.sim(n_pools = n_pools,
                         seedtrt = as.matrix(seedtrt[,4:6]),
                         EG = as.matrix(reduced50EGcover),
                         ERVA = as.matrix(ERVAdens),
                         NF = as.matrix(sumNFcover),
                         aii = alpha_LACO_mean[,5],
                         a1 = alpha_EG_mean[,5],
                         a2 = alpha_ERVA_mean[,5], 
                         a3 = alpha_NF_mean[,5],
                         lambda = lambda_mean[,5],
                         s = s_mean[,5],
                         g = 0.7,
                         glow = 0.2)

reduced75EG_LACO <- bh.sim(n_pools = n_pools,
                           seedtrt = as.matrix(seedtrt[,4:6]),
                           EG = as.matrix(reduced75EGcover),
                           ERVA = as.matrix(ERVAdens),
                           NF = as.matrix(sumNFcover),
                           aii = alpha_LACO_mean[,5],
                           a1 = alpha_EG_mean[,5],
                           a2 = alpha_ERVA_mean[,5], 
                           a3 = alpha_NF_mean[,5],
                           lambda = lambda_mean[,5],
                           s = s_mean[,5],
                           g = 0.7,
                           glow = 0.2)

reduced25EG_LACO <- bh.sim(n_pools = n_pools,
                           seedtrt = as.matrix(seedtrt[,4:6]),
                           EG = as.matrix(reduced25EGcover),
                           ERVA = as.matrix(ERVAdens),
                           NF = as.matrix(sumNFcover),
                           aii = alpha_LACO_mean[,5],
                           a1 = alpha_EG_mean[,5],
                           a2 = alpha_ERVA_mean[,5], 
                           a3 = alpha_NF_mean[,5],
                           lambda = lambda_mean[,5],
                           s = s_mean[,5],
                           g = 0.7,
                           glow = 0.2)
#Combine simulated LACO
colnames(predicted_LACO) <- c("2000", "2001", "2002", "2003", "2004", "2005", "2006",
                              "2007", "2008", "2009", "2010", "2011", "2012", "2013",
                              "2014", "2015", "2016", "2017")
colnames(reduced50EG_LACO) <- c("2000", "2001", "2002", "2003", "2004", "2005", "2006",
                              "2007", "2008", "2009", "2010", "2011", "2012", "2013",
                              "2014", "2015", "2016", "2017")
colnames(reduced75EG_LACO) <- c("2000", "2001", "2002", "2003", "2004", "2005", "2006",
                                "2007", "2008", "2009", "2010", "2011", "2012", "2013",
                                "2014", "2015", "2016", "2017")
colnames(reduced25EG_LACO) <- c("2000", "2001", "2002", "2003", "2004", "2005", "2006",
                                "2007", "2008", "2009", "2010", "2011", "2012", "2013",
                                "2014", "2015", "2016", "2017")
predicted_LACO <- as.data.frame(predicted_LACO) %>% 
  mutate(Pool = row_number()) %>%
  gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, `2007`, `2008`, `2009`, `2010`,
         `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = predicted_LACO)
reduced50EG_LACO <- as.data.frame(reduced50EG_LACO) %>% 
  mutate(Pool = row_number()) %>%
  gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, `2007`, `2008`, `2009`, `2010`,
         `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = reduced50EG_LACO)
reduced75EG_LACO <- as.data.frame(reduced75EG_LACO) %>% 
  mutate(Pool = row_number()) %>%
  gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, `2007`, `2008`, `2009`, `2010`,
         `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = reduced75EG_LACO)
reduced25EG_LACO <- as.data.frame(reduced25EG_LACO) %>% 
  mutate(Pool = row_number()) %>%
  gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, `2007`, `2008`, `2009`, `2010`,
         `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = reduced25EG_LACO)
grass_sim_LACO <- left_join(left_join(left_join(predicted_LACO, reduced50EG_LACO, by = c("Pool", "time")), reduced75EG_LACO, by = c("Pool", "time")), reduced25EG_LACO, by = c("Pool", "time")) %>%
  gather(`predicted_LACO`, `reduced50EG_LACO`, `reduced75EG_LACO`, `reduced25EG_LACO`, key = type, value = LACO) %>%
  mutate(log_LACO = log(LACO)) %>%
  mutate_if(is.numeric, ~replace(., is.infinite(.), 0))

#plot them!
ggplot(grass_sim_LACO, aes(x = time, y = log_LACO, col = type)) +
  geom_jitter()

se<-function(x){
  sd(x)/sqrt(length(x))
} # this is a function for calculating standard error

summary_grass_sim_LACO <- grass_sim_LACO %>%
  group_by(time, type) %>%
  summarise(mean_log_LACO = mean(log_LACO),
            se_log_LACO = se(log_LACO),
            mean_LACO = mean(LACO),
            se_LACO = se(LACO))

ggplot(summary_grass_sim_LACO%>%filter(type != "reduced25EG_LACO"), aes(x = time, y = mean_LACO, col = type)) +
  geom_point() +
  geom_line(aes(x = time, y = mean_LACO, group = type)) +
  geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 1) +
  theme_bw() +
  labs(x = "Year", y = "Mean LACO density") +
  scale_color_discrete(name = "Treatment", labels = c("No grass removal", "50% grass removed", "75% grass removed"))

#Graph mean of exotic grass cover timeseries
EG_summary <- const_dummy_join %>%
  group_by(Year) %>%
  summarise(EG_mean = mean(sum_EG),
            EG_se = se(sum_EG)) %>%
  mutate(type = "constructed")
EG_ref_summary <- ref_com_join %>%
  group_by(Year) %>%
  summarise(EG_mean = mean(sum_EG),
            EG_se = se(sum_EG)) %>%
  mutate(type = "reference")
EG_summary_join <- rbind(EG_summary, EG_ref_summary) 
ggplot(EG_summary_join, aes(x = Year, y = EG_mean, col= type)) +
  geom_point()+
  geom_errorbar(aes(ymin = EG_mean-EG_se, ymax = EG_mean+EG_se)) +
  theme_bw() +
  labs(x = "Year", y = "Mean exotic grass cover (%)") 
