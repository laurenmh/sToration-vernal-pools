#####################################################
#Would adaptive management improve LACO populations?#
#####################################################
#Goal: Simulate exotic grasses (EG) removal to promote LACO persistence

#Step 1. Simulate EG removal
#Step 2. Average the growth rates of LACO over time for all simulation scenarios
#Step 3. Plot modeled abundance and GRWR

# Load packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstan)
library(StanHeaders)

# Remember to set your data pathway first!

# Data
source("data_compiling/compile_composition.R") 
# Run "data_wrangling/prep data before modeling.R"

#Extract parameters for constructed pools
# Run "analysis/complex_belowground_v5.R"
Post <- rstan::extract(BH_fit)
alpha_LACO <- as.matrix(Post$alpha_LACO)
alpha_EG <- as.matrix(Post$alpha_EG)
alpha_ERVA <- as.matrix(Post$alpha_ERVA)
alpha_NF <- as.matrix(Post$alpha_NF)
lambda <- as.matrix(Post$lambda)
s <- as.matrix(Post$survival_LACO)

# Function for standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 

#-----------------------
#Step 1. Simulate EG removal

# Simulation model:
sim_n_pools <- 142 #number of pools
sim_n_years <- 18 #years of data
sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of LACO seed counts
sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of mean LACO seed counts
sim_stem_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) # empty metrix of LACO stem counts

bh.formula <- function(sim_obs_LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g){
  sim_obs_LACO*lambda/(1+sim_obs_LACO*aii+EG*a1+ERVA*a2+NF*a3)+s*(1-g)*sim_obs_LACO/g
} #this is the modified Beverton-Holt model we'll use for LACO stem counts

bh.sim <- function(n_pools, seedtrt, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
  for(i in 1:nrow(sim_mu)){
    for(j in 1:1){
      sim_mu[i,j] <- 100
      sim_stem_LACO[i,j] <- sim_mu[i,j] * g
      sim_obs_LACO[i,j] <- rbinom(1,100,g)
    }
    for(j in 2:3){
      if (EG[i,j-1]> 100){
        g = glow
      }
      else{g = g}
      sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-1],
                                EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                aii = aii[i,j-1], a1 = a1[i,j-1], a2 = a2[i,j-1], a3 = a3[i,j-1],
                                lambda = lambda[i,j-1], s = s[i], g = g)
      sim_stem_LACO[i,j] <- sim_mu[i,j] *g
      if(sim_stem_LACO[i,j] > 0){
        sim_obs_LACO[i,j] <- rpois(1, lambda = (seedtrt[i,j] * g + sim_stem_LACO[i,j]))
      }
      else{
        sim_obs_LACO[i,j] <- rpois(1, lambda = seedtrt[i,j] * g)
      }
    }
    for(j in 4:ncol(sim_mu)){
      if (EG[i,j-1]> 100){
        g = glow
      }
      else{g = g}
      if (sim_obs_LACO[i,j-1] > 0){
        sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-1],
                                  EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                  aii = aii[i,j-1], a1 = a1[i,j-1], a2 = a2[i,j-1], a3 = a3[i,j-1],
                                  lambda = lambda[i,j-1], s = s[i], g = g)
      }
      else {
        sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-2]*lambda[i,j-2]/(1+sim_obs_LACO[i,j-2]*aii[i,j-2]+EG[i,j-2]*a1[i,j-2]+ERVA[i,j-2]*a2[i,j-2]+NF[i,j-2]*a3[i,j-2])+s[i]*(1-g)*sim_obs_LACO[i,j-2]/g,
                                  EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                  aii = aii[i,j-1], a1 = a1[i,j-1], a2 = a2[i,j-1], a3 = a3[i,j-1],
                                  lambda = lambda[i,j-1], s = s[i], g = g)
      }
      sim_stem_LACO[i,j] <- sim_mu[i,j] *g
      if(sim_stem_LACO[i,j] > 0){
        sim_obs_LACO[i,j] <- rpois(1, lambda = sim_mu[i,j]*g)
      }
      else{
        sim_obs_LACO[i,j] = 0
      }
    }
  }
  return(sim_obs_LACO)
}

#Simulate LACO abundance without exotic grass removal
predicted_LACO <- bh.sim(n_pools = n_pools,
                         seedtrt = as.matrix(seedtrt[,4:6]),
                         EG = as.matrix(sumEGcover),
                         ERVA = as.matrix(ERVAdens),
                         NF = as.matrix(sumNFcover),
                         aii = alpha_LACO,
                         a1 = alpha_EG,
                         a2 = alpha_ERVA, 
                         a3 = alpha_NF,
                         lambda = lambda,
                         s = s,
                         g = 0.7,
                         glow = 0.2)

#Remove 25% of EG every year from 2001-2017
mult.25 <- function(x)(x*0.75)
reduced25EGcover <- sumEGcover %>%
  mutate_at(c("2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016", "2017"), mult.25)

#Remove 50% of EG every year from 2001-2017
mult.5 <- function(x)(x*0.5)
reduced50EGcover <- sumEGcover %>%
  mutate_at(c("2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016", "2017"), mult.5)

#Remove 75% of EG every year from 2001-2017
mult.75 <- function(x)(x*0.25)
reduced75EGcover <- sumEGcover %>%
  mutate_at(c("2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016", "2017"), mult.75)

#Simulate LACO abundance with exotic grass removal
reduced50EG_LACO <- bh.sim(n_pools = n_pools,
                         seedtrt = as.matrix(seedtrt[,4:6]),
                         EG = as.matrix(reduced50EGcover),
                         ERVA = as.matrix(ERVAdens),
                         NF = as.matrix(sumNFcover),
                         aii = alpha_LACO,
                         a1 = alpha_EG,
                         a2 = alpha_ERVA, 
                         a3 = alpha_NF,
                         lambda = lambda,
                         s = s,
                         g = 0.7,
                         glow = 0.2)

reduced75EG_LACO <- bh.sim(n_pools = n_pools,
                           seedtrt = as.matrix(seedtrt[,4:6]),
                           EG = as.matrix(reduced75EGcover),
                           ERVA = as.matrix(ERVAdens),
                           NF = as.matrix(sumNFcover),
                           aii = alpha_LACO,
                           a1 = alpha_EG,
                           a2 = alpha_ERVA, 
                           a3 = alpha_NF,
                           lambda = lambda,
                           s = s,
                           g = 0.7,
                           glow = 0.2)

reduced25EG_LACO <- bh.sim(n_pools = n_pools,
                           seedtrt = as.matrix(seedtrt[,4:6]),
                           EG = as.matrix(reduced25EGcover),
                           ERVA = as.matrix(ERVAdens),
                           NF = as.matrix(sumNFcover),
                           aii = alpha_LACO,
                           a1 = alpha_EG,
                           a2 = alpha_ERVA, 
                           a3 = alpha_NF,
                           lambda = lambda,
                           s = s,
                           g = 0.7,
                           glow = 0.2)
#Combine simulated LACO
years <- c("2000", "2001", "2002", "2003", "2004", "2005", "2006",
           "2007", "2008", "2009", "2010", "2011", "2012", "2013",
           "2014", "2015", "2016", "2017")
colnames(predicted_LACO) <- years
colnames(reduced50EG_LACO) <- years
colnames(reduced75EG_LACO) <- years
colnames(reduced25EG_LACO) <- years

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

summary_grass_sim_LACO <- grass_sim_LACO %>%
  group_by(time, type) %>%
  summarise(mean_log_LACO = mean(log_LACO),
            se_log_LACO = se(log_LACO),
            mean_LACO = mean(LACO),
            se_LACO = se(LACO),
            sd_LACO = sd(LACO))

#Step 2. Average the growth rates of LACO over time for all simulation scenarios

#Set up a simpler population model for calculating GRWR
bh.sim.control <- function(LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
  for(i in 1:nrow(sim_LACO)){
    for(j in 1:1){
      sim_LACO[i,j] <- LACO*lambda[i,j]/(1+LACO*aii[i,j]+EG[j]*a1[i,j]+ERVA[j]*a2[i,j]+NF[j]*a3[i,j])+s[i]*(1-g)*LACO/g #this is the modified Beverton-Holt model we'll use for LACO stem counts
    }
    for(j in 2:ncol(sim_LACO)){
      if (EG[j-1]> 100){
        g = glow
      }
      else{g = g}
      sim_LACO[i,j] <- LACO*lambda[i,j]/(1+LACO*aii[i,j]+EG[j]*a1[i,j]+ERVA[j]*a2[i,j]+NF[j]*a3[i,j])+s[i]*(1-g)*LACO/g 
    }
  }
  return(sim_LACO)
}

#Use the control plots in reference pools (no LACO present) to calculate the stable equilibrium frequency of non-LACO species in the model each year. 
    #Non-LACO species in our model:
    #ERVA
    #Exotic grass group - BRHO, HOMA, LOMU
    #Native forb group - PLST, DOCO
const_com_control <- const_com %>% #use constructed pools data
                            filter(Treatment.1999 == "Control") %>% #filter control plots only
                            drop_na() %>% #remove any rows with na
                            filter(LACO <= 0) %>% #remove communities with LACO present
                            mutate(sumEG = BRHO + HOMA + LOMU, 
                                   sumNF = PLST + DOCO) %>% #create a new column for sum of EG and sum of NF
                            group_by(Year)%>% #summarize by year
                            summarize(avg_ERVA = round(mean(ERVA), digits = 0),
                                      avg_sumEG = round(mean(sumEG), digits = 0),
                                      avg_sumNF = round(mean(sumNF), digits = 0)) %>%#take the average freq.
                            filter(Year != "2017") %>% #filter out 2017
                            pivot_longer(-Year)%>%
                            pivot_wider(names_from = Year, values_from = value)

# 0% EG removal
LACO_const <- bh.sim.control(   LACO = 1,
                                EG = as.numeric(const_com_control[2,2:18]),
                                ERVA = as.numeric(const_com_control[1,2:18]),
                                NF = as.numeric(const_com_control[3,2:18]),
                                aii = alpha_LACO,
                                a1 = alpha_EG,
                                a2 = alpha_ERVA,
                                a3 = alpha_NF,
                                lambda = lambda,
                                s = s,
                                g = 0.7,
                                glow = 0.2)
GRWR_LACO_const <- log(LACO_const) #2000-2016

#Remove 25% of EG in all years
sim_LACO <- matrix(nrow = 2000, ncol = 17)
mult.25 <- function(x)(x*0.75)
reduced25EGcover_all <- const_com_control[2,2:18] %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.25)
LACO_25EG_all <- bh.sim.control(  LACO = 1,
                                  EG = as.numeric(reduced25EGcover_all),
                                  ERVA = as.numeric(const_com_control[1,2:18]),
                                  NF = as.numeric(const_com_control[3,2:18]),
                                  aii = alpha_LACO,
                                  a1 = alpha_EG,
                                  a2 = alpha_ERVA,
                                  a3 = alpha_NF,
                                  lambda = lambda,
                                  s = s,
                                  g = 0.7,
                                  glow = 0.2)
GRWR_LACO_25EG_all <- log(LACO_25EG_all)

#Remove 50% of EG in all years
mult.5 <- function(x)(x*0.5)
reduced50EGcover_all <- const_com_control[2,2:18] %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.5)
LACO_50EG_all <- bh.sim.control(  LACO = 1,
                                  EG = as.numeric(reduced50EGcover_all),
                                  ERVA = as.numeric(const_com_control[1,2:18]),
                                  NF = as.numeric(const_com_control[3,2:18]),
                                  aii = alpha_LACO,
                                  a1 = alpha_EG,
                                  a2 = alpha_ERVA,
                                  a3 = alpha_NF,
                                  lambda = lambda,
                                  s = s,
                                  g = 0.7,
                                  glow = 0.2)
GRWR_LACO_50EG_all <- log(LACO_50EG_all)

#Remove 75% of EG in all years
mult.75 <- function(x)(x*0.25)
reduced75EGcover_all <- const_com_control[2,2:18] %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.75)
LACO_75EG_all <- bh.sim.control(  LACO = 1,
                                  EG = as.numeric(reduced75EGcover_all),
                                  ERVA = as.numeric(const_com_control[1,2:18]),
                                  NF = as.numeric(const_com_control[3,2:18]),
                                  aii = alpha_LACO,
                                  a1 = alpha_EG,
                                  a2 = alpha_ERVA,
                                  a3 = alpha_NF,
                                  lambda = lambda,
                                  s = s,
                                  g = 0.7,
                                  glow = 0.2)
GRWR_LACO_75EG_all <- log(LACO_75EG_all)

#Remove 100% of EG in all years
mult.100 <- function(x)(x*0)
reduced100EGcover_all <- const_com_control[2,2:18] %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.100)
LACO_100EG_all <- bh.sim.control( LACO = 1,
                                  EG = as.numeric(reduced100EGcover_all),
                                  ERVA = as.numeric(const_com_control[1,2:18]),
                                  NF = as.numeric(const_com_control[3,2:18]),
                                  aii = alpha_LACO,
                                  a1 = alpha_EG,
                                  a2 = alpha_ERVA,
                                  a3 = alpha_NF,
                                  lambda = lambda,
                                  s = s,
                                  g = 0.7,
                                  glow = 0.2)
GRWR_LACO_100EG_all <- log(LACO_100EG_all)

#Average the growth rates of LACO over time for all simulation scenarios.
GRWR_LACO_const_summary <- as.data.frame(GRWR_LACO_const) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  summarise(mean = mean(GRWR),
            se = se(GRWR)) # -0.2218259
GRWR_LACO_25EG_all_summary <- as.data.frame(GRWR_LACO_25EG_all) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  summarise(mean = mean(GRWR),
            se = se(GRWR))# -0.02672747
GRWR_LACO_50EG_all_summary <- as.data.frame(GRWR_LACO_50EG_all) %>% 
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  summarise(mean = mean(GRWR),
            se = se(GRWR)) # 0.2306799
GRWR_LACO_75EG_all_summary <- as.data.frame(GRWR_LACO_75EG_all) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  summarise(mean = mean(GRWR),
            se = se(GRWR)) # 0.612166
GRWR_simulated_all <- as.data.frame(rbind(GRWR_LACO_const_summary, GRWR_LACO_50EG_all_summary,  GRWR_LACO_75EG_all_summary)) %>%
  mutate(treatment = c("0%", "50%", "75%"))
GRWR_simulated_all$treatment <- ordered(GRWR_simulated_all$treatment, levels = c("0%", "50%", "75%"))

#Step 3. Plot modeled abundance and GRWR

#FIGURE4
sim_timeseries <- ggplot(summary_grass_sim_LACO%>%filter(type != "reduced25EG_LACO"), aes(x = as.numeric(time), y = mean_LACO, group = type)) +
                          geom_point() +
                          geom_line(aes(linetype = type), size = 1.4) +
                          geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 1) +
                          theme(text = element_text(size=18),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                axis.text = element_text(size = 18),
                                axis.line = element_line(colour = "black"),
                                legend.position = c(0.7, 0.8)) +
                          labs(x = "Year", y = bquote(Modeled~italic(L.~conj.)~Density~(stems/m^2))) +
                          scale_linetype_manual(name = "Percent Reduction in Exotic Grasses", 
                                               labels = c("0%", "50%", "75%"),
                                               values = c("solid", "twodash", "dotted"))

sim_GRWR <- ggplot(GRWR_simulated_all , aes(x = treatment, y = mean))+
                          geom_bar(stat = "identity")+
                          theme(text = element_text(size=18),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                axis.text = element_text(size = 18),
                                axis.line = element_line(colour = "black"),
                                legend.position = "bottom")+
                          geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.4, alpha = 0.9, size = 1) +
                          labs(y = "Average Low Density Growth Rate", x = "Percent Reduction in Exotic Grasses")+
                          geom_hline(yintercept = 0)

ggarrange(sim_timeseries, sim_GRWR,  ncol = 2, nrow = 1, 
          labels = c("(a)", "(b)"), widths = c(0.6, 0.4),
          #common.legend = TRUE, legend = "bottom", 
          font.label = list(size = 20))



################################################
#Which year had the greatest effect of removal?#
################################################
# #Calculate the effect size of trt = (mean_LACO_EGreduced-mean_LACO_no_removal)/sd_LACO_no_removal
# removal_eff <- summary_grass_sim_LACO %>%
#   select(time, type, mean_LACO, sd_LACO) %>%
#   pivot_wider(names_from = type, values_from = c(mean_LACO, sd_LACO)) %>%
#   select(-sd_LACO_reduced25EG_LACO, -sd_LACO_reduced50EG_LACO, -sd_LACO_reduced75EG_LACO) %>%
#   mutate(eff_50 = (mean_LACO_reduced50EG_LACO-mean_LACO_predicted_LACO)/sd_LACO_predicted_LACO,
#          eff_75 = (mean_LACO_reduced75EG_LACO-mean_LACO_predicted_LACO)/sd_LACO_predicted_LACO) %>%
#   select(-mean_LACO_predicted_LACO, -mean_LACO_reduced25EG_LACO, -mean_LACO_reduced50EG_LACO, -mean_LACO_reduced75EG_LACO,
#          -sd_LACO_predicted_LACO)
#write.csv(removal_eff, "C:\\Users\\Lina\\Desktop\\Repositories\\sToration-vernal-pools\\data_analysis\\Table2.csv", row.names = FALSE )

###########################################
#Show the relationship between LACO and EG#
###########################################
#Supplemental Figure 2
ggplot(const_dummy_join, aes(x = sum_EG, y= LACOdens))+
  geom_point()+
  scale_y_log10(limits=c(0.1,6000))+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_smooth(method = "lm")+
  labs(x = "Sum of exotic grass cover (%)", y = "LACO density (log)")


