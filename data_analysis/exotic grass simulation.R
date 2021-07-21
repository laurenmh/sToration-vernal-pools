#####################################################
#Would adaptive management improve LACO populations?#
#####################################################

#Simulate LACO populations with and without exotic grass removal

library(tidyverse)
library(ggplot2)
library(ggpubr)
#-----------------------
#Goal: Simulate exotic grasses (EG) removal to promote LACO persistence

#Step 1. Simulate EG removal
#Step 2. Average the growth rates of LACO over time for all simulation scenarios
#Step 3. Plot simulated GRWR

#Step 1. Remove EG in all years

#Remove 25% of EG in all years
sim_LACO <- matrix(nrow = 2000, ncol = 17)
mult.25 <- function(x)(x*0.75)
reduced25EGcover_all <- const_com_control[2,2:18] %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.25)
LACO_25EG_all <- bh.sim.control(
  LACO = 1,
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
LACO_50EG_all <- bh.sim.control(
  LACO = 1,
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
LACO_75EG_all <- bh.sim.control(
  LACO = 1,
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
LACO_100EG_all <- bh.sim.control(
  LACO = 1,
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

#FIGURE4
GRWR_simulated_all <- as.data.frame(rbind(GRWR_LACO_const_summary, GRWR_LACO_50EG_all_summary,  GRWR_LACO_75EG_all_summary)) %>%
  mutate(treatment = c("0%", "50%", "75%"))
GRWR_simulated_all$treatment <- ordered(GRWR_simulated_all$treatment, levels = c("0%", "50%", "75%"))

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

# Simulation model:
sim_obs_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of LACO seed counts
sim_mu <- matrix(nrow = sim_n_pools, ncol = sim_n_years) #empty matrix of mean LACO seed counts
sim_stem_LACO <- matrix(nrow = sim_n_pools, ncol = sim_n_years) # empty metrix of LACO stem counts

bh.formula <- function(sim_obs_LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g){
  sim_obs_LACO*lambda/(1+sim_obs_LACO*aii+EG*a1+ERVA*a2+NF*a3)+s*(1-g)*sim_obs_LACO/g
} #this is the modified Beverton-Holt model we'll use for LACO stem counts

#now simulate LACO stem counts
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
                                aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                lambda = lambda[j-1], s = s, g = g)
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
                                  aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                  lambda = lambda[j-1], s = s, g = g)
      }
      else {
        sim_mu[i,j] <- bh.formula(sim_obs_LACO = sim_obs_LACO[i,j-2]*lambda[j-2]/(1+sim_obs_LACO[i,j-2]*aii[j-2]+EG[i,j-2]*a1[j-2]+ERVA[i,j-2]*a2[j-2]+NF[i,j-2]*a3[j-2])+s*(1-g)*sim_obs_LACO[i,j-2]/g,
                                  EG = EG[i,j-1], ERVA = ERVA[i,j-1], NF = NF[i,j-1],
                                  aii = aii[j-1], a1 = a1[j-1], a2 = a2[j-1], a3 = a3[j-1],
                                  lambda = lambda[j-1], s = s, g = g)
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

#plot them!
ggplot(grass_sim_LACO, aes(x = time, y = log_LACO, col = type)) +
  geom_jitter() 
  

ggplot(grass_sim_LACO%>%filter(type != "reduced25EG_LACO"), aes(x = as.numeric(time), y = LACO, col = type)) +
  geom_smooth(method = "gam",formula = y ~s(x))+
  scale_y_log10()+
  theme_bw()

se<-function(x){
  sd(x)/sqrt(length(x))
} # this is a function for calculating standard error

summary_grass_sim_LACO <- grass_sim_LACO %>%
  group_by(time, type) %>%
  summarise(mean_log_LACO = mean(log_LACO),
            se_log_LACO = se(log_LACO),
            mean_LACO = mean(LACO),
            se_LACO = se(LACO),
            sd_LACO = sd(LACO))

##########
#Figure 4# 
##########
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

summary_grass_sim_LACO %>%
  group_by(type) %>%
  summarize(mean = mean(mean_LACO))

anova(lm(mean_LACO ~ time, data = summary_grass_sim_LACO))

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


###################################
#How grass cover changes over time#
###################################
#Graph mean of exotic grass cover timeseries
EG_summary <- const_dummy_join %>%
  group_by(Year) %>%
  summarise(EG_mean = mean(sum_EG),
            EG_se = se(sum_EG)) %>%
  mutate(type = "constructed")

EG_ref_summary <- ref_com %>%
  mutate(sum_EG = `BRHO`+`HOMA`+ `LOMU`)%>%
  group_by(Year) %>%
  summarise(EG_mean = mean(sum_EG),
            EG_se = se(sum_EG)) %>%
  mutate(type = "reference")
EG_summary_join <- rbind(EG_summary, EG_ref_summary) 

# SUPPLEMENTAL FIGURE 1
ggplot(EG_summary_join, aes(x = Year, y = EG_mean, group= type)) +
  geom_point(aes(col =type))+
  geom_line(aes(col =type))+
  geom_errorbar(aes(ymin = EG_mean-EG_se, ymax = EG_mean+EG_se, col =type)) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom")+
  labs(x = "Time (year)", y = "Mean exotic grass cover (%)") +
  scale_color_manual(name = "", values = c("#000000", "#888888"))

###########################################
#Show the relationship between LACO and EG#
###########################################
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


#################################
#When should grasses be removed?# 
#################################
#Try manipulating the timing of exotic grass removal because doing it every year is not very practical.

# #Remove 75% of EG in 2002
# reduced75EGcover2002 <- sumEGcover %>%
#   mutate_at(c("2002"), mult.75)
# #Remove 75% of EG in 2007
# reduced75EGcover2007 <- sumEGcover %>%
#   mutate_at(c("2007"), mult.75)
# #Remove 75% of EG in 2012
# reduced75EGcover2012 <- sumEGcover %>%
#   mutate_at(c("2012"), mult.75)
# 
# #Simulate with exotic grass removal
# reduced75EG2002_LACO <- bh.sim(n_pools = n_pools,
#                            seedtrt = as.matrix(seedtrt[,4:6]),
#                            EG = as.matrix(reduced75EGcover2002),
#                            ERVA = as.matrix(ERVAdens),
#                            NF = as.matrix(sumNFcover),
#                            aii = alpha_LACO_mean[,5],
#                            a1 = alpha_EG_mean[,5],
#                            a2 = alpha_ERVA_mean[,5], 
#                            a3 = alpha_NF_mean[,5],
#                            lambda = lambda_mean[,5],
#                            s = s_mean[,5],
#                            g = 0.7,
#                            glow = 0.2)
# 
# reduced75EG2007_LACO <- bh.sim(n_pools = n_pools,
#                            seedtrt = as.matrix(seedtrt[,4:6]),
#                            EG = as.matrix(reduced75EGcover2007),
#                            ERVA = as.matrix(ERVAdens),
#                            NF = as.matrix(sumNFcover),
#                            aii = alpha_LACO_mean[,5],
#                            a1 = alpha_EG_mean[,5],
#                            a2 = alpha_ERVA_mean[,5], 
#                            a3 = alpha_NF_mean[,5],
#                            lambda = lambda_mean[,5],
#                            s = s_mean[,5],
#                            g = 0.7,
#                            glow = 0.2)
# 
# reduced75EG2012_LACO <- bh.sim(n_pools = n_pools,
#                            seedtrt = as.matrix(seedtrt[,4:6]),
#                            EG = as.matrix(reduced75EGcover2012),
#                            ERVA = as.matrix(ERVAdens),
#                            NF = as.matrix(sumNFcover),
#                            aii = alpha_LACO_mean[,5],
#                            a1 = alpha_EG_mean[,5],
#                            a2 = alpha_ERVA_mean[,5], 
#                            a3 = alpha_NF_mean[,5],
#                            lambda = lambda_mean[,5],
#                            s = s_mean[,5],
#                            g = 0.7,
#                            glow = 0.2)
# 
# #change column names to years
# colnames(reduced75EG2002_LACO) <- years
# colnames(reduced75EG2007_LACO) <- years
# colnames(reduced75EG2012_LACO) <- years
# 
# #combine the tables
# reduced75EG2002_LACO <- as.data.frame(reduced75EG2002_LACO) %>% 
#   mutate(Pool = row_number()) %>%
#   gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, `2007`, `2008`, `2009`, `2010`,
#          `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = reduced75EG2002)
# reduced75EG2007_LACO <- as.data.frame(reduced75EG2007_LACO) %>% 
#   mutate(Pool = row_number()) %>%
#   gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, `2007`, `2008`, `2009`, `2010`,
#          `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = reduced75EG2007)
# reduced75EG2012_LACO <- as.data.frame(reduced75EG2012_LACO) %>% 
#   mutate(Pool = row_number()) %>%
#   gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`, `2007`, `2008`, `2009`, `2010`,
#          `2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = time, value = reduced75EG2012)
# 
# grass_sim_LACO_timing <- left_join(left_join(left_join(predicted_LACO, reduced75EG2002_LACO, by = c("Pool", "time")), reduced75EG2007_LACO, by = c("Pool", "time")), reduced75EG2012_LACO, by = c("Pool", "time")) %>%
#   gather(`predicted_LACO`, `reduced75EG2002`, `reduced75EG2007`, `reduced75EG2012`, key = type, value = LACO) %>%
#   mutate(log_LACO = log(LACO)) %>%
#   mutate_if(is.numeric, ~replace(., is.infinite(.), 0))
# 
# #plot them
# ggplot(grass_sim_LACO_timing, aes(x = time, y = log_LACO, col = type)) +
#   geom_jitter()
# 
# summary_grass_sim_LACO_timing <- grass_sim_LACO_timing %>%
#   group_by(time, type) %>%
#   summarise(mean_log_LACO = mean(log_LACO),
#             se_log_LACO = se(log_LACO),
#             mean_LACO = mean(LACO),
#             se_LACO = se(LACO),
#             sd_LACO = sd(LACO))
# 
# ggplot(summary_grass_sim_LACO_timing, aes(x = time, y = mean_LACO, col = type)) +
#   geom_point() +
#   geom_line(aes(x = time, y = mean_LACO, group = type)) +
#   geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 1) +
#   theme_bw() +
#   labs(x = "Year", y = "Mean LACO density") +
#   scale_color_discrete(name = "Treatment", labels = c("No grass removal", "2002 removal", "2007 removal", "2012 removal"))