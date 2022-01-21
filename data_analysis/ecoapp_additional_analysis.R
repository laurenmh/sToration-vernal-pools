#Additional analyses from Eco App review 1/20/22

library(ggplot2)
library(ggpubr)
library(HDInterval)
library(tidyverse)
library(dplyr)
library(lme4)

#Function for standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 


#Data: 
#Run "prep data before modeling.R" to get 'const_dummy_join'
#Run "analysis/complex_belowground_v5.R" to get 'BH_fit_AB', 'BH_fit_BA', 'BH_fit_LALA', and 'BH_fit_LANo'


###Compare the effects of pool size and seed addition treatments on restoration outcomes

##LACO density
const_dummy_join$Year <- as.numeric(const_dummy_join$Year)
const_dummy_join <- const_dummy_join %>%
  mutate(Treatment = paste(Treatment.1999, Treatment.2000, sep = "-")) %>%
  filter(Year < 2016) 

const_LACO_size <- const_dummy_join %>%
  select(LACOdens, Year, Size, Pool) %>%
  group_by(Size, Year) %>%
  summarise(mean_LACO = mean(LACOdens), se_LACO = se(LACOdens))

fig_LACO_size <- ggplot(const_LACO_size, aes(x = Year, y = mean_LACO, col = Size))+
                    geom_point()+
                    geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 0.8) +
                    geom_line()+
                    theme(text = element_text(size=16),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"))+
                    ylab(bquote(LACO~Density~(stems/m^2)))+
                    #scale_y_log10()+
                    scale_color_discrete(name = "Pool size", labels = c("large", "medium", "small"))

const_LACO_trt <- const_dummy_join %>%
  select(LACOdens, Year, Treatment, Pool) %>%
  group_by(Treatment, Year) %>%
  summarise(mean_LACO = mean(LACOdens), se_LACO = se(LACOdens))

fig_LACO_trt <- ggplot(const_LACO_trt, aes(x = Year, y = mean_LACO, col = Treatment))+
                    geom_point()+
                    geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 0.8) +
                    geom_line()+
                    theme(text = element_text(size=16),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(colour = "black"))+
                    #scale_y_log10()+
                    ylab(bquote(LACO~Density~(stems/m^2)))



##LACO lambda

#Extract lambda estimates for seeding treatment AB
Post_AB <- rstan::extract(BH_fit_AB) 
CI_lambda_AB <-  as.data.frame(HDInterval::hdi(Post_AB$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_AB <- as.data.frame(Post_AB$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(treatment = "Group A-Group B") %>%
  full_join(., CI_lambda_AB)

#Extract lambda estimates for seeding treatment BA
Post_BA <- rstan::extract(BH_fit_BA) 
CI_lambda_BA <-  as.data.frame(HDInterval::hdi(Post_BA$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_BA <- as.data.frame(Post_BA$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(treatment = "Group B-Group A") %>%
  full_join(., CI_lambda_BA)

#Extract lambda estimates for seeding treatment LALA
Post_LALA <- rstan::extract(BH_fit_LALA) 
CI_lambda_LALA <-  as.data.frame(HDInterval::hdi(Post_LALA$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_LALA <- as.data.frame(Post_LALA$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(treatment = "Lasthenia-Lasthenia") %>%
  full_join(., CI_lambda_LALA)

#Extract lambda estimates for seeding treatment LANo
Post_LANo <- rstan::extract(BH_fit_LANo) 
CI_lambda_LANo <-  as.data.frame(HDInterval::hdi(Post_LANo$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_LANo <- as.data.frame(Post_LANo$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(treatment = "Lasthenia-NO Lasthenia") %>%
  full_join(., CI_lambda_LANo)

lambda_const_trt <- rbind(lambda_const_AB, lambda_const_BA, lambda_const_LALA, lambda_const_LANo)
lambda_const_trt$Year <- as.numeric(lambda_const_trt$Year)

fig_lambda_trt <- ggplot(lambda_const_trt, aes(x = Year, y = mean, col = treatment))+
  geom_point() +
  geom_line(size=1)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        #legend.position = c(0.2, 0.9), 
        axis.title = element_text(size = 14))+
  ylab(bquote(Intrinsic~Growth~Rate ~(lambda[t])))+
  scale_x_continuous(name = NULL, limits = c(1999.5,2015.5))+
  scale_color_discrete(name = "Treatment")


##LACO low-density growth rate

#Calculate the annual growth rate of LACO in restored pools when one LACO is introduced into a stable community.

#Extract parameters for seeding treatment AB
Post_AB <- rstan::extract(BH_fit_AB)
alpha_LACO_AB <- as.matrix(Post_AB$alpha_LACO)
alpha_EG_AB <- as.matrix(Post_AB$alpha_EG)
alpha_ERVA_AB <- as.matrix(Post_AB$alpha_ERVA)
alpha_NF_AB <- as.matrix(Post_AB$alpha_NF)
lambda_AB <- as.matrix(Post_AB$lambda)
s_AB <- as.matrix(Post_AB$survival_LACO)
sim_LACO_AB <- matrix(nrow = 2000, ncol = 15)

#Set up the population model for LACO
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

#Plug in the stable equilibrium freq and parameters in the model
LACO_const_AB <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,3:17]),
  ERVA = as.numeric(const_com_control[1,3:17]),
  NF = as.numeric(const_com_control[3,3:17]),
  aii = alpha_LACO_AB,
  a1 = alpha_EG_AB,
  a2 = alpha_ERVA_AB,
  a3 = alpha_NF_AB,
  lambda = lambda_AB,
  s = s_AB,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_const_all <- log(LACO_const) #2000-2016 #log transform
GRWR_LACO_const <- GRWR_LACO_const_all[,1:15] #2000-2014 truncate the last two points 


###Exotic grass cover time series colored by pool size and treatment 
summary(lmer(sum_EG~Size + (1|Pool), const_dummy_join))
anova(lm(sum_EG~Size + (1|Pool), const_dummy_join))

const_EG_size <- const_dummy_join %>%
  select(sum_EG, Year, Size, Pool) %>%
  group_by(Size, Year) %>%
  summarise(mean_EG = mean(sum_EG), se_EG = se(sum_EG))

ggplot(const_EG_size, aes(x = Year, y = mean_EG, col = Size))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_EG-se_EG, ymax = mean_EG+se_EG), width = 0.4, alpha = 0.9, size = 0.8) +
  geom_line()+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  ylab(bquote(Exotic~Grass~Cover~('%')))+
  scale_color_discrete(name = "Pool size", labels = c("large", "medium", "small"))

summary(lmer(sum_EG~Treatment.2000 + (1|Pool), const_dummy_join))
anova(lm(sum_EG~Treatment.2000 + (1|Pool), const_dummy_join))

const_EG_trt <- const_dummy_join %>%
  select(sum_EG, Year, Treatment, Pool) %>%
  group_by(Treatment, Year) %>%
  summarise(mean_EG = mean(sum_EG), se_EG = se(sum_EG))

ggplot(const_EG_trt, aes(x = Year, y = mean_EG, col = Treatment))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_EG-se_EG, ymax = mean_EG+se_EG), width = 0.4, alpha = 0.9, size = 0.8) +
  geom_line()+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  ylab(bquote(Exotic~Grass~Cover~('%')))



