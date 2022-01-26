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

anova(lm(LACOdens~Size + (1|Pool), const_dummy_join))
anova(lm(LACOdens~Treatment.2000 + (1|Pool), const_dummy_join))

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
                          legend.position = c(0.8, 0.8),
                          axis.line = element_line(colour = "black"))+
                    ylab(bquote(LACO~Density~(stems/m^2)))+
                    scale_x_continuous(name = NULL,
                     limits = c(1999.5,2015.5))+
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
                          legend.position = c(0.8,0.8),
                          axis.line = element_line(colour = "black"))+
                    scale_x_continuous(name = NULL,
                     limits = c(1999.5,2015.5))+
                    #scale_y_log10()+
                    ylab(bquote(LACO~Density~(stems/m^2)))



##LACO lambda

#Extract lambda estimates for seeding treatment AB
Post_AB <- rstan::extract(BH_fit_AB) 
CI_lambda_AB <-  as.data.frame(HDInterval::hdi(Post_AB$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2017)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_AB <- as.data.frame(Post_AB$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(treatment = "Group A-Group B") %>%
  full_join(., CI_lambda_AB)

#Extract lambda estimates for seeding treatment BA
Post_BA <- rstan::extract(BH_fit_BA) 
CI_lambda_BA <-  as.data.frame(HDInterval::hdi(Post_BA$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2017)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_BA <- as.data.frame(Post_BA$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(treatment = "Group B-Group A") %>%
  full_join(., CI_lambda_BA)

#Extract lambda estimates for seeding treatment LALA
Post_LALA <- rstan::extract(BH_fit_LALA) 
CI_lambda_LALA <-  as.data.frame(HDInterval::hdi(Post_LALA$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2017)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_LALA <- as.data.frame(Post_LALA$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(treatment = "Lasthenia-Lasthenia") %>%
  full_join(., CI_lambda_LALA)

#Extract lambda estimates for seeding treatment LANo
Post_LANo <- rstan::extract(BH_fit_LANo) 
CI_lambda_LANo <-  as.data.frame(HDInterval::hdi(Post_LANo$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2017)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_LANo <- as.data.frame(Post_LANo$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(treatment = "Lasthenia-NO Lasthenia") %>%
  full_join(., CI_lambda_LANo)

lambda_const_trt <- rbind(lambda_const_AB, lambda_const_BA, lambda_const_LALA, lambda_const_LANo)
lambda_const_trt$Year <- as.numeric(lambda_const_trt$Year)

fig_lambda_trt <- ggplot(lambda_const_trt%>%filter(Year < 2016), aes(x = Year, y = mean, col = treatment))+
  geom_point() +
  geom_line(size=1)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.title = element_text(size = 14))+
  ylab(bquote(Intrinsic~Growth~Rate ~(lambda[t])))+
  scale_x_continuous(name = NULL, limits = c(1999.5,2015.5))+
  scale_color_discrete(name = "Treatment")

#Extract lambda estimates for small pools
Post_Pool_s <- rstan::extract(BH_fit_Pool_s) 
CI_lambda_Pool_s <-  as.data.frame(HDInterval::hdi(Post_Pool_s$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2017)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_Pool_s <- as.data.frame(Post_Pool_s$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(Size = "small") %>%
  full_join(., CI_lambda_Pool_s)

#Extract lambda estimates for medium pools
Post_Pool_m <- rstan::extract(BH_fit_Pool_m) 
CI_lambda_Pool_m <-  as.data.frame(HDInterval::hdi(Post_Pool_m$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2017)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_Pool_m <- as.data.frame(Post_Pool_m$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(Size = "medium") %>%
  full_join(., CI_lambda_Pool_m)

#Extract lambda estimates for medium pools
Post_Pool_l <- rstan::extract(BH_fit_Pool_l) 
CI_lambda_Pool_l <-  as.data.frame(HDInterval::hdi(Post_Pool_l$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2017)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const_Pool_l <- as.data.frame(Post_Pool_l$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(Size = "large") %>%
  full_join(., CI_lambda_Pool_l)

lambda_const_size <- rbind(lambda_const_Pool_s, lambda_const_Pool_m, lambda_const_Pool_l)
lambda_const_size$Year <- as.numeric(lambda_const_size$Year)

fig_lambda_size <- ggplot(lambda_const_size%>%filter(Year < 2016), aes(x = Year, y = mean, col = Size))+
  geom_point() +
  geom_line(size=1)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.title = element_text(size = 14))+
  ylab(bquote(Intrinsic~Growth~Rate ~(lambda[t])))+
  scale_x_continuous(name = NULL, limits = c(1999.5,2015.5))+
  scale_color_discrete(name = "Pool Size")

##LACO low-density growth rate

#Calculate the annual growth rate of LACO in restored pools when one LACO is introduced into a stable community.

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

#Extract parameters for each seeding treatment
alpha_LACO_AB <- as.matrix(Post_AB$alpha_LACO)
alpha_EG_AB <- as.matrix(Post_AB$alpha_EG)
alpha_ERVA_AB <- as.matrix(Post_AB$alpha_ERVA)
alpha_NF_AB <- as.matrix(Post_AB$alpha_NF)
lambda_AB <- as.matrix(Post_AB$lambda)
s_AB <- as.matrix(Post_AB$survival_LACO)
sim_LACO_AB <- matrix(nrow = 2000, ncol = 17)

alpha_LACO_BA <- as.matrix(Post_BA$alpha_LACO)
alpha_EG_BA <- as.matrix(Post_BA$alpha_EG)
alpha_ERVA_BA <- as.matrix(Post_BA$alpha_ERVA)
alpha_NF_BA <- as.matrix(Post_BA$alpha_NF)
lambda_BA <- as.matrix(Post_BA$lambda)
s_BA <- as.matrix(Post_BA$survival_LACO)
sim_LACO_BA <- matrix(nrow = 2000, ncol = 17)

alpha_LACO_LALA <- as.matrix(Post_LALA$alpha_LACO)
alpha_EG_LALA <- as.matrix(Post_LALA$alpha_EG)
alpha_ERVA_LALA <- as.matrix(Post_LALA$alpha_ERVA)
alpha_NF_LALA <- as.matrix(Post_LALA$alpha_NF)
lambda_LALA <- as.matrix(Post_LALA$lambda)
s_LALA <- as.matrix(Post_LALA$survival_LACO)
sim_LACO_LALA <- matrix(nrow = 2000, ncol = 17)

alpha_LACO_LANo <- as.matrix(Post_LANo$alpha_LACO)
alpha_EG_LANo <- as.matrix(Post_LANo$alpha_EG)
alpha_ERVA_LANo <- as.matrix(Post_LANo$alpha_ERVA)
alpha_NF_LANo <- as.matrix(Post_LANo$alpha_NF)
lambda_LANo <- as.matrix(Post_LANo$lambda)
s_LANo <- as.matrix(Post_LANo$survival_LACO)
sim_LACO_LANo <- matrix(nrow = 2000, ncol = 17)

#Plug in the stable equilibrium freq and parameters in the model
LACO_const_AB <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:18]),
  ERVA = as.numeric(const_com_control[1,2:18]),
  NF = as.numeric(const_com_control[3,2:18]),
  aii = alpha_LACO_AB,
  a1 = alpha_EG_AB,
  a2 = alpha_ERVA_AB,
  a3 = alpha_NF_AB,
  lambda = lambda_AB,
  s = s_AB,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_const_AB_all <- log(LACO_const_AB) #2000-2016 #log transform
GRWR_LACO_const_AB <- GRWR_LACO_const_AB_all[,1:15] #2000-2014 truncate the last two points 

LACO_const_BA <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:18]),
  ERVA = as.numeric(const_com_control[1,2:18]),
  NF = as.numeric(const_com_control[3,2:18]),
  aii = alpha_LACO_BA,
  a1 = alpha_EG_BA,
  a2 = alpha_ERVA_BA,
  a3 = alpha_NF_BA,
  lambda = lambda_BA,
  s = s_BA,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_const_BA_all <- log(LACO_const_BA) #2000-2016 #log transform
GRWR_LACO_const_BA <- GRWR_LACO_const_BA_all[,1:15] #2000-2014 truncate the last two points 

LACO_const_LALA <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:18]),
  ERVA = as.numeric(const_com_control[1,2:18]),
  NF = as.numeric(const_com_control[3,2:18]),
  aii = alpha_LACO_LALA,
  a1 = alpha_EG_LALA,
  a2 = alpha_ERVA_LALA,
  a3 = alpha_NF_LALA,
  lambda = lambda_LALA,
  s = s_LALA,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_const_LALA_all <- log(LACO_const_LALA) #2000-2016 #log transform
GRWR_LACO_const_LALA <- GRWR_LACO_const_LALA_all[,1:15] #2000-2014 truncate the last two points 

LACO_const_LANo <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:18]),
  ERVA = as.numeric(const_com_control[1,2:18]),
  NF = as.numeric(const_com_control[3,2:18]),
  aii = alpha_LACO_LANo,
  a1 = alpha_EG_LANo,
  a2 = alpha_ERVA_LANo,
  a3 = alpha_NF_LANo,
  lambda = lambda_LANo,
  s = s_LANo,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_const_LANo_all <- log(LACO_const_LANo) #2000-2016 #log transform
GRWR_LACO_const_LANo <- GRWR_LACO_const_LANo_all[,1:15] #2000-2014 truncate the last two points 

#calculate 95% CI and mean GRWR
CI_GRWR_AB <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const_AB, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const_AB <- as.data.frame(GRWR_LACO_const_AB) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(treatment = "Group A-Group B") %>%
  full_join(., CI_GRWR_AB)

CI_GRWR_BA <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const_BA, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const_BA <- as.data.frame(GRWR_LACO_const_BA) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(treatment = "Group B-Group A") %>%
  full_join(., CI_GRWR_BA)

CI_GRWR_LALA <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const_LALA, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const_LALA <- as.data.frame(GRWR_LACO_const_LALA) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(treatment = "Lasthenia-Lasthenia") %>%
  full_join(., CI_GRWR_LALA)

CI_GRWR_LANo <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const_LANo, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const_LANo <- as.data.frame(GRWR_LACO_const_LANo) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(treatment = "Lasthenia-NO Lasthenia") %>%
  full_join(., CI_GRWR_LANo)

GRWR_time_trt <- rbind(GRWR_time_const_AB, GRWR_time_const_BA, GRWR_time_const_LALA, GRWR_time_const_LANo) 
GRWR_time_trt$Year <- as.numeric(GRWR_time_trt$Year)

fig_GRWR_trt <- ggplot(GRWR_time_trt, aes(x = Year, y = mean, col = treatment))+
  geom_point() +
  geom_line(size=1)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title = element_text(size = 14))+
  ylab(bquote(Low~Density~Growth~Rate~(italic(r[t]))))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(name = NULL,
                     limits = c(1999.5,2015.5))#+
  #scale_color_manual(name = "", values = c("#000000", "#888888"))

# Average the growth rates of LACO over time for each trt.
GRWR_LACO_const_AB_mean <- as.data.frame(GRWR_LACO_const_AB) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) #mean of GRWR across years 

GRWR_LACO_const_AB_summary <- GRWR_LACO_const_AB_mean %>%
  summarize(Mean = mean(mean_GRWR), CI = hdi(mean_GRWR, credMass = 0.95)) %>% #take the 95% CI across iterations
  mutate(name = c("lowCI", "upCI")) %>% #top values is the low CI and bottom value is the high CI
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for Group A-Group B seeding treatment = -0.903452 (CI: -1.075655, -0.7611101)

GRWR_LACO_const_BA_mean <- as.data.frame(GRWR_LACO_const_BA) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) #mean of GRWR across years 

GRWR_LACO_const_BA_summary <- GRWR_LACO_const_BA_mean %>%
  summarize(Mean = mean(mean_GRWR), CI = hdi(mean_GRWR, credMass = 0.95)) %>% #take the 95% CI across iterations
  mutate(name = c("lowCI", "upCI")) %>% #top values is the low CI and bottom value is the high CI
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for Group B-Group A seeding treatment = -0.4835636 (CI: -0.6493003, -0.3221352)

GRWR_LACO_const_LALA_mean <- as.data.frame(GRWR_LACO_const_LALA) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) #mean of GRWR across years 

GRWR_LACO_const_LALA_summary <- GRWR_LACO_const_LALA_mean %>%
  summarize(Mean = mean(mean_GRWR), CI = hdi(mean_GRWR, credMass = 0.95)) %>% #take the 95% CI across iterations
  mutate(name = c("lowCI", "upCI")) %>% #top values is the low CI and bottom value is the high CI
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for Lasthenia-Lasthenia seeding treatment = -0.6884356 (CI:-0.8465281, -0.52569)

GRWR_LACO_const_LANo_mean <- as.data.frame(GRWR_LACO_const_LANo) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) #mean of GRWR across years 

GRWR_LACO_const_LANo_summary <- GRWR_LACO_const_LANo_mean %>%
  summarize(Mean = mean(mean_GRWR), CI = hdi(mean_GRWR, credMass = 0.95)) %>% #take the 95% CI across iterations
  mutate(name = c("lowCI", "upCI")) %>% #top values is the low CI and bottom value is the high CI
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for Lasthenia-No Lasthenia seeding treatment = -0.5845346 (CI: -0.7934202, -0.3630242)

# Plot LACO density, lambda, and GRWR timeseries by trt
ggarrange(fig_LACO_trt, fig_lambda_trt, fig_GRWR_trt,  ncol = 1, labels = c("(a)", "(b)", "(c)"))

#Extract parameters for each pool size
alpha_LACO_Pool_s <- as.matrix(Post_Pool_s$alpha_LACO)
alpha_EG_Pool_s <- as.matrix(Post_Pool_s$alpha_EG)
alpha_ERVA_Pool_s <- as.matrix(Post_Pool_s$alpha_ERVA)
alpha_NF_Pool_s <- as.matrix(Post_Pool_s$alpha_NF)
lambda_Pool_s <- as.matrix(Post_Pool_s$lambda)
s_Pool_s <- as.matrix(Post_Pool_s$survival_LACO)
sim_LACO_Pool_s <- matrix(nrow = 2000, ncol = 17)

alpha_LACO_Pool_m <- as.matrix(Post_Pool_m$alpha_LACO)
alpha_EG_Pool_m <- as.matrix(Post_Pool_m$alpha_EG)
alpha_ERVA_Pool_m <- as.matrix(Post_Pool_m$alpha_ERVA)
alpha_NF_Pool_m <- as.matrix(Post_Pool_m$alpha_NF)
lambda_Pool_m <- as.matrix(Post_Pool_m$lambda)
s_Pool_m <- as.matrix(Post_Pool_m$survival_LACO)
sim_LACO_Pool_m <- matrix(nrow = 2000, ncol = 17)

alpha_LACO_Pool_l <- as.matrix(Post_Pool_l$alpha_LACO)
alpha_EG_Pool_l <- as.matrix(Post_Pool_l$alpha_EG)
alpha_ERVA_Pool_l <- as.matrix(Post_Pool_l$alpha_ERVA)
alpha_NF_Pool_l <- as.matrix(Post_Pool_l$alpha_NF)
lambda_Pool_l <- as.matrix(Post_Pool_l$lambda)
s_Pool_l <- as.matrix(Post_Pool_l$survival_LACO)
sim_LACO_Pool_l <- matrix(nrow = 2000, ncol = 17)

#Plug in the stable equilibrium freq and parameters in the model
LACO_const_Pool_s <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:18]),
  ERVA = as.numeric(const_com_control[1,2:18]),
  NF = as.numeric(const_com_control[3,2:18]),
  aii = alpha_LACO_Pool_s,
  a1 = alpha_EG_Pool_s,
  a2 = alpha_ERVA_Pool_s,
  a3 = alpha_NF_Pool_s,
  lambda = lambda_Pool_s,
  s = s_Pool_s,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_const_Pool_s_all <- log(LACO_const_Pool_s) #2000-2016 #log transform
GRWR_LACO_const_Pool_s <- GRWR_LACO_const_Pool_s_all[,1:15] #2000-2014 truncate the last two points 

LACO_const_Pool_m <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:18]),
  ERVA = as.numeric(const_com_control[1,2:18]),
  NF = as.numeric(const_com_control[3,2:18]),
  aii = alpha_LACO_Pool_m,
  a1 = alpha_EG_Pool_m,
  a2 = alpha_ERVA_Pool_m,
  a3 = alpha_NF_Pool_m,
  lambda = lambda_Pool_m,
  s = s_Pool_m,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_const_Pool_m_all <- log(LACO_const_Pool_m) #2000-2016 #log transform
GRWR_LACO_const_Pool_m <- GRWR_LACO_const_Pool_m_all[,1:15] #2000-2014 truncate the last two points 

LACO_const_Pool_l <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:18]),
  ERVA = as.numeric(const_com_control[1,2:18]),
  NF = as.numeric(const_com_control[3,2:18]),
  aii = alpha_LACO_Pool_l,
  a1 = alpha_EG_Pool_l,
  a2 = alpha_ERVA_Pool_l,
  a3 = alpha_NF_Pool_l,
  lambda = lambda_Pool_l,
  s = s_Pool_l,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_const_Pool_l_all <- log(LACO_const_Pool_l) #2000-2016 #log transform
GRWR_LACO_const_Pool_l <- GRWR_LACO_const_Pool_l_all[,1:15] #2000-2014 truncate the last two points 

#calculate 95% CI and mean GRWR
CI_GRWR_Pool_s <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const_Pool_s, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const_Pool_s <- as.data.frame(GRWR_LACO_const_Pool_s) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(Size = "small") %>%
  full_join(., CI_GRWR_Pool_s)

CI_GRWR_Pool_m <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const_Pool_m, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const_Pool_m <- as.data.frame(GRWR_LACO_const_Pool_s) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(Size = "medium") %>%
  full_join(., CI_GRWR_Pool_m)

CI_GRWR_Pool_l <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const_Pool_l, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const_Pool_l <- as.data.frame(GRWR_LACO_const_Pool_l) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(Size = "large") %>%
  full_join(., CI_GRWR_Pool_l)

GRWR_time_size <- rbind(GRWR_time_const_Pool_s, GRWR_time_const_Pool_m, GRWR_time_const_Pool_l) 
GRWR_time_size$Year <- as.numeric(GRWR_time_size$Year)

fig_GRWR_size <- ggplot(GRWR_time_size, aes(x = Year, y = mean, col = Size))+
  geom_point() +
  geom_line(size=1)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title = element_text(size = 14))+
  ylab(bquote(Low~Density~Growth~Rate~(italic(r[t]))))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(name = NULL,
                     limits = c(1999.5,2015.5))#+
#scale_color_manual(name = "", values = c("#000000", "#888888"))

# Average the growth rates of LACO over time for each pool size.
GRWR_LACO_const_s_mean <- as.data.frame(GRWR_LACO_const_Pool_s) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) #mean of GRWR across years 

GRWR_LACO_const_s_summary <- GRWR_LACO_const_s_mean %>%
  summarize(Mean = mean(mean_GRWR), CI = hdi(mean_GRWR, credMass = 0.95)) %>% #take the 95% CI across iterations
  mutate(name = c("lowCI", "upCI")) %>% #top values is the low CI and bottom value is the high CI
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for small pools = -0.574 (CI: -0.686, -0.453)

GRWR_LACO_const_m_mean <- as.data.frame(GRWR_LACO_const_Pool_m) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) #mean of GRWR across years 

GRWR_LACO_const_m_summary <- GRWR_LACO_const_m_mean %>%
  summarize(Mean = mean(mean_GRWR), CI = hdi(mean_GRWR, credMass = 0.95)) %>% #take the 95% CI across iterations
  mutate(name = c("lowCI", "upCI")) %>% #top values is the low CI and bottom value is the high CI
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for medium pools = -0.403 (CI: -0.497, -0.324)

GRWR_LACO_const_l_mean <- as.data.frame(GRWR_LACO_const_Pool_l) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) #mean of GRWR across years 

GRWR_LACO_const_l_summary <- GRWR_LACO_const_l_mean %>%
  summarize(Mean = mean(mean_GRWR), CI = hdi(mean_GRWR, credMass = 0.95)) %>% #take the 95% CI across iterations
  mutate(name = c("lowCI", "upCI")) %>% #top values is the low CI and bottom value is the high CI
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for large pools = -0.209 (CI: -0.285, -0.129)

# Plot LACO density, lambda, and GRWR timeseries by trt
ggarrange(fig_LACO_size, fig_lambda_size, fig_GRWR_size,  ncol = 1, labels = c("(a)", "(b)", "(c)"))

###Exotic grass cover time series - were some pools more invaded than others?
const_dummy_join$Year <- as.numeric(const_dummy_join$Year)
summary(lmer(sum_EG~Pool + (1|Year), const_dummy_join)) #EG did not differ by Pool 
summary(lmer(sum_EG~Year + (1|Pool), const_dummy_join)) #EG did differ by Year 

const_EG_initial <- const_dummy_join %>%
  filter(Year == 2000) %>%
  select(Pool, Year, sum_EG) %>%
  mutate(EG_2000 = sum_EG) %>%
  select(Pool, EG_2000)

const_EG_time <- left_join(const_dummy_join, const_EG_initial) %>%
  select(Pool, Year, EG_2000, sum_EG)
const_EG_time$Year <- as.numeric(const_EG_time$Year)

ggplot(const_EG_time%>%filter(Year < 2016), aes(x = Year, y = sum_EG, col = EG_2000))+
  geom_jitter()+
  #geom_line(aes(fill = Pool))+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  ylab(bquote(Exotic~Grass~Cover~('%')))+
  scale_color_continuous(breaks = c(0, 20, 40, 60), low = "grey", high = "brown", name = "Initial EG %") #EG cover increased over the years regardless of initial EG cover

const_EG_trt <- const_dummy_join %>%
    select(sum_EG, Year, Treatment, Pool) %>%
    group_by(Treatment, Year) %>%
    summarise(mean_EG = mean(sum_EG), se_EG = se(sum_EG))

ggplot(const_EG_trt%>%filter(Year < 2016), aes(x = Year, y = mean_EG, col = Treatment))+
  geom_point()+
  geom_line()+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_errorbar(aes(ymin = mean_EG-se_EG, ymax = mean_EG+se_EG), width = 0.4, alpha = 0.9, size = 0.8) + 
  ylab(bquote(Exotic~Grass~Cover~('%')))

const_EG_size <- const_dummy_join %>%
    select(sum_EG, Year, Size, Pool) %>%
    group_by(Size, Year) %>%
    summarise(mean_EG = mean(sum_EG), se_EG = se(sum_EG))

ggplot(const_EG_size%>%filter(Year < 2016), aes(x = Year, y = mean_EG, col = Size))+
  geom_point()+
  geom_line()+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_errorbar(aes(ymin = mean_EG-se_EG, ymax = mean_EG+se_EG), width = 0.4, alpha = 0.9, size = 0.8) + 
  ylab(bquote(Exotic~Grass~Cover~('%')))+
  scale_color_discrete(name = "Pool size", labels = c("large", "medium", "small")) 

#Effect of exotic grass on density, lambda, GRWR 
const_EG_mean <- const_dummy_join %>%
  select(Year, Pool, LACOdens, sum_EG) %>%
  group_by(Year) %>%
  summarise(LACO = mean(LACOdens), 
            EG = mean(sum_EG))
Post <- rstan::extract(BH_fit)
lambda_const_mean <- as.data.frame(Post$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(lambda = mean(lambda))
GRWR_const_mean <- as.data.frame(GRWR_LACO_const) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(GRWR = mean(GRWR))

join_means <- left_join(const_EG_mean, lambda_const_mean) %>%
  left_join(., GRWR_const_mean) %>%
  filter(Year < 2016) %>%
  filter(Year > 2000)

join_means$Year <- as.numeric(join_means$Year)
fig_LACO_EG <-ggplot(join_means, aes(x = EG, y = LACO))+
  geom_point()+
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  ylab(bquote(Mean~Annual~LACO~Density~(stems/m^2)))+
  xlab(bquote(Mean~Annual~Exotic~Grass~Cover~('%')))

fig_lambda_EG <-ggplot(join_means, aes(x = EG, y = lambda))+
  geom_point()+
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  ylab(bquote(Mean~Annual~Intrinsic~Growth~Rate~(lambda[t])))+
  xlab(bquote(Mean~Annual~Exotic~Grass~Cover~('%')))

fig_GRWR_EG <-ggplot(join_means, aes(x = EG, y = GRWR))+
  geom_point()+
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  ylab(bquote(Mean~Annual~Low~Density~Growth~Rate~(italic(r[t]))))+
  xlab(bquote(Mean~Annual~Exotic~Grass~Cover~('%')))

ggarrange(fig_LACO_EG, fig_lambda_EG, fig_GRWR_EG, ncol = 1, align = "v")


#Effect of pool depth on each parameter
aLACO_const_mean <- as.data.frame(Post$alpha_LACO) %>% #Calculate mean alpha_LACO each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "aLACO")) %>%
  group_by(Year) %>%
  summarise(aLACO = mean(aLACO))
aEG_const_mean <- as.data.frame(Post$alpha_EG) %>% #Calculate mean alpha_EG each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "aEG")) %>%
  group_by(Year) %>%
  summarise(aEG = mean(aEG))
aNF_const_mean <- as.data.frame(Post$alpha_NF) %>% #Calculate mean alpha_NF each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "aNF")) %>%
  group_by(Year) %>%
  summarise(aNF = mean(aNF))
aERVA_const_mean <- as.data.frame(Post$alpha_ERVA) %>% #Calculate mean alpha_ERVA each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "aERVA")) %>%
  group_by(Year) %>%
  summarise(aERVA = mean(aERVA))
