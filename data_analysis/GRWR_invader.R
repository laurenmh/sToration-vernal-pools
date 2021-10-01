#PARTS:
#I. Long-term growth rate when rare (r_invader) of LACO
#II. r_invader partitioning

# Remember to set your data pathway first!

# Data
source("data_compiling/compile_composition.R") 

# Packages
library(rstan)
library(StanHeaders)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(stargazer)
library(ggpubr)
library(HDInterval)

#--------------------------
#Goal: Calculate the long-term growth rate when rare (r_invader) of LACO 

  #Step 1. Calculate the stable equilibrium frequency of non-LACO species in the model each year. 
  #Step 2. Calculate the annual growth rate of LACO in restored pools when one LACO is introduced into a stable community. 
  #Step 3. Do the same step for the reference pools.
  #Step 4. Average the growth rates of LACO over time for restored and reference pools.
  

#Step 1. Calculate the stable equilibrium frequency of non-LACO species in the model each year. 

#Use just the control plots in reference pools (no LACO present). Average the frequency of non-LACO species across space each year. 
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

#Step 2. Calculate the annual growth rate of LACO in restored pools when one LACO is introduced into a stable community.

#Extract parameters for restored pools. Run "analysis/complex_belowground_v5.R".
Post <- rstan::extract(BH_fit)
alpha_LACO <- as.matrix(Post$alpha_LACO)
alpha_EG <- as.matrix(Post$alpha_EG)
alpha_ERVA <- as.matrix(Post$alpha_ERVA)
alpha_NF <- as.matrix(Post$alpha_NF)
lambda <- as.matrix(Post$lambda)
s <- as.matrix(Post$survival_LACO)
sim_LACO <- matrix(nrow = 2000, ncol = 17)
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
LACO_const <- bh.sim.control(
                      LACO = 1,
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
GRWR_LACO_const_all <- log(LACO_const) #2000-2016 #log transform
GRWR_LACO_const <- GRWR_LACO_const_all[,1:15] #2000-2014 truncate the last two points 

#Step 3. Do the same step for the reference pools.

#Extract parameters for reference pools. Run "analysis/reference_pool_model.R".
Post_ref <- rstan::extract(BH_ref_fit)
refalpha_LACO <- as.matrix(Post_ref$alpha_LACO)
refalpha_EG <- as.matrix(Post_ref$alpha_EG)
refalpha_ERVA <- as.matrix(Post_ref$alpha_ERVA)
refalpha_NF <- as.matrix(Post_ref$alpha_NF)
reflambda <- as.matrix(Post_ref$lambda)
refs <- as.matrix(Post_ref$survival_LACO)
sim_LACO <- matrix(nrow = 2000, ncol = 13)
#Plug in the stable equilibrium freq and parameters in the model.
LACO_ref <- bh.sim.control(
                    LACO = 1,
                    EG = as.numeric(const_com_control[2,4:16]),
                    ERVA = as.numeric(const_com_control[1,4:16]),
                    NF = as.numeric(const_com_control[3,4:16]),
                    aii = refalpha_LACO,
                    a1 = refalpha_EG,
                    a2 = refalpha_ERVA,
                    a3 = refalpha_NF,
                    lambda = reflambda,
                    s = refs,
                    g = 0.7,
                    glow = 0.2)
GRWR_LACO_ref <- log(LACO_ref) #2002-2014 #log transform

#Step 4. Average the growth rates of LACO over time for restored and reference pools.
GRWR_LACO_const_mean <- as.data.frame(GRWR_LACO_const) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) #mean of GRWR across years 

GRWR_LACO_const_summary <- GRWR_LACO_const_mean %>%
  summarize(Mean = mean(mean_GRWR), CI = hdi(mean_GRWR, credMass = 0.95)) %>% #take the 95% CI across iterations
  mutate(name = c("lowCI", "upCI")) %>% #top values is the low CI and bottom value is the high CI
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for constructed pools = -0.4274569 (CI: -0.4631997, -0.3918765)

GRWR_LACO_ref_mean <- as.data.frame(GRWR_LACO_ref) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) 

GRWR_LACO_ref_summary <- GRWR_LACO_ref_mean %>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Average GRWR LACO for reference pools = 0.3650499 (CI: 0.1973319, 0.5323367)

#-----------------------
#Goal: Partition r_LACO by Ellner et al. (2019) method

  #Step 1. Contribution to the overall GRWR when all variation is removed
  #Step 2. Vary lambda while keeping everything else constant
  #Step 3. Vary alphas while keeping everything else constant
  #Step 4. Vary belowground parameters (survival and germination) while keeping everything else constant (removed in the manuscript)
  #Step 5. Interactive effect from simultaneous variation in lambda, alpha, and belowground after accounting for each main effect
  #Step 6. Make sure all the epsilon estimates add up to overall GRWR

#Step 1. Contribution to the overall GRWR when all variation is removed

#calculate mean of each parameter over the whole time series
alpha_LACO_knot <- matrix(rep(rowMeans(alpha_LACO[,1:15]), 15), ncol = 15, nrow = 2000)
alpha_EG_knot <- matrix(rep(rowMeans(alpha_EG[,1:15]), 15), ncol = 15, nrow = 2000)
alpha_ERVA_knot <- matrix(rep(rowMeans(alpha_ERVA[,1:15]), 15), ncol = 15, nrow = 2000)
alpha_NF_knot <- matrix(rep(rowMeans(alpha_NF[,1:15]), 15), ncol = 15, nrow = 2000)
lambda_mean_knot <- matrix(rep(rowMeans(lambda[,1:15]), 15), ncol = 15, nrow = 2000)
g_knot <- (0.7+0.2)/2

refalpha_LACO_knot <- matrix(rep(rowMeans(refalpha_LACO), 13), ncol = 13, nrow = 2000)
refalpha_EG_knot <- matrix(rep(rowMeans(refalpha_EG), 13), ncol = 13, nrow = 2000)
refalpha_ERVA_knot <- matrix(rep(rowMeans(refalpha_ERVA), 13), ncol = 13, nrow = 2000)
refalpha_NF_knot <- matrix(rep(rowMeans(refalpha_NF), 13), ncol = 13, nrow = 2000)
reflambda_knot <- matrix(rep(rowMeans(reflambda), 13), ncol = 13, nrow = 2000)
refg_knot <- (0.7+0.2)/2

#calculate mean of each plant group over the whole time series
EG_knot <- as.numeric(rep(rowMeans(const_com_control[2,2:16]), 15))
ERVA_knot <- as.numeric(rep(rowMeans(const_com_control[1,2:16]), 15))
NF_knot <- as.numeric(rep(rowMeans(const_com_control[3,2:16]), 15))

#modify the model with constant g
bh.sim.control.knot <- function(LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
  for(i in 1:nrow(sim_LACO)){
    for(j in 1:ncol(sim_LACO)){
      sim_LACO[i,j] <- LACO*lambda[i,j]/(1+LACO*aii[i,j]+EG[j]*a1[i,j]+ERVA[j]*a2[i,j]+NF[j]*a3[i,j])+s[i]*(1-g)*LACO/g 
    }
  }
  return(sim_LACO)
}

#run the model with constant parameters and constant equilibrium community
sim_LACO <- matrix(nrow = 2000, ncol = 15)
LACO_const_knot <- bh.sim.control.knot(
                        LACO = 1,
                        EG = EG_knot,
                        ERVA = ERVA_knot,
                        NF = NF_knot,
                        aii = alpha_LACO_knot,
                        a1 = alpha_EG_knot,
                        a2 = alpha_ERVA_knot,
                        a3 = alpha_NF_knot,
                        lambda = lambda_mean_knot,
                        s = s,
                        g = g_knot)
GRWR_LACO_const_knot <- log(LACO_const_knot) #log transform

sim_LACO <- matrix(nrow = 2000, ncol = 13)
LACO_ref_knot <- bh.sim.control.knot(
                        LACO = 1,
                        EG = EG_knot[1:13],
                        ERVA = ERVA_knot[1:13],
                        NF = NF_knot[1:13],
                        aii = refalpha_LACO_knot,
                        a1 = refalpha_EG_knot,
                        a2 = refalpha_ERVA_knot,
                        a3 = refalpha_NF_knot,
                        lambda = reflambda_knot,
                        s = refs,
                        g = refg_knot)
GRWR_LACO_ref_knot <- log(LACO_ref_knot) #log transform

#calculate mean and 95% CI:

LACO_const_epsilon_0_mean <- as.data.frame(GRWR_LACO_const_knot) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR))

LACO_const_epsilon_0_summary <- LACO_const_epsilon_0_mean  %>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Epsilon knot for constructed pools is -0.3137382 (CI: -0.4318092, -0.1885801)

LACO_ref_epsilon_0_mean <- as.data.frame(GRWR_LACO_ref_knot) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR))

LACO_ref_epsilon_0_summary <- LACO_ref_epsilon_0_mean  %>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Epsilon knot for reference pools is 0.4966552 (CI:0.3233172, 0.6778781)

#Step 2. Vary lambda while keeping everything else constant
sim_LACO <- matrix(nrow = 2000, ncol = 15)
LACO_const_lambda <- bh.sim.control.knot(
                        LACO = 1,
                        EG = EG_knot,
                        ERVA = ERVA_knot,
                        NF = NF_knot,
                        aii = alpha_LACO_knot,
                        a1 = alpha_EG_knot,
                        a2 = alpha_ERVA_knot,
                        a3 = alpha_NF_knot,
                        lambda = lambda[,1:15],
                        s = s,
                        g = g_knot)

GRWR_LACO_const_lambda <- log(LACO_const_lambda) #log transform
LACO_const_epsilon_lambda <- GRWR_LACO_const_lambda - GRWR_LACO_const_knot #subtract GRWR knot (from above) from GRWR with varying lambda

LACO_const_epsilon_lambda_mean <- as.data.frame(LACO_const_epsilon_lambda) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR))
  
LACO_const_epsilon_lambda_summary <- LACO_const_epsilon_lambda_mean %>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Epsilon lambda for constructed pools is -0.9676916 (CI: -1.067676, -0.8599437)

sim_LACO <- matrix(nrow = 2000, ncol = 13)
LACO_ref_lambda <- bh.sim.control.knot(
                        LACO = 1,
                        EG = EG_knot[1:13],
                        ERVA = ERVA_knot[1:13],
                        NF = NF_knot[1:13],
                        aii = refalpha_LACO_knot,
                        a1 = refalpha_EG_knot,
                        a2 = refalpha_ERVA_knot,
                        a3 = refalpha_NF_knot,
                        lambda = reflambda,
                        s = refs,
                        g = refg_knot)

GRWR_LACO_ref_lambda <- log(LACO_ref_lambda) 
LACO_ref_epsilon_lambda <- GRWR_LACO_ref_lambda - GRWR_LACO_ref_knot 

LACO_ref_epsilon_lambda_mean <- as.data.frame(LACO_ref_epsilon_lambda) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR))

LACO_ref_epsilon_lambda_summary <- LACO_ref_epsilon_lambda_mean %>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Epsilon lambda for reference pools is -0.3247265 (CI:-0.377565, -0.2748186)

#Step 3. Vary alphas while keeping everything else constant
sim_LACO <- matrix(nrow = 2000, ncol = 15)
LACO_const_alpha <- bh.sim.control.knot(
                        LACO = 1,
                        EG = as.numeric(const_com_control[2,2:16]),
                        ERVA = as.numeric(const_com_control[1,2:16]),
                        NF = as.numeric(const_com_control[3,2:16]),
                        aii = alpha_LACO[,1:15],
                        a1 = alpha_EG[,1:15],
                        a2 = alpha_ERVA[,1:15],
                        a3 = alpha_NF[,1:15],
                        lambda = lambda_mean_knot,
                        s = s,
                        g = g_knot)

GRWR_LACO_const_alpha <- log(LACO_const_alpha) #log transform
LACO_const_epsilon_alpha <- GRWR_LACO_const_alpha - GRWR_LACO_const_knot #subtract GRWR knot (from above) from GRWR with varying alpha

LACO_const_epsilon_alpha_mean <- as.data.frame(LACO_const_epsilon_alpha) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR))

LACO_const_epsilon_alpha_summary <- LACO_const_epsilon_alpha_mean %>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Epsilon alpha for constructed pools is 0.8543256 (CI: 0.8070232, 0.8935593)
  
sim_LACO <- matrix(nrow = 2000, ncol = 13)
LACO_ref_alpha <- bh.sim.control.knot(
                        LACO = 1,
                        EG = as.numeric(const_com_control[2,4:16]),
                        ERVA = as.numeric(const_com_control[1,4:16]),
                        NF = as.numeric(const_com_control[3,4:16]),
                        aii = refalpha_LACO,
                        a1 = refalpha_EG,
                        a2 = refalpha_ERVA,
                        a3 = refalpha_NF,
                        lambda = reflambda_knot,
                        s = refs,
                        g = refg_knot)

GRWR_LACO_ref_alpha <- log(LACO_ref_alpha)
LACO_ref_epsilon_alpha <- GRWR_LACO_ref_alpha - GRWR_LACO_ref_knot

LACO_ref_epsilon_alpha_mean <- as.data.frame(LACO_ref_epsilon_alpha) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) 

LACO_ref_epsilon_alpha_summary <- LACO_ref_epsilon_alpha_mean %>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Epsilon alpha for reference pools is 0.197456 (CI: 0.08089213, 0.3079596)

#Step 4. Vary belowground parameters (survival and germination) while keeping everything else constant
# sim_LACO <- matrix(nrow = 17, ncol = 1)
# LACO_const_bg <- bh.sim.control(
#                         LACO = 1,
#                         EG = const_com_control$avg_sumEG,
#                         ERVA = const_com_control$avg_ERVA,
#                         NF = const_com_control$avg_sumNF,
#                         aii = alpha_LACO_knot,
#                         a1 = alpha_EG_knot,
#                         a2 = alpha_ERVA_knot,
#                         a3 = alpha_NF_knot,
#                         lambda = lambda_mean_knot,
#                         s = s_mean[,5],
#                         g = 0.7,
#                         glow = 0.2)
# GRWR_LACO_const_bg <- log(LACO_const_bg)
# LACO_const_epsilon_bg <- mean(GRWR_LACO_const_bg) - LACO_const_epsilon_0 # 0.1116214
# sim_LACO <- matrix(nrow = 13, ncol = 1)
# LACO_ref_bg <- bh.sim.control(
#                         LACO = 1,
#                         EG = const_com_control$avg_sumEG,
#                         ERVA = const_com_control$avg_ERVA,
#                         NF = const_com_control$avg_sumNF,
#                         aii = refalpha_LACO_knot,
#                         a1 = refalpha_EG_knot,
#                         a2 = refalpha_ERVA_knot,
#                         a3 = refalpha_NF_knot,
#                         lambda = reflambda_knot,
#                         s = refs_mean[,5],
#                         g = 0.7,
#                         glow = 0.2)
# GRWR_LACO_ref_bg <- log(LACO_ref_bg)
# LACO_ref_epsilon_bg <- mean(GRWR_LACO_ref_bg) - LACO_ref_epsilon_0 # -0.002154302

#Step 5. Interactive effect from simultaneous variation in lambda, alpha, and belowground after accounting for each main effect
GRWR_LACO_const_interaction <- GRWR_LACO_const - (GRWR_LACO_const_knot + LACO_const_epsilon_alpha + LACO_const_epsilon_lambda) 

LACO_const_epsilon_interaction_mean <- as.data.frame(GRWR_LACO_const_interaction) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR))

LACO_const_epsilon_interaction_summary <- LACO_const_epsilon_interaction_mean%>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Epsilon interaction for constructed pools is -0.0003528306 (CI: -0.001945592, 0.0001624778)

GRWR_LACO_ref_interaction <- GRWR_LACO_ref - (GRWR_LACO_ref_knot + LACO_ref_epsilon_alpha + LACO_ref_epsilon_lambda)

LACO_ref_epsilon_interaction_mean <- as.data.frame(GRWR_LACO_ref_interaction) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  rownames_to_column(., var = "iteration") %>%
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR)) 

LACO_ref_epsilon_interaction_summary <- LACO_ref_epsilon_interaction_mean %>%
  summarise(Mean = mean(mean_GRWR),
            CI = hdi(mean_GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Epsilon interaction for reference pools is -0.004334681 (CI: -0.01677704, -1.046061e-09)

#Step 6: Make sure all the epsilon estimates add up to overall GRWR
# Constructed overall GRWR: -0.4274569
-0.0003528306 + 0.8543256 -0.9676916 -0.3137382 
# Reference overall GRWR: 0.3650499
-0.004334681 + 0.197456 -0.3247265 +0.4966552

#FIGURE 3
Partitioning_GRWR_LACO_const <- as.data.frame(rbind(GRWR_LACO_const_summary, LACO_const_epsilon_0_summary, LACO_const_epsilon_alpha_summary, LACO_const_epsilon_lambda_summary, LACO_const_epsilon_interaction_summary)) %>%
  mutate(mechanism = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Partitioning_GRWR_LACO_const$mechanism <- ordered(Partitioning_GRWR_LACO_const$mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
xlabels <- c("r_overall" = expression(bar("r")[i]^" "), 
             "epsilon_0" = expression(epsilon[i]^0), 
             "epsilon_alpha" = expression(epsilon[i]^alpha),
             "epsilon_lambda" = expression(epsilon[i]^lambda),
             "epsilon_int" = expression(epsilon[i]^{alpha*lambda}))
Part_const <- ggplot(Partitioning_GRWR_LACO_const, aes(x = mechanism, y = Mean, fill = mechanism))+
            geom_bar(stat = "identity")+
            geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 0.8) +
            theme_classic(base_size = 20)+
            theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
            geom_hline(yintercept = 0)+
            scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
            ylim(-1.2, 1.2)+
            annotate("text", x = 2, y = 1.1, label = "Constructed", size = 6)+
            scale_x_discrete(labels = xlabels)

Partitioning_GRWR_LACO_ref <- as.data.frame(rbind(GRWR_LACO_ref_summary, LACO_ref_epsilon_0_summary, LACO_ref_epsilon_alpha_summary, LACO_ref_epsilon_lambda_summary, LACO_ref_epsilon_interaction_summary)) %>%
  mutate(mechanism = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Partitioning_GRWR_LACO_ref$mechanism <- ordered(Partitioning_GRWR_LACO_ref$mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Part_ref <- ggplot(Partitioning_GRWR_LACO_ref, aes(x = mechanism, y = Mean, fill = mechanism))+
            geom_bar(stat = "identity")+
            geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 0.8) +
            theme_classic(base_size = 20)+
            theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
            geom_hline(yintercept = 0)+
            scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
            ylim(-1.2, 1.2)+
            annotate("text", x = 2, y = 1.1, label = "Reference", size = 6)+
            scale_x_discrete(labels = xlabels)
figure_partitioning <- ggarrange(Part_ref, Part_const, ncol = 2, nrow = 1, legend = "none", 
                                  labels = c("(a)", "(b)"), font.label = list(size = 20))
annotate_figure(figure_partitioning, bottom = text_grob("Mechanisms", size = 20),
                left = text_grob("Partitioning of Low Density Growth Rate", size = 18, rot = 90))

