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
GRWR_LACO_const <- log(LACO_const) #2000-2016

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
GRWR_LACO_ref <- log(LACO_ref) #2002-2014

#Step 4. Average the growth rates of LACO over time for restored and reference pools.
GRWR_LACO_const_summary <- as.data.frame(GRWR_LACO_const) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  summarise(mean = mean(GRWR),
            CI = hdi(GRWR, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI)#Average GRWR LACO for constructed pools = -0.2204419

GRWR_LACO_ref_summary <- as.data.frame(GRWR_LACO_ref) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  summarise(mean = mean(GRWR),
            CI = hdi(GRWR, credMass = 0.95))%>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) #Average GRWR LACO for reference pools = 0.2770058

#-----------------------
#Goal: Partition r_LACO by Ellner et al. (2019) method

  #Step 1. Contribution to the overall GRWR when all variation is removed
  #Step 2. Vary lambda while keeping everything else constant
  #Step 3. Vary alphas while keeping everything else constant
  #Step 4. Vary belowground parameters (survival and germination) while keeping everything else constant (removed in the manuscript)
  #Step 5. Interactive effect from simultaneous variation in lambda, alpha, and belowground after accounting for each main effect

#Step 1. Contribution to the overall GRWR when all variation is removed

#calculate mean of each parameter over the whole time series
alpha_LACO_knot <- matrix(rep(rowMeans(alpha_LACO), 17), ncol = 17, nrow = 2000)
alpha_EG_knot <- matrix(rep(rowMeans(alpha_EG), 17), ncol = 17, nrow = 2000)
alpha_ERVA_knot <- matrix(rep(rowMeans(alpha_ERVA), 17), ncol = 17, nrow = 2000)
alpha_NF_knot <- matrix(rep(rowMeans(alpha_NF), 17), ncol = 17, nrow = 2000)
lambda_mean_knot <- matrix(rep(rowMeans(lambda), 17), ncol = 17, nrow = 2000)
g_knot <- (0.7+0.2)/2

refalpha_LACO_knot <- matrix(rep(rowMeans(refalpha_LACO), 13), ncol = 13, nrow = 2000)
refalpha_EG_knot <- matrix(rep(rowMeans(refalpha_EG), 13), ncol = 13, nrow = 2000)
refalpha_ERVA_knot <- matrix(rep(rowMeans(refalpha_ERVA), 13), ncol = 13, nrow = 2000)
refalpha_NF_knot <- matrix(rep(rowMeans(refalpha_NF), 13), ncol = 13, nrow = 2000)
reflambda_knot <- matrix(rep(rowMeans(reflambda), 13), ncol = 13, nrow = 2000)
refg_knot <- (0.7+0.2)/2

#calculate mean of each plant group over the whole time series
EG_knot <- as.numeric(rep(rowMeans(const_com_control[2,2:18]), 17))
ERVA_knot <- as.numeric(rep(rowMeans(const_com_control[1,2:18]), 17))
NF_knot <- as.numeric(rep(rowMeans(const_com_control[3,2:18]), 17))

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
sim_LACO <- matrix(nrow = 2000, ncol = 17)
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
GRWR_LACO_const_knot <- log(LACO_const_knot) 
LACO_const_epsilon_0_summary <- as.data.frame(GRWR_LACO_const_knot) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "epsilon")) %>%
  summarise(mean = mean(epsilon),
            CI = hdi(epsilon, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI)

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
GRWR_LACO_ref_knot <- log(LACO_ref_knot) 
LACO_ref_epsilon_0_summary <- as.data.frame(GRWR_LACO_ref_knot) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "epsilon")) %>%
  summarise(mean = mean(epsilon),
            CI = hdi(epsilon, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI)

#Step 2. Vary lambda while keeping everything else constant
sim_LACO <- matrix(nrow = 2000, ncol = 17)
LACO_const_lambda <- bh.sim.control.knot(
                        LACO = 1,
                        EG = EG_knot,
                        ERVA = ERVA_knot,
                        NF = NF_knot,
                        aii = alpha_LACO_knot,
                        a1 = alpha_EG_knot,
                        a2 = alpha_ERVA_knot,
                        a3 = alpha_NF_knot,
                        lambda = lambda,
                        s = s,
                        g = g_knot)
GRWR_LACO_const_lambda <- log(LACO_const_lambda) 
LACO_const_epsilon_lambda <- GRWR_LACO_const_lambda - GRWR_LACO_const_knot 
LACO_const_epsilon_lambda_summary <- as.data.frame(LACO_const_epsilon_lambda) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "epsilon")) %>%
  summarise(mean = mean(epsilon),
            CI = hdi(epsilon, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI) 

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
LACO_ref_epsilon_lambda_summary <- as.data.frame(LACO_ref_epsilon_lambda) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "epsilon")) %>%
  summarise(mean = mean(epsilon),
            CI = hdi(epsilon, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI)

#Step 3. Vary alphas while keeping everything else constant
sim_LACO <- matrix(nrow = 2000, ncol = 17)
LACO_const_alpha <- bh.sim.control.knot(
                        LACO = 1,
                        EG = as.numeric(const_com_control[2,2:18]),
                        ERVA = as.numeric(const_com_control[1,2:18]),
                        NF = as.numeric(const_com_control[3,2:18]),
                        aii = alpha_LACO,
                        a1 = alpha_EG,
                        a2 = alpha_ERVA,
                        a3 = alpha_NF,
                        lambda = lambda_mean_knot,
                        s = s,
                        g = g_knot)
GRWR_LACO_const_alpha <- log(LACO_const_alpha) 
LACO_const_epsilon_alpha <- GRWR_LACO_const_alpha - GRWR_LACO_const_knot 
LACO_const_epsilon_alpha_summary <- as.data.frame(LACO_const_epsilon_alpha) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "epsilon")) %>%
  summarise(mean = mean(epsilon),
            CI = hdi(epsilon, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI)
  
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
LACO_ref_epsilon_alpha_summary <- as.data.frame(LACO_ref_epsilon_alpha) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "epsilon")) %>%
  summarise(mean = mean(epsilon),
            CI = hdi(epsilon, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI)  

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
LACO_const_interaction_summary <- as.data.frame(GRWR_LACO_const_interaction) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "epsilon")) %>%
  summarise(mean = mean(epsilon),
            CI = hdi(epsilon, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI)  

GRWR_LACO_ref_interaction <- GRWR_LACO_ref - (GRWR_LACO_ref_knot + LACO_ref_epsilon_alpha + LACO_ref_epsilon_lambda)
LACO_ref_interaction_summary <- as.data.frame(GRWR_LACO_ref_interaction) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "epsilon")) %>%
  summarise(mean = mean(epsilon),
            CI = hdi(epsilon, credMass = 0.95)) %>%
  mutate(name = c("lowCI", "upCI")) %>%
  pivot_wider(names_from = name, values_from = CI)

#FIGURE 3
Partitioning_GRWR_LACO_const <- as.data.frame(rbind(GRWR_LACO_const_summary, LACO_const_epsilon_0_summary, LACO_const_epsilon_alpha_summary, LACO_const_epsilon_lambda_summary, LACO_const_interaction_summary)) %>%
  mutate(mechanism = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Partitioning_GRWR_LACO_const$mechanism <- ordered(Partitioning_GRWR_LACO_const$mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
xlabels <- c("r_overall" = expression(bar("r")[i]^" "), 
             "epsilon_0" = expression(epsilon[i]^0), 
             "epsilon_alpha" = expression(epsilon[i]^alpha),
             "epsilon_lambda" = expression(epsilon[i]^lambda),
             "epsilon_int" = expression(epsilon[i]^{alpha*lambda}))
Part_const <- ggplot(Partitioning_GRWR_LACO_const, aes(x = mechanism, y = mean, fill = mechanism))+
            geom_bar(stat = "identity")+
            geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 0.8) +
            theme_classic(base_size = 20)+
            theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
            geom_hline(yintercept = 0)+
            scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
            #ylim(-1.2, 1.2)+
            annotate("text", x = 2, y = 1.0, label = "Constructed", size = 6)+
            scale_x_discrete(labels = xlabels)

Partitioning_GRWR_LACO_ref <- as.data.frame(rbind(GRWR_LACO_ref_summary, LACO_ref_epsilon_0_summary, LACO_ref_epsilon_alpha_summary, LACO_ref_epsilon_lambda_summary, LACO_ref_interaction_summary)) %>%
  mutate(mechanism = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Partitioning_GRWR_LACO_ref$mechanism <- ordered(Partitioning_GRWR_LACO_ref$mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Part_ref <- ggplot(Partitioning_GRWR_LACO_ref, aes(x = mechanism, y = mean, fill = mechanism))+
            geom_bar(stat = "identity")+
            geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 0.8) +
            theme_classic(base_size = 20)+
            theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
            geom_hline(yintercept = 0)+
            scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
            #ylim(-1.2, 1.2)+
            annotate("text", x = 2, y = 1.0, label = "Reference", size = 6)+
            scale_x_discrete(labels = xlabels)
figure_partitioning <- ggarrange(Part_ref, Part_const, ncol = 2, nrow = 1, legend = "none", 
                                  labels = c("(a)", "(b)"), font.label = list(size = 20))
annotate_figure(figure_partitioning, bottom = text_grob("Mechanisms", size = 20),
                left = text_grob("Partitioning of Low Density Growth Rate", size = 18, rot = 90))

