#PARTS:
#I. Long-term growth rate when rare (r_invader) of LACO
#II. r_invader partitioning
#III. EG removal simulation

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
                    EG = as.numeric(const_com_control[2,2:18]),
                    ERVA = as.numeric(const_com_control[1,2:18]),
                    NF = as.numeric(const_com_control[3,2:18]),
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
mean(GRWR_LACO_const) #Average GRWR LACO for constructed pools = -0.2203983
mean(GRWR_LACO_ref) #Average GRWR LACO for reference pools = 0.3754354 

#-----------------------
#Goal: Partition r_LACO by Ellner et al. (2019) method

  #Step 1. Contribution to the overall GRWR when all variation is removed
  #Step 2. Vary lambda while keeping everything else constant
  #Step 3. Vary alphas while keeping everything else constant
  #Step 4. Vary belowground parameters (survival and germination) while keeping everything else constant (removed in the manuscript)
  #Step 5. Interactive effect from simultaneous variation in lambda, alpha, and belowground after accounting for each main effect

#Step 1. Contribution to the overall GRWR when all variation is removed

#calculate mean of each parameter over the whole time series
alpha_LACO_knot <- rep(mean(alpha_LACO), 17)
alpha_EG_knot <- rep(mean(alpha_EG), 17)
alpha_ERVA_knot <- rep(mean(alpha_ERVA), 17)
alpha_NF_knot <- rep(mean(alpha_NF), 17)
lambda_mean_knot <- rep(mean(lambda), 17)
s_mean_knot <- mean(s[,5])
g_knot <- (0.7+0.2)/2

refalpha_LACO_knot <- rep(mean(refalpha_LACO), 13)
refalpha_EG_knot <- rep(mean(refalpha_EG), 13)
refalpha_ERVA_knot <- rep(mean(refalpha_ERVA), 13)
refalpha_NF_knot <- rep(mean(refalpha_NF), 13)
reflambda_knot <- rep(mean(reflambda), 13)
refs_knot <- mean(refs)
refg_knot <- (0.7+0.2)/2

#modify the model with constant g
bh.sim.control.knot <- function(LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
  for(i in 1:nrow(sim_LACO)){
    sim_LACO[i] <- LACO*lambda[i]/(1+LACO*aii[i]+EG[i]*a1[i]+ERVA[i]*a2[i]+NF[i]*a3[i])+s*(1-g)*LACO/g #this is the modified Beverton-Holt model we'll use for LACO stem counts
  }
  return(sim_LACO)
}

#run the model with constant parameters
sim_LACO <- matrix(nrow = 17, ncol = 1)
LACO_const_knot <- bh.sim.control.knot(
                        LACO = 1,
                        EG = const_com_control$avg_sumEG,
                        ERVA = const_com_control$avg_ERVA,
                        NF = const_com_control$avg_sumNF,
                        aii = alpha_LACO_knot,
                        a1 = alpha_EG_knot,
                        a2 = alpha_ERVA_knot,
                        a3 = alpha_NF_knot,
                        lambda = lambda_mean_knot,
                        s = s_mean_knot,
                        g = g_knot)
GRWR_LACO_const_knot <- log(LACO_const_knot) 
LACO_const_epsilon_0 <- mean(GRWR_LACO_const_knot) # -0.04333787
sim_LACO <- matrix(nrow = 13, ncol = 1)
LACO_ref_knot <- bh.sim.control.knot(
                        LACO = 1,
                        EG = const_com_control$avg_sumEG,
                        ERVA = const_com_control$avg_ERVA,
                        NF = const_com_control$avg_sumNF,
                        aii = refalpha_LACO_knot,
                        a1 = refalpha_EG_knot,
                        a2 = refalpha_ERVA_knot,
                        a3 = refalpha_NF_knot,
                        lambda = reflambda_knot,
                        s = refs_knot,
                        g = refg_knot)
GRWR_LACO_ref_knot <- log(LACO_ref_knot) 
LACO_ref_epsilon_0 <- mean(GRWR_LACO_ref_knot) # 0.2524456

#Step 2. Vary lambda while keeping everything else constant
sim_LACO <- matrix(nrow = 17, ncol = 1)
LACO_const_lambda <- bh.sim.control.knot(
                        LACO = 1,
                        EG = const_com_control$avg_sumEG,
                        ERVA = const_com_control$avg_ERVA,
                        NF = const_com_control$avg_sumNF,
                        aii = alpha_LACO_knot,
                        a1 = alpha_EG_knot,
                        a2 = alpha_ERVA_knot,
                        a3 = alpha_NF_knot,
                        lambda = lambda,
                        s = s_mean_knot,
                        g = g_knot)
GRWR_LACO_const_lambda <- log(LACO_const_lambda) 
LACO_const_epsilon_lambda <- mean(GRWR_LACO_const_lambda) - LACO_const_epsilon_0 # -1.194625
sim_LACO <- matrix(nrow = 13, ncol = 1)
LACO_ref_lambda <- bh.sim.control.knot(
                        LACO = 1,
                        EG = const_com_control$avg_sumEG,
                        ERVA = const_com_control$avg_ERVA,
                        NF = const_com_control$avg_sumNF,
                        aii = refalpha_LACO_knot,
                        a1 = refalpha_EG_knot,
                        a2 = refalpha_ERVA_knot,
                        a3 = refalpha_NF_knot,
                        lambda = reflambda,
                        s = refs_knot,
                        g = refg_knot)
GRWR_LACO_ref_lambda <- log(LACO_ref_lambda) 
LACO_ref_epsilon_lambda <- mean(GRWR_LACO_ref_lambda) - LACO_ref_epsilon_0 #-0.3529654

#Step 3. Vary alphas while keeping everything else constant
sim_LACO <- matrix(nrow = 17, ncol = 1)
LACO_const_alpha <- bh.sim.control.knot(
                        LACO = 1,
                        EG = const_com_control$avg_sumEG,
                        ERVA = const_com_control$avg_ERVA,
                        NF = const_com_control$avg_sumNF,
                        aii = alpha_LACO_mean[,5],
                        a1 = alpha_EG_mean[,5],
                        a2 = alpha_ERVA_mean[,5],
                        a3 = alpha_NF_mean[,5],
                        lambda = lambda_mean_knot,
                        s = s_mean_knot,
                        g = g_knot)
GRWR_LACO_const_alpha <- log(LACO_const_alpha) 
LACO_const_epsilon_alpha <- mean(GRWR_LACO_const_alpha) - LACO_const_epsilon_0 #0.6620671
sim_LACO <- matrix(nrow = 13, ncol = 1)
LACO_ref_alpha <- bh.sim.control.knot(
                        LACO = 1,
                        EG = const_com_control$avg_sumEG,
                        ERVA = const_com_control$avg_ERVA,
                        NF = const_com_control$avg_sumNF,
                        aii = refalpha_LACO_mean[,5],
                        a1 = refalpha_EG_mean[,5],
                        a2 = refalpha_ERVA_mean[,5],
                        a3 = refalpha_NF_mean[,5],
                        lambda = reflambda_knot,
                        s = refs_knot,
                        g = refg_knot)
GRWR_LACO_ref_alpha <- log(LACO_ref_alpha)
LACO_ref_epsilon_alpha <- mean(GRWR_LACO_ref_alpha) - LACO_ref_epsilon_0 #0.128182

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
GRWR_LACO_const_interaction <- mean(GRWR_LACO_const) - (LACO_const_epsilon_0 + LACO_const_epsilon_alpha + LACO_const_epsilon_lambda)
GRWR_LACO_ref_interaction <- mean(GRWR_LACO_ref) - (LACO_ref_epsilon_0 + LACO_ref_epsilon_alpha + LACO_ref_epsilon_lambda)

#plot partitioned GRWR
Partitioning_GRWR_LACO_const <- as.data.frame(c(mean(GRWR_LACO_const), LACO_const_epsilon_0, LACO_const_epsilon_alpha, LACO_const_epsilon_lambda, GRWR_LACO_const_interaction)) %>%
  mutate(mechanism = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
colnames(Partitioning_GRWR_LACO_const) <- c("Partitioned_GRWR_LACO_const","Mechanism") 
Partitioning_GRWR_LACO_const$Mechanism <- ordered(Partitioning_GRWR_LACO_const$Mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
xlabels <- c("r_overall" = expression(bar("r")[i]^" "), 
             "epsilon_0" = expression(epsilon[i]^0), 
             "epsilon_alpha" = expression(epsilon[i]^alpha),
             "epsilon_lambda" = expression(epsilon[i]^lambda),
             "epsilon_int" = expression(epsilon[i]^{alpha*lambda}))
Part_const <- ggplot(Partitioning_GRWR_LACO_const, aes(x = Mechanism, y = Partitioned_GRWR_LACO_const, fill = Mechanism))+
            geom_bar(stat = "identity")+
            theme_classic(base_size = 14)+
            theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
            geom_hline(yintercept = 0)+
            scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
            ylim(-1.4, 0.8)+
            scale_x_discrete(labels = xlabels)

Partitioning_GRWR_LACO_ref <- as.data.frame(c(mean(GRWR_LACO_ref), LACO_ref_epsilon_0, LACO_ref_epsilon_alpha, LACO_ref_epsilon_lambda, GRWR_LACO_ref_interaction)) %>%
  mutate(mechanism = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
colnames(Partitioning_GRWR_LACO_ref) <- c("Partitioned_GRWR_LACO_ref","Mechanism") 
Partitioning_GRWR_LACO_ref$Mechanism <- ordered(Partitioning_GRWR_LACO_ref$Mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Part_ref <- ggplot(Partitioning_GRWR_LACO_ref, aes(x = Mechanism, y = Partitioned_GRWR_LACO_ref, fill = Mechanism))+
            geom_bar(stat = "identity")+
            theme_classic(base_size = 14)+
            theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
            geom_hline(yintercept = 0)+
            scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
            ylim(-1.4, 0.8)+
            scale_x_discrete(labels = xlabels)
figure_partitioning <- ggarrange(Part_ref, Part_const, ncol = 2, nrow = 1, legend = "none", 
                                  labels = c("(a) Reference", "(b) Constructed"), font.label = list(size = 11))
annotate_figure(figure_partitioning, bottom = "Mechanisms", left = "Partitioning of average population growth rate")

#-----------------------
#Goal: Simulate exotic grasses (EG) removal to promote LACO persistence

  #Step 1. Simulate EG removal
  #Step 2. Average the growth rates of LACO over time for all simulation scenarios
  #Step 3. Plot simulated GRWR

#Step 1. Remove EG in all years

#Remove 25% of EG in all years
mult.25 <- function(x)(x*0.75)
reduced25EGcover_all <- sumEG %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.25)
LACO_25EG_all <- bh.sim.control(
                      LACO = 1,
                      EG = as.matrix(reduced25EGcover_all),
                      ERVA = const_com_control$avg_ERVA,
                      NF = const_com_control$avg_sumNF,
                      aii = alpha_LACO_mean[,5],
                      a1 = alpha_EG_mean[,5],
                      a2 = alpha_ERVA_mean[,5],
                      a3 = alpha_NF_mean[,5],
                      lambda = lambda_mean[,5],
                      s = s_mean[,5],
                      g = 0.7,
                      glow = 0.2)
GRWR_LACO_25EG_all <- log(LACO_25EG_all)

#Remove 50% of EG in all years
mult.5 <- function(x)(x*0.5)
reduced50EGcover_all <- sumEG %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.5)
LACO_50EG_all <- bh.sim.control(
                      LACO = 1,
                      EG = as.matrix(reduced50EGcover_all),
                      ERVA = const_com_control$avg_ERVA,
                      NF = const_com_control$avg_sumNF,
                      aii = alpha_LACO_mean[,5],
                      a1 = alpha_EG_mean[,5],
                      a2 = alpha_ERVA_mean[,5],
                      a3 = alpha_NF_mean[,5],
                      lambda = lambda_mean[,5],
                      s = s_mean[,5],
                      g = 0.7,
                      glow = 0.2)
GRWR_LACO_50EG_all <- log(LACO_50EG_all)

#Remove 75% of EG in all years
mult.75 <- function(x)(x*0.25)
reduced75EGcover_all <- sumEG %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.75)
LACO_75EG_all <- bh.sim.control(
                      LACO = 1,
                      EG = as.matrix(reduced75EGcover_all),
                      ERVA = const_com_control$avg_ERVA,
                      NF = const_com_control$avg_sumNF,
                      aii = alpha_LACO_mean[,5],
                      a1 = alpha_EG_mean[,5],
                      a2 = alpha_ERVA_mean[,5],
                      a3 = alpha_NF_mean[,5],
                      lambda = lambda_mean[,5],
                      s = s_mean[,5],
                      g = 0.7,
                      glow = 0.2)
GRWR_LACO_75EG_all <- log(LACO_75EG_all)

#Remove 100% of EG in all years
mult.100 <- function(x)(x*0)
reduced100EGcover_all <- sumEG %>%
  mutate_at(c("2000", "2001", "2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010", "2011", 
              "2012", "2013", "2014", "2015", "2016"), mult.100)
LACO_100EG_all <- bh.sim.control(
                      LACO = 1,
                      EG = as.matrix(reduced100EGcover_all),
                      ERVA = const_com_control$avg_ERVA,
                      NF = const_com_control$avg_sumNF,
                      aii = alpha_LACO_mean[,5],
                      a1 = alpha_EG_mean[,5],
                      a2 = alpha_ERVA_mean[,5],
                      a3 = alpha_NF_mean[,5],
                      lambda = lambda_mean[,5],
                      s = s_mean[,5],
                      g = 0.7,
                      glow = 0.2)
GRWR_LACO_100EG_all <- log(LACO_100EG_all)

#Average the growth rates of LACO over time for all simulation scenarios.
mean(GRWR_LACO_25EG_all) # -0.3718841
mean(GRWR_LACO_50EG_all) # -0.1032695
mean(GRWR_LACO_75EG_all) # 0.2953273
mean(GRWR_LACO_100EG_all) # 1.313785

#Plot simulated GRWR
GRWR_simulated_all <- cbind(GRWR_LACO, GRWR_LACO_25EG_all, GRWR_LACO_50EG_all, GRWR_LACO_75EG_all, GRWR_LACO_100EG_all) 
colnames(GRWR_simulated_all) <- c("0% removed", "Year", "reference", "25% removed", "50% removed", "75% removed", "100% removed")
GRWR_simulated_all <- GRWR_simulated_all %>% 
  gather(key = "treatment", "GRWR", -Year)

GRWR_simulated_all$treatment <- ordered(GRWR_simulated_all$treatment, levels = c("reference", "0% removed", "25% removed", "50% removed", "75% removed", "100% removed"))


ggplot(GRWR_simulated_all%>%filter(!treatment %in% c("25% removed", "100% removed")), aes(x = Year, y = GRWR, group = treatment))+
  geom_rect(aes(xmin = c(2001.5), xmax = c(2003.5), ymin = -4, ymax = 4), fill = "#f2f2f2")+
  geom_rect(aes(xmin = c(2004.5), xmax = c(2006.5), ymin = -4, ymax = 4), fill = "#f2f2f2")+
  geom_rect(aes(xmin = c(2010.5), xmax = c(2011.5), ymin = -4, ymax = 4), fill = "#f2f2f2")+
  geom_rect(aes(xmin = c(2012.5), xmax = c(2013.5), ymin = -4, ymax = 4), fill = "#f2f2f2")+
  geom_rect(aes(xmin = c(2014.5), xmax = c(2015.5), ymin = -4, ymax = 4), fill = "#f2f2f2")+
  geom_point(aes(color = treatment))+
  geom_line(size=0.8, aes(color = treatment, linetype = treatment))+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom")+
  labs(x = "Time (year)", y = "Growth rate when rare") +
  scale_color_manual( values = c("#888888", "#000000", "#000000",  "#000000"))+
  scale_linetype_manual(values = c("solid", "solid",  "dashed", "dotted"))

AvgGRWR_all <- GRWR_simulated_all %>%
  filter(!is.na(GRWR))%>%
  group_by(treatment) %>%
  summarize(overall_GRWR = mean(GRWR))

sim_GRWR <- ggplot(AvgGRWR_all %>% filter(!treatment %in% c("reference",
                    "25% removed", "100% removed")), aes(x = treatment, y = overall_GRWR))+
                    geom_bar(stat = "identity")+
                    theme(text = element_text(size=16),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.text = element_text(size = 18),
                          axis.line = element_line(colour = "black"),
                          legend.position = "bottom")+
                    labs(y = "Average GRWR", x = "Grass Removal Treatment")+
                    geom_hline(yintercept = 0) +
                    annotate("text", x = 1, y = 0.05, label = "Does not persist", size = 5)+
                    annotate("text", x = 2, y = 0.05, label = "Does not persist", size = 5)+
                    annotate("text", x = 3, y = 0.35, label = "Persists", size = 5)
                           