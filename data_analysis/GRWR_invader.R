#Goal: Calculate the long-term growth rate when rare (r_invader) of LACO 

#Step 1. Calculate the stable equilibrium frequency of non-LACO species in the model each year. 
#Step 2. Calculate the annual growth rate of LACO in restored pools when one LACO is introduced into a stable community. 
#Step 3. Do the same step for the reference pools.
#Step 4. Average the growth rates of LACO over time for restored and reference pools.

#--------------------------

# Load data and package
# Remember to set your data pathway first!
source("data_compiling/compile_composition.R") 

#--------------------------

#Step 1. Calculate the stable equilibrium frequency of non-LACO species in the model each year. 

#Filter out just the control plots in reference pools (no LACO present). Average the frequency of non-LACO species across space each year. 
#Non-LACO species in our model:
#ERVA
#Exotic grass group - BRHO, HOMA, LOMU
#Native forb group - PLST, DOCO

const_com_control <- const_com %>% #constructed pools data
  filter(Treatment.1999 == "Control") %>% #filter control plots only
  drop_na() %>% #remove any rows with na
  filter(LACO <= 0) %>% #remove communities with LACO present
    mutate(sumEG = BRHO + HOMA + LOMU, 
         sumNF = PLST + DOCO) %>% #create a new column for sum of EG and sum of NF
  group_by(Year)%>% #summarize by year
  summarize(avg_ERVA = round(mean(ERVA), digits = 0),
            avg_sumEG = round(mean(sumEG), digits = 0),
            avg_sumNF = round(mean(sumNF), digits = 0)) %>%#take the average freq.
  filter(Year != "2017")

#Step 2. Calculate the annual growth rate of LACO in restored pools when one LACO is introduced into a stable community.

#Extract parameters for restored pools. Run "complex_belowground_v5.R".
alpha_LACO_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_LACO")))
alpha_EG_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_EG")))
alpha_ERVA_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_ERVA")))
alpha_NF_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_NF")))
lambda_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("lambda")))
s_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("survival_LACO")))
sim_LACO <- matrix(nrow = 17, ncol = 1)
#Set up the population model for LACO
bh.sim.control <- function(LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
  for(i in 1:1){
    sim_LACO[i] <- LACO*lambda[i]/(1+LACO*aii[i]+EG[i]*a1[i]+ERVA[i]*a2[i]+NF[i]*a3[i])+s*(1-g)*LACO/g #this is the modified Beverton-Holt model we'll use for LACO stem counts
  }
  for(i in 2:nrow(sim_LACO)){
    if (EG[i-1]> 100){
      g = glow
    }
    else{g = g}
    sim_LACO[i] <- LACO*lambda[i]/(1+LACO*aii[i]+EG[i]*a1[i]+ERVA[i]*a2[i]+NF[i]*a3[i])+s*(1-g)*LACO/g 
  }
 return(sim_LACO)
}
#Plug in the stable equilibrium freq and parameters in the model
LACO_const <- bh.sim.control(
                      LACO = 1,
                      EG = const_com_control$avg_sumEG,
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
GRWR_LACO_const <- log(LACO_const) #2000-2016

#Step 3. Do the same step for the reference pools.

#Extract parameters for reference pools. Run "reference_pool_model.R".
refalpha_LACO_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("alpha_LACO")))
refalpha_EG_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("alpha_EG")))
refalpha_ERVA_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("alpha_ERVA")))
refalpha_NF_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("alpha_NF")))
reflambda_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("lambda")))
refs_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("survival_LACO")))
sim_LACO <- matrix(nrow = 13, ncol = 1)
#Plug in the stable equilibrium freq and parameters in the model.
LACO_ref <- bh.sim.control(
                    LACO = 1,
                    EG = const_com_control$avg_sumEG,
                    ERVA = const_com_control$avg_ERVA,
                    NF = const_com_control$avg_sumNF,
                    aii = refalpha_LACO_mean[,5],
                    a1 = refalpha_EG_mean[,5],
                    a2 = refalpha_ERVA_mean[,5],
                    a3 = refalpha_NF_mean[,5],
                    lambda = reflambda_mean[,5],
                    s = refs_mean[,5],
                    g = 0.7,
                    glow = 0.2)
GRWR_LACO_ref <- log(LACO_ref) #2002-2014

#Step 4. Average the growth rates of LACO over time for restored and reference pools.
mean(GRWR_LACO_const) #Average GRWR LACO for constructed pools = -0.5818442
mean(GRWR_LACO_ref) #Average GRWR LACO for reference pools = 0.02065058 

#-----------------------
#Goal: simulate removing exotic grasses (EG) to promote LACO persistence
library(dplyr)
library(ggplot2)
library(stargazer)
library(ggpubr)

#Step 1. Figure out when to remove EG
#Join GRWR from restored and reference pools
GRWR_LACO_const <- as.data.frame(GRWR_LACO_const) %>%
  mutate(Year = c(2000:2016))
colnames(GRWR_LACO_const) <- c("GRWR_const", "Year")
GRWR_LACO_ref <- as.data.frame(GRWR_LACO_ref) %>%
  mutate(Year = c(2002:2014))
colnames(GRWR_LACO_ref) <- c("GRWR_ref", "Year")
GRWR_LACO <- left_join(GRWR_LACO_const, GRWR_LACO_ref)

#Load precipitation data
PPT <- read_csv(paste(datpath, "Monthly precip averages/Fairfield_precip.csv", sep="")) 

#Summarize PPT data 2001-2017
PPT_0117 <- PPT %>%
  filter(!Year %in% c("1999", "2017", "2018" , "2019")) %>%
  select(Year, Oct_Dec_cm, Jan_March_cm, Total_ppt_cm)

mean(PPT_0117$Total_ppt_cm)
mean(PPT_0117$Oct_Dec_cm)
mean(PPT_0117$Jan_March_cm)

#Combine GRWR_LACO and PPT data
GRWR_PPT <- left_join(GRWR_LACO, PPT_0117)

#EG in restored and reference pools
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

#Plot

ggplot(GRWR_PPT, aes(x = Oct_Dec_cm, y = GRWR_const))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(GRWR_PPT, aes(x = Oct_Dec_cm, y = GRWR_ref))+
  geom_point()+
  geom_smooth(method = "lm")

a <- ggplot(GRWR_PPT, aes(x = Year, y = Total_ppt_cm))+
        geom_hline(yintercept = 57.19, color = "red")+
        geom_bar(stat= "identity")+
        theme(axis.title.x = element_blank())

b <- ggplot(GRWR_PPT, aes(x = Year, y = Oct_Dec_cm))+
        geom_hline(yintercept = 23.86, color = "red")+
        geom_bar(stat="identity")+
        theme(axis.title.x = element_blank())

c <- ggplot(GRWR_PPT, aes(x = Year, y = Jan_March_cm))+
        geom_hline(yintercept = 27.14, color = "red")+      
        geom_bar(stat = "identity")+
        theme(axis.title.x = element_blank())

d <- ggplot(GRWR_PPT, aes(x = Year, y = GRWR_const))+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_hline(yintercept = -0.5818442, color = "red")+
        geom_point()+
        geom_line()+
        theme(axis.title.x = element_blank())

e <- ggplot(GRWR_PPT, aes(x = Year, y = GRWR_ref))+
        geom_hline(yintercept = 0, linetype = "dashed")+
        geom_hline(yintercept = 0.02065058, color = "red")+
        geom_point()+
        geom_line()+
        theme(axis.title.x = element_blank())

f <- ggplot(EG_summary_join%>%filter(Year != "2017"), aes(x = Year, y = EG_mean, group= type)) +
  geom_point(aes(col =type))+
  geom_line(aes(col =type))+
  geom_errorbar(aes(ymin = EG_mean-EG_se, ymax = EG_mean+EG_se, col =type)) +
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom")+
  labs(x = "Time (year)", y = "Mean exotic grass cover (%)") +
  scale_color_manual(name = "", values = c("#000000", "#888888"))

ggarrange(e, d, ncol = 1, nrow =2, labels = c("a)", "b"))

ggarrange(b, c, f,  ncol = 1, nrow = 3, 
          labels = c("a)", "b)",
                     "c)"), 
          #common.legend = TRUE, legend = "bottom", 
          font.label = list(size = 10),
          heights = c(1,1,1.5))
#Remove EG in years with above average Oct-Dec rain (2002, 2003, 2005, 2006, 2011, 2013, 2015)

#Step 2. Simulate EG removal
sumEG <- const_com_control %>%
  select(Year, avg_sumEG) %>%
  spread(key = "Year", value = "avg_sumEG")

#Remove 25% of EG in wet early-season year
mult.25 <- function(x)(x*0.75) 
reduced25EGcover <- sumEG %>%
  mutate_at(c("2002","2003", "2005","2006", "2011", "2013", "2015"), mult.25)

LACO_25EG <- bh.sim.control(
                      LACO = 1,
                      EG = as.matrix(reduced25EGcover),
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
GRWR_LACO_25EG <- log(LACO_25EG)

#Remove 50% of EG in wet early-season year
mult.5 <- function(x)(x*0.5)
reduced50EGcover <- sumEG %>%
  mutate_at(c("2002","2003", "2005","2006", "2011", "2013", "2015"), mult.5)

LACO_50EG <- bh.sim.control(
                      LACO = 1,
                      EG = as.matrix(reduced50EGcover),
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
GRWR_LACO_50EG <- log(LACO_50EG)

#Remove 75% of EG in wet early-season year
mult.75 <- function(x)(x*0.25)
reduced75EGcover <- sumEG %>%
  mutate_at(c("2002","2003", "2005","2006", "2011", "2013", "2015"), mult.75)

LACO_75EG <- bh.sim.control(
                      LACO = 1,
                      EG = as.matrix(reduced75EGcover),
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
GRWR_LACO_75EG <- log(LACO_75EG)

#Remove 100% of EG in wet early-season year
mult.100 <- function(x)(x*0)
reduced100EGcover <- sumEG %>%
  mutate_at(c("2002","2003", "2005","2006", "2011", "2013", "2015"), mult.100)

LACO_100EG <- bh.sim.control(
                      LACO = 1,
                      EG = as.matrix(reduced100EGcover),
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
GRWR_LACO_100EG <- log(LACO_100EG)


#Step 3. Average the growth rates of LACO over time for all simulation scenarios.
mean(GRWR_LACO_25EG) #Average GRWR when 25% EG cover removed = -0.5008274
mean(GRWR_LACO_50EG) #Average GRWR when 50% EG cover removed = -0.3951378
mean(GRWR_LACO_75EG) #Average GRWR when 75% EG cover removed = -0.2377705
mean(GRWR_LACO_100EG) #Average GRWR when 100% EG cover removed = 0.1683965

#Step 4. Plot simulated GRWR
GRWR_simulated <- cbind(GRWR_LACO, GRWR_LACO_25EG, GRWR_LACO_50EG, GRWR_LACO_75EG, GRWR_LACO_100EG) 
colnames(GRWR_simulated) <- c("constructed", "Year", "reference", "25% EG removed", "50% EG removed", "75% EG removed", "100% EG removed")
GRWR_simulated <- GRWR_simulated %>% 
  gather(key = "treatment", "GRWR", -Year)

GRWR_simulated$treatment <- ordered(GRWR_simulated$treatment, levels = c("reference", "constructed", "25% EG removed", "50% EG removed", "75% EG removed", "100% EG removed"))

ggplot(GRWR_simulated%>%filter(treatment != "25% EG removed"), aes(x = Year, y = GRWR, col = treatment))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom")+
  labs(x = "Year", y = "GRWR of LACO") +
  scale_color_manual(name = "", values = c("#ff4c4c", "#000000", "#0818A8", "#3F00FF",  "#40E0D0"))
