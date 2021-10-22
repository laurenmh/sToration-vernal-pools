# Script for NEW FIGURE 1
# Timeseries of LACO abundance, lambda, and GRWR

# Set your data pathway first!

# Data
source("data_wrangling/prep data before modeling.R")
# Go to "analysis/complex_belowground_v5.R" for lambda values in restored pools
# Go to "analysis/reference_pool_model.R" for lambda values in reference pools 
# Go to "analysis/GRWR_invader.R" for GRWR values

# Packages
library(ggplot2)
library(ggpubr)
library(HDInterval)
library(tidyverse)
library(dplyr)

# Function for standard error
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 

# REFERENCE pool LACO cover data
ref_LACO <- ref_com_LACO %>% #9 reference pools that have consecutive data on LACO abundance from 2002 to 2015
  group_by(Pool) %>%
  gather(key = "Year", value = "LACO" , - Pool) 
ref_LACO$Pool <- as.character(ref_LACO$Pool)

ref_LACOden_edge <- ref_com %>% #any reference LACO abundance data in 2000, 2001, 2016, 2017
  select(Year, Pool, LACO) %>%
  group_by(Year, Pool) %>%
  filter(Year %in% c(2000, 2001, 2016, 2017)) %>% 
  summarise_each(funs(mean))
ref_LACOden_edge$LACO <- as.integer(ref_LACOden_edge$LACO)

ref_LACOdens <- full_join(ref_LACOden_edge, ref_LACO)%>% #join reference LACO abundance data 2000-2017
  mutate(LACOdens = round(exp(-0.42)+LACO^1.36)) %>% # Convert LACO frequency to density 
  mutate(type = "reference")
ref_LACOdens$Year <- as.numeric(ref_LACOdens$Year)

# CONSTRUCTED pool LACO density data 
const_LACOden <- const_com_noNA %>% #72 pools
  select(- c(Size)) %>%
  group_by(Pool) %>%
  gather(key = "Year", value = "LACOdens", -Pool) %>%
  mutate(type = "constructed")
const_LACO <- left_join(const_LACOden, (const_com %>% select (Year, Pool, LACO, Treatment.1999, Treatment.2000)))

const_LACO$Year <- as.numeric(const_LACO$Year)
const_LACO$LACO <- as.integer(const_LACO$LACO)
const_LACO$Pool <- as.character(const_LACO$Pool)

# Join tables
join_LACO <- full_join(ref_LACOdens, const_LACO, all = TRUE) %>%
  mutate(treatment = paste(Treatment.1999, Treatment.2000)) %>%
  mutate(log_LACOdens = log(LACOdens)) %>%
  mutate_if(is.numeric, ~replace(., is.infinite(.), 0)) 

# Visualize timeseries of observed LACOdens
se <- function(x){
  sd(x)/sqrt(length(x))# this is a function for calculating standard error
} 

mean_join_LACO <- join_LACO %>%
  group_by(Year, type) %>%
  summarise(mean_LACOdens = mean(LACOdens),
            se_LACOdens = se(LACOdens))
mean_join_LACO$Year <- as.numeric(mean_join_LACO$Year)

fabundance <- ggplot(mean_join_LACO%>%filter(Year %in% c(2000:2015)), aes(x = Year, y = mean_LACOdens, col = type)) +
  geom_point()+
  geom_line(size=1)+
  scale_y_log10()+
  geom_errorbar(aes(ymin = mean_LACOdens-se_LACOdens, ymax = mean_LACOdens+se_LACOdens), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(.2, .3),
        axis.title = element_text(size = 14))+
  ylab(bquote(Density~(stems/m^2))) +
  scale_x_continuous(name = NULL, limits = c(1999.5,2015.5))+
  scale_color_manual(name = "", values = c("#000000", "#888888"))

# Visualize timeseries of predicted LACO lambda 
Post_ref_00_01 <- rstan::extract(BH_ref_fit_00_15) #Run the model with 2000-2015 data (7 pools)
Post_ref_02_14 <- rstan::extract(BH_ref_fit) #Run the model with 2002-2015 data (9 pools)
CI_lambda_ref <- as.data.frame(HDInterval::hdi(Post_ref_00_01$lambda[,1:2], credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2002)) %>%
  cbind(., as.data.frame(HDInterval::hdi(Post_ref_02_14$lambda, credMass = 0.95))) %>% #piece together 2000-2001 data with 2002-2014 data
  magrittr::set_colnames(c(2001:2015))%>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)

lambda_ref <- as.data.frame(Post_ref_00_01$lambda[,1:2])%>% #2000-2001 lambda data 
  magrittr::set_colnames(c(2001:2002))%>%
  cbind(., as.data.frame(Post_ref_02_14$lambda)) %>%#piece together 2000-2001 data with 2002-2014 data
  magrittr::set_colnames(c(2001:2015))%>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(type = "reference") %>%
  full_join(., CI_lambda_ref)

Post <- rstan::extract(BH_fit) #Extract parameter estimates from 'complex_belowground_v5.R'
CI_lambda <-  as.data.frame(HDInterval::hdi(Post$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2001:2017)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const <- as.data.frame(Post$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2001:2017)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(type = "constructed") %>%
  full_join(., CI_lambda)

lambda_const_ref <- rbind(lambda_ref, lambda_const) %>% #combine lambda tables
  filter(Year < 2016) # cut the last two points in constructed pools to match reference pools
lambda_const_ref$Year <- as.numeric(lambda_const_ref$Year)

flambda <- ggplot(lambda_const_ref, aes(x = Year, y = mean, col = type))+
  geom_point() +
  geom_line(size=1)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.2, 0.9), 
        axis.title = element_text(size = 14))+
  ylab(bquote(Intrinsic~Growth~Rate ~(lambda[t])))+
  scale_x_continuous(name = NULL, limits = c(1999.5,2015.5))+
  ylim(0,150)+
  scale_color_manual(name = "", values = c("#000000", "#888888"))

#Visualize timeseries of GRWR (see 'GRWR_invader.R')
CI_GRWR_ref <-  as.data.frame(HDInterval::hdi(GRWR_LACO_ref_00_15[,1:2], credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2002)) %>%
  cbind(., as.data.frame(HDInterval::hdi(GRWR_LACO_ref, credMass = 0.95))) %>% #piece together 2000-2001 data with 2002-2014 data
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_ref <- as.data.frame(GRWR_LACO_ref_00_15[,1:2]) %>%
  cbind(., as.data.frame(GRWR_LACO_ref))%>% #piece together 2000-2001 data with 2002-2014 data
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(type = "reference") %>%
  full_join(., CI_GRWR_ref)

CI_GRWR <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2001:2015)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const <- as.data.frame(GRWR_LACO_const) %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(type = "constructed") %>%
  full_join(., CI_GRWR)

GRWR_time <- rbind(GRWR_time_ref, GRWR_time_const) %>% #combine lambda tables
  filter(Year < 2016) # cut the last two points in constructed pools to match reference pools
GRWR_time$Year <- as.numeric(GRWR_time$Year)

fGRWR <- ggplot(GRWR_time, aes(x = Year, y = mean, col = type))+
  geom_point() +
  geom_line(size=1)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.2,0.2),
        axis.title = element_text(size = 14))+
  ylab(bquote(Low~Density~Growth~Rate~(italic(r[t]))))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(name = NULL,
                     limits = c(1999.5,2015.5))+
  scale_color_manual(name = "", values = c("#000000", "#888888"))

#FIGURE 1
Fig1_v2 <- ggarrange(fabundance, flambda, fGRWR, ncol = 1, nrow = 3, align = "v", 
                  labels = c("(a)", "(b)", "(c)"), font.label = list(size = 16))
