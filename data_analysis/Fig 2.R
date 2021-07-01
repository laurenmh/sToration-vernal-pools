# Script for FIGURE 2

# Set your data pathway first!

# Data 
# Go to "analysis/complex_belowground_v5.R" for lambda values in restored pools
# Go to "analysis/reference_pool_model.R" for lambda values in reference pools 
# Go to "analysis/GRWR_invader.R" for GRWR values

# Packages
library(ggplot2)
library(ggpubr)
library(HDInterval)

# Visualize timeseries of predicted LACO lambda 
Post_ref_00_01 <- rstan::extract(BH_ref_fit) #Run the model with 2000-2015 data (7 pools)
Post_ref_02_14 <- rstan::extract(BH_ref_fit) #Run the model with 2002-2015 data (9 pools)
CI_lambda_ref <- as.data.frame(HDInterval::hdi(Post_ref_00_01$lambda[,1:2], credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2000:2001)) %>%
  cbind(., as.data.frame(HDInterval::hdi(Post_ref_02_14$lambda, credMass = 0.95))) %>%
  magrittr::set_colnames(c(2000:2014))%>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)

lambda_ref <- as.data.frame(Post_ref_00_01$lambda[,1:2])%>% #2000-2001 lambda data 
  magrittr::set_colnames(c(2000:2001))%>%
  cbind(., as.data.frame(Post_ref_02_14$lambda)) %>%#piece together 2000-2001 data with 2002-2014 data
  magrittr::set_colnames(c(2000:2014))%>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(type = "reference") %>%
  full_join(., CI_lambda_ref)

Post <- rstan::extract(BH_fit) #Extract parameter estimates from 'complex_belowground_v5.R'
CI_lambda <-  as.data.frame(HDInterval::hdi(Post$lambda, credMass = 0.95)) %>% #Calculate 95% credible interval of lambda estimates
  magrittr::set_colnames(c(2000:2016)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
lambda_const <- as.data.frame(Post$lambda) %>% #Calculate mean lambda each year
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "lambda")) %>%
  group_by(Year) %>%
  summarise(mean = mean(lambda)) %>%
  mutate(type = "constructed") %>%
  full_join(., CI_lambda)

lambda_const_ref <- rbind(lambda_ref, lambda_const) #combine lambda tables
lambda_const_ref$Year <- as.numeric(lambda_const_ref$Year)

flambda <- ggplot(lambda_const_ref, aes(x = Year, y = mean, col = type))+
  geom_point() +
  geom_line(size=2)+
  geom_errorbar(aes(ymin = mean-lowCI, ymax = mean+upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.title = element_text(size = 14))+
  ylab(bquote(Intrinsic~Growth~Rate ~(lambda[t])))+
  scale_x_continuous(name = NULL,
                     limits = c(1999.5,2016.5))+
  scale_color_manual(name = "", values = c("#000000", "#888888"))

#Visualize timeseries of GRWR (see 'GRWR_invader.R')
CI_GRWR_ref <-  as.data.frame(HDInterval::hdi(GRWR_LACO_ref, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2002:2014)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_ref <- as.data.frame(GRWR_LACO_ref) %>%
  magrittr::set_colnames(c(2002:2014)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(type = "reference") %>%
  full_join(., CI_GRWR_ref)

CI_GRWR <-  as.data.frame(HDInterval::hdi(GRWR_LACO_const, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2000:2016)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_const <- as.data.frame(GRWR_LACO_const) %>%
  magrittr::set_colnames(c(2000:2016)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(type = "constructed") %>%
  full_join(., CI_GRWR)

GRWR_time <- rbind(GRWR_time_ref, GRWR_time_const) #combine lambda tables
GRWR_time$Year <- as.numeric(GRWR_time$Year)
  
fGRWR <- ggplot(GRWR_time, aes(x = Year, y = mean, col = type))+
  geom_point() +
  geom_line(size=2)+
  geom_errorbar(aes(ymin = mean-lowCI, ymax = mean+upCI), width = 0.4, alpha = 0.9, size = 1) +
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
                     limits = c(1999.5,2016.5))+
  scale_color_manual(name = "", values = c("#000000", "#888888"))

#FIGURE 2
Fig2 <- ggarrange(flambda, fGRWR, ncol = 1, nrow = 2, align = "v", 
                  labels = c("(a)", "(b)"), font.label = list(size = 14))
