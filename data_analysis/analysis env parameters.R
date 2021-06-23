################################################################################
#What is the relationship between env and model parameters (alphas and lambdas)#
################################################################################
#Packages
library(tidyverse)
library(ggplot2)
library(stargazer)
library(ggpubr)

#Load data
PPT <- read_csv(paste(datpath, "Monthly precip averages/Fairfield_precip.csv", sep="")) 
source("data_compiling/compile_constructed_depth.R")

#Summarize PPT data 2001-2017
PPT_0117 <- PPT %>%
  filter(!Year %in% c("1999", "2017", "2018" , "2019")) %>%
  select(Year, Oct_Dec_cm, Jan_March_cm, Total_ppt_cm)

#Summarize inundation data 
inundation <- const_depth %>%
  filter(Location == "SeedPlot", Treatment.1999 != "Control") %>%
  group_by(Year, Pool, Distance, Size) %>%
  summarise(max_depth = as.numeric(max(Depth)),
            duration_wk = as.numeric(max(Duration.weeks))) %>% #inundation data available for 2000, 2002, 2009, 2010, 2011, 2012
  filter(max_depth <= 10) %>%
  filter(duration_wk >= 0) %>%
  group_by(Year) %>%
  summarise(mean_max_depth = mean(max_depth),
            mean_duration_wk = mean(duration_wk))

# Make a lambda table  2000-2016
lambda_ref_full <- as.data.frame(reflambda_mean[,5]) %>%
  mutate(Year = c(2002:2014))
colnames(lambda_ref_full) <- c("lambda_ref", "Year")  
lambda_const_full <- as.data.frame(lambda_mean[,5]) %>%
  mutate(Year = c(2000:2016))
colnames(lambda_const_full) <- c("lambda_const", "Year")
lambda_join <- left_join(lambda_const_full, lambda_ref_full)
lambda_long <- lambda_join %>%
  gather(key = "type", value = "lambda", lambda_const, lambda_ref) %>%
  mutate(type = ifelse(type == "lambda_const", "constructed", "reference"))

# Make an alpha_LACO table 2000-2016
alpha_LACO_ref <- as.data.frame(refalpha_LACO_mean[,5]) %>%
  mutate(Year = c(2002:2014))
colnames(alpha_LACO_ref) <- c("alpha_LACO_ref", "Year")
alpha_LACO_const <- as.data.frame(alpha_LACO_mean[,5]) %>%
  mutate(Year = c(2000:2016))
colnames(alpha_LACO_const) <- c("alpha_LACO_const", "Year")
alpha_LACO_join <- left_join(alpha_LACO_const, alpha_LACO_ref)
alpha_LACO_long <- alpha_LACO_join %>%
  gather(key = "type", value = "alpha_LACO", alpha_LACO_const, alpha_LACO_ref) %>%
  mutate(type = ifelse(type == "alpha_LACO_const", "constructed", "reference"))

# Make an alpha_EG table 2000-2016
alpha_EG_ref <- as.data.frame(refalpha_EG_mean[,5]) %>%
  mutate(Year = c(2002:2014))
colnames(alpha_EG_ref) <- c("alpha_EG_ref", "Year")
alpha_EG_const <- as.data.frame(alpha_EG_mean[,5]) %>%
  mutate(Year = c(2000:2016))
colnames(alpha_EG_const) <- c("alpha_EG_const", "Year")
alpha_EG_join <- left_join(alpha_EG_const, alpha_EG_ref)
alpha_EG_long <- alpha_EG_join %>%
  gather(key = "type", value = "alpha_EG", alpha_EG_const, alpha_EG_ref) %>%
  mutate(type = ifelse(type == "alpha_EG_const", "constructed", "reference"))

# Make an alpha_ERVA table 2000-2016
alpha_ERVA_ref <- as.data.frame(refalpha_ERVA_mean[,5]) %>%
  mutate(Year = c(2002:2014))
colnames(alpha_ERVA_ref) <- c("alpha_ERVA_ref", "Year")
alpha_ERVA_const <- as.data.frame(alpha_ERVA_mean[,5]) %>%
  mutate(Year = c(2000:2016))
colnames(alpha_ERVA_const) <- c("alpha_ERVA_const", "Year")
alpha_ERVA_join <- left_join(alpha_ERVA_const, alpha_ERVA_ref)
alpha_ERVA_long <- alpha_ERVA_join %>%
  gather(key = "type", value = "alpha_ERVA", alpha_ERVA_const, alpha_ERVA_ref) %>%
  mutate(type = ifelse(type == "alpha_ERVA_const", "constructed", "reference"))

# Make an alpha_NF table 2000-2016
alpha_NF_ref <- as.data.frame(refalpha_NF_mean[,5]) %>%
  mutate(Year = c(2002:2014))
colnames(alpha_NF_ref) <- c("alpha_NF_ref", "Year")
alpha_NF_const <- as.data.frame(alpha_NF_mean[,5]) %>%
  mutate(Year = c(2000:2016))
colnames(alpha_NF_const) <- c("alpha_NF_const", "Year")
alpha_NF_join <- left_join(alpha_NF_const, alpha_NF_ref)
alpha_NF_long <- alpha_NF_join %>%
  gather(key = "type", value = "alpha_NF", alpha_NF_const, alpha_NF_ref) %>%
  mutate(type = ifelse(type == "alpha_NF_const", "constructed", "reference"))

# Join all parameters and PPT
PPT_parameters <- left_join(PPT_0117, lambda_join) %>%
  left_join(., alpha_LACO_join ) %>%
  left_join(., alpha_EG_join) %>%
  left_join(., alpha_ERVA_join) %>%
  left_join(., alpha_NF_join) %>%
  left_join(., inundation)

# Join just the parameters
Parameters <- left_join(lambda_long, alpha_LACO_long) %>%
  left_join(., alpha_EG_long) %>%
  left_join(., alpha_ERVA_long) %>%
  left_join(., alpha_NF_long)

Parameters$type <- factor(Parameters$type, levels = c("reference", "constructed"))

f1 <- ggplot(Parameters, aes(x = lambda, y = alpha_LACO)) +
  geom_point() +
  theme_bw()+
  facet_wrap(~type)+
  labs(x = NULL, y = "alpha LACO") +
  ylim(0, 1.00)

f2 <- ggplot(Parameters, aes(x = lambda, y = alpha_EG)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~type) +
  labs(x = NULL, y = "alpha EG") +
  ylim(0, 1.00)

f3 <- ggplot(Parameters, aes(x = lambda, y = alpha_ERVA)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~type) +
  labs(x = NULL, y = "alpha ERVA") +
  ylim(0, 1.00)

f4 <- ggplot(Parameters, aes(x = lambda, y = alpha_NF)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~type) +
  labs(x = "lambda", y = "alpha NF") +
  ylim(0, 1.00)

ggarrange(f1, f2, f3, f4,  ncol = 1, nrow = 4, 
          labels = c("a)", "b)",
                     "c)", "d)"), 
          common.legend = TRUE, legend = "bottom", 
          font.label = list(size = 10),
          heights = c(1,1,1,1.15))

pairs(~lambda_const + alpha_LACO_const + alpha_EG_const + alpha_ERVA_const + Oct_Dec_cm + Jan_March_cm + Total_ppt_cm +  mean_max_depth, PPT_parameters)
pairs(~lambda_ref + alpha_LACO_ref + alpha_EG_ref + alpha_ERVA_ref + Oct_Dec_cm + Jan_March_cm + Total_ppt_cm + mean_max_depth, PPT_parameters)
pairs(~lambda_const + alpha_LACO_const + alpha_EG_const + alpha_ERVA_const, PPT_parameters)
pairs(~lambda_ref + alpha_LACO_ref + alpha_EG_ref + alpha_ERVA_ref, PPT_parameters)

# Make a regression table
# dt <- matrix(ncol=5, nrow = 5)
# rownames(dt) <- c("early ppt", "late ppt", "total ppt", "max depth", "max duration")
# colnames(dt) <- c("lambda", "alpha_LACO", "alpha_EG", "alpha_ERVA", "alpha_NF")
# 
# parenthesis <- function(x){
#     x[i] <- paste("(" , x[i]) %>%
#             paste(., ")")
# }
summary(lm(lambda_const~Oct_Dec_cm, PPT_parameters))
summary(lm(lambda_const~Jan_March_cm, PPT_parameters))
summary(lm(lambda_const~Total_ppt_cm, PPT_parameters))
summary(lm(lambda_const~mean_max_depth, PPT_parameters))
summary(lm(lambda_const~mean_duration_wk, PPT_parameters))

summary(lm(alpha_LACO_const~Oct_Dec_cm, PPT_parameters))
summary(lm(alpha_LACO_const~Jan_March_cm, PPT_parameters))
summary(lm(alpha_LACO_const~Total_ppt_cm, PPT_parameters))
summary(lm(alpha_LACO_const~mean_max_depth, PPT_parameters))
summary(lm(alpha_LACO_const~mean_duration_wk, PPT_parameters))

summary(lm(alpha_EG_const~Oct_Dec_cm, PPT_parameters))
summary(lm(alpha_EG_const~Jan_March_cm, PPT_parameters))
summary(lm(alpha_EG_const~Total_ppt_cm, PPT_parameters))
summary(lm(alpha_EG_const~mean_max_depth, PPT_parameters))
summary(lm(alpha_EG_const~mean_duration_wk, PPT_parameters))


summary(lm(alpha_ERVA_const~Oct_Dec_cm, PPT_parameters))
summary(lm(alpha_ERVA_const~Jan_March_cm, PPT_parameters))
summary(lm(alpha_ERVA_const~Total_ppt_cm, PPT_parameters))
summary(lm(alpha_ERVA_const~mean_max_depth, PPT_parameters))
summary(lm(alpha_ERVA_const~mean_duration_wk, PPT_parameters))

summary(lm(alpha_NF_const~Oct_Dec_cm, PPT_parameters))
summary(lm(alpha_NF_const~Jan_March_cm, PPT_parameters))
summary(lm(alpha_NF_const~Total_ppt_cm, PPT_parameters))
summary(lm(alpha_NF_const~mean_max_depth, PPT_parameters))
summary(lm(alpha_NF_const~mean_duration_wk, PPT_parameters))

summary(lm(lambda_ref~Oct_Dec_cm, PPT_parameters))
summary(lm(alpha_LACO_ref~Oct_Dec_cm, PPT_parameters))
summary(lm(alpha_EG_ref~Oct_Dec_cm, PPT_parameters))
summary(lm(alpha_ERVA_ref~Oct_Dec_cm, PPT_parameters))
summary(lm(alpha_NF_ref~Oct_Dec_cm, PPT_parameters))

summary(lm(lambda_ref~Jan_March_cm, PPT_parameters))
summary(lm(alpha_LACO_ref~Jan_March_cm, PPT_parameters))
summary(lm(alpha_EG_ref~Jan_March_cm, PPT_parameters))
summary(lm(alpha_ERVA_ref~Jan_March_cm, PPT_parameters))
summary(lm(alpha_NF_ref~Jan_March_cm, PPT_parameters))

summary(lm(lambda_ref~Total_ppt_cm, PPT_parameters))
summary(lm(alpha_LACO_ref~Total_ppt_cm, PPT_parameters))
summary(lm(alpha_EG_ref~Total_ppt_cm, PPT_parameters))
summary(lm(alpha_ERVA_ref~Total_ppt_cm, PPT_parameters))
summary(lm(alpha_NF_ref~Total_ppt_cm, PPT_parameters))

summary(lm(lambda_ref~mean_max_depth, PPT_parameters))
summary(lm(alpha_LACO_ref~mean_max_depth, PPT_parameters))
summary(lm(alpha_EG_ref~mean_max_depth, PPT_parameters))
summary(lm(alpha_ERVA_ref~mean_max_depth, PPT_parameters))
summary(lm(alpha_NF_ref~mean_max_depth, PPT_parameters))

summary(lm(lambda_ref~mean_duration_wk, PPT_parameters))
summary(lm(alpha_LACO_ref~mean_duration_wk, PPT_parameters))
summary(lm(alpha_EG_ref~mean_duration_wk, PPT_parameters))
summary(lm(alpha_ERVA_ref~mean_duration_wk, PPT_parameters))
summary(lm(alpha_NF_ref~mean_duration_wk, PPT_parameters))

#stargazer(dt, type = "text")

