# Run "prep data before modeling.R" first.
# Run "complex_belowground_v5.R" and "reference_pool_model.R" with real data.
# Run "timeseries visuals.R"

# Question:
## Does lambda provide a better regression than abundance measures?

# Outline:
## I. Plot LACO dens constructed vs. reference
## II. Plot Lambda constructed vs. reference
## III. Plot LACO Nt/Nt-1 constructed vs. reference

# Set up:
library(ggplot2)
library(tidyverse)
library(ggpubr)

##############################################
# I. Plot LACO dens constructed vs. reference#
##############################################
#join_LACO from "timeseries visuals.R"
mean_LACOdens <-join_LACO %>%
  filter(Year !=c("2000","2001")) %>%
  group_by(type, Year) %>%
  summarise(mean_LACO = mean(LACOdens),
            mean_logLACO = mean(log_LACOdens)) %>%
  pivot_wider(names_from = type, values_from = c(mean_LACO, mean_logLACO))

summary(lm(mean_LACO_constructed ~ mean_LACO_reference, mean_LACOdens))
ggplot(mean_LACOdens, aes(y = mean_LACO_constructed, x = mean_LACO_reference)) +
  geom_point() +
  labs(y = "Constructed mean LACO density",
       x = "Reference mean LACO density") +
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method = "lm") +
  annotate("text", label = "R^2 = 0.232
           adjusted R^2 = 0.168
           p-value = 0.081", x = 60, y = 400) +
  theme_bw() +
  geom_text(aes(label = Year), hjust = 0, vjust = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted")

#Log-log analysis of LACO density
ggplot(mean_LACOdens, aes(y = mean_logLACO_constructed, x = mean_logLACO_reference)) +
  geom_point()+
  theme_bw()+
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  ylim(0,6.2)+
  xlim(0,6.2)+
  geom_text(aes(label = Year), hjust = 0, vjust = 0) +
  labs(y = "Constructed mean LACO density (log)",
       x = "Reference mean LACO density (log)")

#Differentiate by seeding treatment 
mean_LACOdens_const <- join_LACO %>%
  filter(Year !=c("2000","2001")) %>%
  filter(type == "constructed") %>%
  group_by(Year, treatment) %>%
  summarise(mean_LACO_const = mean(LACOdens),
            mean_logLACO_const = mean(log_LACOdens)) 

mean_LACOdens_ref <- join_LACO %>%
  filter(Year !=c("2000", "2001")) %>%
  group_by(Year) %>%
  summarise(mean_LACO_ref = mean(LACOdens),
            mean_logLACO_ref = mean(log_LACOdens)) 

mean_LACOdens_trt <- merge(mean_LACOdens_const, mean_LACOdens_ref)

summary(lm(mean_LACO_const ~ mean_LACO_ref, mean_LACOdens_trt))
ggplot(mean_LACOdens_trt, aes(y = mean_LACO_const, x = mean_LACO_ref)) +
  geom_point(aes(col = treatment)) +
  labs(y = "Constructed mean LACO density",
       x = "Reference mean LACO density",  col = "Seeding treatment") +
  geom_smooth(aes(col = treatment), method = "lm", se = FALSE) + 
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()

#############################################
# II. Plot lambda constructed vs. reference #
#############################################
# combine lambda estimates from constructed and reference models
const_lambda_trim <- as.data.frame(lambda_mean[-c(1,2,16,17),5]) #trim 2001, 2002, 2016, and 2017 estimates
ref_lambda_trim <- as.data.frame(reflambda_mean[,5]) 
join_lambda_trim <- cbind(const_lambda_trim, ref_lambda_trim)
row.names(join_lambda_trim) <- c(2003:2015)
colnames(join_lambda_trim) <- c("constructed", "reference")

summary(lm(constructed ~ reference, join_lambda_trim))
ggplot(join_lambda_trim, aes(y = constructed, x = reference)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO lambda", x = "Reference LACO lambda") +
  annotate("text", label = "R^2 = 0.1791
           adjusted R^2 = 0.1045
           p-value = 0.1496", x = 20, y = 50) +
  theme_bw() +
  geom_text(aes(label = row.names(join_lambda_trim)), hjust = 0, vjust =0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted")+
  xlim(0, 80) +
  ylim(-5, 70)

####################################################
# III. Plot LACO Nt/Nt-1 constructed vs. reference #
####################################################
#observed Nt/Nt-1 constructed vs. reference
const_relative_GR <- const_com_noNA %>% 
  mutate_if(is.numeric, ~replace(., . == 0, 1)) %>% #to avoid NAN and infinite numbers, we'll replace all zeros with one.
  group_by(Pool) %>%
  summarize('2003G' = `2003`/`2002`,
            '2004G' = `2004`/`2003`,
            '2005G' = `2005`/`2004`,
            '2006G' = `2006`/`2005`,
            '2007G' = `2007`/`2006`,
            '2008G' = `2008`/`2007`,
            '2009G' = `2009`/`2008`,
            '2010G' = `2010`/`2009`,
            '2011G' = `2011`/`2010`,
            '2012G' = `2012`/`2011`,
            '2013G' = `2013`/`2012`,
            '2014G' = `2014`/`2013`,
            '2015G' = `2015`/`2014`) %>% 
  gather(`2003G`, `2004G`, `2005G`, `2006G`, `2007G`, `2008G`, 
         `2009G`, `2010G`, `2011G`, `2012G`, `2013G`, `2014G`, `2015G`, key = Year, value = const_GR) %>%
  group_by(Year) %>%
  summarize(mean_const_GR = mean(const_GR),
            log_mean_const_GR = log(mean_const_GR)) %>%
  mutate(Year = gsub("G", "", Year) )

ref_relative_GR <- ref_com_LACO %>%
  mutate_if(is.numeric, ~replace(., . == 0, 1)) %>% #to avoid NAN and infinite numbers, we'll replace all zeros with one.
  group_by(Pool) %>%
  summarize('2003G' = `2003`/`2002`,
            '2004G' = `2004`/`2003`,
            '2005G' = `2005`/`2004`,
            '2006G' = `2006`/`2005`,
            '2007G' = `2007`/`2006`,
            '2008G' = `2008`/`2007`,
            '2009G' = `2009`/`2008`,
            '2010G' = `2010`/`2009`,
            '2011G' = `2011`/`2010`,
            '2012G' = `2012`/`2011`,
            '2013G' = `2013`/`2012`,
            '2014G' = `2014`/`2013`,
            '2015G' = `2015`/`2014`) %>%
  gather(`2003G`, `2004G`, `2005G`, `2006G`, `2007G`, `2008G`, 
         `2009G`, `2010G`, `2011G`, `2012G`, `2013G`, `2014G`, `2015G`, key = Year, value = ref_GR) %>%
  group_by(Year) %>%
  summarize(mean_ref_GR = mean(ref_GR), 
            log_mean_ref_GR = log(mean_ref_GR)) %>%
  mutate(Year = gsub("G", "", Year) )

join_relative_GR <- left_join(const_relative_GR, ref_relative_GR, by = "Year") #join two tables

summary(lm(log_mean_const_GR ~ log_mean_ref_GR, join_relative_GR))
ggplot(join_relative_GR, aes(y = log_mean_const_GR, x = log_mean_ref_GR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO observed growth rate (log)", x = "Reference LACO observed growth rate (log)") +
  theme_bw() +
  geom_text(aes(label = row.names(join_lambda_trim)), hjust = 0, vjust =0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  annotate("text", label = "R^2 = 0.614
           adjusted R^2 = 0.579
           p-value = 0.002", x = 0, y = 4) +
  xlim(-1, 3.5)

#predicted Nt/Nt-1 constructed vs. reference
predicted_const_GR <- predicted_LACO %>%
  spread(time, predicted_LACO) %>%
  mutate_if(is.numeric, ~replace(., . == 0, 1)) %>%
  group_by(Pool) %>%
  summarize('2003G' = `2003`/`2002`,
            '2004G' = `2004`/`2003`,
            '2005G' = `2005`/`2004`,
            '2006G' = `2006`/`2005`,
            '2007G' = `2007`/`2006`,
            '2008G' = `2008`/`2007`,
            '2009G' = `2009`/`2008`,
            '2010G' = `2010`/`2009`,
            '2011G' = `2011`/`2010`,
            '2012G' = `2012`/`2011`,
            '2013G' = `2013`/`2012`,
            '2014G' = `2014`/`2013`,
            '2015G' = `2015`/`2014`) %>%
  gather(`2003G`, `2004G`, `2005G`, `2006G`, `2007G`, `2008G`, 
         `2009G`, `2010G`, `2011G`, `2012G`, `2013G`, `2014G`, `2015G`, key = Year, value = const_GR) %>%
  group_by(Year) %>%
  summarize(mean_const_GR = mean(const_GR), 
            log_mean_const_GR = log(mean_const_GR)) %>%
  mutate(Year = gsub("G", "", Year) ) 

predicted_ref_GR <- ref_predicted_LACO %>%
  spread(time, predicted_LACO) %>%
  mutate_if(is.numeric, ~replace(., . == 0, 1)) %>%
  group_by(Pool) %>%
  summarize('2003G' = `2003`/`2002`,
            '2004G' = `2004`/`2003`,
            '2005G' = `2005`/`2004`,
            '2006G' = `2006`/`2005`,
            '2007G' = `2007`/`2006`,
            '2008G' = `2008`/`2007`,
            '2009G' = `2009`/`2008`,
            '2010G' = `2010`/`2009`,
            '2011G' = `2011`/`2010`,
            '2012G' = `2012`/`2011`,
            '2013G' = `2013`/`2012`,
            '2014G' = `2014`/`2013`,
            '2015G' = `2015`/`2014`) %>%
  gather(`2003G`, `2004G`, `2005G`, `2006G`, `2007G`, `2008G`, 
         `2009G`, `2010G`, `2011G`, `2012G`, `2013G`, `2014G`, `2015G`, key = Year, value = ref_GR) %>%
  group_by(Year) %>%
  summarize(mean_ref_GR = mean(ref_GR), 
            log_mean_ref_GR = log(mean_ref_GR)) %>%
  mutate(Year = gsub("G", "", Year) ) 

join_predicted_GR <- left_join(predicted_const_GR, predicted_ref_GR, by = "Year") #join two tables
summary(lm(log_mean_const_GR ~ log_mean_ref_GR, join_predicted_GR))
ggplot(join_predicted_GR, aes(y = log_mean_const_GR, x = log_mean_ref_GR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO predicted growth rate (log)", x = "Reference LACO predicted growth rate (log)") +
  theme_bw() +
  geom_text(aes(label = row.names(join_lambda_trim)), hjust = 0, vjust =0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  annotate("text", label = "R^2 = 0.604
           adjusted R^2 = 0.568
           p-value = 0.002", x = -0.5, y = 2) +
  xlim(-1.5, 2) 


### Are lambdas (intrinsic growth rate) correlated with LACO actual growth rate?
join_lambda_GR <- cbind(const_relative_GR, lambda_mean[-c(1,2,16,17),5])
colnames(join_lambda_GR) <- c("Year", "mean_GR", "log_mean_GR", "lambda")
ggplot(join_lambda_GR, aes(x = lambda, y = log_mean_GR)) +
  geom_point()+
  geom_text(aes(label = join_lambda_GR$Year), hjust = 0, vjust =0) +
  xlim(0, 75)+
  theme_bw()

############
# Figure 2 #
############

f2a <- ggplot(mean_LACOdens, aes(y = mean_LACO_constructed, x = mean_LACO_reference)) +
  geom_point() +
  labs(y = "Constructed LACO density",
       x = "Reference LACO density") +
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method = "lm") +
  geom_text(aes(x = 100, y = 1000, label = "R^2 == 0.232~~~p == 0.081"), parse = TRUE) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = '#F39C12', size = 1.2)+
  geom_text(aes(label = Year), hjust = +0.5, vjust = -0.5)

f2b <- ggplot(join_relative_GR, aes(y = mean_const_GR, x = mean_ref_GR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO observed growth rate", x = "Reference LACO observed growth rate") +
  theme_bw() +
  geom_text(aes(label = row.names(join_lambda_trim)), hjust = +0.5, vjust = -0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = '#F39C12', size = 1.2) +
  geom_text(aes(x =3, y = 65, label = "R^2 == 0.614~~~p == 0.002"), parse = TRUE) +
  scale_x_log10()+
  scale_y_log10()

f2c <- ggplot(join_lambda_trim, aes(y = constructed, x = reference)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO lambda", x = "Reference LACO lambda") +
  geom_text(aes(x = 35, y = 60, label = "R^2 == 0.179~~~p == 0.149"), parse = TRUE) +
  theme_bw() +
  geom_text(aes(label = row.names(join_lambda_trim)), hjust = +0.5, vjust = -0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = '#F39C12', size = 1.2)+
  xlim(0, 80) +
  ylim(-5, 70)

f2d <- ggplot(join_predicted_GR, aes(y = mean_const_GR, x = mean_ref_GR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO predicted growth rate", x = "Reference LACO predicted growth rate") +
  theme_bw() +
  geom_text(aes(label = row.names(join_lambda_trim)), hjust = +0.5, vjust = -0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = '#F39C12', size = 1.2) +
  geom_text(aes(x =1, y = 12, label = "R^2 == 0.604~~~p == 0.002"), parse = TRUE) +
  scale_x_log10()+
  scale_y_log10()

ggarrange(f2a, f2b, f2c, f2d,  ncol = 2, nrow = 2, 
          labels = c("a)", "b)",
                     "c)", "d)"),
          font.label = list(size = 12))
