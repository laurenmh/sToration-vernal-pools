# Run "prep data before modeling.R" first.
# Run "complex_belowground_v5.R" and "reference_pool_model.R" with real data.
# Run "timeseries visuals.R"

# Question:
## Does lambda provide a better regression than abundance measures?

# Outline:
## I. Plot LACO dens constructed vs. reference
## II. Plot Lambda constructed vs. reference
## III. Plot Relative growth rate constructed vs. reference

# Set up:
library(ggplot2)
library(tidyverse)

##############################################
# I. Plot LACO dens constructed vs. reference#
##############################################
spread_join_LACO <- spread(join_LACO, type, LACOdens) #join_LACO from "timeseries visuals.R"
mean_LACOdens <- spread_join_LACO %>%
  select(Year, constructed, reference) %>%
  filter(Year != c("2000", "2001")) %>%
  group_by(Year) %>%
  summarise_at(c("constructed", "reference"), mean, na.rm = TRUE)

summary(lm(constructed ~ reference, mean_LACOdens))
ggplot(mean_LACOdens, aes(y = constructed, x = reference)) +
  geom_point() +
  ylim(0, 500) +
  xlim(0, 250) +
  labs(y = "Constructed mean LACO density",
       x = "Reference mean LACO density") +
  geom_smooth(method = "lm") +
  annotate("text", label = "R^2 = 0.232
           adjusted R^2 = 0.168
           p-value = 0.081", x = 50, y = 400) +
  theme_bw() +
  geom_text(aes(label = Year), hjust = 0, vjust = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted")

#Differentiate by seeding treatment 
mean_ref_LACOdens <- spread_join_LACO %>%
  select(Year, reference) %>%
  filter(Year != c("2000", "2001")) %>%
  group_by(Year) %>%
  summarise_at(c("reference"), mean, na.rm = TRUE)
mean_const_LACOdens <- spread_join_LACO %>%
  select(Year, constructed, treatment) %>%
  filter(Year != c("2000", "2001")) %>%
  filter(treatment != c("NA NA")) %>%
  group_by(Year, treatment) %>%
  summarise_at(c("constructed"), mean, na.rm = TRUE)
mean_LACOdens_trt <- merge(mean_const_LACOdens, mean_ref_LACOdens, by = "Year")

summary(lm(constructed ~ reference, mean_LACOdens_trt))
ggplot(mean_LACOdens_trt, aes(y = constructed, x = reference)) +
  geom_point(aes(col = treatment)) +
  ylim(0, 500) +
  xlim(0, 250) +
  geom_smooth(method = "lm") +
  labs(y = "Constructed mean LACO density",
       x = "Reference mean LACO density") +
  annotate("text", label = "R^2 = 0.1822 
          adjusted R^2 = 0.167
          p-value = 0.001", x = 50, y = 400) + 
  theme_bw()

ggplot(mean_LACOdens_trt, aes(y = constructed, x = reference, col = treatment)) +
  geom_point() +
  ylim(0, 500) +
  xlim(0, 250) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = "Constructed mean LACO density",
       x = "Reference mean LACO density", col = "Seeding treatment") +
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
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO lambda", x = "Reference LACO lambda") +
  annotate("text", label = "R^2 = 0.1791
           adjusted R^2 = 0.1045
           p-value = 0.1496", x = 20, y = 50) +
  theme_bw() +
  geom_text(aes(label = row.names(join_lambda_trim)), hjust = 0, vjust =0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted")

############################################################
# III. Plot Relative Growth Rate constructed vs. reference #
############################################################
# relative growth rate of constructed 2002-2015
const_GR <- long_relative_GR %>% 
  select(Year, Pool, GR) %>% 
  mutate_if(is.numeric, ~replace(., is.infinite(.), 1)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  group_by(Year) %>%
  summarise(constructed = mean(GR)) %>%
  filter(Year != 2016) %>%
  filter(Year != 2017) %>%
  filter(Year != 2001) %>%
  filter(Year != 2002)

ref_GR <- ref_com_LACO %>%
  group_by(Pool) %>%
  summarise('2003' = `2003`/`2002`,
            '2004' = `2004`/`2003`,
            '2005' = `2005`/`2004`,
            '2006' = `2006`/`2005`,
            '2007' = `2007`/`2006`,
            '2008' = `2008`/`2007`,
            '2009' = `2009`/`2008`,
            '2010' = `2010`/`2009`,
            '2011' = `2011`/`2010`,
            '2012' = `2012`/`2011`,
            '2013' = `2013`/`2012`,
            '2014' = `2014`/`2013`,
            '2015' = `2015`/`2014`) %>%
  gather(`2003`, `2003`, `2004`, `2005`, `2006`, 
         `2007`, `2008`, `2009`, `2010`, `2011`, `2012`, `2013`,
         `2014`, `2015`,
         key = Year, value = GR) %>%
  mutate_if(is.numeric, ~replace(., is.infinite(.), 1)) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  group_by(Year) %>%
  summarise(reference = mean(GR))

join_GR <- left_join(ref_GR, const_GR, by = 'Year')

ggplot(join_GR, aes(x = reference, y = constructed)) +
  geom_point() +
  theme_bw() +
  geom_text(aes(label = Year), hjust = .5, vjust = 0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  labs(y = "Constructed mean LACO growth rate",
       x = "Reference mean LACO growth rate") 
