# Run complex_belowground_v2 with multiple year seeding addition
# Question:
## How does the environment, namely inundation, change LACO growth rate and competition interactions?
# Outline:
## I. Inundation vs. lambda
## II. Inundation vs. alpha
library(ggplot2)

### I. INUNDATION VS. LAMBDA ###

#1.Extract lambda values from the model

# Run the model first
se<-function(x){
  sd(x)/sqrt(length(x))
} # this is a function for calculating standard error

lambda_data <- as.data.frame(get_posterior_mean(BH_fit, pars = c("lambda")))

lambda_gather <- lambda_data %>%
  mutate(Year = c(2001:2017)) %>%
  select(Year, `mean-chain:1`, `mean-chain:2`, `mean-chain:3`, `mean-chain:4`) %>%
  gather(`mean-chain:1`, `mean-chain:2`, `mean-chain:3`, `mean-chain:4`, key = chain, value = lambda)
  
lambda_summary <- lambda_gather %>%
  group_by(Year) %>%
  summarize(mean_lambda = mean(lambda), se_lambda = se(lambda)) #this is used for plotting

#2.Regression of max water depth and precipitation to predict inundation of all years

PPT <- read_csv(paste(datpath, "Monthly precip averages/Fairfield_precip.csv", sep="")) 
source("data_compiling/compile_constructed_depth.R")

max_depth <- const_depth %>%
  filter(Location == "SeedPlot", Treatment.1999 != "Control") %>%
  group_by(Year, Pool, Distance, Size) %>%
  summarize(max_depth = as.numeric(max(Depth)))

summary(lm(max_depth ~ Size, max_depth)) #no effect of size on max_depth
summary(lm(max_depth ~ Distance, max_depth)) #no effect of distance on max_depth

ggplot(max_depth, aes(x = as.factor(Year), y = max_depth, color = Size)) +
  geom_point()
ggplot(max_depth, aes(x = as.factor(Year), y = max_depth, color = Distance)) +
  geom_point()

PPT_subset <- PPT %>%
  select(Year, Jan_March_cm) %>%
  filter(Year %in% c(2000, 2002, 2009, 2010, 2011, 2012)) #subset the precip data by specific year

PPT_depth <- left_join(max_depth, PPT_subset, by = "Year") %>%
  filter(max_depth != 'NA') #combine ppt and water depth data

summary(lm(max_depth ~ Jan_March_cm, PPT_depth)) #max water depth is significantly correlated with Jan - March precipitation
# Equation: max_depth = 5.09207 + 0.03201 * PPT_Jan_March_cm #R-squared is 0.007

ggplot() +
  geom_bar(aes(x = PPT_subset$Year, y = PPT_subset$Jan_March_cm), stat= "identity", fill = "lightblue") +
  geom_point(aes(x = PPT_depth$Year, y = PPT_depth$max_depth)) +
  scale_y_continuous(name = "Max water depth (cm)",
                     sec.axis = sec_axis(~./1, name = "Precip Jan-Mar (cm)")) +
  theme_bw()

pred_max_depth <- PPT %>%
  select(Year, Jan_March_cm) %>%
  mutate(pred_max_depth = 5.09207 + 0.03201 * Jan_March_cm) %>%
  filter(Year %in% c(2001:2017))

#3.Join lambda and max water depth by year

lambda_real_depth <- left_join(max_depth_mean, lambda_gather, by = "Year") %>%
  filter(Year != "2000")

lambda_depth <- left_join(lambda_gather, pred_max_depth, by = "Year")

#4.Plot lambda and max water depth

ggplot(lambda_real_depth, aes(x = mean_max_depth, y = lambda)) +
  geom_point()

ggplot(lambda_depth, aes(x = pred_max_depth, y = lambda)) +
  geom_point()

#5.What's the fit?

summary(lm(pred_max_depth ~ lambda, lambda_depth)) 

### II. INUNDATION VS. ALPHA ###