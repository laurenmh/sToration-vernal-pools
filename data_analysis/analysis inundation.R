# Run "prep data before modeling.R" then "complex_belowground_v3.R" first.

# Question:
## How does the environment, namely inundation, change LACO growth rate and competition interactions?

# Outline:
## I. Env vs. LACO counts
## II. Inundation vs. lambda
## III. Inundation vs. alpha

# Set up:
library(ggplot2)

se<-function(x){
  sd(x)/sqrt(length(x))
} # this is a function for calculating standard error

PPT <- read_csv(paste(datpath, "Monthly precip averages/Fairfield_precip.csv", sep="")) 
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R") #I haven't use this yet

##########################
# I. Env vs. LACO counts #
##########################

#Does the pond size matter?
#-> what do the timeseries of LACO look like by pond size?
ggplot(const_com, aes(x = Year, y = log(LACOdens), col = Size)) +
  geom_jitter() #no obvious difference of and LACO counts among pond sizes
summary_LACO_size <- const_com %>%
  drop_na() %>%
  group_by(Year, Size) %>%
  summarise(mean_LACO = mean(LACOdens),
            se_LACO = se(LACOdens))
ggplot(summary_LACO_size, aes(x = Year, col = Size)) +
  geom_point(aes(y = mean_LACO))+
  geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 1)

#Are relative growth rates of LACO responding similarly to env fluctuations across all pond sizes?
#-> what do the timeseries of relative growth rates look like by pond size?
relative_GR <- const_com_noNA %>%
  group_by(Pool, Size) %>%
  summarize('2001' = `2001`/`2000`,
            '2002' = `2002`/`2001`,
            '2003' = `2003`/`2002`,
            '2004' = `2004`/`2003`,
            '2005' = `2005`/`2004`,
            '2006' = `2006`/`2005`,
            '2007' = `2007`/`2006`,
            '2008' = `2008`/`2007`)
long_relative_GR <- relative_GR %>% gather(`2001`, `2002`, `2003`, `2003`, `2004`, `2005`, `2006`, `2007`, `2008`, key = Year, value = GR)
ggplot(long_relative_GR, aes(x = Year, y = log(GR), col = Size)) +
  geom_jitter() #no difference in LACO growth rates among pond sizes

#Which environmental variable is important for LACO?
#-> what does the correlation matrix of LACO stem counts, early ppt, max depth, inundation duration look like?
inundation_summary <- const_depth %>%
  filter(Location == "SeedPlot", Treatment.1999 != "Control") %>%
  group_by(Year, Pool, Distance, Size) %>%
  summarize(max_depth = as.numeric(max(Depth)),
            duration_wk = as.numeric(max(Duration.weeks))) %>% #inundation data available for 2000, 2002, 2009, 2010, 2011, 2012
  filter(max_depth <= 10)

env_summary <- left_join(inundation_summary, PPT%>% select(Year, Jan_March_cm), by = "Year") 
env_summary$Year <- as.character(env_summary$Year)
LACO_env <- left_join(env_summary, const_com%>%select(Year, Pool, LACOdens), by = c("Year", "Pool")) %>%
  drop_na() #use this combined data

pairs(~log(LACOdens)+max_depth+duration_wk+Jan_March_cm, data = LACO_env) #this makes a pairwise comparison  
LACO_env_cor <- cor(LACO_env[,5:8], method = c("pearson")) #this makes a correlation matrix
library(corrplot)
corrplot(LACO_env_cor, type = "upper", order = "hclust") #this makes a correlogram
library(Hmisc)
rcorr(as.matrix(LACO_env[,5:8])) #this makes a correlation matrix with significance levels
    #LACO count looks positively correlated with max depth (increases up to 10cm)
    #max depth and duration are strongly positively correlated
    #early PPT is lightly correlated with max depth
    #early PPT is lightly correlated with LACO count

ggplot(LACO_env, aes(x = max_depth, y = log(LACOdens))) +
  geom_point() +
  geom_smooth(method = "lm")+
  annotate(geom = "text", x = 2, y = 7, label = "R2 = 0.09, y = 0.197x + 1.378") +
  theme_bw()#LACO count pos correlates with max depth 
LACO_env$logLACOdens <- log(LACO_env$LACOdens)
summary(lm(logLACOdens ~ max_depth, LACO_env %>% filter(logLACOdens>=0))) 

ggplot(LACO_env, aes(x = Jan_March_cm, y = max_depth)) +
  geom_jitter() #max depth pos correlates with early PPT

#-> What is the overall effect of ppt on LACO?
PPT$Year <- as.character(PPT$Year)
LACO_ppt <- left_join(const_com %>% select(Year, Size, Pool, LACOdens), PPT %>% select(Year, Jan_March_cm, Jan, Feb, March), by = "Year") %>%
  drop_na()
summary(lm(LACOdens ~ Jan_March_cm, LACO_ppt)) #ppt has a significant effect on LACO counts
ggplot(LACO_ppt, aes(x = Jan_March_cm, y = log(LACOdens), col = Size))+
  geom_point()

#-> what are the effects of pool size and distance from runway on pool depth?
summary(lm(max_depth ~ Size, inundation_summary)) #no effect of size on max_depth
summary(lm(max_depth ~ Distance, inundation_summary)) #no effect of distance on max_depth

ggplot(inundation_summary, aes(x = as.factor(Year), y = max_depth, color = Size)) +
  geom_jitter()
ggplot(inundation_summary, aes(x = as.factor(Year), y = max_depth, color = Distance)) +
  geom_jitter()


#############################
# II. INUNDATION VS. LAMBDA #
#############################

#1.Extract lambda values from the model

# Run the model first
lambda_data <- as.data.frame(get_posterior_mean(BH_fit, pars = c("lambda")))

lambda_gather <- lambda_data %>%
  mutate(Year = c(2001:2017)) %>%
  select(Year, `mean-chain:1`, `mean-chain:2`, `mean-chain:3`, `mean-chain:4`) %>%
  gather(`mean-chain:1`, `mean-chain:2`, `mean-chain:3`, `mean-chain:4`, key = chain, value = lambda)
  
lambda_summary <- lambda_gather %>%
  group_by(Year) %>%
  summarize(mean_lambda = mean(lambda), se_lambda = se(lambda)) #this is used for plotting

#2.Regression of max water depth and precipitation to predict inundation of all years
summary(lm(max_depth ~ Jan_March_cm, env_summary)) #max water depth is significantly correlated with Jan - March precipitation
# Equation: max_depth = 5.09207 + 0.03201 * PPT_Jan_March_cm #R-squared is 0.007

ggplot() +
  geom_bar(aes(x = PPT$Year, y = PPT$Jan_March_cm), stat= "identity", fill = "lightblue") +
  geom_point(aes(x = LACO_env$Year, y = LACO_env$max_depth)) +
  scale_y_continuous(name = "Max water depth (cm)",
                     sec.axis = sec_axis(~./1, name = "Precip Jan-Mar (cm)")) +
  theme_bw()

pred_max_depth <- PPT %>%
  select(Year, Jan_March_cm) %>%
  mutate(pred_max_depth = 5.09207 + 0.03201 * Jan_March_cm) %>%
  filter(Year %in% c(2001:2017))

#3.Join lambda and max water depth by year
env_summary$Year <- as.integer(env_summary$Year)
lambda_real_depth <- left_join(env_summary, lambda_gather, by = "Year") %>%
  filter(Year != "2000")

pred_max_depth$Year <- as.integer(pred_max_depth$Year)
lambda_depth <- left_join(lambda_gather, pred_max_depth, by = "Year")

#4.Plot lambda and max water depth

ggplot(lambda_real_depth, aes(x = max_depth, y = lambda)) +
  geom_point()

ggplot(lambda_depth, aes(x = pred_max_depth, y = lambda)) +
  geom_jitter()+
  geom_smooth(method= "lm") +
  annotate(geom = "text", label = "R2 = 0.398", x = 5.75, y = 40)

#5.What's the fit?

summary(lm(pred_max_depth ~ lambda, lambda_depth)) 

#############################
# III. INUNDATION VS. ALPHA #
#############################
# INTRASPECIFIC COMPETITION
alphaLACO_data <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_LACO")))

alphaLACO_gather <- alphaLACO_data %>%
  mutate(Year = c(2001:2006)) %>%
  select(Year, `mean-chain:1`, `mean-chain:2`, `mean-chain:3`, `mean-chain:4`) %>%
  gather(`mean-chain:1`, `mean-chain:2`, `mean-chain:3`, `mean-chain:4`, key = chain, value = alphaLACO)

alphaLACO_summary <- alphaLACO_gather %>%
  group_by(Year) %>%
  summarize(mean_alphaLACO = mean(alphaLACO), se_alphaLACO = se(alphaLACO)) #this is used for plotting

alphaLACO_depth <- left_join(alphaLACO_gather, pred_max_depth, by = "Year")

ggplot(alphaLACO_depth, aes(x = pred_max_depth, y = alphaLACO)) +
  geom_jitter()+
  geom_smooth(method= "lm") +
  annotate(geom = "text", label = "R2 = 0.4773", x = 5.75, y = 0.4)

summary(lm(pred_max_depth ~ alphaLACO, alphaLACO_depth)) 

# INTERSPECIFIC COMPETITION
alphaEG_data <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_EG")))

alphaEG_gather <- alphaEG_data %>%
  mutate(Year = c(2001:2006)) %>%
  select(Year, `mean-chain:1`, `mean-chain:2`, `mean-chain:3`, `mean-chain:4`) %>%
  gather(`mean-chain:1`, `mean-chain:2`, `mean-chain:3`, `mean-chain:4`, key = chain, value = alphaEG)

alphaEG_summary <- alphaEG_gather %>%
  group_by(Year) %>%
  summarize(mean_alphaEG = mean(alphaEG), se_alphaEG = se(alphaEG)) #this is used for plotting

alphaEG_depth <- left_join(alphaEG_gather, pred_max_depth, by = "Year")

ggplot(alphaEG_depth, aes(x = pred_max_depth, y = alphaEG)) +
  geom_jitter()+
  geom_smooth(method= "lm") +
  annotate(geom = "text", label = "R2 = 0.546", x = 5.75, y = 0.4)

summary(lm(pred_max_depth ~ alphaEG, alphaEG_depth)) 

