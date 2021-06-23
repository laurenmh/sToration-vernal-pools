# Script for FIGURE 1 and 2

# Set your data pathway first!
# Go to "data_wrangling/prep data before modeling.R" to get data on LACO abundance and exotic grass cover
PPT <- read_csv(paste(datpath, "Monthly precip averages/Fairfield_precip.csv", sep="")) #precip data
# Go to "analysis/complex_belowground_v5.R" for lambda values
# Go to "analysis/GRWR_invader" for GRWR values

library(ggplot2)
library(ggpubr)


# REFERENCE pool LACO cover data (see 'prep data before modeling.R')
ref_LACO <- ref_com_LACO %>% #9 reference pools that have consecutive data on LACO abundance from 2002 to 2015
  group_by(Pool) %>%
  gather(key = "Year", value = "LACO" , - Pool) 

ref_LACOden_edge <- ref_com %>% #any reference LACO abundance data in 2000, 2001, 2016, 2017
  select(Year, Pool, LACO) %>%
  group_by(Year, Pool) %>%
  filter(Year %in% c(2000, 2001, 2016, 2017)) %>% 
  summarise_each(funs(mean)) %>%
  mutate_each(funs(as.integer(.)))

ref_LACOdens <- full_join(ref_LACOden_edge, ref_LACO)%>% #join reference LACO abundance data 2000-2017
  mutate(LACOdens = round(exp(-0.42)+LACO^1.36)) %>% # Convert LACO frequency to density 
  mutate(type = "reference")
ref_LACOdens$Year <- as.numeric(ref_LACOdens$Year)

# CONSTRUCTED pool LACO density data (see 'prep data before modeling.R')
const_LACOden <- const_com_noNA %>% #72 pools
  select(- c(Size)) %>%
  group_by(Pool) %>%
  gather(key = "Year", value = "LACOdens", -Pool) %>%
  mutate(type = "constructed")
const_LACO <- left_join(const_LACOden, (const_com %>% select (Year, Pool, LACO, Treatment.1999, Treatment.2000)))

const_LACO$Year <- as.numeric(const_LACO$Year)
const_LACO$LACO <- as.numeric(const_LACO$LACO)

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

fabundance <- ggplot(mean_join_LACO, aes(x = Year, y = mean_LACOdens, col = type)) +
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
        ylab(bquote(~italic(L.~conj.)~Density~(stems/m^2))) +
        scale_x_continuous(name = NULL, limits = c(1999.5,2017.5))+
        scale_color_manual(name = "", values = c("#000000", "#888888"))

# Visualize timeseries of exotic grass cover (see 'prep data before modeling.R'))
ref_EG <- ref_com %>% 
  filter(Pool %in% c(9, 20, 34, 38, 52, 63, 77, 31, 27)) %>%
  select(-Quadrat) %>%
  group_by(Year, Pool) %>%
  summarise_each(funs(mean)) %>%
  group_by(Year, Pool) %>%
  summarise(sum_EG = sum(BRHO, HOMA, LOMU)) %>%
  mutate(type = "reference") %>%
  group_by(Year, type) %>%
  summarise(mean_EG = mean(sum_EG),
            se_EG = se(sum_EG))
const_EG <- const_dummy_join %>% 
  select(Year, Pool, sum_EG) %>% 
  mutate(type = "constructed") %>%
  group_by(Year, type) %>%
  summarise(mean_EG = mean(sum_EG),
            se_EG = se(sum_EG))
EG_all <- rbind(ref_EG, const_EG) 
EG_all$Year <- as.numeric(EG_all$Year)

fexoticgrass <- ggplot(EG_all, aes(x = Year, y = mean_EG, col = type)) +
                      geom_point() +
                      geom_line(size = 1) +
                      geom_errorbar(aes(ymin = mean_EG-se_EG, ymax = mean_EG+se_EG), width = 0.4, alpha = 0.9, size = 1) +
                      ylab(bquote(Exotic~Grass~Cover~('%'))) +
                      scale_x_continuous(name = NULL, limits = c(1999.5,2017.5))+
                      scale_color_manual(name = "", values = c("#000000", "#888888")) +
                      theme(text = element_text(size=16),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            legend.position = c(0.2, 0.9), 
                            axis.title = element_text(size = 14))

# Visualize timeseries of precipitation (PPT object)
PPT_long <- PPT %>% select(Year, Jan_March_cm, Oct_Dec_cm) %>% 
  pivot_longer(cols = c('Jan_March_cm', 'Oct_Dec_cm'),names_to = "season", values_to= "rain") #Oct_Dec_cm rain is calculated from t-1, while Jan_March_cm is calculated from t.
PPT_long$season <- as.factor(PPT_long$season)
frain <- ggplot(PPT_long %>%filter(Year  %in%  c(2000:2017)), aes(fill = season, y = rain, x = Year)) +
                geom_bar(position = "stack", stat = "identity") +
                ylab(bquote(Precipitation~(cm)))+
                scale_x_continuous(name = NULL, limits = c(1999.5,2017.5))+
                scale_fill_manual(name = "", labels = c("January-March", "October-December"), values = c("#000000", "#888888")) +
                theme(text = element_text(size=16),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position = c(0.7, 0.85), 
                      axis.title = element_text(size = 14))

# Visualize timeseries of predicted LACO lambda (see 'complex_belowground_v5.R')
lambda_ref_00_01 <- as.data.frame(reflambda_mean[1:2,5])%>% #2000-2001 lambda data 
  mutate(Year = c(2000:2001))%>%
  mutate(type = "reference")
lambda_ref_02_14 <- as.data.frame(reflambda_mean[,5]) %>% #2002-2014 lambda data
  mutate(Year = c(2002:2014))%>%
  mutate(type = "reference")
colnames(lambda_ref_00_01) <- c("lambda", "Year", "type")  
colnames(lambda_ref_02_14) <- c("lambda", "Year", "type")
lambda_ref <- full_join(lambda_ref_00_01, lambda_ref_02_14) #piece together 2000-2001 data with 2002-2014 data

lambda_const <- as.data.frame(lambda_mean[,5]) %>%
  mutate(Year = c(2000:2016))%>%
  mutate(type = "constructed")
colnames(lambda_const) <- c("lambda", "Year", "type")
lambda_const_ref <- rbind(lambda_ref, lambda_const) 

flambda <- ggplot(lambda_const_ref, aes(x = Year, y = lambda, col = type))+
          geom_point() +
          geom_line(size=1)+
          theme(text = element_text(size=16),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position = "none", 
                axis.title = element_text(size = 14))+
          ylab(bquote(Intrinsic~Growth~Rate~lambda))+
          scale_x_continuous(name = NULL,
                     limits = c(1999.5,2017.5))+
          scale_color_manual(name = "", values = c("#000000", "#888888"))

#Visualize timeseries of GRWR (see 'GRWR_invader.R')
GRWR_time <- GRWR_LACO %>%
  gather(key = "type", "GRWR", -Year)

fGRWR <- ggplot(GRWR_time, aes(x = Year, y = GRWR, col = type))+
            geom_point() +
            geom_line(size=1)+
            theme(text = element_text(size=16),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.position = "none",
                  axis.title = element_text(size = 14))+
            ylab(bquote(Growth~rate~when~rare))+
            geom_hline(yintercept = 0, linetype = "dashed")+
            scale_x_continuous(name = NULL,
                               limits = c(1999.5,2017.5))+
            scale_color_manual(name = "", values = c("#000000", "#888888"))

#FIGURE 1
Fig1 <- ggarrange(fabundance, fexoticgrass, frain,  ncol=1, nrow=3, align = "v", 
                  font.label = list(size = 14), hjust = 0.9)
#annotate_figure(Fig1, bottom = text_grob("Time (year)", size = 14))

#FIGURE 2
Fig2 <- ggarrange(flambda, fGRWR, ncol = 1, nrow = 2, align = "v", 
                  labels = c("(a)", "(b)"), font.label = list(size = 14))

# Visualize timeseries of observed LACO growth rate
#observed Nt/Nt-1 constructed vs. reference
#const_GR <- const_com_noNA %>% 
#  mutate_if(is.numeric, ~replace(., . == 0, 1)) %>% #to avoid NAN and infinite numbers, we'll replace all zeros with one.
#  group_by(Pool) %>%
#  summarize('2003G' = `2003`/`2002`,
#            '2004G' = `2004`/`2003`,
#            '2005G' = `2005`/`2004`,
#            '2006G' = `2006`/`2005`,
#            '2007G' = `2007`/`2006`,
#            '2008G' = `2008`/`2007`,
#            '2009G' = `2009`/`2008`,
#            '2010G' = `2010`/`2009`,
#            '2011G' = `2011`/`2010`,
#            '2012G' = `2012`/`2011`,
#            '2013G' = `2013`/`2012`,
#            '2014G' = `2014`/`2013`,
#            '2015G' = `2015`/`2014`) %>% 
#  gather(`2003G`, `2004G`, `2005G`, `2006G`, `2007G`, `2008G`, 
#         `2009G`, `2010G`, `2011G`, `2012G`, `2013G`, `2014G`, `2015G`, key = Year, value = const_GR) %>%
#  group_by(Year) %>%
#  summarize(mean_GR = mean(const_GR)) %>%
#  mutate(Year = gsub("G", "", Year)) %>%
#  mutate(type = "constructed")

# ref_GR <- ref_com_LACO %>%
#   mutate_if(is.numeric, ~replace(., . == 0, 1)) %>% #to avoid NAN and infinite numbers, we'll replace all zeros with one.
#   group_by(Pool) %>%
#   summarize('2003G' = `2003`/`2002`,
#             '2004G' = `2004`/`2003`,
#             '2005G' = `2005`/`2004`,
#             '2006G' = `2006`/`2005`,
#             '2007G' = `2007`/`2006`,
#             '2008G' = `2008`/`2007`,
#             '2009G' = `2009`/`2008`,
#             '2010G' = `2010`/`2009`,
#             '2011G' = `2011`/`2010`,
#             '2012G' = `2012`/`2011`,
#             '2013G' = `2013`/`2012`,
#             '2014G' = `2014`/`2013`,
#             '2015G' = `2015`/`2014`) %>%
#   gather(`2003G`, `2004G`, `2005G`, `2006G`, `2007G`, `2008G`, 
#          `2009G`, `2010G`, `2011G`, `2012G`, `2013G`, `2014G`, `2015G`, key = Year, value = ref_GR) %>%
#   group_by(Year) %>%
#   summarize(mean_GR = mean(ref_GR)) %>%
#   mutate(Year = gsub("G", "", Year)) %>%
#   mutate(type = "reference")

# join_GR <- rbind(const_GR, ref_GR) #join two tables
# 
# ggplot(join_GR, aes(x = as.numeric(Year), y = mean_GR, col = type))+
#   geom_point() +
#   geom_line()+
#   theme_classic()+
#   ylab("Observed LACO growth rate") +
#   scale_x_continuous(name = NULL,
#                      limits = c(2000,2015))


# Visualize predicted LACO growth rate
#predicted Nt/Nt-1 constructed vs. reference
# predicted_const_GR <- predicted_LACO %>%
#   spread(time, predicted_LACO) %>%
#   mutate_if(is.numeric, ~replace(., . == 0, 1)) %>%
#   group_by(Pool) %>%
#   summarize('2003G' = `2003`/`2002`,
#             '2004G' = `2004`/`2003`,
#             '2005G' = `2005`/`2004`,
#             '2006G' = `2006`/`2005`,
#             '2007G' = `2007`/`2006`,
#             '2008G' = `2008`/`2007`,
#             '2009G' = `2009`/`2008`,
#             '2010G' = `2010`/`2009`,
#             '2011G' = `2011`/`2010`,
#             '2012G' = `2012`/`2011`,
#             '2013G' = `2013`/`2012`,
#             '2014G' = `2014`/`2013`,
#             '2015G' = `2015`/`2014`) %>%
#   gather(`2003G`, `2004G`, `2005G`, `2006G`, `2007G`, `2008G`, 
#          `2009G`, `2010G`, `2011G`, `2012G`, `2013G`, `2014G`, `2015G`, key = Year, value = const_GR) %>%
#   group_by(Year) %>%
#   summarize(mean_GR = mean(const_GR)) %>%
#   mutate(Year = gsub("G", "", Year)) %>%
#   mutate(type = "constructed")
# 
# predicted_ref_GR <- ref_predicted_LACO %>%
#   spread(time, predicted_LACO) %>%
#   mutate_if(is.numeric, ~replace(., . == 0, 1)) %>%
#   group_by(Pool) %>%
#   summarize('2003G' = `2003`/`2002`,
#             '2004G' = `2004`/`2003`,
#             '2005G' = `2005`/`2004`,
#             '2006G' = `2006`/`2005`,
#             '2007G' = `2007`/`2006`,
#             '2008G' = `2008`/`2007`,
#             '2009G' = `2009`/`2008`,
#             '2010G' = `2010`/`2009`,
#             '2011G' = `2011`/`2010`,
#             '2012G' = `2012`/`2011`,
#             '2013G' = `2013`/`2012`,
#             '2014G' = `2014`/`2013`,
#             '2015G' = `2015`/`2014`) %>%
#   gather(`2003G`, `2004G`, `2005G`, `2006G`, `2007G`, `2008G`, 
#          `2009G`, `2010G`, `2011G`, `2012G`, `2013G`, `2014G`, `2015G`, key = Year, value = ref_GR) %>%
#   group_by(Year) %>%
#   summarize(mean_GR = mean(ref_GR)) %>%
#   mutate(Year = gsub("G", "", Year)) %>%
#   mutate(type = "reference")
# 
# join_predicted_GR <- rbind(predicted_const_GR, predicted_ref_GR) #join two tables
# 
# ggplot(join_predicted_GR, aes(x = as.numeric(Year), y = mean_GR, col = type))+
#   geom_point() +
#   geom_line()+
#   theme_classic()+
#   ylab("Predicted LACO growth rate") +
#   scale_x_continuous(name = NULL,
#                      limits = c(2000,2015))

#EXTRA PLOTS:

# Visualize LACOdens separately
#ref <- ggplot(ref_LACO, aes(x = Year, y = LACOdens)) +
#  geom_jitter(col = "blue") + 
#  geom_smooth(method = "loess") +
#  ylab("Reference pools LACO density") +
#  theme_bw()
#const <- ggplot(const_LACO, aes(x = Year, y = LACOdens)) +
#  geom_jitter(col = "pink") + 
#  geom_smooth(method = "loess", col = "red") +
#  ylab("Constructed pools LACO density") +
#  theme_bw()
#ggarrange(ref, const, ncol = 1, nrow = 2)

# Visualize LACOdens together
#ggplot(join_LACO, aes(x = Year, y = LACOdens, col = type)) +
#  geom_jitter() +
#  scale_y_log10()+
#  geom_smooth() +
#  theme_bw() +
#  ylab("LACO density (log)") 

#ggplot(join_LACO, aes(x = Year, y = log_LACOdens, col = treatment)) +
#  geom_jitter() +
#  geom_smooth(method = "loess", se = FALSE) +
#  theme_bw() +
#  ylab("LACO density (log)")

# Visualize LACO cover separately
#refcov <- ggplot(ref_LACO, aes(x = Year, y = LACO)) +
#  geom_jitter(col = "blue") +
#  geom_smooth(method = "loess") +
#  ylab("Reference pools LACO frequency") +
#  theme_bw()
#constcov <- ggplot(const_LACO, aes(x = Year, y = LACO)) +
#  geom_jitter(col = "pink") +
#  geom_smooth(method = "loess", col = "red") +
#  ylab("Constructed pool LACO frequency") +
#  theme_bw()
#ggarrange(refcov, constcov, ncol = 1, nrow  =2)

# Visualize LACO cover together
#ggplot(join_LACO, aes(x = Year, y = LACO, col = type)) + 
#  geom_jitter() +
#  geom_smooth() +
#  theme_bw() +
#  ylab("LACO frequency (%)")

#ggplot(join_LACO, aes(x = Year, y = LACO, col = treatment)) + 
#  geom_jitter() +
#  geom_smooth(se = FALSE) +
#  theme_bw() +
#  ylab("LACO frequency (%)")
