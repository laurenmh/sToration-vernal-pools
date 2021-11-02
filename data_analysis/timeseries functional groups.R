source("data_compiling/compile_composition.R")

library(ggplot2)
library(tidyverse)

#calculate standard error
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#timeseries of exotic grasses, native forbs, ERVA and LACO
FunctionalGroup_timeseries <- const_com %>%
  filter(Treatment.1999 != "Control") %>%
  #mutate(treatment = paste(Treatment.1999, Treatment.2000, sep = "-")) %>%
  filter(!is.na(BRHO), !is.na(HOMA), !is.na(LOMU), !is.na(PLST), !is.na(DOCO), !is.na(ERVA), !is.na(LACO)) %>%
  select(Year, BRHO, HOMA, LOMU, PLST, DOCO, ERVA, LACO) %>%
  pivot_longer(!Year, names_to = "Species", values_to = "frequency") %>%
  group_by(Year, Species) %>%
  summarize(mean=mean(frequency), se=calcSE(frequency)) %>%
  filter(Year < 2016)

FunctionalGroup_timeseries$Species <- ordered(as.factor(FunctionalGroup_timeseries$Species), levels = c("LACO", "BRHO","HOMA","LOMU", "ERVA",
                                                                                                        "PLST", "DOCO"))

ggplot(FunctionalGroup_timeseries, aes(x = as.numeric(Year), y = mean))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
  facet_wrap(~Species, ncol= 1)+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.1, 0.9), 
        axis.title = element_text(size = 14))+
  ylab("Mean Frequency (%)")+
  xlab("Year")

# What percentage of the total frequency is the six species we included in the model?
per_freq <- ref_com %>%
  mutate(total_freq = select(., ACMO:VUBR) %>% rowSums(na.rm = TRUE)) %>%
  mutate(frac_freq = BRHO+ERVA+LOMU+HOMA+PLST+DOCO) %>%
  mutate(ratio = frac_freq/total_freq)

sum(per_freq$frac_freq)/sum(per_freq$total_freq)

# #Native forb frequency timeseries
# NF_timeseries <- const_com %>%
#   filter(Treatment.1999 != "Control") %>%
#   #mutate(treatment = paste(Treatment.1999, Treatment.2000, sep = "-")) %>%
#   filter(!is.na(PLST), !is.na(DOCO)) %>%
#   select(Year, PLST, DOCO) %>%
#   pivot_longer(!Year, names_to = "Species", values_to = "frequency") %>%
#   group_by(Year, Species) %>%
#   summarize(mean=mean(frequency), se=calcSE(frequency)) 
# 
# ggplot(NF_timeseries, aes(x = as.numeric(Year), y = mean))+
#   geom_point()+
#   geom_line()+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
#   facet_wrap(~Species, ncol= 1)+
#   theme(text = element_text(size=16),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = c(0.1, 0.9), 
#         axis.title = element_text(size = 14))+
#   ylab("Mean Frequency (%)")+
#   xlab("Year")
# 
# #ERVA frequency timeseries
# ERVA_timeseries <- const_com %>%
#   filter(Treatment.1999 != "Control") %>%
#   #mutate(treatment = paste(Treatment.1999, Treatment.2000, sep = "-")) %>%
#   filter(!is.na(ERVA)) %>%
#   select(Year, ERVA) %>%
#   pivot_longer(!Year, names_to = "Species", values_to = "frequency") %>%
#   group_by(Year, Species) %>%
#   summarize(mean=mean(frequency), se=calcSE(frequency)) 
# 
# ggplot(ERVA_timeseries, aes(x = as.numeric(Year), y = mean))+
#   geom_point()+
#   geom_line()+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
#   #facet_wrap(~Species, ncol= 1)+
#   theme(text = element_text(size=16),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = c(0.1, 0.9), 
#         axis.title = element_text(size = 14))+
#   ylab("Mean Frequency (%)")+
#   xlab("Year")
# 
# #LACO frequency timeseries
# LACO_timeseries <- const_com %>%
#   filter(Treatment.1999 != "Control") %>%
#   #mutate(treatment = paste(Treatment.1999, Treatment.2000, sep = "-")) %>%
#   filter(!is.na(LACO)) %>%
#   select(Year, LACO) %>%
#   pivot_longer(!Year, names_to = "Species", values_to = "frequency") %>%
#   group_by(Year, Species) %>%
#   summarize(mean=mean(frequency), se=calcSE(frequency)) 
# 
# ggplot(LACO_timeseries, aes(x = as.numeric(Year), y = mean))+
#   geom_point()+
#   geom_line()+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
#   #facet_wrap(~Species, ncol= 1)+
#   theme(text = element_text(size=16),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = c(0.1, 0.9), 
#         axis.title = element_text(size = 14))+
#   ylab("Mean Frequency (%)")+
#   xlab("Year")
