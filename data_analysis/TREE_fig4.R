library(tidyverse)
library(ggplot2)
#4 csv files located in sToration-vernal-pools > data_analysis
#"const_LACO_depth.csv"
#"GRWR_depth_all.csv"
#"const_EG_depth.csv"
#"GRWR_simulated_all.csv"

#LACO by depth
ggplot(const_LACO_depth, aes(x = Year, y = mean_LACO, col = Depth))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 0.8) +
  geom_line(size=1)+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.8,0.8),
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = NULL,
                     limits = c(1999.5,2015.5))+
  #scale_y_log10()+
  ylab(bquote(LACO~Density~(stems/m^2)))

#GRWR by depth
ggplot(GRWR_depth_all , aes(x = Depth, y = Mean))+
  geom_bar(stat = "identity", col = "grey27")+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 18),
        axis.line = element_line(colour = "black"),
        legend.position = "none")+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  labs(y = "Average Low Density Growth Rate", x = "Pool Depth")+
  geom_hline(yintercept = 0)

#EG cover by depth
ggplot(const_EG_depth%>%filter(Year < 2016), aes(x = Year, y = mean_EG, col = Depth))+
  geom_point()+
  geom_line()+
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.8,0.2))+
  geom_errorbar(aes(ymin = mean_EG-se_EG, ymax = mean_EG+se_EG), width = 0.4, alpha = 0.9, size = 0.8) + 
  ylab(bquote(Exotic~Grass~Cover~('%')))+
  scale_color_discrete(name = "Pool depth", labels = c("deep", "medium", "shallow"))

#GRWR by EG reduction
ggplot(GRWR_simulated_all , aes(x = treatment, y = Mean))+
  geom_bar(stat = "identity", col = "grey27")+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 18),
        axis.line = element_line(colour = "black"),
        legend.position = "none")+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  labs(y = "Average Low Density Growth Rate", x = "Reduction in Exotic Grasses")+
  geom_hline(yintercept = 0)
