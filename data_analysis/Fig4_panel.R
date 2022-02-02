#RUN exotic grass simulation.R
#RUN ecoapp_additional_analysis.R

library(tidyverse)
library(ggplot2)
library(ggpubr)

#reorder levels in Depth column
const_LACO_depth$Depth <- ordered(as.factor(const_LACO_depth$Depth), levels = c("shallow", "medium", "deep"))
GRWR_depth_all$Depth <- ordered(as.factor(GRWR_depth_all$Depth), levels = c("shallow", "medium", "deep"))

#FIGURE4
fig_LACO_depth <- ggplot(const_LACO_depth, aes(x = Year, y = mean_LACO))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 0.8) +
  geom_line(size=1, aes(linetype = Depth))+
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 18),
        legend.position = c(0.8,0.9),
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(name = "Year",
                     limits = c(1999.5,2015.5))+
  scale_y_log10()+
  ylab(bquote(Observed~LACO~Density~(stems/m^2)))+
  scale_linetype_manual(values=c("solid", "twodash", "dotted"),
                        name = "Pool depth",
                        labels = c("shallow (< 3.8 cm)", "medium (3.8-8 cm)", "deep (> 8 cm)"))

figbar_GRWR_depth <- ggplot(GRWR_depth_all , aes(x = Depth, y = Mean))+
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

sim_timeseries <- ggplot(summary_grass_sim_LACO%>%filter(time %in% c(2000:2015)), aes(x = as.numeric(time), y = mean_LACO)) +
  geom_point() +
  geom_line(size = 1, aes(linetype = type)) +
  geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 18),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.4, 0.2)) +
  scale_y_log10()+
  labs(y = bquote(Simulated~LACO~Density~(stems/m^2))) +
  scale_x_continuous(name = "Year", limits = c(1999.5,2015.5))+
  scale_linetype_manual(values=c("solid", "twodash", "dotted"),
                        name = "Reduction in Exotic Grasses",
                        labels = c("0%", "50%", "75%"))

sim_GRWR <- ggplot(GRWR_simulated_all , aes(x = treatment, y = Mean))+
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

ggarrange(fig_LACO_depth, figbar_GRWR_depth, sim_timeseries, sim_GRWR,  ncol = 2, nrow = 2, 
          labels = c("(a)", "(b)", "(c)", "(d)"), widths = c(0.6, 0.4),
          #common.legend = TRUE, legend = "bottom", 
          font.label = list(size = 20), align = "v")
