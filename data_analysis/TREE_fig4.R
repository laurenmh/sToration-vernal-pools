library(tidyverse)
library(ggplot2)
library(ggpubr)

###DATA

#csv files located in sToration-vernal-pools > data_analysis

#data for LACO density (C)
#"mean_join_LACO.csv"

#data for LACO GRWR (D)
#"Partitioning_GRWR_LACO_const.csv"
#"Partitioning_GRWR_LACO_ref.csv"

#data for LACO by depth (F)
#"const_LACO_depth.csv" 

#distribution of GRWR by depth for new (G) - rows = iterations, columns = years
#"GRWR_LACO_const_deep.csv"
#"GRWR_LACO_const_medium.csv"
#"GRWR_LACO_const_shallow.csv"

#data for avg EG timeseries (I)
#"EG_all.csv"

#distribution of GRWR by EG reduction for new (J) - rows = iterations, columns = years
#"GRWR_LACO_const.csv" - 0% EG reduction
#"GRWR_LACO_50EG.csv" - 50% EG reduction
#"GRWR_LACO_75EG.csv" - 75% EG reduction

#extras
#"GRWR_depth_all.csv" - old (G) summarized GRWR by depth
#"const_EG_depth.csv" - old (I) EG by depth
#"GRWR_simulated_all.csv"- old (J) summarized GRWR by EG reduction

###PLOTS

#(C) LACO density
ggplot(mean_join_LACO%>%filter(Year %in% c(2000:2015)), aes(x = Year, y = mean_LACOdens, col = type)) +
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
  ylab(bquote(Density~(stems/m^2))) +
  scale_x_continuous(name = NULL, limits = c(1999.5,2015.5))+
  scale_color_manual(name = "", values = c("#000000", "#888888"))

#(D) GRWR partitioning
#FIGURE 3
Partitioning_GRWR_LACO_const$mechanism <- ordered(Partitioning_GRWR_LACO_const$mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
xlabels <- c("r_overall" = expression(bar("r")[i]^" "), 
             "epsilon_0" = expression(epsilon[i]^0), 
             "epsilon_alpha" = expression(epsilon[i]^alpha),
             "epsilon_lambda" = expression(epsilon[i]^lambda),
             "epsilon_int" = expression(epsilon[i]^{alpha*lambda}))
Part_const <- ggplot(Partitioning_GRWR_LACO_const, aes(x = mechanism, y = Mean, fill = mechanism))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 0.8) +
  theme_classic(base_size = 20)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
  ylim(-1.2, 1.2)+
  annotate("text", x = 2, y = 1.1, label = "Constructed", size = 6)+
  scale_x_discrete(labels = xlabels)

Partitioning_GRWR_LACO_ref$mechanism <- ordered(Partitioning_GRWR_LACO_ref$mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Part_ref <- ggplot(Partitioning_GRWR_LACO_ref, aes(x = mechanism, y = Mean, fill = mechanism))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 0.8) +
  theme_classic(base_size = 20)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
  ylim(-1.2, 1.2)+
  annotate("text", x = 2, y = 1.1, label = "Reference", size = 6)+
  scale_x_discrete(labels = xlabels)
figure_partitioning <- ggarrange(Part_ref, Part_const, ncol = 2, nrow = 1, legend = "none", 
                                 labels = c("(a)", "(b)"), font.label = list(size = 20))
annotate_figure(figure_partitioning, bottom = text_grob("Mechanisms", size = 20),
                left = text_grob("Partitioning of Low Density Growth Rate", size = 18, rot = 90))

#(F) LACO by depth
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

#old (G) GRWR by depth
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

#old (I) EG cover by depth
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

#old (J) GRWR by EG reduction
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
