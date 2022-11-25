library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)



# Depth

const_LACO_depth <- read_csv("const_LACO_depth.csv") %>%
  mutate(Depth = recode(Depth, deep = "Deep", medium = "Medium", shallow = "Shallow")) %>%
  mutate(`Pool depth` = factor(Depth, levels = c("Deep", "Medium", "Shallow" )))



f <- ggplot(const_LACO_depth, aes(x = Year, y = mean_LACO, color = `Pool depth`)) + 
  geom_point() + geom_line() +
  labs(x = "Year", y =  bquote(`Lasthenia density`~(stems/m^2))) +
  geom_errorbar(aes(ymin = mean_LACO - se_LACO, ymax = mean_LACO + se_LACO)) +
  scale_y_log10() + 
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.1, 0.2), legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.background = element_blank(), 
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_color_manual(values = c("#023858", "#3690c0", "#d0d1e6"))

# Depth simulations


deep <- read_csv("GRWR_LACO_const_deep_mean.csv") %>%
  mutate(Depth = "Deep")
medium <- read_csv("GRWR_LACO_const_medium_mean.csv") %>%
  mutate(Depth = "Medium")
shallow <- read_csv("GRWR_LACO_const_shallow_mean.csv") %>%
  mutate(Depth = "Shallow")

depthtog <- rbind(deep, medium, shallow) %>%
  mutate(`Pool depth` = factor(Depth, levels = c("Deep", "Medium", "Shallow" )))


g <- ggplot(depthtog, aes(x=mean_GRWR, color = `Pool depth`, fill=`Pool depth`)) + 
  geom_density(alpha = .2, bw = 0.04) +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.2, 0.8), legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.background = element_blank(), 
        legend.box.background = element_blank(),
        legend.key = element_blank()) + 
  xlim(-2, 1.2) +
  labs(x = "Lasthenia low-density growth rate", y = "Probability density") +
  geom_vline(xintercept = 0, alpha = 0.5, linetype = "dashed") +
  geom_point(data = . %>%
               group_by(`Pool depth`) %>%
               summarize(mean = mean(mean_GRWR)) %>%
               mutate(y = -0.2),
             aes(x = mean, y = y, colour = `Pool depth`), size = 2) +
  scale_color_manual(values = c("#023858", "#3690c0", "#d0d1e6")) +
  scale_fill_manual(values = c("#023858", "#3690c0", "#d0d1e6"))
  

# invasive grass patterns 
const_EG_depth <- read_csv("const_EG_depth.csv") %>%
  mutate(Depth = recode(Depth, deep = "Deep", medium = "Medium", shallow = "Shallow")) %>%
  mutate(`Pool depth` = factor(Depth, levels = c("Deep", "Medium", "Shallow" )))



i <- ggplot(const_EG_depth, aes(x = Year, y = mean_EG, color = `Pool depth`)) + 
  geom_point() + geom_line() +
  labs(x = "Year", y =  "Annual grass density (% cover)") +
  geom_errorbar(aes(ymin = mean_EG - se_EG, ymax = mean_EG + se_EG)) +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.7, 0.15), legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.background = element_blank(), legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_color_manual(values = c("#023858", "#3690c0", "#d0d1e6")) 
  


## INVasive grass simulations

GRWR_LACO_50EG_mean <- read_csv("GRWR_LACO_50EG_mean.csv") %>%
  mutate(Scenario = "50% grass removal")
GRWR_LACO_75EG_mean <- read_csv("GRWR_LACO_75EG_mean.csv") %>%
  mutate(Scenario = "75% grass removal")
GRWR_LACO_const_mean <- read_csv("GRWR_LACO_const_mean.csv") %>%
  mutate(Scenario = "No grass removal")

simtog <- rbind(GRWR_LACO_50EG_mean, GRWR_LACO_75EG_mean, GRWR_LACO_const_mean) %>%
  mutate(Scenario = factor(Scenario, levels = c("No grass removal", "50% grass removal", "75% grass removal" )))


j <- ggplot(simtog, aes(x=mean_GRWR, color = Scenario, fill=Scenario)) + 
  geom_density(alpha = .2, bw = 0.04) +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.35, 0.8), legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.background = element_blank(), 
        legend.box.background = element_blank(),
        legend.key = element_blank()) + 
  xlim(-.6, .6) +
  labs(x = "Lasthenia low-density growth rate", y = "Probability density") +
  geom_vline(xintercept = 0, alpha = 0.5, linetype = "dashed") + 
  geom_point(data = . %>%
               group_by(Scenario) %>%
               summarize(mean = mean(mean_GRWR)) %>%
               mutate(y = -0.2),
             aes(x = mean, y = y, colour = Scenario), size = 2) +
  scale_color_manual(values = c("#662506", "#ec7014", "#fee391")) +
  scale_fill_manual(values = c("#662506", "#ec7014", "#fee391"))

plot_grid(f,g,i,j)

ggsave(filename = "fig4a.pdf", width = 9, height = 6)

# PUT THEM TOGETHER


mean_join_LACO <- read_csv("mean_join_LACO.csv") %>%
  mutate(type = recode(type, constructed = "Constructed", reference = "Reference"),
         `Pool type` = factor(type, levels = c("Reference", "Constructed")))
#(C) LACO density
c <- ggplot(mean_join_LACO%>%filter(Year %in% c(2000:2015)), 
       aes(x = Year, y = mean_LACOdens, col = `Pool type`)) +
  geom_point()+
  geom_line()+
  scale_y_log10()+
  geom_errorbar(aes(ymin = mean_LACOdens-se_LACOdens, ymax = mean_LACOdens+se_LACOdens), 
                width = 0.4, alpha = 0.9) +
  # theme(text = element_text(size=16),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       panel.background = element_blank(),
  #       axis.line = element_line(colour = "black"),
  #       legend.position = c(.2, .3),
  #       axis.title = element_text(size = 14))+
  labs(y = bquote(`Lasthenia density`~(stems/m^2)), x = "Year") +
  scale_x_continuous(limits = c(1999.5,2015.5))+
  scale_color_manual( values = c("#000000", "#888888")) +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.3, 0.3), legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.background = element_blank(), legend.box.background = element_blank(),
        legend.key = element_blank()) 

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
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9) +
  theme_classic(base_size = 8) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
  ylim(-1.2, 1.2)+
  annotate("text", x = 2, y = 1.1, label = "Constructed", size = 2)+
  scale_x_discrete(labels = xlabels)

Partitioning_GRWR_LACO_ref$mechanism <- ordered(Partitioning_GRWR_LACO_ref$mechanism, levels = c("r_overall", "epsilon_0", "epsilon_alpha", "epsilon_lambda", "epsilon_int"))
Part_ref <- ggplot(Partitioning_GRWR_LACO_ref, aes(x = mechanism, y = Mean, fill = mechanism))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9) +
  theme_classic(base_size = 8) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c("grey27", "grey60", "grey60", "grey60", "grey60"))+
  ylim(-1.2, 1.2)+
  annotate("text", x = 2, y = 1.1, label = "Reference", size = 2)+
  scale_x_discrete(labels = xlabels)

figure_partitioning <- ggarrange(Part_ref, Part_const, ncol = 2, nrow = 1, legend = "none")

d <- annotate_figure(figure_partitioning, bottom = text_grob("Mechanisms", size = 8),
                left = text_grob("Partitioning of low-density growth rate", size = 8, rot = 90))
ggsave(filename = "fig4a2.pdf", width = 3, height = 3)

plot_grid(c,d)

ggsave(filename = "fig4b.pdf", width = 9, height = 3)

c
ggsave(filename = "fig4a1.pdf", width = 3, height = 3)




# 
# #(F) LACO by depth
# f <- ggplot(const_LACO_depth,aes(x = Year, y = mean_LACO, color = Depth)) + 
#   geom_point() + geom_line() +
#   labs(x = "Year", y =  bquote(Density~(stems/m^2))) +
#   geom_errorbar(aes(ymin = mean_LACO - se_LACO, ymax = mean_LACO + se_LACO)) +
#   scale_y_log10() + 
#   theme_classic(base_size = 8) +
#   theme(legend.position = c(0.8, 0.8), legend.margin = margin(0, 0, 0, 0, "cm"),
#         legend.background = element_blank(), legend.box.background = element_blank(),
#         legend.key = element_blank()) 

# ggplot(const_LACO_depth, aes(x = Year, y = mean_LACO, col = Depth))+
#   geom_point()+
#   geom_errorbar(aes(ymin = mean_LACO-se_LACO, ymax = mean_LACO+se_LACO), width = 0.4, alpha = 0.9, size = 0.8) +
#   geom_line(size=1)+
#   theme(text = element_text(size=16),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         legend.position = c(0.8,0.8),
#         axis.line = element_line(colour = "black"))+
#   scale_x_continuous(name = NULL,
#                      limits = c(1999.5,2015.5))+
#   #scale_y_log10()+
#   ylab(bquote(LACO~Density~(stems/m^2)))

# #old (G) GRWR by depth


# ggplot(GRWR_depth_all , aes(x = Depth, y = Mean))+
#   geom_bar(stat = "identity", col = "grey27")+
#   theme(text = element_text(size=18),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_text(size = 18),
#         axis.line = element_line(colour = "black"),
#         legend.position = "none")+
#   geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
#   labs(y = "Average Low Density Growth Rate", x = "Pool Depth")+
#   geom_hline(yintercept = 0)

# #old (I) EG cover by depth
# ggplot(const_EG_depth%>%filter(Year < 2016), aes(x = Year, y = mean_EG, col = Depth))+
#   geom_point()+
#   geom_line()+
#   theme(text = element_text(size=16),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position = c(0.8,0.2))+
#   geom_errorbar(aes(ymin = mean_EG-se_EG, ymax = mean_EG+se_EG), width = 0.4, alpha = 0.9, size = 0.8) + 
#   ylab(bquote(Exotic~Grass~Cover~('%')))+
#   scale_color_discrete(name = "Pool depth", labels = c("deep", "medium", "shallow"))

# #old (J) GRWR by EG reduction
# ggplot(GRWR_simulated_all , aes(x = treatment, y = Mean))+
#   geom_bar(stat = "identity", col = "grey27")+
#   theme(text = element_text(size=18),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_text(size = 18),
#         axis.line = element_line(colour = "black"),
#         legend.position = "none")+
#   geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
#   labs(y = "Average Low Density Growth Rate", x = "Reduction in Exotic Grasses")+
#   geom_hline(yintercept = 0)
