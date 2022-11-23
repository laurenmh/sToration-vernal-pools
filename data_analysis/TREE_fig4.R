library(tidyverse)
library(ggplot2)
library(ggpubr)

###DATA

#csv files located in sToration-vernal-pools > data_analysis

#data for LACO density (C)
#"mean_join_LACO.csv"
mean_join_LACO <- read_csv("mean_join_LACO.csv")


#data for LACO GRWR (D)
#"Partitioning_GRWR_LACO_const.csv"
#"Partitioning_GRWR_LACO_ref.csv"
Partitioning_GRWR_LACO_const <- read_csv("Partitioning_GRWR_LACO_const.csv")
Partitioning_GRWR_LACO_ref <- read_csv("Partitioning_GRWR_LACO_ref.csv")


#data for LACO by depth (F)
#"const_LACO_depth.csv" 
const_LACO_depth <- read_csv("const_LACO_depth.csv")


#distribution of GRWR by depth for new (G) - rows = iterations, columns = years
#"GRWR_LACO_const_deep.csv"
#"GRWR_LACO_const_medium.csv"
#"GRWR_LACO_const_shallow.csv"

GRWR_LACO_const_deep <- read_csv("GRWR_LACO_const_deep.csv") %>%
  select(-1) %>%
  rowwise() %>%
  mutate(meandeep = mean(c(V1:V15))) 

GRWR_LACO_const_medium <- read_csv("GRWR_LACO_const_medium.csv") %>%
  select(-1) %>%
  rowwise() %>%
  mutate(meanmedium = mean(c(V1:V15)))

GRWR_LACO_const_shallow <- read_csv("GRWR_LACO_const_shallow.csv") %>%
  select(-1) %>%
  rowwise() %>%
  mutate(meanshallow = mean(c(V1:V15)))


depthdist <- as_tibble(cbind(GRWR_LACO_const_deep$meandeep, GRWR_LACO_const_medium$meanmedium, GRWR_LACO_const_shallow$meanshallow))
names(depthdist) = c("Deep", "Medium", "Shallow")

depthdist2 <- pivot_longer(depthdist, 1:3, "Depth")

g <- ggplot(depthdist2, aes(x=value, color = Depth, fill = Depth)) + 
  geom_density(alpha = .2, bw = 0.04) +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.2, 0.8), legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.background = element_blank(), legend.box.background = element_blank(),
        legend.key = element_blank()) + 
  labs(x = "Low-density growth rate", y = "Probability density") +
  geom_vline(xintercept = 0, alpha = 0.5, linetype = "dashed") 


#data for avg EG timeseries (I)
#"EG_all.csv"
EG_all <- read_csv("EG_all.csv")

# #distribution of GRWR by EG reduction for new (J) - rows = iterations, columns = years
# #"GRWR_LACO_const.csv" - 0% EG reduction
# GRWR_LACO_const <- read_csv("GRWR_LACO_const.csv") %>%
#   select(-1) %>%
#   rowwise() %>%
#   mutate(meanconst = mean(c(V1:V15)))
# 
# GRWR_LACO_const2 <- read_csv("GRWR_LACO_const.csv") %>%
#   select(-1) %>%
#   pivot_longer(1:15, "year")
#   
# mean(GRWR_LACO_const2$value)
# ggplot(GRWR_LACO_const2, aes(x=value)) + geom_density()

#"GRWR_LACO_50EG.csv" - 50% EG reduction
GRWR_LACO_50EG <- read_csv("GRWR_LACO_50EG.csv") %>%
  select(-1) %>%
  rowwise() %>%
  mutate(mean50 = mean(c(V1:V15)))

GRWR_LACO_50EG_mean <- read_csv("GRWR_LACO_50EG.csv") %>%
  magrittr::set_colnames(c(2001:2015)) %>%
  rownames_to_column(., var = "iteration") %>% #2000 iterations from Bayesian modeling
  pivot_longer(!iteration, names_to = "Year", values_to = "GRWR") %>%
  group_by(iteration) %>%
  summarize(mean_GRWR = mean(GRWR))  #mean of GRWR across years

#"GRWR_LACO_75EG.csv" - 75% EG reduction
GRWR_LACO_75EG <- read_csv("GRWR_LACO_75EG.csv") %>%
  select(-1) %>%
  rowwise() %>%
  mutate(mean75 = mean(c(V1:V15)))

grasssim <- as_tibble(cbind(GRWR_LACO_const$meanconst, GRWR_LACO_50EG$mean50, GRWR_LACO_75EG$mean75))
names(grasssim) = c("No reduction", "50% reduction", "75% reduction")

grasssim2 <- pivot_longer(grasssim, 1:3, "Scenario")

ggplot(grasssim2, aes(x=value, color = Scenario, fill = Scenario)) +
  geom_density(alpha = .2, bw = 0.04) +
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.2, 0.8), legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.background = element_blank(), legend.box.background = element_blank(),
        legend.key = element_blank()) + 
  labs(x = "Low-density growth rate", y = "Probability density") +
  geom_vline(xintercept = 0, alpha = 0.5, linetype = "dashed") 

  

## alternative 

GRWR_LACO_const2 <- read_csv("GRWR_LACO_const.csv") %>%
  select(-1) %>%
  summarize_all(.funs="mean") %>%
  pivot_longer(1:15, "year")

ggplot(GRWR_LACO_const2, aes(x=value)) +
  geom_density()


meanconst <- GRWR_LACO_const$meanconst

ggplot(GRWR_LACO_const, aes(x=meanconst)) + geom_density()

#"GRWR_LACO_50EG.csv" - 50% EG reduction
GRWR_LACO_50EG <- read_csv("GRWR_LACO_50EG.csv") %>%
  select(-1) %>%
  rowwise() %>%
  mutate(mean50 = mean(c(V1:V15)))

#"GRWR_LACO_75EG.csv" - 75% EG reduction
GRWR_LACO_75EG <- read_csv("GRWR_LACO_75EG.csv") %>%
  select(-1) %>%
  rowwise() %>%
  mutate(mean75 = mean(c(V1:V15)))

grasssim <- as_tibble(cbind(GRWR_LACO_const$meanconst, GRWR_LACO_50EG$mean50, GRWR_LACO_75EG$mean75))
names(grasssim) = c("No reduction", "50% reduction", "75% reduction")

grasssim2 <- pivot_longer(grasssim, 1:3, "Scenario")

ggplot(grasssim2, aes(x=value, color = Scenario)) +
  geom_density()

#extras
#"GRWR_depth_all.csv" - old (G) summarized GRWR by depth
#"const_EG_depth.csv" - old (I) EG by depth
#"GRWR_simulated_all.csv"- old (J) summarized GRWR by EG reduction

###PLOTS

#(C) LACO density
c <- ggplot(mean_join_LACO%>%filter(Year %in% c(2000:2015)), 
       aes(x = Year, y = mean_LACOdens, col = type)) +
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
  ylab(bquote(Density~(stems/m^2))) +
  xlab("Year") +
  scale_x_continuous(name = NULL, limits = c(1999.5,2015.5))+
  scale_color_manual(name = "", values = c("#000000", "#888888")) +
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
  annotate("text", x = 2, y = 1.1, label = "Constructed", size = 6)+
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
  annotate("text", x = 2, y = 1.1, label = "Reference", size = 6)+
  scale_x_discrete(labels = xlabels)
figure_partitioning <- ggarrange(Part_ref, Part_const, ncol = 2, nrow = 1, legend = "none", 
                                 labels = c("(a)", "(b)"), font.label = list(size = 20))
annotate_figure(figure_partitioning, bottom = text_grob("Mechanisms", size = 20),
                left = text_grob("Partitioning of Low Density Growth Rate", size = 18, rot = 90))
d <- figure_partitioning

#(F) LACO by depth
f <- ggplot(const_LACO_depth,aes(x = Year, y = mean_LACO, color = Depth)) + 
  geom_point() + geom_line() +
  labs(x = "Year", y =  bquote(Density~(stems/m^2))) +
  geom_errorbar(aes(ymin = mean_LACO - se_LACO, ymax = mean_LACO + se_LACO)) +
  scale_y_log10() + 
  theme_classic(base_size = 8) +
  theme(legend.position = c(0.8, 0.8), legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.background = element_blank(), legend.box.background = element_blank(),
        legend.key = element_blank()) 

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
