Post_ref_00_15 <- rstan::extract(BH_ref_fit_00_15)
refalpha_LACO_00_15 <- as.matrix(Post_ref_00_15$alpha_LACO)
refalpha_EG_00_15 <- as.matrix(Post_ref_00_15$alpha_EG)
refalpha_ERVA_00_15 <- as.matrix(Post_ref_00_15$alpha_ERVA)
refalpha_NF_00_15 <- as.matrix(Post_ref_00_15$alpha_NF)
reflambda_00_15 <- as.matrix(Post_ref_00_15$lambda)
refs_00_15 <- as.matrix(Post_ref_00_15$survival_LACO)
sim_LACO <- matrix(nrow = 2000, ncol = 15)
#Plug in the stable equilibrium freq and parameters in the model.
LACO_ref_00_15 <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:16]),
  ERVA = as.numeric(const_com_control[1,2:16]),
  NF = as.numeric(const_com_control[3,2:16]),
  aii = refalpha_LACO_00_15,
  a1 = refalpha_EG_00_15,
  a2 = refalpha_ERVA_00_15,
  a3 = refalpha_NF_00_15,
  lambda = reflambda_00_15,
  s = refs_00_15,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_ref_00_15 <- log(LACO_ref_00_15) #2002-2014 #log transform

ggplot(GRWR_time_ref_00_15, aes(x = Year, y = mean))+
  geom_point() +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1)

CI_GRWR_ref_00_15 <-  as.data.frame(HDInterval::hdi(GRWR_LACO_ref_00_15, credMass = 0.95)) %>% #Calculate 95% credible interval of GRWR
  magrittr::set_colnames(c(2000:2014)) %>%
  mutate(CI = c("lowCI", "upCI")) %>%
  pivot_longer(!CI, names_to = "Year", values_to = "CI_values") %>%
  pivot_wider( names_from = CI, values_from = CI_values)
GRWR_time_ref_00_15 <- as.data.frame(GRWR_LACO_ref_00_15) %>%
  magrittr::set_colnames(c(2000:2014)) %>%
  pivot_longer(cols = everything()) %>%
  magrittr::set_colnames(c("Year", "GRWR")) %>%
  group_by(Year) %>%
  summarise(mean = mean(GRWR)) %>%
  mutate(type = "reference") %>%
  full_join(., CI_GRWR_ref_00_15)


GRWR_time_00_15 <- rbind(GRWR_time_ref_00_15, GRWR_time_const) %>% #combine lambda tables
  filter(Year < 2015) # cut the last two points in constructed pools to match reference pools
GRWR_time_00_15$Year <- as.numeric(GRWR_time_00_15$Year)

fGRWR <- ggplot(GRWR_time_00_15, aes(x = Year, y = mean, col = type))+
  geom_point() +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.1,0.2),
        axis.title = element_text(size = 14))+
  ylab(bquote(Low~Density~Growth~Rate~(italic(r[t]))))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(name = NULL,
                     limits = c(1999.5,2015.5))+
  scale_color_manual(name = "", values = c("#000000", "#888888"))

#FIGURE 2
Fig2 <- ggarrange(flambda, fGRWR, ncol = 1, nrow = 2, align = "v", 
                  labels = c("(a)", "(b)"), font.label = list(size = 14))
