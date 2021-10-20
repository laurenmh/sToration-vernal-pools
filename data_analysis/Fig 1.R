# Script for FIGURE 1

# Set your data pathway first!

# Data
source("data_wrangling/prep data before modeling.R") #Get data on LACO abundance and exotic grass cover
PPT <- read_csv(paste(datpath, "Monthly precip averages/Fairfield_precip.csv", sep="")) #Get precip data

# Packages
library(ggplot2)
library(ggpubr)

# REFERENCE pool LACO cover data
ref_LACO <- ref_com_LACO %>% #9 reference pools that have consecutive data on LACO abundance from 2002 to 2015
  group_by(Pool) %>%
  gather(key = "Year", value = "LACO" , - Pool) 
ref_LACO$Pool <- as.character(ref_LACO$Pool)

ref_LACOden_edge <- ref_com %>% #any reference LACO abundance data in 2000, 2001, 2016, 2017
  select(Year, Pool, LACO) %>%
  group_by(Year, Pool) %>%
  filter(Year %in% c(2000, 2001, 2016, 2017)) %>% 
  summarise_each(funs(mean))
ref_LACOden_edge$LACO <- as.integer(ref_LACOden_edge$LACO)

ref_LACOdens <- full_join(ref_LACOden_edge, ref_LACO)%>% #join reference LACO abundance data 2000-2017
  mutate(LACOdens = round(exp(-0.42)+LACO^1.36)) %>% # Convert LACO frequency to density 
  mutate(type = "reference")
ref_LACOdens$Year <- as.numeric(ref_LACOdens$Year)

# CONSTRUCTED pool LACO density data 
const_LACOden <- const_com_noNA %>% #72 pools
  select(- c(Size)) %>%
  group_by(Pool) %>%
  gather(key = "Year", value = "LACOdens", -Pool) %>%
  mutate(type = "constructed")
const_LACO <- left_join(const_LACOden, (const_com %>% select (Year, Pool, LACO, Treatment.1999, Treatment.2000)))

const_LACO$Year <- as.numeric(const_LACO$Year)
const_LACO$LACO <- as.integer(const_LACO$LACO)
const_LACO$Pool <- as.character(const_LACO$Pool)

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
  geom_line(size=1.5)+
  scale_y_log10()+
  geom_errorbar(aes(ymin = mean_LACOdens-se_LACOdens, ymax = mean_LACOdens+se_LACOdens), width = 0.4, alpha = 0.9, size = 1) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(.2, .3),
        axis.title = element_text(size = 14))+
  ylab(bquote(~italic(L.conj.)~Density~(stems/m^2))) +
  scale_x_continuous(name = NULL, limits = c(1999.5,2017.5))+
  scale_color_manual(name = "", values = c("#000000", "#888888"))

# Visualize timeseries of exotic grass cover
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
  geom_line(size = 1.5) +
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
  scale_fill_manual(name = "", labels = c("January-March", "October-December"), values = c("#888888", "#000000")) +
  theme(text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.7, 0.85), 
        axis.title = element_text(size = 14)) +
  guides(fill = guide_legend(reverse = TRUE))

#FIGURE 1
SuppFig_EG_PPT <- ggarrange(fexoticgrass, frain,  ncol=1, nrow=2, align = "v", 
                  font.label = list(size = 14),  labels = c("(a)", "(b)"))
#annotate_figure(Fig1, bottom = text_grob("Time (year)", size = 14))
