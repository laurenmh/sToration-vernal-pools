# Make a timeseries figure of LACO density in REFERENCE and CONSTRUCTED pools from 2000-2015.

# Load data and packages
# Set your data pathway first!
# Go to "data_wrangling/prep data before modeling.R"
# Run "const_com_noNA" in option 1 and "ref_com_LACO"
library(ggplot2)
library(ggpubr)


# REFERENCE pool LACO cover data
ref_LACO <- ref_com_LACO %>% #9 pools that have consecutive data from 2002 to 2015
  group_by(Pool) %>%
  gather(key = "Year", value = "LACO" , - Pool) 

ref_LACOden_edge <- ref_com %>%
  select(Year, Pool, LACO) %>%
  group_by(Year, Pool) %>%
  filter(Year %in% c(2000, 2001, 2016, 2017)) %>% #pull out any reference data in 2000, 2001, 2016, 2017
  summarise_each(funs(mean)) %>%
  mutate_each(funs(as.integer(.)))

ref_LACOdens <- full_join(ref_LACOden_edge, ref_LACO)%>%
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

f1a <- ggplot(mean_join_LACO, aes(x = Year, y = mean_LACOdens, col = type)) +
        geom_point()+
        geom_line(size=1)+
        scale_y_log10()+
        geom_errorbar(aes(ymin = mean_LACOdens-se_LACOdens, ymax = mean_LACOdens+se_LACOdens), width = 0.4, alpha = 0.9, size = 1) +
        theme(text = element_text(size=16),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"))+
        ylab(bquote('Count (individuals/ '~ m^2~')')) +
        scale_x_continuous(name = NULL, limits = c(1999.5,2017.5))+
        scale_color_manual(name = "", values = c("#000000", "#888888"))

# Visualize timeseries of predicted LACO lambda
lambda_ref <- as.data.frame(reflambda_mean[,5]) %>%
  mutate(Year = c(2002:2014))%>%
  mutate(type = "reference")
colnames(lambda_ref) <- c("lambda", "Year", "type")  
lambda_const <- as.data.frame(lambda_mean[,5]) %>%
  mutate(Year = c(2000:2016))%>%
  mutate(type = "constructed")
colnames(lambda_const) <- c("lambda", "Year", "type")
lambda_const_ref <- rbind(lambda_ref, lambda_const) 

f1b <- ggplot(lambda_const_ref, aes(x = Year, y = lambda, col = type))+
          geom_point() +
          geom_line(size=1)+
          theme(text = element_text(size=16),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"))+
          ylab("Intrinsic growth rate")+
          scale_x_continuous(name = "Time (year)",
                     limits = c(1999.5,2017.5))+
          scale_color_manual(name = "", values = c("#000000", "#888888"))

#FIGURE 1
ggarrange(f1a, f1b, ncol=1, nrow=2, common.legend = TRUE, legend = "bottom", align = "v", labels = c("a)", "b)"))

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
