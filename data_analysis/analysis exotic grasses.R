# Analyze the relationship between germination of LACO and exotic grasses
# According to field observation from Akasha, LACO germination is very low when there is high litter cover.
# We can incorporate this dynamic in the model by changing the germination rate depending on previous year's exotic grass cover.


source("data_wrangling/prep data before modeling.R")
library(ggplot2)
library(gridExtra)

#EG group - LOMU, BRHO, HOMA
ggplot(const_dummy_join, aes(x = sum_EG, y = log(LACOdens))) +
  geom_point()+ 
  geom_smooth(se= FALSE, method = "gam", formula = y ~ poly(x, 2))

plot1 <- ggplot(const_dummy_join, aes(x = as.factor(Year))) +
            geom_boxplot(aes(y = LACOdens)) + 
            ylim(0, 120)
plot2 <- ggplot(const_dummy_join, aes(x = as.factor(Year))) +
            geom_boxplot(aes(y = sum_EG)) +
            ylim(0, 120)
grid.arrange(plot1, plot2, nrow=2, ncol=1)

EG <- lm(LACOdens ~ poly(sum_EG, 2, raw = TRUE), const_dummy_join)
summary(EG)

