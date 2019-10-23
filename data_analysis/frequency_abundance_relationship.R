# Exploring the relationship between frequency and abundance
# This will allow reference data (frequency) to be transformed to abundance
# The abundance estimates will be used to parameterize a Beverton-Holt population model
# This addresses Q3 and Q4 in the sToration group

library(tidyverse)
library(gridExtra)

### Source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

### Get model for LACO
# Isolate LACO as species of interest
lasth <- const_com %>%
  select(Year, Pool, Complex, 
         Distance, Size, Treatment.1999,
         Treatment.2000, LACOdens, LACO) %>%
  filter(!is.na(LACOdens))

# Data visualization
p1 <- ggplot(lasth, aes(x = LACO, y = LACOdens)) +
  geom_point()  +
  annotate("text", x = 10, y = 2000,
           label = "Natural log") +
  scale_y_continuous(trans = "log") +
  xlab("LACO frequency") +
  ylab("LACO density")

p2 <- ggplot(lasth, aes(x = LACO, y = LACOdens)) +
  geom_point() +
  annotate("text", x = 10, y = 2000,
           label = "Log2") +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(trans = "log2") +
  xlab("LACO frequency") +
  ylab("LACO density")

p3 <- ggplot(lasth, aes(x = LACO, y = LACOdens)) +
  geom_point() +
  annotate("text", x = 10, y = 2000,
           label = "Log10") +
  scale_y_continuous(trans = "log10") +
  xlab("LACO frequency") +
  ylab("LACO density")

p4 <- ggplot(lasth, aes(x = LACO, y = LACOdens)) +
  geom_point() +
  annotate("text", x = 10, y = 2000,
           label = "Squareroot") +
  scale_y_continuous(trans = "sqrt") +
  xlab("LACO frequency") +
  ylab("LACO density")

grid.arrange(p1, p2, p3, p4, ncol = 2)

# Square root looks most linear; fit lm
m1 <- lm(sqrt(LACOdens) ~ LACO, data = lasth)
# Intercept is negative

# Try plotting
lasth$predict <- predict(m1)^2
lasth$se <- predict(m1, se.fit = TRUE)$se.fit^2

ggplot(lasth) +
  geom_point(aes(x = LACO, y = LACOdens)) +
  geom_ribbon(aes(x = LACO, ymin = predict - se,
                  ymax = predict + se)) +
  geom_line(aes(x = LACO, y = predict),
            color = "red")

# Log-log looks even better!
lasth_pres <- lasth %>%
  filter(LACO > 0 & LACOdens > 0)

m2 <- lm(log(LACOdens) ~ log(LACO), data = lasth_pres)

# Try plotting
lasth_pres$predict <- predict(m2)
lasth_pres$se <- predict(m2, se.fit = TRUE)$se.fit

ggplot(lasth_pres) +
  geom_point(aes(x = LACO, y = LACOdens)) +
  geom_ribbon(aes(x = LACO, ymin = exp(predict - se),
                  ymax = exp(predict + se))) +
  geom_line(aes(x = LACO, y = exp(predict)),
            color = "red") +
  scale_y_continuous(trans = "log") +
  scale_x_continuous(trans = "log")
  
# Dude, that's pretty awesome.
a <- as.numeric(coefficients(m2)[1])
b <- as.numeric(coefficients(m2)[2])

### Predicting abundance for reference data
ref_com$LACOdens <- exp(a) * ref_com$LACO ^ b
write.csv(ref_com, "data_analysis//Predicted_Ref_Com.csv")
