# Run "prep data before modeling.R" then "complex_belowground_v5.R" first.
# Run "timeseries visuals.R"

# Question:
## Does lambda provide a better regression than abundance measures?

# Outline:
## I. Plot LACO dens constructed vs. reference
## II. Plot Lambda constructed vs. reference

# Set up:
library(ggplot2)

##############################################
# I. Plot LACO dens constructed vs. reference#
##############################################
spread_join_LACO <- spread(join_LACO, type, LACOdens) #join_LACO from "timeseries visuals.R"
mean_LACOdens <- spread_join_LACO %>%
  select(Year, constructed, reference) %>%
  filter(Year != c("2000", "2001")) %>%
  group_by(Year) %>%
  summarise_at(c("constructed", "reference"), mean, na.rm = TRUE)

summary(lm(reference ~ constructed, mean_LACOdens))
ggplot(mean_LACOdens, aes(x = constructed, y = reference)) +
  geom_point() +
  ylim(0, 300) +
  xlim(0, 500) +
  labs(x = "Constructed mean LACO density",
       y = "Reference mean LACO density") +
  geom_smooth(method = "lm") +
  annotate("text", label = "R^2 = 0.168", x = 100, y = 250) +
  theme_bw()

#############################################
# II. Plot lambda constructed vs. reference #
#############################################


