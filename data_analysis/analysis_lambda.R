# Run "prep data before modeling.R" first.
# Run "complex_belowground_v5.R" and "reference_pool_model.R" with real data.
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

summary(lm(constructed ~ reference, mean_LACOdens))
ggplot(mean_LACOdens, aes(y = constructed, x = reference)) +
  geom_point() +
  ylim(0, 500) +
  xlim(0, 250) +
  labs(y = "Constructed mean LACO density",
       x = "Reference mean LACO density") +
  geom_smooth(method = "lm") +
  annotate("text", label = "R^2 = 0.168", x = 100, y = 400) +
  theme_bw()

#############################################
# II. Plot lambda constructed vs. reference #
#############################################
# combine lambda estimates from constructed and reference models
const_lambda_trim <- as.data.frame(lambda_mean[-c(1,2,16,17),5]) #trim 2001, 2002, 2016, and 2017 estimates
ref_lambda_trim <- as.data.frame(reflambda_mean[,5]) 
join_lambda <- cbind(const_lambda_trim, ref_lambda_trim)
row.names(join_lambda) <- c(2003:2015)
colnames(join_lambda) <- c("constructed", "reference")

summary(lm(constructed ~ reference, join_lambda))
ggplot(join_lambda, aes(y = constructed, x = reference)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Constructed LACO lambda", x = "Reference LACO lambda") +
  annotate("text", label = "R^2 = 0.1045", x = 20, y = 50) +
  theme_bw()
