# Make a timeseries figure of LACO density in REFERENCE and CONSTRUCTED pools from 2000-2015.

# Load data and packages
# Set your data pathway first!
# Go to "data_wrangling/prep data before modeling.R"
# Run "const_com_noNA" in option 1 and "ref_com_LACO"
library(ggplot2)
library(ggpubr)

# Convert LACO frequency to density in REFERENCE pools
ref_LACOden <- ref_com_LACO %>%
  group_by(Pool) %>%
  gather(key = "Year", value = "LACOcover" , - Pool) %>%
  mutate(refLACOdens = round(exp(-0.42)+LACOcover^1.36))
ref_LACOden$Year <- as.numeric(ref_LACOden$Year)

# Truncate CONSTRUCTED pool data to 2000-2015
const_LACOden <- const_com_noNA %>%
  select(- c(Size, `2016`, `2017`)) %>%
  group_by(Pool) %>%
  gather(key = "Year", value = "constLACOdens", -Pool) 
const_LACOden$Pool <- as.character(const_LACOden$Pool)
const_LACOden$Year <- as.numeric(const_LACOden$Year)

# Join tables
join_LACOden <- full_join(ref_LACOden, const_LACOden, all = TRUE)

# Visualize separately
ref <- ggplot(ref_LACOden, aes(x = Year, y = refLACOdens)) +
  geom_jitter(col = "blue") + 
  geom_smooth(method = "loess") +
  ylab("Reference pools LACO density") +
  theme_bw()
const <- ggplot(const_LACOden, aes(x = Year, y = constLACOdens)) +
  geom_jitter(col = "pink") + 
  geom_smooth(method = "loess", col = "red") +
  ylab("Constructed pools LACO density") +
  theme_bw()
ggarrange(ref, const, ncol = 1, nrow = 2)

# Visualize together
ggplot(join_LACOden, aes(x = Year)) +
  geom_jitter(aes(y = log(refLACOdens)), col = "lightblue") +
  geom_jitter(aes(y = log(constLACOdens)), col = "pink") + 
  geom_smooth(aes(y = log(refLACOdens)), method = "loess", na.rm = TRUE) +
  geom_smooth(aes(y = log(constLACOdens)), method = "loess", na.rm = TRUE, col = "red") +
  theme_bw()
