# Make a timeseries figure of LACO density in REFERENCE and CONSTRUCTED pools from 2000-2015.

# Load data and packages
# Set your data pathway first!
# Go to "data_wrangling/prep data before modeling.R"
# Run "const_com_noNA" in option 1 and "ref_com_LACO"
library(ggplot2)
library(ggpubr)

# Convert LACO frequency to density in REFERENCE pools
ref_LACO <- ref_com_LACO %>% #7pools
  group_by(Pool) %>%
  gather(key = "Year", value = "LACO" , - Pool) %>%
  mutate(LACOdens = round(exp(-0.42)+LACO^1.36)) %>%
  mutate(type = "reference")
ref_LACO$Year <- as.numeric(ref_LACO$Year)

# Truncate CONSTRUCTED pool LACO data to 2000-2015
const_LACOden <- const_com_noNA %>% #72 pools
  select(- c(Size, `2016`, `2017`)) %>%
  group_by(Pool) %>%
  gather(key = "Year", value = "LACOdens", -Pool) %>%
  mutate(type = "constructed")
const_LACO <- left_join(const_LACOden, (const_com %>% select (Year, Pool, LACO)))
const_LACO$Pool <- as.character(const_LACO$Pool)
const_LACO$Year <- as.numeric(const_LACO$Year)
const_LACO$LACO <- as.numeric(const_LACO$LACO)

# Join tables
join_LACO <- full_join(ref_LACO, const_LACO, all = TRUE)

# Visualize LACOdens separately
ref <- ggplot(ref_LACO, aes(x = Year, y = LACOdens)) +
  geom_jitter(col = "blue") + 
  geom_smooth(method = "loess") +
  ylab("Reference pools LACO density") +
  theme_bw()
const <- ggplot(const_LACO, aes(x = Year, y = LACOdens)) +
  geom_jitter(col = "pink") + 
  geom_smooth(method = "loess", col = "red") +
  ylab("Constructed pools LACO density") +
  theme_bw()
ggarrange(ref, const, ncol = 1, nrow = 2)

# Visualize LACOdens together
ggplot(join_LACO, aes(x = Year, y = LACOdens, col = type)) +
  geom_jitter() +
  geom_smooth() +
  theme_bw() +
  ylab("LACO density")

# Visualize LACO cover seperately
refcov <- ggplot(ref_LACO, aes(x = Year, y = LACO)) +
  geom_jitter(col = "blue") +
  geom_smooth(method = "loess") +
  ylab("Reference pools LACO frequency") +
  theme_bw()
constcov <- ggplot(const_LACO, aes(x = Year, y = LACO)) +
  geom_jitter(col = "pink") +
  geom_smooth(method = "loess", col = "red") +
  ylab("Constructed pool LACO frequency") +
  theme_bw()
ggarrange(refcov, constcov, ncol = 1, nrow  =2)

# Visualize LACO cover together
ggplot(join_LACO, aes(x = Year, y = LACO, col = type)) + 
  geom_jitter() +
  geom_smooth() +
  theme_bw() +
  ylab("LACO frequency (%)")