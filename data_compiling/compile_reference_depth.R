## Make sure to manually set your data pathway in set_pathway.R first!!

library(tidyverse)

## Read in data from the reference vernal pools
ref_center_depth_2000 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data reference pools/2000 Center depth data REFERENCE pools.csv", sep="")) %>%
  gather(Date, Depth, `14-Jan`:`30-Mar`) %>%
  rename(Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2000, Location = "PoolCenter")


ref_center_depth_2002 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data reference pools/2002 Center depth data REFERENCE pools.csv", sep="")) %>%
  gather(Date, Depth, `20-Dec`:`9-Apr`) %>%
  rename(Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2002, Location = "PoolCenter")

ref_seed_depth_2002 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data reference pools/2002 Plot depth data REFERENCE pools.csv", sep="")) %>%
  gather(Date, Depth, `20-Dec`:`9-Apr`) %>%
  rename(Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2002, Location = "SeedPlot")

ref_center_depth_2009 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data reference pools/2009 Center depth data REFERENCE pools.csv", sep="")) %>%
  gather(Date, Depth, `19-Feb`:`19-Mar`) %>%
  rename(Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2009, Location = "PoolCenter")


ref_seed_depth_2009 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data reference pools/2009 Plot depth data REFERENCE pools.csv", sep="")) %>%
  gather(Date, Depth, `19-Feb`:`19-Mar`) %>%
  rename(Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2009, Location = "SeedPlot")


ref_seed_depth_2010 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data reference pools/2010 Plot depth data REFERENCE pools.csv", sep="")) %>%
  gather(Date, Depth, `24-Jan`:`16-Apr`) %>%
  rename(Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2010, Location = "SeedPlot")

ref_seed_depth_2011 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data reference pools/2011 Plot depth data REFERENCE pools.csv", sep="")) %>%
  gather(Date, Depth, `21-Dec`:`15-Apr`) %>%
  rename(Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2011, Location = "SeedPlot") 

ref_seed_depth_2012 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data reference pools/2012 Plot depth data REFERENCE pools.csv", sep="")) %>%
  select(-X17) %>%
  gather(Date, Depth, `11-Jan`:`22-Apr`) %>%
  rename(Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2012, Location = "SeedPlot") 

# Combine all the depth data
ref_depth <- rbind(ref_center_depth_2000 %>% mutate(Quadrat = NA), ref_center_depth_2002 %>% mutate(Quadrat = NA), ref_seed_depth_2002, ref_center_depth_2009 %>% mutate(Quadrat = NA),
                   ref_seed_depth_2009,ref_seed_depth_2010, ref_seed_depth_2011, ref_seed_depth_2012)            
  
# Distill down to unique duration data
ref_duration <- ref_depth %>%
  select(Pool:Duration.days, Year, Location, Quadrat) %>%
  unique()

rm(ref_center_depth_2000, ref_center_depth_2002 , ref_seed_depth_2002, ref_center_depth_2009,
   ref_seed_depth_2009,ref_seed_depth_2010, ref_seed_depth_2011, ref_seed_depth_2012)
