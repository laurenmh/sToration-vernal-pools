## Make sure to manually set your data pathway in set_pathway.R first!!

library(tidyverse)

## Read in data from the constructed vernal pools
construct_center_depth_2000 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2000 Center depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `14-Jan`:`30-Mar`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2000, Location = "PoolCenter")

construct_seed_depth_2000 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2000 Seed plot depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `14-Jan`:`30-Mar`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2000, Location = "SeedPlot")

construct_center_depth_2002 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2002 Center depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `20-Dec`:`11-Apr`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2002, Location = "PoolCenter")

construct_seed_depth_2002 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2002 Seed plot depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `20-Dec`:`11-Apr`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration(days)`) %>%
  mutate(Year = 2002, Location = "SeedPlot")

construct_center_depth_2009 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2009 Center depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `26-Jan`:`19-Mar`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2009, Location = "PoolCenter")


construct_seed_depth_2009 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2009 Seed plot depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `26-Jan`:`19-Mar`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2009, Location = "SeedPlot")

construct_seed_depth_2010 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2010 Seed plot depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `24-Jan`:`16-Apr`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2010, Location = "SeedPlot")

construct_seed_depth_2011 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2011 Seed plot depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `21-Dec`:`15-Apr`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2011, Location = "SeedPlot") 

construct_seed_depth_2012 <- read_csv(paste(datpath, "Environmental weekly water depth data/Env water depth data constructed pools/2012 Seed plot depth data CONSTRUCTED pools.csv", sep="")) %>%
  gather(Date, Depth, `11-Jan`:`22-Apr`) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`,
         Duration.weeks = `Duration (weeks)`,
         Duration.days = `Duration (days)`) %>%
  mutate(Year = 2012, Location = "SeedPlot") 

# Compile all the depth data
const_depth <- rbind(construct_center_depth_2000, construct_seed_depth_2000, construct_center_depth_2002, construct_seed_depth_2002,
                           construct_center_depth_2009, construct_seed_depth_2009, construct_seed_depth_2010, construct_seed_depth_2011,
                           construct_seed_depth_2012) %>%
  mutate(Size = tolower(Size))

# Distill down to unique duration data
const_duration <- constructed_depth %>%
  select(Pool:Duration.days, Year, Location) %>%
  unique()

rm(construct_center_depth_2000, construct_seed_depth_2000, construct_center_depth_2002, construct_seed_depth_2002,
   construct_center_depth_2009, construct_seed_depth_2009, construct_seed_depth_2010, construct_seed_depth_2011,
   construct_seed_depth_2012)
