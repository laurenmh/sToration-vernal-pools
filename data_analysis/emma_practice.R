library(tidyverse)

## Read in data from the constructed vernal pools
construct_com <- read_csv(paste(datpath, "Aboveground annual vegetation data/dat_csv_constructed_pools.csv", sep="")) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`) 

## Read in data from the reference vernal pools
ref_com <- read_csv(paste(datpath, "Aboveground annual vegetation data/dat_csv_reference_pools.csv", sep="")) 


View(ref_com)

colnames(ref_com)

ref_com_long<-gather(ref_com,sp,pres,ACMO:VUBR)

View(ref_com_long)
