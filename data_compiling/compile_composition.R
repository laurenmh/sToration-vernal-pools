## Make sure to manually set your data pathway in set_pathway.R first!!

library(tidyverse)

## Read in data from the constructed vernal pools
const_com <- read_csv(paste(datpath, "Aboveground annual vegetation data/dat_csv_constructed_pools.csv", sep="")) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`) %>%
  mutate(Treatment.1999 = fct_recode(Treatment.1999, "Group A" = "Group 1", "Group B" = "Group 2")) %>%
  mutate(Treatment.2000 = fct_recode(Treatment.2000, "Group A" = "Group 1", "Group B" = "Group 2"))


## Read in data from the reference vernal pools
ref_com <- read_csv(paste(datpath, "Aboveground annual vegetation data/dat_csv_reference_pools.csv", sep="")) 
## Read in the species list and traits - constructed
const_spp <- read_csv(paste(datpath, "Aboveground annual vegetation data/Constructed pool metadata.csv", sep=""), skip = 29) %>%
  filter(!is.na(Genus))

## Read in the species list and traits - reference
ref_spp <- read_excel(paste(datpath, "Aboveground annual vegetation data/Reference pool metadata.xlsx", sep=""), skip = 26) %>%
  filter(!is.na(Genus)) %>%
  select(1:7)

