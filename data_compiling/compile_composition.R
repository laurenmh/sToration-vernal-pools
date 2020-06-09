## Make sure to manually set your data pathway in set_pathway.R first!!

library(tidyverse)

## Read in data from the constructed vernal pools
const_com <- read_csv(paste(datpath, "Aboveground annual vegetation data/dat_csv_constructed_pools.csv", sep="")) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`) %>%
  mutate(Treatment.1999 = fct_recode(Treatment.1999, "Group A" = "Group 1", "Group B" = "Group 2")) %>%
  mutate(Treatment.2000 = fct_recode(Treatment.2000, "Group A" = "Group 1", "Group B" = "Group 2"))

const_com_2000 <- read_csv(paste(datpath, "Aboveground annual vegetation data/2000 veg data CONSTRUCTED pools.csv", sep = "")) %>%
  rename(Treatment.1999 = `Treatment 1999`,
         Treatment.2000 = `Treatment 2000`) %>%
  mutate(Year = "2000")

const_com <-  rbind(const_com, const_com_2000)

## Read in data from the reference vernal pools
ref_com <- read_csv(paste(datpath, "Aboveground annual vegetation data/dat_csv_reference_pools.csv", sep="")) %>%
  filter(Year != "0")

ref_com_2000 <- read_csv(paste(datpath, "Aboveground annual vegetation data/2000 veg data REFERENCE pools.csv", sep = "")) %>%
  mutate(Year = "2000")

ref_com <- rbind(ref_com, ref_com_2000)

## Read in the species list and traits - constructed
const_spp <- read_csv(paste(datpath, "Aboveground annual vegetation data/Constructed pool metadata.csv", sep=""), skip = 29) %>%
  filter(!is.na(Genus))

## Read in the species list and traits - reference
ref_spp <- read_excel(paste(datpath, "Aboveground annual vegetation data/Reference pool metadata.xlsx", sep=""), skip = 26) %>%
  filter(!is.na(Genus)) %>%
  select(1:7)

