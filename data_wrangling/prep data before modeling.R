# This file preps the data to fit it in the model
# Outline of 3 options:
##1a. subset complete data 2000-2017 - we'll use this to plot timeseries of real observed data
##1b. subset complete data 2000-2006 - we'll use this to parameterize the model
##1c. fill in the missing data - we may also use this to plot timeseries of real observed data
##2. count the number of pools 
##3. count the number of years
##4. sum frequency data
##5. spread the dataset
##6. create a matrix of seeds added each year

# Remember to set your data pathway!

# Load data and package
source("data_compiling/compile_composition.R")
# View(const_com)
library(tidyverse)

#############################################
#How many pools have complete LACOdens data?#
#############################################
const_com_LACO <- const_com %>%
  filter(Treatment.1999 != "Control") %>% #remove control plots
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) #191 pools with seeding treatment #lots of missing data in 2007
const_com_noNA <- const_com_LACO[complete.cases(const_com_LACO),] #only 72 pools have complete data

#1a. subset the complete data
const_com_subset <- inner_join(const_com, const_com_noNA, by.y = "Pool")

#2a. count the number of pools
n_pools <- length(unique(const_com_subset$Pool))

#3a. count the number of years
n_years <- length(unique(const_com_subset$Year))

#4a. sum frequency data
#We'll use observed abundance data for LACO and ERVA.
#We'll sum the percent cover of species in each group 

#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
const_com_subset$sum_EG <- rowSums(cbind(const_com_subset$BRHO, const_com_subset$HOMA, const_com_subset$LOMU))

#Native forb (NF) group contains PLST and DOCO:
const_com_subset$sum_NF <- rowSums(cbind(const_com_subset$PLST, const_com_subset$DOCO))

#5a. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
LACOdens <- const_com_subset %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens <- const_com_subset %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGdens <- const_com_subset %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFdens <- const_com_subset %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

#6a. create a matrix of seeds added each year
seedtrt <- const_com_subset %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(const_com_subset$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

########################################## 
#Subset complete LACO data from 2000-2006#
##########################################
#1b. subset complete LACO data 2000-2006
const_com_LACO_short <- const_com_LACO %>%
  select(Pool, `2000`, `2001`, `2002`, `2003`, `2004`, `2005`, `2006`) 
const_com_noNA_short <- const_com_LACO_short[complete.cases(const_com_LACO_short),] #152 pools 

const_com_short <- inner_join(const_com, const_com_noNA_short, by.y = "Pool") %>%
  filter(Year %in% c(2000, 2001, 2002, 2003, 2004, 2005, 2006))

#2b. count the number of pools
n_pools <- length(unique(const_com_short$Pool))

#3b. count the number of years
n_years <- length(unique(const_com_short$Year))

#4b. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
const_com_short$sum_EG <- rowSums(cbind(const_com_short$BRHO, const_com_short$HOMA, const_com_short$LOMU))

#Native forb (NF) group contains PLST and DOCO:
const_com_short$sum_NF <- rowSums(cbind(const_com_short$PLST, const_com_short$DOCO))

#5b. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
LACOdens <- const_com_short %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens <- const_com_short %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGdens <- const_com_short %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFdens <- const_com_short %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

#6b. create a matrix of seeds added each year
seedtrt <- const_com_short %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(const_com_short$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

############################################################################## 
#Rreplace missing values with the mean value of pre and post NA#
##############################################################################

#remove rows that have more than two consecutive NAs 
const_com_revised <- const_com_LACO[-c(2,18,32,96,101,113,142,156,158,159,176,178,190),]

#fill in NA in 2007 with zeros if the population remains zero for the rest of the timeseries
const_com_revised$`2007` <- ifelse(const_com_revised$Pool %in% c(319,320,324,348,361,367,369,482,525,527,542,545), 0, const_com_revised$`2007`)
  
#1c.fill in the rest of NAs with the mean value of before and after NAs
#We'll use 2001 data for 2000 and 2016 data for 2017, otherwise take the average of before and after
const_com_revised$`2000` <- ifelse(is.na(const_com_revised$`2000`), const_com_revised$`2001`, const_com_revised$`2000`)
const_com_revised$`2002` <- ifelse(is.na(const_com_revised$`2002`), as.integer((const_com_revised$`2001`+const_com_revised$`2003`)/2), const_com_revised$`2002`)
const_com_revised$`2004` <- ifelse(is.na(const_com_revised$`2004`), as.integer((const_com_revised$`2003`+const_com_revised$`2005`)/2), const_com_revised$`2004`)
const_com_revised$`2007` <- ifelse(is.na(const_com_revised$`2007`), as.integer((const_com_revised$`2006`+const_com_revised$`2008`)/2), const_com_revised$`2007`)
const_com_revised$`2010` <- ifelse(is.na(const_com_revised$`2010`), as.integer((const_com_revised$`2009`+const_com_revised$`2011`)/2), const_com_revised$`2010`)
const_com_revised$`2013` <- ifelse(is.na(const_com_revised$`2013`), as.integer((const_com_revised$`2012`+const_com_revised$`2014`)/2), const_com_revised$`2013`)
const_com_revised$`2016` <- ifelse(is.na(const_com_revised$`2016`), as.integer((const_com_revised$`2015`+const_com_revised$`2017`)/2), const_com_revised$`2016`)
const_com_revised$`2017` <- ifelse(is.na(const_com_revised$`2017`), const_com_revised$`2016`, const_com_revised$`2017`)

#rearrange the table
const_com_revised <- const_com_revised %>%
  gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`,`2007`,`2008`,`2009`,`2010`,`2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = "Year", value = "LACOdens")
const_com_nafilled <- left_join(const_com_revised, const_com_completeveg, by = c("Pool", "Year"))

#2c. count the number of pools
n_pools <- length(unique(const_com_nafilled$Pool))

#3c. count the number of years
n_years <- length(unique(const_com_nafilled$Year))

#4c. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
const_com_nafilled$sum_EG <- rowSums(cbind(const_com_nafilled$BRHO, const_com_nafilled$HOMA, const_com_nafilled$LOMU))

#Native forb (NF) group contains PLST and DOCO:
const_com_nafilled$sum_NF <- rowSums(cbind(const_com_nafilled$PLST, const_com_nafilled$DOCO))

#5c. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
LACOdens <- const_com_nafilled %>%
  select(Year, Pool, LACOdens.x) %>%
  spread(key = Year, value = LACOdens.x) %>%
  select(-Pool) #no NAs
ERVAdens <- const_com_nafilled %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
ERVAdens$`2000` <- ifelse(is.na(ERVAdens$`2000`), ERVAdens$`2001`, ERVAdens$`2000`)
ERVAdens$`2002` <- ifelse(is.na(ERVAdens$`2002`), floor((ERVAdens$`2001`+ERVAdens$`2003`)/2), ERVAdens$`2002`)
ERVAdens$`2004` <- ifelse(is.na(ERVAdens$`2004`), floor((ERVAdens$`2003`+ERVAdens$`2005`)/2), ERVAdens$`2004`)
ERVAdens$`2007` <- ifelse(is.na(ERVAdens$`2007`), floor((ERVAdens$`2006`+ERVAdens$`2008`)/2), ERVAdens$`2007`)
ERVAdens$`2010` <- ifelse(is.na(ERVAdens$`2010`), floor((ERVAdens$`2009`+ERVAdens$`2011`)/2), ERVAdens$`2010`)
ERVAdens$`2013` <- ifelse(is.na(ERVAdens$`2013`), floor((ERVAdens$`2012`+ERVAdens$`2014`)/2), ERVAdens$`2013`)
ERVAdens$`2016` <- ifelse(is.na(ERVAdens$`2016`), floor((ERVAdens$`2015`+ERVAdens$`2017`)/2), ERVAdens$`2016`)
ERVAdens$`2017` <- ifelse(is.na(ERVAdens$`2017`), ERVAdens$`2016`, ERVAdens$`2017`)
sumEGdens <- const_com_nafilled %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumEGdens$`2000` <- ifelse(is.na(sumEGdens$`2000`), sumEGdens$`2001`, sumEGdens$`2000`)
sumEGdens$`2002` <- ifelse(is.na(sumEGdens$`2002`), floor((sumEGdens$`2001`+sumEGdens$`2003`)/2), sumEGdens$`2002`)
sumEGdens$`2004` <- ifelse(is.na(sumEGdens$`2004`), floor((sumEGdens$`2003`+sumEGdens$`2005`)/2), sumEGdens$`2004`)
sumEGdens$`2007` <- ifelse(is.na(sumEGdens$`2007`), floor((sumEGdens$`2006`+sumEGdens$`2008`)/2), sumEGdens$`2007`)
sumEGdens$`2010` <- ifelse(is.na(sumEGdens$`2010`), floor((sumEGdens$`2009`+sumEGdens$`2011`)/2), sumEGdens$`2010`)
sumEGdens$`2013` <- ifelse(is.na(sumEGdens$`2013`), floor((sumEGdens$`2012`+sumEGdens$`2014`)/2), sumEGdens$`2013`)
sumEGdens$`2016` <- ifelse(is.na(sumEGdens$`2016`), floor((sumEGdens$`2015`+sumEGdens$`2017`)/2), sumEGdens$`2016`)
sumEGdens$`2017` <- ifelse(is.na(sumEGdens$`2017`), sumEGdens$`2016`, sumEGdens$`2017`)
sumNFdens <- const_com_nafilled %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)
sumNFdens$`2000` <- ifelse(is.na(sumNFdens$`2000`), sumNFdens$`2001`, sumNFdens$`2000`)
sumNFdens$`2002` <- ifelse(is.na(sumNFdens$`2002`), floor((sumNFdens$`2001`+sumNFdens$`2003`)/2), sumNFdens$`2002`)
sumNFdens$`2004` <- ifelse(is.na(sumNFdens$`2004`), floor((sumNFdens$`2003`+sumNFdens$`2005`)/2), sumNFdens$`2004`)
sumNFdens$`2007` <- ifelse(is.na(sumNFdens$`2007`), floor((sumNFdens$`2006`+sumNFdens$`2008`)/2), sumNFdens$`2007`)
sumNFdens$`2010` <- ifelse(is.na(sumNFdens$`2010`), floor((sumNFdens$`2009`+sumNFdens$`2011`)/2), sumNFdens$`2010`)
sumNFdens$`2013` <- ifelse(is.na(sumNFdens$`2013`), floor((sumNFdens$`2012`+sumNFdens$`2014`)/2), sumNFdens$`2013`)
sumNFdens$`2016` <- ifelse(is.na(sumNFdens$`2016`), floor((sumNFdens$`2015`+sumNFdens$`2017`)/2), sumNFdens$`2016`)
sumNFdens$`2017` <- ifelse(is.na(sumNFdens$`2017`), sumNFdens$`2016`, sumNFdens$`2017`)

#6c. create a matrix of seeds added each year
seedtrt <- const_com_nafilled %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  na.omit(seedtrt$Treatment.1999) %>%
  unique(seedtrt$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 100, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0)) 

