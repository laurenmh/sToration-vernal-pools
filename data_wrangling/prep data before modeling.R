# This file preps the veg data to fit it in the model

# REFERENCE POOLS
# subset complete data 2000-2015
# 1. subset the data
# 2. count the number of pools
# 3. count the number of years
# 4. sum frequency data
# 5. spread the table

# CONSTRUCTED POOLS
# Option a. subset complete data 2000-2017 - use this to plot timeseries
# Option b. subset complete data 2000-2006 
# Option c. fill in the missing data 
# Option d. fill in 2007 data and subset complete data 2000-2017 - use this to parameterize
# Option e. follow option d and subset by seeding treatment
# Option f. follow option d and subset by pool size
# Option g. follow option d and subset by average pool depth
# For each option, I am doing the following:
# 1. subset the data
# 2. count the number of pools 
# 3. count the number of years
# 4. sum frequency data
# 5. spread the table
# 6. create a matrix of seeds added each year

#------------------------------
# Load data and package
# Remember to set your data pathway first!
source("data_compiling/compile_composition.R") 
# View(const_com)
# View(ref_com)
library(tidyverse)

#----------------------
# REFERENCE POOLS:
#Nine pools with complete LACO data 2002-2015 
#1. Take the mean of sub-sampling plots and organize by year, pool, and LACO
ref_com_mean <- ref_com %>%
  select(-Quadrat) %>%
  group_by(Year, Pool) %>%
  filter(Year %in% c(2002:2015)) %>%
  summarise_each(funs(mean)) %>%
  mutate_each(funs(as.integer(.)))

ref_com_LACO <- ref_com_mean %>%
  select(Year, Pool, LACO) %>%
  spread(key = Year, value = LACO) %>%  
  filter(Pool %in% c(9, 20, 34, 38, 52, 63, 77, 31, 27)) %>%
  mutate_each(funs(as.integer(.)))

ref_com_join <- inner_join(ref_com_mean, ref_com_LACO, by.y = "Pool") %>%
  filter(Year %in% c(2002:2015))

#2. count the number of pools
ref_n_pools <- length(unique(ref_com_join$Pool))

#3. count the number of years
ref_n_years <- length(unique(ref_com_join$Year))

#4. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
ref_com_join$sum_EG <- as.integer(rowSums(cbind(ref_com_join$BRHO, ref_com_join$HOMA, ref_com_join$LOMU)))

#Native forb (NF) group contains PLST and DOCO:
ref_com_join$sum_NF <- as.integer(rowSums(cbind(ref_com_join$PLST, ref_com_join$DOCO)))

#5b. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
ref_LACOcover <- ref_com_join %>%
  select(Year, Pool, LACO) %>%
  spread(key = Year, value = LACO) %>%
  select(-Pool)
ref_ERVAcover <-  ref_com_join  %>%
  select(Year, Pool, ERVA) %>%
  spread(key = Year, value = ERVA) %>%
  select(-Pool)
ref_sumEGcover <-  ref_com_join  %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
ref_sumNFcover <-  ref_com_join  %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

# #Seven pools with complete data 2000-2015.
ref_com_mean_00_15 <- ref_com %>%
  select(-Quadrat) %>%
  group_by(Year, Pool) %>%
  filter(Year %in% c(2000:2015)) %>%
  summarise_each(funs(mean)) %>%
  mutate_each(funs(as.integer(.)))

ref_com_LACO_00_15 <- ref_com_mean_00_15 %>%
  select(Year, Pool, LACO) %>%
  spread(key = Year, value = LACO) %>%
  filter(Pool %in% c(9, 20, 34, 38, 52, 63, 77)) %>%
  mutate_each(funs(as.integer(.)))

ref_com_join_00_15 <- inner_join(ref_com_mean_00_15, ref_com_LACO_00_15, by.y = "Pool") %>%
  filter(Year %in% c(2000:2015))

#2. count the number of pools
ref_n_pools_00_15 <- length(unique(ref_com_join_00_15$Pool))

#3. count the number of years
ref_n_years_00_15 <- length(unique(ref_com_join_00_15$Year))

#4. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
ref_com_join_00_15$sum_EG <- as.integer(rowSums(cbind(ref_com_join_00_15$BRHO, ref_com_join_00_15$HOMA, ref_com_join_00_15$LOMU)))

#Native forb (NF) group contains PLST and DOCO:
ref_com_join_00_15$sum_NF <- as.integer(rowSums(cbind(ref_com_join_00_15$PLST, ref_com_join_00_15$DOCO)))

#5b. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
ref_LACOcover_00_15 <- ref_com_join_00_15 %>%
  select(Year, Pool, LACO) %>%
  spread(key = Year, value = LACO) %>%
  select(-Pool)
ref_ERVAcover_00_15 <-  ref_com_join_00_15  %>%
  select(Year, Pool, ERVA) %>%
  spread(key = Year, value = ERVA) %>%
  select(-Pool)
ref_sumEGcover_00_15 <-  ref_com_join_00_15  %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
ref_sumNFcover_00_15 <-  ref_com_join_00_15  %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

#------------------------------
# CONSTRUCTED POOLS:
# Remove control plots and organize by year, pool, and LACOdens
const_com_LACO <- const_com %>%
  filter(Treatment.1999 != "Control") %>% #remove control plots
  select(Year, Pool, LACOdens, Size) %>%
  spread(key = Year, value = LACOdens) #202 pools with seeding treatment #lots of missing data in 2007

#######################################################
#Option a. How many pools have complete LACOdens data?#
#######################################################
const_com_noNA <- const_com_LACO[complete.cases(const_com_LACO),] #only 72 pools have complete data

#1a. subset the complete data
#const_com_subset <- inner_join(const_com, const_com_noNA, by.y = "Pool")

#2a. count the number of pools
#n_pools <- length(unique(const_com_subset$Pool))

#3a. count the number of years
#n_years <- length(unique(const_com_subset$Year))

#4a. sum frequency data
#We'll use observed abundance data for LACO and ERVA.
#We'll sum the percent cover of species in each group 

#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
#const_com_subset$sum_EG <- rowSums(cbind(const_com_subset$BRHO, const_com_subset$HOMA, const_com_subset$LOMU))

#Native forb (NF) group contains PLST and DOCO:
#const_com_subset$sum_NF <- rowSums(cbind(const_com_subset$PLST, const_com_subset$DOCO))

#5a. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
#LACOdens <- const_com_subset %>%
#  select(Year, Pool, LACOdens) %>%
#  spread(key = Year, value = LACOdens) %>%
#  select(-Pool)
#ERVAdens <- const_com_subset %>%
#  select(Year, Pool, ERVAdens) %>%
#  spread(key = Year, value = ERVAdens) %>%
#  select(-Pool)
#sumEGcover <- const_com_subset %>%
#  select(Year, Pool, sum_EG) %>%
#  spread(key = Year, value = sum_EG) %>%
#  select(-Pool)
#sumNFcover <- const_com_subset %>%
#  select(Year, Pool, sum_NF) %>%
#  spread(key = Year, value = sum_NF) %>%
#  select(-Pool)

#6a. create a matrix of seeds added each year
#seedtrt <- const_com_subset %>%
#  select(Pool, Treatment.1999, Treatment.2000) %>%
#  unique(const_com_subset$Pool, incomparables = FALSE) %>%
#  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
#  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
#  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

####################################################
#Option b. Subset complete LACO data from 2000-2006#
####################################################

#1b. subset complete LACO data 2000-2006
#const_com_LACO_short <- const_com_LACO %>%
#  select(Pool, `2000`, `2001`, `2002`, `2003`, `2004`, `2005`, `2006`) 
#const_com_noNA_short <- const_com_LACO_short[complete.cases(const_com_LACO_short),] #152 pools 

#const_com_short <- inner_join(const_com, const_com_noNA_short, by.y = "Pool") %>%
#  filter(Year %in% c(2000, 2001, 2002, 2003, 2004, 2005, 2006))

#2b. count the number of pools
#n_pools <- length(unique(const_com_short$Pool))

#3b. count the number of years
#n_years <- length(unique(const_com_short$Year))

#4b. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
#const_com_short$sum_EG <- rowSums(cbind(const_com_short$BRHO, const_com_short$HOMA, const_com_short$LOMU))

#Native forb (NF) group contains PLST and DOCO:
#const_com_short$sum_NF <- rowSums(cbind(const_com_short$PLST, const_com_short$DOCO))

#5b. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
#LACOdens <- const_com_short %>%
#  select(Year, Pool, LACOdens) %>%
#  spread(key = Year, value = LACOdens) %>%
#  select(-Pool)
#ERVAdens <- const_com_short %>%
#  select(Year, Pool, ERVAdens) %>%
#  spread(key = Year, value = ERVAdens) %>%
#  select(-Pool)
#sumEGcover <- const_com_short %>%
#  select(Year, Pool, sum_EG) %>%
#  spread(key = Year, value = sum_EG) %>%
#  select(-Pool)
#sumNFcover <- const_com_short %>%
#  select(Year, Pool, sum_NF) %>%
#  spread(key = Year, value = sum_NF) %>%
#  select(-Pool)

#6b. create a matrix of seeds added each year
#seedtrt <- const_com_short %>%
#  select(Pool, Treatment.1999, Treatment.2000) %>%
#  unique(const_com_short$Pool, incomparables = FALSE) %>%
#  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
#  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
#  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

########################################################################## 
#Option c. Replace missing values with the mean value of pre and post NA#
##########################################################################

#remove rows that have more than two consecutive NAs 
#const_com_revised <- const_com_LACO[-c(2,18,32,96,101,113,142,156,158,159,176,178,190),]

#fill in NA in 2007 with zeros if the population remains zero for the rest of the timeseries
#const_com_revised$`2007` <- ifelse(const_com_revised$Pool %in% c(319,320,324,348,361,367,369,482,525,527,542,545), 0, const_com_revised$`2007`)
  
#1c.fill in the rest of NAs with the mean value of before and after NAs
#We'll use 2001 data for 2000 and 2016 data for 2017, otherwise take the average of before and after
#const_com_revised$`2000` <- ifelse(is.na(const_com_revised$`2000`), const_com_revised$`2001`, const_com_revised$`2000`)
#const_com_revised$`2002` <- ifelse(is.na(const_com_revised$`2002`), as.integer((const_com_revised$`2001`+const_com_revised$`2003`)/2), const_com_revised$`2002`)
#const_com_revised$`2004` <- ifelse(is.na(const_com_revised$`2004`), as.integer((const_com_revised$`2003`+const_com_revised$`2005`)/2), const_com_revised$`2004`)
#const_com_revised$`2007` <- ifelse(is.na(const_com_revised$`2007`), as.integer((const_com_revised$`2006`+const_com_revised$`2008`)/2), const_com_revised$`2007`)
#const_com_revised$`2010` <- ifelse(is.na(const_com_revised$`2010`), as.integer((const_com_revised$`2009`+const_com_revised$`2011`)/2), const_com_revised$`2010`)
#const_com_revised$`2013` <- ifelse(is.na(const_com_revised$`2013`), as.integer((const_com_revised$`2012`+const_com_revised$`2014`)/2), const_com_revised$`2013`)
#const_com_revised$`2016` <- ifelse(is.na(const_com_revised$`2016`), as.integer((const_com_revised$`2015`+const_com_revised$`2017`)/2), const_com_revised$`2016`)
#const_com_revised$`2017` <- ifelse(is.na(const_com_revised$`2017`), const_com_revised$`2016`, const_com_revised$`2017`)

#rearrange the table
#const_com_revised <- const_com_revised %>%
# gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`,`2007`,`2008`,`2009`,`2010`,`2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = "Year", value = "LACOdens")
#const_com_nafilled <- left_join(const_com_revised, const_com_completeveg, by = c("Pool", "Year"))

#2c. count the number of pools
#n_pools <- length(unique(const_com_nafilled$Pool))

#3c. count the number of years
#n_years <- length(unique(const_com_nafilled$Year))

#4c. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
#const_com_nafilled$sum_EG <- rowSums(cbind(const_com_nafilled$BRHO, const_com_nafilled$HOMA, const_com_nafilled$LOMU))

#Native forb (NF) group contains PLST and DOCO:
#const_com_nafilled$sum_NF <- rowSums(cbind(const_com_nafilled$PLST, const_com_nafilled$DOCO))

#5c. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
#LACOdens <- const_com_nafilled %>%
# select(Year, Pool, LACOdens.x) %>%
# spread(key = Year, value = LACOdens.x) %>%
#  select(-Pool) #no NAs
#ERVAdens <- const_com_nafilled %>%
#  select(Year, Pool, ERVAdens) %>%
#  spread(key = Year, value = ERVAdens) %>%
#  select(-Pool)
#ERVAdens$`2000` <- ifelse(is.na(ERVAdens$`2000`), ERVAdens$`2001`, ERVAdens$`2000`)
#ERVAdens$`2002` <- ifelse(is.na(ERVAdens$`2002`), floor((ERVAdens$`2001`+ERVAdens$`2003`)/2), ERVAdens$`2002`)
#ERVAdens$`2004` <- ifelse(is.na(ERVAdens$`2004`), floor((ERVAdens$`2003`+ERVAdens$`2005`)/2), ERVAdens$`2004`)
#ERVAdens$`2007` <- ifelse(is.na(ERVAdens$`2007`), floor((ERVAdens$`2006`+ERVAdens$`2008`)/2), ERVAdens$`2007`)
#ERVAdens$`2010` <- ifelse(is.na(ERVAdens$`2010`), floor((ERVAdens$`2009`+ERVAdens$`2011`)/2), ERVAdens$`2010`)
#ERVAdens$`2013` <- ifelse(is.na(ERVAdens$`2013`), floor((ERVAdens$`2012`+ERVAdens$`2014`)/2), ERVAdens$`2013`)
#ERVAdens$`2016` <- ifelse(is.na(ERVAdens$`2016`), floor((ERVAdens$`2015`+ERVAdens$`2017`)/2), ERVAdens$`2016`)
#ERVAdens$`2017` <- ifelse(is.na(ERVAdens$`2017`), ERVAdens$`2016`, ERVAdens$`2017`)
#sumEGcover <- const_com_nafilled %>%
#  select(Year, Pool, sum_EG) %>%
#  spread(key = Year, value = sum_EG) %>%
#  select(-Pool)
#sumEGcover$`2000` <- ifelse(is.na(sumEGcover$`2000`), sumEGcover$`2001`, sumEGcover$`2000`)
#sumEGcover$`2002` <- ifelse(is.na(sumEGcover$`2002`), floor((sumEGcover$`2001`+sumEGcover$`2003`)/2), sumEGcover$`2002`)
#sumEGcover$`2004` <- ifelse(is.na(sumEGcover$`2004`), floor((sumEGcover$`2003`+sumEGcover$`2005`)/2), sumEGcover$`2004`)
#sumEGcover$`2007` <- ifelse(is.na(sumEGcover$`2007`), floor((sumEGcover$`2006`+sumEGcover$`2008`)/2), sumEGcover$`2007`)
#sumEGcover$`2010` <- ifelse(is.na(sumEGcover$`2010`), floor((sumEGcover$`2009`+sumEGcover$`2011`)/2), sumEGcover$`2010`)
#sumEGcover$`2013` <- ifelse(is.na(sumEGcover$`2013`), floor((sumEGcover$`2012`+sumEGcover$`2014`)/2), sumEGcover$`2013`)
#sumEGcover$`2016` <- ifelse(is.na(sumEGcover$`2016`), floor((sumEGcover$`2015`+sumEGcover$`2017`)/2), sumEGcover$`2016`)
#sumEGcover$`2017` <- ifelse(is.na(sumEGcover$`2017`), sumEGcover$`2016`, sumEGcover$`2017`)
#sumNFcover <- const_com_nafilled %>%
#  select(Year, Pool, sum_NF) %>%
#  spread(key = Year, value = sum_NF) %>%
#  select(-Pool)
#sumNFcover$`2000` <- ifelse(is.na(sumEGcover$`2000`), sumEGcover$`2001`, sumEGcover$`2000`)
#sumNFcover$`2002` <- ifelse(is.na(sumEGcover$`2002`), floor((sumEGcover$`2001`+sumEGcover$`2003`)/2), sumEGcover$`2002`)
#sumNFcover$`2004` <- ifelse(is.na(sumEGcover$`2004`), floor((sumEGcover$`2003`+sumEGcover$`2005`)/2), sumEGcover$`2004`)
#sumNFcover$`2007` <- ifelse(is.na(sumEGcover$`2007`), floor((sumEGcover$`2006`+sumEGcover$`2008`)/2), sumEGcover$`2007`)
#sumNFcover$`2010` <- ifelse(is.na(sumEGcover$`2010`), floor((sumEGcover$`2009`+sumEGcover$`2011`)/2), sumEGcover$`2010`)
#sumNFcover$`2013` <- ifelse(is.na(sumEGcover$`2013`), floor((sumEGcover$`2012`+sumEGcover$`2014`)/2), sumEGcover$`2013`)
#sumNFcover$`2016` <- ifelse(is.na(sumEGcover$`2016`), floor((sumEGcover$`2015`+sumEGcover$`2017`)/2), sumEGcover$`2016`)
#sumNFcover$`2017` <- ifelse(is.na(sumEGcover$`2017`), sumEGcover$`2016`, sumEGcover$`2017`)

#6c. create a matrix of seeds added each year
#seedtrt <- const_com_nafilled %>%
#  select(Pool, Treatment.1999, Treatment.2000) %>%
#  na.omit(seedtrt$Treatment.1999) %>%
#  unique(seedtrt$Pool, incomparables = FALSE) %>%
#  mutate(Y1 = ifelse(Treatment.1999 == "Control", 100, 100)) %>%
#  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
#  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0)) 

################################################################################ 
#Option d. Fill in 2007 data with dummy data and subset complete data 2000-2017#
################################################################################
const_com_dummy <- const_com_LACO 
const_com_dummy$`2007` <- 1
const_dummy_sub <- const_com_dummy[complete.cases(const_com_dummy),] #142 pools

#1d. subset the complete data and fill in 2007 data with zeros
const_dummy_join <- inner_join(const_com, const_dummy_sub, by.y = "Pool") %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0))

#2d. count the number of pools
n_pools <- length(unique(const_dummy_join$Pool))

#3d. count the number of years
n_years <- length(unique(const_dummy_join$Year))

#4b. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
const_dummy_join$sum_EG <- rowSums(cbind(const_dummy_join$BRHO, const_dummy_join$HOMA, const_dummy_join$LOMU))

#Native forb (NF) group contains PLST and DOCO:
const_dummy_join$sum_NF <- rowSums(cbind(const_dummy_join$PLST, const_dummy_join$DOCO))

#5b. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
LACOdens <- const_dummy_join %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens <- const_dummy_join %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover <- const_dummy_join %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover <- const_dummy_join %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

#6b. create a matrix of seeds added each year
seedtrt <- const_dummy_join %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(const_dummy_join$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

###################################################################
#Option e. Follow Option d. step 1 and subset by seeding treatment#
###################################################################
#1e. subset by seeding treatment
GrAB <- const_dummy_join %>%
  mutate(treatment = paste(Treatment.1999, Treatment.2000)) %>%
  filter(treatment == "Group A Group B")
GrBA <- const_dummy_join %>%
  mutate(treatment = paste(Treatment.1999, Treatment.2000)) %>%
  filter(treatment == "Group B Group A")
LALA <- const_dummy_join %>%
  mutate(treatment = paste(Treatment.1999, Treatment.2000)) %>%
  filter(treatment == "Lasthenia Lasthenia")
LANo <- const_dummy_join %>%
  mutate(treatment = paste(Treatment.1999, Treatment.2000)) %>%
  filter(treatment == "Lasthenia NO Lasthenia")

#2e. count the number of pools
n_pools_AB <- length(unique(GrAB$Pool))
n_pools_BA <- length(unique(GrBA$Pool))
n_pools_LALA <- length(unique(LALA$Pool))
n_pools_LANo <- length(unique(LANo$Pool))

#3e. count the number of years
n_years_AB <- length(unique(GrAB$Year))
n_years_BA <- length(unique(GrBA$Year))
n_years_LALA <- length(unique(LALA$Year))
n_years_LANo <- length(unique(LANo$Year))

#4e. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
GrAB$sum_EG <- rowSums(cbind(GrAB$BRHO, GrAB$HOMA, GrAB$LOMU))
GrBA$sum_EG <- rowSums(cbind(GrBA$BRHO, GrBA$HOMA, GrBA$LOMU))
LALA$sum_EG <- rowSums(cbind(LALA$BRHO, LALA$HOMA, LALA$LOMU))
LANo$sum_EG <- rowSums(cbind(LANo$BRHO, LANo$HOMA, LANo$LOMU))

#Native forb (NF) group contains PLST and DOCO:
GrAB$sum_NF <- rowSums(cbind(GrAB$PLST, GrAB$DOCO))
GrBA$sum_NF <- rowSums(cbind(GrBA$PLST, GrBA$DOCO))
LALA$sum_NF <- rowSums(cbind(LALA$PLST, LALA$DOCO))
LANo$sum_NF <- rowSums(cbind(LANo$PLST, LANo$DOCO))

#5e. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
LACOdens_AB <- GrAB %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_AB <- GrAB %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_AB <- GrAB %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_AB <- GrAB %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

LACOdens_BA <- GrBA %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_BA <- GrBA %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_BA <- GrBA %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_BA <- GrBA %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

LACOdens_LALA <- LALA %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_LALA <- LALA %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_LALA <- LALA %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_LALA <- LALA %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

LACOdens_LANo <- LANo %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_LANo <- LANo %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_LANo <- LANo %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_LANo <- LANo %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

#6e. create a matrix of seeds added each year
seedtrt_AB <- GrAB %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(GrAB$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

seedtrt_BA <- GrBA %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(GrBA$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

seedtrt_LALA <- LALA %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(LALA$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

seedtrt_LANo <- LANo %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(LANo$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

############################################################ 
#Option f. Follow Option d. step 1 and subset by pool size#
############################################################
#1f. subset by pool size
Pool_s <- const_dummy_join %>%
  filter(Size == "s")
Pool_m <- const_dummy_join %>%
  filter(Size == "m")
Pool_l <- const_dummy_join %>%
  filter(Size == "l")

#2f. count the number of pools
n_pools_s <- length(unique(Pool_s$Pool))
n_pools_m <- length(unique(Pool_m$Pool))
n_pools_l <- length(unique(Pool_l$Pool))

#3f. count the number of years
n_years_s <- length(unique(Pool_s$Year))
n_years_m <- length(unique(Pool_m$Year))
n_years_l <- length(unique(Pool_l$Year))

#4f. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
Pool_s$sum_EG <- rowSums(cbind(Pool_s$BRHO, Pool_s$HOMA, Pool_s$LOMU))
Pool_m$sum_EG <- rowSums(cbind(Pool_m$BRHO, Pool_m$HOMA, Pool_m$LOMU))
Pool_l$sum_EG <- rowSums(cbind(Pool_l$BRHO, Pool_l$HOMA, Pool_l$LOMU))

#Native forb (NF) group contains PLST and DOCO:
Pool_s$sum_NF <- rowSums(cbind(Pool_s$PLST, Pool_s$DOCO))
Pool_m$sum_NF <- rowSums(cbind(Pool_m$PLST, Pool_m$DOCO))
Pool_l$sum_NF <- rowSums(cbind(Pool_l$PLST, Pool_l$DOCO))

#5f. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
LACOdens_s <- Pool_s %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_s <- Pool_s %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_s <- Pool_s %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_s <- Pool_s %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

LACOdens_m <- Pool_m %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_m <- Pool_m %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_m <- Pool_m %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_m <- Pool_m %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

LACOdens_l <- Pool_l %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_l <- Pool_l %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_l <- Pool_l %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_l <- Pool_l %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

#6f. create a matrix of seeds added each year
seedtrt_s <- Pool_s %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(Pool_s$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))
seedtrt_m <- Pool_m %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(Pool_s$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))
seedtrt_l <- Pool_l %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(Pool_s$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

############################################################
#Option g. Follow Option d. step 1 and subset by pool depth#
############################################################
source("data_compiling/compile_constructed_depth.R")
pool_depth_annual <- const_depth %>%
  drop_na() %>%
  filter(Location == "SeedPlot", Treatment.1999 != "Control") %>%
  group_by(Year, Pool) %>%
  summarize(max_depth = as.numeric(max(Depth))) #inundation data available for 2000, 2002, 2009, 2010, 2011, 2012
pool_depth_average <- pool_depth_annual %>%
  group_by(Pool) %>%
  summarize(avg_max_depth = mean(max_depth)) #191 pools
#25th, 50th, and 75th percentiles of the data
quantile(pool_depth_average$avg_max_depth, c(0.25, 0.50, 0.75)) 
#     25%      50%      75% 
#3.875000 6.800000 8.008333 

#1g. subset by pool depth class
depth_s <- left_join(const_dummy_join, pool_depth_average) %>%
  filter(avg_max_depth <= 3.875000)
depth_m <- left_join(const_dummy_join, pool_depth_average) %>%
  filter(avg_max_depth >= 3.875000 & avg_max_depth < 8.008333)
depth_d <- left_join(const_dummy_join, pool_depth_average) %>%
  filter(avg_max_depth >= 8.008333)

#2g. count the number of pools
n_depth_s <- length(unique(depth_s$Pool))
n_depth_m <- length(unique(depth_m$Pool))
n_depth_d <- length(unique(depth_d$Pool))

#3g. count the number of years
n_years_depth_s <- length(unique(depth_s$Year))
n_years_depth_m <- length(unique(depth_m$Year))
n_years_depth_d <- length(unique(depth_d$Year))

#4g. sum frequency data
#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
depth_s$sum_EG <- rowSums(cbind(depth_s$BRHO, depth_s$HOMA, depth_s$LOMU))
depth_m$sum_EG <- rowSums(cbind(depth_m$BRHO, depth_m$HOMA, depth_m$LOMU))
depth_d$sum_EG <- rowSums(cbind(depth_d$BRHO, depth_d$HOMA, depth_d$LOMU))

#Native forb (NF) group contains PLST and DOCO:
depth_s$sum_NF <- rowSums(cbind(depth_s$PLST, depth_s$DOCO))
depth_m$sum_NF <- rowSums(cbind(depth_m$PLST, depth_m$DOCO))
depth_d$sum_NF <- rowSums(cbind(depth_d$PLST, depth_d$DOCO))

#5g. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
LACOdens_depth_s <- depth_s %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_depth_s <- depth_s %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_depth_s <- depth_s %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_depth_s <- depth_s %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

LACOdens_depth_m <- depth_m %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_depth_m <- depth_m %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_depth_m <- depth_m %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_depth_m <- depth_m %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

LACOdens_depth_d <- depth_d %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens) %>%
  select(-Pool)
ERVAdens_depth_d <- depth_d %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGcover_depth_d <- depth_d %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFcover_depth_d <- depth_d %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

#6g. create a matrix of seeds added each year
seedtrt_depth_s <- depth_s %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(Pool_s$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))
seedtrt_depth_m <- depth_m %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(Pool_s$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))
seedtrt_depth_d <- depth_d %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(Pool_s$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))
