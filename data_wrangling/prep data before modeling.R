#Prep the data to fit it in the model
#1a. subset complete data
#1b. fill in the missing data
#2. count the number of pools 
#3. count the number of years
#4. sum frequency data
#5. spread the dataset
#6. create a matrix of seeds added each year

#set data pathway!

#load data
source("data_compiling/compile_composition.R")
#View(const_com)

########################################
#How many pools have complete data?#
########################################
const_com_LACO <- const_com %>%
  select(Year, Pool, LACOdens) %>%
  spread(key = Year, value = LACOdens)
const_com_noNA <- const_com_LACO[complete.cases(const_com_LACO),] #only 108 pools have complete data

#1a. subset the complete data
const_com_subset <- inner_join(const_com, const_com_noNA, by.y = "Pool")

#2a. count the number of pools
n_pools <- length(unique(const_com_subset$Pool))

#3a. count the number of years
n_years <- length(unique(const_com_subset$Year))

#4a. convert frequency to abundance data
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
  unique(seedtrt$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))

##################################################################################  
#1b. ALTERNATIVELY, replace missing values with the mean value of pre and post NA#
##################################################################################
#remove rows that have more than two consecutive NAs
const_com_revised <- const_com_LACO[-c(4,20,21,25,26,43,47,49,62,68,70,130,137,159,183,194,211,213,214,226,228,236,240,243,246,254),]

const_com_revised$`2000` <- ifelse(is.na(const_com_revised$`2000`), const_com_revised$`2001`, const_com_revised$`2000`)
const_com_revised$`2002` <- ifelse(is.na(const_com_revised$`2002`), as.integer((const_com_revised$`2001`+const_com_revised$`2003`)/2), const_com_revised$`2002`)
const_com_revised$`2004` <- ifelse(is.na(const_com_revised$`2004`), as.integer((const_com_revised$`2003`+const_com_revised$`2005`)/2), const_com_revised$`2004`)
const_com_revised$`2007` <- ifelse(is.na(const_com_revised$`2007`), as.integer((const_com_revised$`2006`+const_com_revised$`2008`)/2), const_com_revised$`2007`)
const_com_revised$`2010` <- ifelse(is.na(const_com_revised$`2010`), as.integer((const_com_revised$`2009`+const_com_revised$`2011`)/2), const_com_revised$`2010`)
const_com_revised$`2013` <- ifelse(is.na(const_com_revised$`2013`), as.integer((const_com_revised$`2012`+const_com_revised$`2014`)/2), const_com_revised$`2013`)
const_com_revised$`2016` <- ifelse(is.na(const_com_revised$`2016`), as.integer((const_com_revised$`2015`+const_com_revised$`2017`)/2), const_com_revised$`2016`)
const_com_revised$`2017` <- ifelse(is.na(const_com_revised$`2017`), const_com_revised$`2016`, const_com_revised$`2017`)

const_com_revised <- const_com_revised %>%
  gather(`2000`,`2001`,`2002`,`2003`,`2004`,`2005`,`2006`,`2007`,`2008`,`2009`,`2010`,`2011`,`2012`,`2013`,`2014`,`2015`,`2016`,`2017`, key = "Year", value = "LACOdens")

const_com_nafilled <- left_join(const_com_revised, const_com, by = c("Pool", "Year"))

#2b. count the number of pools
n_pools <- length(unique(const_com_nafilled$Pool))

#3b. count the number of years
n_years <- length(unique(const_com_nafilled$Year))

#4b. convert frequency to abundance data
#We'll use observed abundance data for LACO and ERVA.
#We'll sum the percent cover of species in each group 

#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
const_com_nafilled$sum_EG <- rowSums(cbind(const_com_nafilled$BRHO, const_com_nafilled$HOMA, const_com_nafilled$LOMU))

#Native forb (NF) group contains PLST and DOCO:
const_com_nafilled$sum_NF <- rowSums(cbind(const_com_nafilled$PLST, const_com_nafilled$DOCO))

#5b. spread the dataset
#each matrix should have a [n_pools x n_years] dimension
LACOdens <- const_com_nafilled %>%
  select(Year, Pool, LACOdens.x) %>%
  spread(key = Year, value = LACOdens.x) %>%
  select(-Pool)
ERVAdens <- const_com_nafilled %>%
  select(Year, Pool, ERVAdens) %>%
  spread(key = Year, value = ERVAdens) %>%
  select(-Pool)
sumEGdens <- const_com_nafilled %>%
  select(Year, Pool, sum_EG) %>%
  spread(key = Year, value = sum_EG) %>%
  select(-Pool)
sumNFdens <- const_com_nafilled %>%
  select(Year, Pool, sum_NF) %>%
  spread(key = Year, value = sum_NF) %>%
  select(-Pool)

#6b. create a matrix of seeds added each year
seedtrt <- const_com_nafilled %>%
  select(Pool, Treatment.1999, Treatment.2000) %>%
  unique(seedtrt$Pool, incomparables = FALSE) %>%
  mutate(Y1 = ifelse(Treatment.1999 == "Control", 0, 100)) %>%
  mutate(Y2 = ifelse(Treatment.2000 %in% c("Control", "NO Lasthenia"), 0, 100)) %>%
  mutate(Y3 = ifelse(Treatment.2000 == "Lasthenia", 100, 0))
