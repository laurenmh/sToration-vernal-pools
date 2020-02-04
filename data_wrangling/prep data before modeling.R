#Prep the data to fit it in the model
#1. count the number of pools
#2. count the number of years
#3. fill in the missing data
#4. convert frequency to abundance data

#set data pathway!

#load data
source("data_compiling/compile_composition.R")
View(const_com)

#1. count the number of pools
n_pools <- length(unique(const_com$Pool))

#2. count the number of years
n_years <- length(unique(const_com$Year))

#3. fill in the missing data
#For now, I will just replace missing values with the mean value within each year
const_com_mean <- const_com %>%
  group_by(Year) %>%
  summarize(mean_LACOdens = mean(LACOdens, na.rm = TRUE),
         mean_ERVAdens = mean(ERVAdens, na.rm = TRUE),
         mean_BRHO = mean(BRHO, na.rm = TRUE),
         mean_HOMA = mean(HOMA, na.rm = TRUE),
         mean_LOMU = mean(LOMU, na.rm = TRUE),
         mean_PLST = mean(PLST, na.rm = TRUE),
         mean_DOCO = mean(DOCO, na.rm = TRUE))

const_com_revised <- const_com %>%
  mutate(LACOdens_filled = ifelse(is.na(LACOdens), const_com_mean$mean_LACOdens[match(const_com$Year, const_com_mean$Year)], LACOdens),
         ERVAdens_filled = ifelse(is.na(ERVAdens), const_com_mean$mean_ERVAdens[match(const_com$Year, const_com_mean$Year)], ERVAdens),
         BRHO_filled = ifelse(is.na(BRHO), const_com_mean$mean_BRHO[match(const_com$Year, const_com_mean$Year)], BRHO),
         HOMA_filled = ifelse(is.na(HOMA), const_com_mean$mean_HOMA[match(const_com$Year, const_com_mean$Year)], HOMA),
         LOMU_filled = ifelse(is.na(LOMU), const_com_mean$mean_LOMU[match(const_com$Year, const_com_mean$Year)], LOMU),
         PLST_filled = ifelse(is.na(PLST), const_com_mean$mean_PLST[match(const_com$Year, const_com_mean$Year)], PLST),
         DOCO_filled = ifelse(is.na(DOCO), const_com_mean$mean_DOCO[match(const_com$Year, const_com_mean$Year)], DOCO))

#4. convert frequency to abundance data
#We'll use observed abundance data for LACO and ERVA.
#Nancy identified that log-log is the best transformation:
#density ~ exp -0.415 + frequency ^ 1.356
#We'll do this transformation after summing the percent cover of all species in each group 

#Exotic grass (EG) group contains BRHO, HOMA, and LOMU:
const_com_revised$sum_EG <- rowSums(cbind(const_com_revised$BRHO_filled, const_com_revised$HOMA_filled, const_com_revised$LOMU_filled))
const_com_revised$sum_EGdens <- exp(-0.415) * const_com_revised$sum_EG ^ 1.356

#Native forb (NF) group contains PLST and DOCO:
const_com_revised$sum_NF <- rowSums(cbind(const_com_revised$PLST_filled, const_com_revised$DOCO_filled))
const_com_revised$sum_NFdens <- exp(-0.415) * const_com_revised$sum_NF ^ 1.356
