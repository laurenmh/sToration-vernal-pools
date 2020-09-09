# Covariance of lambda and competition term

# Run "prep data before modeling.R" first.
# Run "complex_belowground_v5.R" and "reference_pool_model.R" with real data.

# Load packages
library(ggplot2)
library(tidyverse)

# Calculate mean of LACO, EG, ERVA, and NF each year
# Constructed pools
LACO_mean <- const_dummy_join %>%
  filter(Year != 2000) %>%
  group_by(Year) %>%
  summarise(LACO_mean = mean(LACOdens)) %>%
  as.matrix()
EG_mean <- const_dummy_join %>%
  filter(Year != 2000) %>%
  group_by(Year) %>%
  summarise(EG_mean = mean(sum_EG)) %>%
  as.matrix()
ERVA_mean <- const_dummy_join %>%
  filter(Year != 2000) %>%
  group_by(Year) %>%
  summarise(ERVA_mean = mean(ERVAdens)) %>%
  as.matrix()
NF_mean <- const_dummy_join %>%
  filter(Year != 2000) %>%
  group_by(Year) %>%
  summarise(NF_mean = mean(sum_NF)) %>%
  as.matrix()
#Reference pools
refLACO_mean <- ref_com_join %>%
  group_by(Year) %>%
  summarise(refLACO_mean = mean(LACO)) %>%
  as.matrix()
refEG_mean <- ref_com_join %>%
  group_by(Year) %>%
  summarise(refEG_mean = mean(sum_EG)) %>%
  as.matrix()
refERVA_mean <- ref_com_join %>%
  group_by(Year) %>%
  summarise(refERVA_mean = mean(ERVA)) %>%
  as.matrix()
refNF_mean <- ref_com_join %>%
  group_by(Year) %>%
  summarise(refNG_mean = mean(sum_NF)) %>%
  as.matrix()

# Calculate the competition term - inter and intra-specific competition
Cx <- c()
comp.fxn <- function(aii, a1, a2, a3, LACO, EG, ERVA, NF){
  for(i in 1:17){
    Cx[i] = 1-(1/(1+(aii[i]*LACO[i])+(a1[i]*EG[i])+(a2[i]*ERVA[i])+(a3[i]*NF[i])))
  }  
return(Cx)
}

comp.ref.fxn <- function(aii, a1, a2, a3, LACO, EG, ERVA, NF){
  for(i in 1:13){
    Cx[i] = 1-(1/(1+(aii[i]*LACO[i])+(a1[i]*EG[i])+(a2[i]*ERVA[i])+(a3[i]*NF[i])))
  }  
  return(Cx)
}
comp <- comp.fxn(aii = alpha_LACO_mean[,5],
                 a1 =  alpha_EG_mean[,5],
                 a2 = alpha_ERVA_mean[,5],
                 a3 = alpha_NF_mean[,5],
                 LACO = as.numeric(LACO_mean[,2]),
                 EG = as.numeric(EG_mean[,2]),
                 ERVA = as.numeric(ERVA_mean[,2]),
                 NF = as.numeric(NF_mean[,2])) %>%
        as.data.frame() # for constructed pools

comp_ref <- comp.ref.fxn(aii = refalpha_LACO_mean[,5],
                     a1 =  refalpha_EG_mean[,5],
                     a2 = refalpha_ERVA_mean[,5],
                     a3 = refalpha_NF_mean[,5],
                     LACO = as.numeric(refLACO_mean[,2]),
                     EG = as.numeric(refEG_mean[,2]),
                     ERVA = as.numeric(refERVA_mean[,2]),
                     NF = as.numeric(refNF_mean[,2])) %>%
            as.data.frame()# for reference pools

# Make a lambda table  2003-2015
lambda_ref <- as.data.frame(reflambda_mean[,5]) %>%
  mutate(Year = c(2003:2015))
colnames(lambda_ref) <- c("reference", "Year")  
lambda_const <- as.data.frame(lambda_mean[-c(1,2,16,17),5]) %>%
  mutate(Year = c(2003:2015))
colnames(lambda_const) <- c("constructed", "Year")
lambda_const_ref <- merge(lambda_ref, lambda_const) %>%
  gather(`reference`, `constructed`, key = type, value = lambda)

# Make a competition term table 
comp_ref <- as.data.frame(comp_ref[-c(14,15,16,17),])%>%
  mutate(Year = c(2003:2015))
colnames(comp_ref) <- c("reference", "Year")  
comp <- as.data.frame(comp[-c(1,2,16,17),]) %>%
  mutate(Year = c(2003:2015))
colnames(comp) <- c("constructed", "Year")
comp_const_ref <- merge(comp_ref, comp) %>%
  gather(`reference`, `constructed`, key = type, value = comp)

# Join lambda and competition tables
join_comp_lambda <- left_join(lambda_const_ref, comp_const_ref, by = c("Year", "type"))

# Plot lambda and competition terms
join_comp_lambda$type <- factor(join_comp_lambda$type, levels = c("reference", "constructed"))
ggplot(join_comp_lambda , aes(x = lambda, y = comp, col = type)) +
  geom_point()+
  theme_bw()+
  geom_text(aes(label = join_comp_lambda$Year), hjust = 0, vjust =0)+
  labs(y = "Competition term", x = "Lambda") +
  xlim(0, 80)+
  facet_wrap(~type)

###############################
#competition term with just EG#
###############################
Cx <- c()
comp.EG.fxn <- function(a1, EG){
  for(i in 1:17){
    Cx[i] = 1-(1/(1+(a1[i]*EG[i])))
  }  
  return(Cx)
}

compEG <- comp.EG.fxn(a1 =  alpha_EG_mean[,5],
                      EG = as.numeric(EG_mean[,2])) %>%
  as.data.frame() # for constructed pools

comp_refEG <- comp.EG.fxn(a1 =  refalpha_EG_mean[,5],
                          EG = as.numeric(refEG_mean[,2])) %>%
  as.data.frame()# for reference pools

# Make a competition term table 
comp_refEG <- as.data.frame(comp_refEG[-c(14,15,16,17),])%>%
  mutate(Year = c(2003:2015))
colnames(comp_refEG) <- c("reference", "Year")  
compEG <- as.data.frame(compEG[-c(1,2,16,17),]) %>%
  mutate(Year = c(2003:2015))
colnames(compEG) <- c("constructed", "Year")
comp_const_refEG <- merge(comp_refEG, compEG) %>%
  gather(`reference`, `constructed`, key = type, value = comp)

# Join lambda and competition tables
join_comp_lambdaEG <- left_join(lambda_const_ref, comp_const_refEG, by = c("Year", "type"))

join_comp_lambdaEG$type <- factor(join_comp_lambdaEG$type, levels = c("reference", "constructed"))
ggplot(join_comp_lambdaEG , aes(x = lambda, y = comp)) +
  geom_point()+
  theme_bw()+
  geom_text(aes(label = join_comp_lambdaEG$Year), hjust = 0, vjust =0)+
  labs(y = "Competition term EG", x = "Lambda") +
  xlim(0, 80)+
  facet_wrap(~type)

#################################
#competition term with just LACO#
#################################
Cx <- c()
comp.LACO.fxn <- function(aii, LACO){
  for(i in 1:17){
    Cx[i] = 1-(1/(1+(aii[i]*LACO[i])))
  }  
  return(Cx)
}

compLACO <- comp.LACO.fxn(aii =  alpha_LACO_mean[,5],
                          LACO = as.numeric(LACO_mean[,2])) %>%
  as.data.frame() # for constructed pools

comp_refLACO <- comp.LACO.fxn(aii =  refalpha_LACO_mean[,5],
                              LACO = as.numeric(refLACO_mean[,2])) %>%
  as.data.frame()# for reference pools

# Make a competition term table 
comp_refLACO <- as.data.frame(comp_refLACO[-c(14,15,16,17),])%>%
  mutate(Year = c(2003:2015))
colnames(comp_refLACO) <- c("reference", "Year")  
compLACO <- as.data.frame(compLACO[-c(1,2,16,17),]) %>%
  mutate(Year = c(2003:2015))
colnames(compLACO) <- c("constructed", "Year")
comp_const_refLACO <- merge(comp_refLACO, compLACO) %>%
  gather(`reference`, `constructed`, key = type, value = comp)

# Join lambda and competition tables
join_comp_lambdaLACO <- left_join(lambda_const_ref, comp_const_refLACO, by = c("Year", "type"))

join_comp_lambdaLACO$type <- factor(join_comp_lambdaLACO$type, levels = c("reference", "constructed"))
ggplot(join_comp_lambdaLACO , aes(x = lambda, y = comp)) +
  geom_point()+
  theme_bw()+
  geom_text(aes(label = join_comp_lambdaEG$Year), hjust = 0, vjust =0)+
  labs(y = "Competition term LACO", x = "Lambda") +
  xlim(0, 80)+
  facet_wrap(~type)


