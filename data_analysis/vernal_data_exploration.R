source("data_compiling/compile_composition.R")

#calculate standard error
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

# graph density of 5 species over time by treatment
const_com_tr <- const_com %>%
  mutate(treatment = paste(Treatment.1999, Treatment.2000, sep = "-")) 
const_com_tr_2 <- const_com_tr %>%
  filter(!is.na(LACOdens), !is.na(DEDAdens), !is.na(ERVAdens), !is.na(LACHdens), !is.na(PLSTdens)) %>%
  group_by(Year, treatment) %>%
  summarize(meanLACOdens=mean(LACOdens), se_LACOdens=calcSE(LACOdens),
            meanERVAdens=mean(ERVAdens), se_ERVAdens=calcSE(ERVAdens),
            meanDEDAdens=mean(DEDAdens), se_DEDAdens=calcSE(DEDAdens), 
            meanLACHdens=mean(LACHdens), se_LACHdens=calcSE(LACHdens),
            meanPLSTdens=mean(PLSTdens), se_PLSTdens=calcSE(PLSTdens))

## raw data
ggplot(const_com_tr, aes(x = Year, y = LACOdens)) + geom_point(aes(col = as.factor(treatment)))

## box plot
library(ggpubr)
f1 <- ggplot(const_com_tr_2, aes(x = Year, y = meanLACOdens)) +  
        geom_line(aes(color=treatment)) + 
        geom_point(aes(color=treatment)) +
        geom_errorbar(aes(ymin=meanLACOdens-se_LACOdens, ymax=meanLACOdens+se_LACOdens, color=treatment), width=.2) +
        ggtitle("")
f2 <- ggplot(const_com_tr_2, aes(x = Year, y = meanERVAdens)) +  
        geom_line(aes(color=treatment)) + 
        geom_point(aes(color=treatment)) +
        geom_errorbar(aes(ymin=meanERVAdens-se_ERVAdens, ymax=meanERVAdens+se_ERVAdens, color=treatment), width=.2) +
        ggtitle("")
f3 <- ggplot(const_com_tr_2, aes(x = Year, y = meanDEDAdens)) +  
        geom_line(aes(color=treatment)) + 
        geom_point(aes(color=treatment)) +
        geom_errorbar(aes(ymin=meanDEDAdens-se_DEDAdens, ymax=meanDEDAdens+se_DEDAdens, color=treatment), width=.2) +
        ggtitle("")
f4 <- ggplot(const_com_tr_2, aes(x = Year, y = meanLACHdens)) +  
        geom_line(aes(color=treatment)) + 
        geom_point(aes(color=treatment)) +
        geom_errorbar(aes(ymin=meanLACHdens-se_LACHdens, ymax=meanLACHdens+se_LACHdens, color=treatment), width=.2) +
        ggtitle("")
f5 <- ggplot(const_com_tr_2, aes(x = Year, y = meanPLSTdens)) +  
        geom_line(aes(color=treatment)) + 
        geom_point(aes(color=treatment)) +
        geom_errorbar(aes(ymin=meanPLSTdens-se_PLSTdens, ymax=meanPLSTdens+se_PLSTdens, color=treatment), width=.2) +
        ggtitle("")
ggarrange(f1, f2, f3, f4, f5, ncol = 1, nrow = 5, 
          labels = c("a) LACO", "b) ERVA",
                     "c) DEDA", "d) LACH", "e) PLST"),
          common.legend = TRUE, legend = "right", 
          font.label = list(size = 10))

# LACO vs ERVA timeseries
ggplot(const_com_tr_2, aes(x = Year, y = meanLACOdens)) +
  facet_wrap(.~treatment) +
  geom_line(col = "blue") + 
  geom_point(col = "blue") +
  geom_errorbar(aes(ymin=meanLACOdens-se_LACOdens, ymax=meanLACOdens+se_LACOdens), col = "blue", width=.2) +
  geom_line(data = const_com_tr_2, aes(x = Year, y = meanERVAdens), col = "orange") +
  geom_point(data = const_com_tr_2, aes(x = Year, y = meanERVAdens), col = "orange") +
  geom_errorbar(aes(ymin=meanERVAdens-se_ERVAdens, ymax=meanERVAdens+se_ERVAdens), col = "orange", width = 0.2)

# calculate freq of LACO / non-native grass ratio
const_com_ig <- const_com_tr %>%
        mutate(ig = AVFA + BRDI + BRHO + AICA + ALCA + BRMI + HOMA + LOMU + PHMI + POMA + TACA + VUBR) %>%
        filter(!is.na(LACO), !is.na(ig))
ggplot(const_com_ig, aes(x = Year, y = LACO/ig)) +
        geom_point(aes(col = as.factor(treatment))) 

# Non-native grass variation over time
const_com_tr_3 <- const_com_tr %>%
  filter(!is.na(AVFA), !is.na(BRDI), !is.na(BRHO), !is.na(AICA), !is.na(ALCA), !is.na(BRMI), !is.na(HOMA), 
         !is.na(LOMU), !is.na(PHMI), !is.na(POMA), !is.na(TACA), !is.na(VUBR)) %>%
  group_by(Pool, treatment) %>%
  summarise(meanAVFA = mean(AVFA), se_AVFA = calcSE(AVFA),
            meanBRDI = mean(BRDI), se_BRDI = calcSE(BRDI),
            meanBRHO = mean(BRHO), se_BRHO = calcSE(BRHO),
            meanAICA = mean(AICA), se_AICA = calcSE(AICA),
            meanALCA = mean(ALCA), se_ALCA = calcSE(ALCA),
            meanBRMI = mean(BRMI), se_BRMI = calcSE(BRMI),
            meanHOMA = mean(HOMA), se_HOMA = calcSE(HOMA),
            meanLOMU = mean(LOMU), se_LOMU = calcSE(LOMU),
            meanPHMI = mean(PHMI), se_PHMI = calcSE(PHMI),
            meanPOMA = mean(POMA), se_POMA = calcSE(POMA),
            meanTACA = mean(TACA), se_TACA = calcSE(TACA),
            meanVUBR = mean(VUBR), se_VUBR = calcSE(VUBR))
ggplot(const_com_tr_3, aes(x = Pool, y = meanAVFA, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanAVFA-se_AVFA, ymax=meanAVFA+se_AVFA)) +
  facet_wrap(.~treatment)
ggplot(const_com_tr_3, aes(x = Pool, y = meanBRDI, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanBRDI-se_BRDI, ymax=meanBRDI+se_BRDI)) +
  facet_wrap(.~treatment)
# BRHO variable over time
ggplot(const_com_tr_3, aes(x = Pool, y = meanBRHO, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanBRHO-se_BRHO, ymax=meanBRHO+se_BRHO)) +
  facet_wrap(.~treatment)
ggplot(const_com_tr_3, aes(x = Pool, y = meanAICA, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanAICA-se_AICA, ymax=meanAICA+se_AICA)) +
  facet_wrap(.~treatment)
ggplot(const_com_tr_3, aes(x = Pool, y = meanALCA, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanALCA-se_ALCA, ymax=meanALCA+se_ALCA)) +
  facet_wrap(.~treatment)                                                                                                                         
ggplot(const_com_tr_3, aes(x = Pool, y = meanBRMI, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanBRMI-se_BRMI, ymax=meanBRMI+se_BRMI)) +
  facet_wrap(.~treatment)
# HOMA variable over time
ggplot(const_com_tr_3, aes(x = Pool, y = meanHOMA, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanHOMA-se_HOMA, ymax=meanHOMA+se_HOMA)) +
  facet_wrap(.~treatment)
# LOMU variable over time
ggplot(const_com_tr_3, aes(x = Pool, y = meanLOMU, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanLOMU-se_LOMU, ymax=meanLOMU+se_LOMU)) +
  facet_wrap(.~treatment)
ggplot(const_com_tr_3, aes(x = Pool, y = meanPHMI, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanPHMI-se_PHMI, ymax=meanPHMI+se_PHMI)) +
  facet_wrap(.~treatment)
ggplot(const_com_tr_3, aes(x = Pool, y = meanPOMA, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanPOMA-se_POMA, ymax=meanPOMA+se_POMA)) +
  facet_wrap(.~treatment)
# TACA sort of variable over time
ggplot(const_com_tr_3, aes(x = Pool, y = meanTACA, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanTACA-se_TACA, ymax=meanTACA+se_TACA)) +
  facet_wrap(.~treatment)
# VUBR sort of variable over time
ggplot(const_com_tr_3, aes(x = Pool, y = meanVUBR, col = as.factor(treatment))) +
  geom_point() +
  geom_errorbar(aes(ymin=meanVUBR-se_VUBR, ymax=meanVUBR+se_VUBR)) +
  facet_wrap(.~treatment)
