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


# calculate freq of LACO / non-native ratio
const_com_ig <- const_com_tr %>%
        mutate(ig = AVFA + BRDI + BRHO + AICA + ALCA + BRMI + HOMA + LOMU + PHMI + POMA + TACA + VUBR)
ggplot(const_com_ig, aes(x = Year, y = LACO)) +
        geom_point(aes(col = as.factor(treatment)))

                                                                                                                                             