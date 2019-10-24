library(corrplot)

source("data_compiling/compile_composition.R")


#calculate standard error
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

# only keep pools that never have NA
const_com_noNA <- const_com %>%
  # remove the controls
  filter(Treatment.1999 != "Control") %>%
  # identify if a pool ever had an NA
  group_by(Pool) %>%
  mutate(tocheck = ifelse(is.na(LACOdens), 1, 0),
         tocheck2 = max(tocheck)) %>%
  # only keep pools that never had NAs
  filter(tocheck2 == 0) %>%
  tbl_df() %>%
  # remove the hack columns that got us there
  select(-tocheck, -tocheck2)


# check the balance of those pools - quite well balanced!
const_com_noNA %>%
  group_by(Treatment.2000, Year) %>%
  summarize(count = n()) %>%
  filter(Year == 2002)

# do a quick check that Sharon and Akasha are right about dominants (they are!)
spvec <- names(const_com_noNA[,13:94])
abundvec <- colSums(const_com_noNA[,13:94])

abundcheck <- as.data.frame(cbind(spvec, abundvec)) %>%
  tbl_df() %>%
  mutate(abundvec = as.numeric(as.character(abundvec))) %>%
  arrange(-abundvec)

## look at frequency that a focal species was present over time
com_focal_PA <- const_com_noNA %>%
  mutate(ig = BRHO + LOMU + HOMA,
         f = PLST + DOCO + LYHY) %>%
  select(Year:Distance, Size, BRHO, LOMU, HOMA, LYHY, PLST, DOCO, ERVA, LACO, ig, f) %>%
  gather(species, cover, BRHO:f) %>%
  mutate(PA = ifelse(cover == 0, 0, 1)) %>%
  group_by(species, Pool) %>%
  mutate(totyr = n(),
         sumPA = sum(PA), 
         propP = sumPA/totyr,
         tempmean = mean(cover),
         tempdev = sd(cover),
         cv = tempmean/tempdev) %>%
  select(species, Pool, Size, sumPA, propP, tempmean, tempdev, cv) %>%
  unique()

# boxplot of frequency, temporal mean and cv for each species across pools
ggplot(com_focal_PA, aes(x=Size, y=propP)) + geom_boxplot() + 
  facet_wrap(~species)
# ggsave("figures_prelim/focalspp_propPresent.pdf", width = 8, height = 6)

ggplot(com_focal_PA, aes(x=Size, y=tempmean)) + geom_boxplot() + 
  facet_wrap(~species)
# ggsave("figures_prelim/focalspp_temporalmean.pdf", width = 8, height = 6)

ggplot(com_focal_PA, aes(x=Size, y=cv)) + geom_boxplot() + 
  facet_wrap(~species)
ggplot(com_focal_PA, aes(x=Size, y=tempdev)) + geom_boxplot() + 
  facet_wrap(~species)

## look at average abundance of focal species over time
com_focal_mean <- const_com_noNA %>%
  mutate(ig = BRHO + LOMU + HOMA,
         f = PLST + DOCO + LYHY) %>%
  select(Year:Distance, Size,  BRHO, LOMU, HOMA, LYHY, PLST, DOCO, ERVA, ig, LACO, f) %>%
  gather(species, cover, BRHO:f) %>%
  group_by(Year, species, Size) %>%
  summarize(meancover = mean(cover), secover = calcSE(cover))

ggplot(com_focal_mean, aes(x=Year, y=meancover, color = Size)) + 
  geom_line() + facet_wrap(~species)
# ggsave("figures_prelim/focalspp_temporaltrend.pdf", width = 8, height = 6)

## now compare species against each other
com_focal_spread <- com_focal_mean %>%
  filter(species != "ig" & species != "f") %>%
  select(-secover) %>%
  spread(species, meancover) 

cormat <- cor(com_focal_spread[,3:10])

# pdf("figures_prelim/focalspp_corrplot.pdf", width = 6, height = 6)
corrplot(cormat,method="color",type="upper", tl.col="black", tl.cex = .8, diag=F)
# dev.off()

## Does LACO hace a rescue effect?

const_com_rescue <- const_com_noNA %>%
  select(Year:Treatment.2000, LACO) %>%
  arrange(Pool, Year) %>%
  group_by(Pool, LACO) %>%
  mutate(absent = ifelse(LACO == 0, 1, 0)) %>%
  mutate(runningabsent = absent,
         runningabsent = ifelse( absent ==1, lag(runningabsent) + 1 , runningabsent),
         newrunningabsent = ifelse(absent == 1 & lag(runningabsent) == 2, 3, runningabsent)) %>%
  group_by(Pool, LACO) %>%
  mutate(maxabsent = max(newrunningabsent)) 

const_com_rescue2 <-  const_com_rescue %>%
  group_by(Year, absent) %>%
  summarize(count = n())

# a lot go to 0 and stay there
ggplot(subset(const_com_rescue2, absent == 0), aes(x=Year, y= count)) + geom_point()


const_com_rescue3 <- const_com_noNA %>%
  select(Year:Treatment.2000, LACO) %>%
  mutate(present = ifelse(LACO == 0, 0, 1)) 

ggplot(const_com_rescue3, aes(x=Year, y=present)) + geom_point() + 
  geom_line() + facet_wrap(~Pool)
# ggsave("figures_prelim/LACO_PA.pdf", width = 12, height = 8)
