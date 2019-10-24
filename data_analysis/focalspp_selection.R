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
         nf = PLST + DOCO) %>%
  select(Year:Distance, Size, BRHO, LOMU, HOMA, LYHY, PLST, DOCO, ERVA, LACO, ig, nf) %>%
  gather(species, cover, BRHO:nf) %>%
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
ggplot(com_focal_PA, aes(x=Size, y=tempmean)) + geom_boxplot() + 
  facet_wrap(~species)
ggplot(com_focal_PA, aes(x=Size, y=cv)) + geom_boxplot() + 
  facet_wrap(~species)
ggplot(com_focal_PA, aes(x=Size, y=tempdev)) + geom_boxplot() + 
  facet_wrap(~species)

## look at average abundance of focal species over time
com_focal_mean <- const_com_noNA %>%
  mutate(ig = BRHO + LOMU + HOMA,
         nf = PLST + DOCO) %>%
  select(Year:Distance, Size,  BRHO, LOMU, HOMA, LYHY, PLST, DOCO, ERVA, ig, LACO, nf) %>%
  gather(species, cover, BRHO:nf) %>%
  group_by(Year, species, Size) %>%
  summarize(meancover = mean(cover), secover = calcSE(cover))

ggplot(com_focal_mean, aes(x=Year, y=meancover, color = Size)) + 
  geom_line() + facet_wrap(~species)

## now compare species against each other
com_focal_spread <- com_focal_mean %>%
  filter(species != "ig" & species != "nf") %>%
  select(-secover) %>%
  spread(species, meancover) 

cormat <- cor(com_focal_spread[,3:10])

corrplot(cormat,method="color",type="upper", tl.col="black", tl.cex = .8, diag=F)


