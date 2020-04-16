
# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

require(dplyr)
# graph consted abundance of LACO 

ggplot(const_com, aes(x = Year, y = LACOdens,group=Size,col=Size)) + geom_smooth() +
  facet_wrap(~Treatment.1999,nrow=2)

seeded<-subset(const_com,Treatment.1999!="Control") #removing those artificial pools 
  #that were not seeded
seeded$Size<-as.factor(seeded$Size)

seeded<-subset(seeded,select=c(Year,Pool,Complex,Distance,Size,Treatment.1999,Treatment.2000,LACOdens,LACO))

pools<-unique(seeded$Pool);length(pools)

min.0<-max.pres<-rep(NA,length(pools))
for (i in 1:length(pools)){
      tmp<-subset(seeded,Pool==pools[i])
      tmp<-tmp[order(tmp$Year),];tmp
      tmp2<-subset(tmp,LACOdens==0)
      if (dim(tmp2)[[1]]>0){min.0[i]<-tmp2$Year[1]}
      tmp3<-subset(tmp,LACOdens>0)
      if (dim(tmp3)[[1]]>0){max.pres[i]<-max(tmp3$Year)}
}
sum.dat<-data.frame(pools,min.0,max.pres)
sum.dat

sum.dat$emerge<-ifelse(sum.dat$max.pres>sum.dat$min.0,1,0)
summary(sum.dat)

