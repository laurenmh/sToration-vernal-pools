
# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# graph constructed abundance of ERVA
head(construct_com)
ggplot(construct_com, aes(x = Year, y = ERVAdens,group=Size,col=Size)) + geom_smooth()


#now looking at frequency
ggplot(construct_com, aes(x = Year, y = ERVA,group=Size,col=Size)) + geom_smooth()


#same process but for LACO
ggplot(construct_com, aes(x = Year, y = LACOdens,group=Size,col=Size)) + geom_smooth()

ggplot(construct_com, aes(x = Year, y = LACO,group=Size,col=Size)) + geom_smooth()

ggplot(construct_com, aes(x = Year, y = LACO,group=Size,col=Size)) + geom_smooth() +
  facet_wrap(~Treatment.1999,nrow=2)

ggplot(construct_com, aes(x = Year, y = LACOdens,group=Size,col=Size)) + geom_smooth() +
  facet_wrap(~Treatment.1999,nrow=2)

ggplot(construct_com, aes(x = Year, y = LACOdens,group=Size,col=Size)) + geom_point() +
  facet_wrap(~Treatment.1999,nrow=2)
subset(construct_com,Treatment.1999=="Group 1" )
