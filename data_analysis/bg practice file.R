
# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# graph consted abundance of ERVA
head(const_com)
ggplot(const_com, aes(x = Year, y = ERVAdens,group=Size,col=Size)) + geom_smooth()


#now looking at frequency
ggplot(const_com, aes(x = Year, y = ERVA,group=Size,col=Size)) + geom_smooth()


#same process but for LACO
ggplot(const_com, aes(x = Year, y = LACOdens,group=Size,col=Size)) + geom_smooth()

ggplot(const_com, aes(x = Year, y = LACO,group=Size,col=Size)) + geom_smooth()

ggplot(const_com, aes(x = Year, y = LACO,group=Size,col=Size)) + geom_smooth() +
  facet_wrap(~Treatment.1999,nrow=2)

ggplot(const_com, aes(x = Year, y = LACOdens,group=Size,col=Size)) + geom_smooth() +
  facet_wrap(~Treatment.1999,nrow=2)

ggplot(const_com, aes(x = Year, y = LACOdens,group=Size,col=Size)) + geom_point() +
   scale_y_continuous(trans='log2') + facet_wrap(~Treatment.1999,nrow=2)

