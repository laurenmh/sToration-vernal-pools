# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# graph constructed pond duration over time
ggplot(constructed_depth, aes(x = Year, y = Duration.days)) + geom_point()+ facet_wrap(~Size)

head(constructed_depth)
