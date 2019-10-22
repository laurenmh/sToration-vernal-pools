#here it is

source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# graph constructed pond duration over time
ggplot(construct_com, aes(x = VIVI, y = VISA)) + geom_point()


aou