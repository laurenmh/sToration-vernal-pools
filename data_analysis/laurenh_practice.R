## Make sure to manually set your data pathway in set_pathway.R first!!

# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# graph constructed pond duration over time
ggplot(ref_duration, aes(x = Year, y = Duration.weeks)) + geom_point()
