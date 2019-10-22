# Example script for sourcing data the begin analysis

# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# graph constructed pond duration over time
ggplot(constructed_duration, aes(x = Year, y = Duration.weeks)) + geom_point()
#just practicing making changes blah blah blah