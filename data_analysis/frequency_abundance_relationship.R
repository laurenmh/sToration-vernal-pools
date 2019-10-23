# Exploring the relationship between frequency and abundance
# This will allow reference data (frequency) to be transformed to abundance
# The abundance estimates will be used to parameterize a Beverton-Holt population model
# This addresses Q3 and Q4 in the sToration group

library(tidyverse)

# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# Data visualization
