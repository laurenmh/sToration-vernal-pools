# Example script for sourcing data the begin analysis

# source the data
source("data_compiling/compile_composition.R")
source("data_compiling/compile_constructed_depth.R")
source("data_compiling/compile_reference_depth.R")

# graph constructed pond duration over time
ggplot(constructed_duration, aes(x = Year, y = Duration.weeks)) + geom_point()

# plot LACO absolute abundance
subset(ref_com, Year != 0) %>% 
  select(Year,Pool,LACO) %>% 
  group_by(Year) %>% 
  summarise(LACO = sum(LACO) ) %>% 
  ggplot( aes(x = Year, y = LACO)) + 
  geom_line() +
  theme_classic()

# create population growth rates
laco_t0 <- subset(ref_com, Year != 0) %>% 
             select(Year,Pool,Quadrat,LACO) %>% 
             mutate( Year = Year + 1 ) %>% 
             rename( laco_t0 = LACO )

laco_t1 <- subset(ref_com, Year != 0) %>% 
            select(Year,Pool,Quadrat,LACO) %>% 
            rename( laco_t1 = LACO )


laco_t  <- inner_join( laco_t0, laco_t1 ) %>% 
            mutate( Year = Year - 1 )
  

# stan model
dat_stan <- list(
  n      = nrow(laco_t),
  year_n = laco_t$Year %>% unique %>% length, 
  pool_n = laco_t$Pool %>% unique %>% length,
  year_i = laco_t$Year %>% unique %>% as.factor %>% as.numeric, 
  pool_i = laco_t$Pool %>% unique %>% as.factor %>% as.numeric,
  y      = laco_t1,
  x      = laco_t0
)

"
data{
  int<lower=0> n;         // N. of years. Same for all vital rates.
  int<lower=0> year_n;         // N. of years. Same for all vital rates.
  int<lower=0> pool_n;         // N. of years. Same for all vital rates.

  int<lower=0> year_i[year_n]; // Index for years
  int<lower=0> pool_i[pool_n]; // Index for years
  
  int<lower=0> x[n];     // abundance at time t0
  int<lower=0> y[n];     // abundance at time t1
}

parameters{
  
  real r[year_n];
  real a;

}

transformed parameters{

  abun  
  
}

model{
  
  

}
"