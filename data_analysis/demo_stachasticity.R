### CREATE A STAN MODEL ###

BH_model_block <- "
data{
    int n_pools; // number of pools
    int n_years; // number of years
    int obs_LACO [n_pools, n_years] ; // observed LACO density // this is an array of integers because the model is a discrete poisson model (see below)
    matrix[n_pools, n_years] obs_EG; // exotic grass density
    matrix[n_pools, n_years] obs_ERVA; // ERVA density
    matrix[n_pools, n_years] obs_NF; // non-native forb density
    int seeds_added [n_pools, 3]; // number of seeds added in 1999-2001
    real low_germ_LACO; // germination rate of LACO in normal years
    real high_germ_LACO; // germination rate of LACO in high litter years
}
parameters{
    vector<lower = 0>[n_years-1] lambda; // max growth rate of LACO in absence of competition
    vector<lower = 0, upper = 1>[n_years-1] alpha_LACO; // competition term for LACO-LACO
    vector<lower = 0, upper = 1>[n_years-1] alpha_EG; // competition term for LACO-exotic grass
    vector<lower = 0, upper = 1>[n_years-1] alpha_ERVA;// competition term for LACO-ERVA
    vector<lower = 0, upper = 1>[n_years-1] alpha_NF; // competition term for LACO-non-native forb
    real <lower = 0, upper = 1> survival_LACO; // survival rate of LACO seeds in the seedbank
}
transformed parameters{
    matrix [n_pools, n_years-1] mu_LACO;// mean expected value of LACO at time t from a discrete BH model
    matrix [n_pools, n_years-3] int_LACO;// intermediate matrix of LACO at time t-1 estimated from values at t-2
    matrix [n_pools, n_years] seedbank;// expected number of LACO seeds in the seedbank at time t
    
    for(i in 1:n_pools){  
        for(j in 1:1){
            mu_LACO[i,j] = 100;
            seedbank[i,j] = 0;
        }
        for(j in 2:2){
            real germ_LACO;
            if (obs_EG[i,j-1] > 100)
                germ_LACO = low_germ_LACO;
            else
                germ_LACO = high_germ_LACO;
            mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
                            obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
                            seedbank[i,j] ./ germ_LACO; // modified Beverton-Holt model
        }
        for(j in 3:n_years-1){
            real germ_LACO;
            if (obs_EG[i,j-1] > 100)
                germ_LACO = low_germ_LACO;
            else
                germ_LACO = high_germ_LACO;
            if (obs_LACO[i,j-1] > 0){
                int_LACO[i,j-2] = (obs_LACO[i,j-2] * lambda[j-2])./(1 + obs_LACO[i,j-2] * alpha_LACO[j-2] + 
                                  obs_EG[i,j-2] * alpha_EG[j-2] + obs_ERVA[i,j-2] * alpha_ERVA[j-2] + obs_NF[i,j-2] * alpha_NF[j-2]) +
                                  seedbank[i,j] ./ germ_LACO;
                mu_LACO[i,j] = (obs_LACO[i,j-1] * lambda[j-1])./(1 + obs_LACO[i,j-1] * alpha_LACO[j-1] + 
                                obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
                                seedbank[i,j] ./ germ_LACO; 
            }
            else{
                int_LACO[i,j-2] = (obs_LACO[i,j-2] * lambda[j-2])./(1 + obs_LACO[i,j-2] * alpha_LACO[j-2] + 
                                  obs_EG[i,j-2] * alpha_EG[j-2] + obs_ERVA[i,j-2] * alpha_ERVA[j-2] + obs_NF[i,j-2] * alpha_NF[j-2]) +
                                  seedbank[i,j] ./ germ_LACO; // LACO at t-2
                mu_LACO[i,j] = (int_LACO[i,j-2] * lambda[j-1])./(1 + int_LACO[i,j-2] * alpha_LACO[j-1] + 
                                obs_EG[i,j-1] * alpha_EG[j-1] + obs_ERVA[i,j-1] * alpha_ERVA[j-1] + obs_NF[i,j-1] * alpha_NF[j-1]) +
                                seedbank[i,j] ./ germ_LACO; // plugging in LACO t-2 for obs LACO t-1
            }
        }
    }
}
model{
    for(a in 1:n_pools){
        for(b in 1:1){
            real germ_LACO;
            germ_LACO = high_germ_LACO;
            obs_LACO[a,b] ~ binomial(seeds_added[a,b], germ_LACO); //the first year's obs_LACO is the initial germination of seeds added in 1999
        }
        for(b in 2:3){
            real germ_LACO;
            if (obs_EG[a,b-1] > 100)
                germ_LACO = low_germ_LACO;
            else
                germ_LACO = high_germ_LACO;
            seeds_germinated ~ binomial((mu_LACO[a,b]), germ_LACO);
            seedbank ~ binomial((mu_LACO[a,b] - seeds_germinated), survival_LACO);
            obs_LACO[a,b] ~ poisson(mu_LACO[a,b] + seeds_added[a,b] * germ_LACO); //the second and third year's obs_LACO is the sum of germination of seeds added in 2000 and 2001 and previous year's population. 
        }
        for(b in 4:(n_years-1)){
            seeds_germinated ~ binomial((mu_LACO[a,b]), germ_LACO);
            seedbank ~ binomial((mu_LACO[a,b] - seeds_germinated), survival_LACO);
            obs_LACO[a,b] ~ poisson(mu_LACO[a,b]); //the rest of the year's obs_LACO is from a poisson distribution of mu_LACO. 
        }
    }
    lambda ~ normal(60,20); //get partially-informed priors from lit
    alpha_LACO ~ normal(0,1);
    alpha_EG ~ normal(0,1);
    alpha_ERVA ~ normal(0,1);
    alpha_NF ~ normal(0,1);
    survival_LACO ~ beta(0.5,0.5);
}"