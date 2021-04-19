#Goal: Calculate the long-term growth rate when rare (r_invader) of LACO 

#Step 1. Calculate the stable equilibrium frequency of non-LACO species in the model each year. 
#Step 2. Calculate the annual growth rate of LACO in restored pools when one LACO is introduced into a stable community. 
#Step 3. Do the same step for the reference pools.
#Step 4. Average the growth rates of LACO over time for restored and reference pools.

#--------------------------

# Load data and package
# Remember to set your data pathway first!
source("data_compiling/compile_composition.R") 

#--------------------------

#Step 1. Calculate the stable equilibrium frequency of non-LACO species in the model each year. 

#Filter out just the control plots in reference pools (no LACO present). Average the frequency of non-LACO species across space each year. 
#Non-LACO species in our model:
#ERVA
#Exotic grass group - BRHO, HOMA, LOMU
#Native forb group - PLST, DOCO

const_com_control <- const_com %>% #constructed pools data
  filter(Treatment.1999 == "Control") %>% #filter control plots only
  drop_na() %>% #remove any rows with na
  filter(LACO <= 0) %>% #remove communities with LACO present
    mutate(sumEG = BRHO + HOMA + LOMU, 
         sumNF = PLST + DOCO) %>% #create a new column for sum of EG and sum of NF
  group_by(Year)%>% #summarize by year
  summarize(avg_ERVA = round(mean(ERVA), digits = 0),
            avg_sumEG = round(mean(sumEG), digits = 0),
            avg_sumNF = round(mean(sumNF), digits = 0)) %>%#take the average freq.
  filter(Year != "2017")

#Step 2. Calculate the annual growth rate of LACO in restored pools when one LACO is introduced into a stable community.

#Extract parameters for restored pools. Run "complex_belowground_v5.R".
alpha_LACO_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_LACO")))
alpha_EG_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_EG")))
alpha_ERVA_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_ERVA")))
alpha_NF_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("alpha_NF")))
lambda_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("lambda")))
s_mean <- as.data.frame(get_posterior_mean(BH_fit, pars = c("survival_LACO")))
sim_LACO <- matrix(nrow = 17, ncol = 1)
#Set up the population model for LACO
bh.sim.control <- function(LACO, EG, ERVA, NF, aii, a1, a2, a3, lambda, s, g, glow){
  for(i in 1:1){
    sim_LACO[i] <- LACO*lambda[i]/(1+LACO*aii[i]+EG[i]*a1[i]+ERVA[i]*a2[i]+NF[i]*a3[i])+s*(1-g)*LACO/g #this is the modified Beverton-Holt model we'll use for LACO stem counts
  }
  for(i in 2:nrow(sim_LACO)){
    if (EG[i-1]> 100){
      g = glow
    }
    else{g = g}
    sim_LACO[i] <- LACO*lambda[i]/(1+LACO*aii[i]+EG[i]*a1[i]+ERVA[i]*a2[i]+NF[i]*a3[i])+s*(1-g)*LACO/g 
  }
 return(sim_LACO)
}
#Plug in the stable equilibrium freq and parameters in the model
LACO_const <- bh.sim.control(
                      LACO = 1,
                      EG = const_com_control$avg_sumEG,
                      ERVA = const_com_control$avg_ERVA,
                      NF = const_com_control$avg_sumNF,
                      aii = alpha_LACO_mean[,5],
                      a1 = alpha_EG_mean[,5],
                      a2 = alpha_ERVA_mean[,5],
                      a3 = alpha_NF_mean[,5],
                      lambda = lambda_mean[,5],
                      s = s_mean[,5],
                      g = 0.7,
                      glow = 0.2)
GRWR_LACO_const <- log(LACO_const) #2000-2016

#Step 3. Do the same step for the reference pools.

#Extract parameters for reference pools. Run "reference_pool_model.R".
refalpha_LACO_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("alpha_LACO")))
refalpha_EG_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("alpha_EG")))
refalpha_ERVA_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("alpha_ERVA")))
refalpha_NF_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("alpha_NF")))
reflambda_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("lambda")))
refs_mean <- as.data.frame(get_posterior_mean(BH_ref_fit, pars = c("survival_LACO")))
sim_LACO <- matrix(nrow = 13, ncol = 1)
#Plug in the stable equilibrium freq and parameters in the model.
LACO_ref <- bh.sim.control(
                    LACO = 1,
                    EG = const_com_control$avg_sumEG,
                    ERVA = const_com_control$avg_ERVA,
                    NF = const_com_control$avg_sumNF,
                    aii = refalpha_LACO_mean[,5],
                    a1 = refalpha_EG_mean[,5],
                    a2 = refalpha_ERVA_mean[,5],
                    a3 = refalpha_NF_mean[,5],
                    lambda = reflambda_mean[,5],
                    s = refs_mean[,5],
                    g = 0.7,
                    glow = 0.2)
GRWR_LACO_ref <- log(LACO_ref) #2002-2014

#Step 4. Average the growth rates of LACO over time for restored and reference pools.
mean(GRWR_LACO_const) #Average GRWR LACO for constructed pools
mean(GRWR_LACO_ref) #Average GRWR LACO for reference pools
