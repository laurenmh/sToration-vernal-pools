Post_ref_00_15 <- rstan::extract(BH_ref_fit_00_15)
refalpha_LACO_00_15 <- as.matrix(Post_ref_00_15$alpha_LACO)
refalpha_EG_00_15 <- as.matrix(Post_ref_00_15$alpha_EG)
refalpha_ERVA_00_15 <- as.matrix(Post_ref_00_15$alpha_ERVA)
refalpha_NF_00_15 <- as.matrix(Post_ref_00_15$alpha_NF)
reflambda_00_15 <- as.matrix(Post_ref_00_15$lambda)
refs_00_15 <- as.matrix(Post_ref_00_15$survival_LACO)
sim_LACO_00_15 <- matrix(nrow = 2000, ncol = 15)
#Plug in the stable equilibrium freq and parameters in the model.
LACO_ref_00_15 <- bh.sim.control(
  LACO = 1,
  EG = as.numeric(const_com_control[2,2:16]),
  ERVA = as.numeric(const_com_control[1,2:16]),
  NF = as.numeric(const_com_control[3,2:16]),
  aii = refalpha_LACO,
  a1 = refalpha_EG,
  a2 = refalpha_ERVA,
  a3 = refalpha_NF,
  lambda = reflambda,
  s = refs,
  g = 0.7,
  glow = 0.2)
GRWR_LACO_ref_00_15 <- log(LACO_ref) #2002-2014 #log transform

ggplot(GRWR_time_ref_00_15, aes(x = Year, y = mean))+
  geom_point() +
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin = lowCI, ymax = upCI), width = 0.4, alpha = 0.9, size = 1)
