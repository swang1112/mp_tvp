rm(list = ls())
source('./rfun/simu_DGP.R')

# Get parameters of VAR(2) ------------------------------------------------
# Autoregressive parameter: B1, B2 in our paper
Par$A_1 
Par$A_2

# structural impact multiplier: A0^-1 in our paper
Par$B
# A0 in our paper
solve(Par$B)

write.csv(cbind(Par$A_1, Par$A_2) , file = '../dat/simu/simu_Bs.csv')
write.csv(solve(Par$B), file = '../dat/simu/simu_A0.csv')

# Set parameters of stochastic vola equation ------------------------------
Pi = matrix(c( 1.24,  0.09, -0.16,
               0.12,  0.96,  0.06,
               0.24,  0.30,  0.8), nrow = 3, ncol = 3, byrow = T)
Q = c(2.4, 1.2, 0.6)

write.csv(Pi, file = '../dat/simu/simu_Gamma.csv')
write.csv(diag(Q), file = '../dat/simu/simu_Q.csv')

# Set parameters of MP-IV equation ----------------------------------------
# MP shock ist the third one
phi = 1
eta_sd = 0.5
(phi/eta_sd)^2


# Make data ---------------------------------------------------------------
Tob = 480 # sample size
set.seed(1234)
dat = .DGP1(Par, Tob, phi, eta_sd, stv_ar = Pi, stv_sd = sqrt(Q), burn = 1000)
write.csv(dat$Y, file = '../dat/simu/simu_f.csv')
write.csv(dat$shock, file = '../dat/simu/simu_eps.csv')
write.csv(dat$iv, file = '../dat/simu/simu_z.csv')
