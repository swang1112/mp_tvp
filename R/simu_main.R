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

# y_t = A1 y_t-1 + A2 y_t-2 + u_t
#     = A1 y_t-1 + A2 y_t-2 + B eps_t
# in our paper
# y_t = A1 y_t-1 + A2 y_t-2 + inv(A0) eps_t
# a1' = [-0.1996225, -0.71816621, 1.00020814]

write.csv(cbind(Par$A_1, Par$A_2) , file = '../dat/simu/simu_Bs.csv')
write.csv(solve(Par$B), file = '../dat/simu/simu_A0.csv')

# Set parameters of stochastic vola equation ------------------------------

Pi = matrix(c(0.72,  0.08,
              0.2,  0.9), nrow = 2, ncol = 2, byrow = T)
Q = c(0.3, 1.2)

# foo_t = (h1_t, lamb_t)
# foo_t = Pi foo_t-1 + zeta_t, zeta_t ~ N(0, diag(Q))


write.csv(Pi, file = '../dat/simu/simu_Gamma.csv')
write.csv(diag(Q), file = '../dat/simu/simu_Q.csv')

# Set parameters of MP-IV equation ----------------------------------------
# MP shock is the third one
phi = 1
eta_sd = 0.2
(phi/eta_sd)^2

# eps3_t = inv(B)[3,] u_t = a1' u_t
# z_t = phi eps3_t + eta_t, eta_t ~ N(0, eta_sd^2)


# Make data ---------------------------------------------------------------
Tob = 4800 # sample size
set.seed(1234)
dat = .DGP1(Par, Tob, phi, eta_sd, stv_ar = Pi, stv_sd = sqrt(Q), burn = 1000)
write.csv(dat$Y, file = '../dat/simu/simu_f.csv')
write.csv(dat$shock, file = '../dat/simu/simu_eps.csv')
write.csv(dat$iv, file = '../dat/simu/simu_z.csv')
write.csv(dat$h, file = '../dat/simu/simu_h.csv')
write.csv(dat$lambda, file = '../dat/simu/simu_lambda.csv')
