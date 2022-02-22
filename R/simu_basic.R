rm(list = ls())
library(magrittr)
source('./rfun/simu_DGP.R')

# Get parameters of VAR(2) ------------------------------------------------
# Autoregressive parameter: B1, B2 in our paper
Par$A_1 
Par$A_2

# structural impact multiplier: A0^-1 in our paper
Par$B
# A0 in our paper
solve(Par$B)


# Make data ---------------------------------------------------------------
Tob = 480 # sample size
set.seed(1234)
dat = .DGP0(Par, Tob, phi = 1, eta_sd = 0.5, burn = 1000)
y = dat$Y


require(vars)
var_est = VAR(y, p = 2, type = 'none')
u = var_est %>% resid
z = dat$iv[-c(1:2)]

# priors and training sample ----------------------------------------------
train_ob = 60
train_u  = u[1:train_ob,]
train_z  = z[1:train_ob]


train_phi2= (crossprod(train_z, train_u)/train_ob) %*% solve(crossprod(train_u)/train_ob) %*% (crossprod(train_u, train_z)/train_ob)
M_phi = sqrt(train_phi2)
V_phi = 0.5

tau_eta = 2
theta_eta = 1


train_omic = 1/.25

V_a1 = solve(crossprod(train_u) * train_omic)
M_a1 = V_a1 %*% (crossprod(train_u, train_z)* train_omic)

# eta = train_z - 1 * train_u %*% M_a1
# var(eta)



# MCMC --------------------------------------------------------------------
burn = 10000
nMC  = burn + 2000

Teff = Tob - train_ob - 2
zeff = z[-(1:train_ob)]
ueff = u[-(1:train_ob),]

phi_star = sqrt(train_phi2)
a1_star = M_a1


# memory reservation --
MC_omic     = matrix(NA, nrow = nMC - burn, ncol = 1)
MC_phi      = matrix(NA, nrow = nMC - burn, ncol = 1)
MC_eta_sig2 = matrix(NA, nrow = nMC - burn, ncol = 1)
MC_a1       = matrix(NA, nrow = nMC - burn, ncol = 3)
MC_eps      = matrix(NA, nrow = nMC - burn, ncol = Teff)

for (i in 1:nMC) {
  
  eps_star = ueff %*% a1_star
  eta = zeff - as.numeric(phi_star) * eps_star
  
  # draw eta_sig
  tau_eta_star    = tau_eta + Teff;
  theta_eta_star  = theta_eta + crossprod(eta);
  foo0            = matrix(rnorm(tau_eta_star), nrow = tau_eta_star, ncol = 1)
  sigma_eta2_star = theta_eta_star / as.numeric(crossprod(foo0))
  
  # draw phi
  V_phi_star = 1/(crossprod(eps_star)/sigma_eta2_star + 1/V_phi) 
  M_phi_star = V_phi_star * ( t(eps_star) %*% zeff /sigma_eta2_star + M_phi/V_phi)
  phi_star   = rnorm(1, mean = as.numeric(M_phi_star), sd = sqrt(as.numeric(V_phi_star)))
  
  omic_star = phi_star^2/sigma_eta2_star
  
  # draw a1
  V_a1_star = solve(crossprod(ueff) * as.numeric(omic_star) + solve(V_a1))
  M_a1_star = V_a1_star %*% (crossprod(ueff, zeff)* as.numeric(omic_star) + solve(V_a1) %*% M_a1 )
  a1_star   = M_a1_star + t(chol(V_a1_star)) %*% matrix(rnorm(3), nrow = 3, ncol = 1)
  
  # store
  if (i > burn){
    MC_omic[i-burn,]      = omic_star
    MC_phi[i-burn,]       = phi_star
    MC_eta_sig2[i-burn,]  = sigma_eta2_star
    MC_a1[i-burn,]        = a1_star
    MC_eps[i-burn,]       = eps_star
  }
  
}

par(mfrow = c(3,2), mar = c(1,4,1,1))
plot(MC_omic, type = 'l'); abline(h = 1/.25, col = 2)
plot(MC_phi, type = 'l'); abline(h = 1, col = 2)
plot(MC_eta_sig2, type = 'l'); abline(h = .25, col = 2)
plot(MC_a1[,1], type = 'l'); abline(h = -0.1996225, col = 2)
plot(MC_a1[,2], type = 'l'); abline(h = -0.71816621, col = 2)
plot(MC_a1[,3], type = 'l'); abline(h = 1.00020814, col = 2)
par(mfrow = c(1,1))


par(mfrow = c(3,2), mar = c(4.5,1,1,1))
hist(MC_omic, breaks = 100, main = round(median(MC_omic),4)); abline( v = 1/.25, col = 2)
hist(MC_phi , breaks = 100, main = round(median(MC_phi),4)); abline( v = 1, col = 2)
hist(MC_eta_sig2, breaks = 100, main = round(median(MC_eta_sig2),4)); abline( v = .25, col = 2)
hist(MC_a1[,1], breaks = 100, main = round(median(MC_a1[,1]),4)); abline( v = -0.1996225, col = 2)
hist(MC_a1[,2], breaks = 100, main = round(median(MC_a1[,2]),4)); abline( v = -0.71816621, col = 2)
hist(MC_a1[,3], breaks = 100, main = round(median(MC_a1[,3]),4)); abline( v = 1.00020814, col = 2)
par(mfrow = c(1,1))



