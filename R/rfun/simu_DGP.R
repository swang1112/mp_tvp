require(magrittr)

#' Solve a 3 dimensional DSGE model and convert it into VAR(2)
#'
#' @param alpha
#' @param beta
#' @param k
#' @param gamma
#' @param delta_x
#' @param tau_x
#' @param tau_pi
#' @param tau_r
#' @param rho_x
#' @param rho_pi
#' @param rho_r
#'
#' @return A list containing VAR and structural parameters.
#'
#' @import dplyr
#' @importFrom nleqslv nleqslv
#' @export DSGE2VAR.DGP1
#'
#' @examples NULL
DSGE2VAR.DGP1 <- function(alpha, beta, k, gamma, delta_x, tau_x, tau_pi, tau_r, rho_x, rho_pi, rho_r){
  # parameter initialization
  Par <- rep(0, 18)
  # function to solve
  optim.fct1 <- function(Par, alpha, beta, k, gamma, delta_x, tau_x, tau_pi, tau_r, rho_x, rho_pi, rho_r){
    # DSGE Model: Gamma_0 %*% z_t = Cmat %*% E_t(z_t+1) + Gamma_1 %*% z_t-1 + Hmat %*% w_t
    Gamma_0 <- matrix(c(1, 0, delta_x,
                        -k, 1, 0,
                        (tau_r - 1)*tau_x, (tau_r - 1)*tau_pi, 1), nrow = 3, ncol = 3, byrow = T)
    Cmat    <- matrix(c(gamma, delta_x, 0,
                        0,  beta*(1 + alpha * beta)^-1, 0,
                        0, 0, 0), nrow = 3, ncol = 3,  byrow = T)
    Gamma_1 <- matrix(c(1-gamma, 0,0,
                        0, alpha * (1 + alpha * beta)^-1, 0,
                        0, 0, tau_r), nrow = 3, ncol = 3,  byrow = T)
    Hmat    <- diag(3)
    
    # VAR Model: z_t = PHImat %*% z_t+1 + Bmat %*% w_t
    # E_t(z_t+1) = PHImat %*% z_t + Bmat %*% Fmat %*% w_t
    Fmat   <-  diag(c(rho_x, rho_pi, rho_r))
    PHImat <- matrix(Par[1 : 9], nrow = 3, ncol = 3)
    Bmat   <- matrix(Par[10:18], nrow = 3, ncol = 3)
    
    obj.eq1 <- c(Gamma_1) - kronecker(diag(3), Gamma_0) %*% c(PHImat) + kronecker(t(PHImat) %*% t(PHImat), diag(3)) %*% c(Cmat)
    obj.eq2 <- kronecker(diag(3), Gamma_0) %*% c(Bmat) - kronecker(t(Bmat) %*% t(PHImat), diag(3)) %*%c(Cmat) - kronecker(t(Fmat) %*% t(Bmat), diag(3)) %*% c(Cmat) - c(diag(3))
    erg     <- rbind(obj.eq1, obj.eq2) %>% as.numeric()
    
    return(erg)
  }
  Par_erg <- nleqslv::nleqslv(Par, optim.fct1,
                              alpha = alpha,
                              beta = beta,
                              k = k,
                              gamma = gamma,
                              delta_x = delta_x,
                              tau_x = tau_x,
                              tau_pi = tau_pi,
                              tau_r = tau_r,
                              rho_x = rho_x,
                              rho_pi = rho_pi,
                              rho_r = rho_r)$x %>% round(2)
  # VAR parameters
  PHImat <- matrix(Par_erg[1 : 9], nrow = 3, ncol = 3)
  Bmat   <- matrix(Par_erg[10:18], nrow = 3, ncol = 3)
  Fmat   <- diag(c(rho_x, rho_pi, rho_r))
  
  # VAR(2) representation of DSGE
  A_1       <-   PHImat + Bmat %*% Fmat %*% solve(Bmat)
  A_2       <- - Bmat %*% Fmat %*% solve(Bmat) %*% PHImat
  
  erg          <- list()
  erg[["Phi"]] <- PHImat
  erg[["B"]]   <- Bmat
  erg[["A_1"]] <- A_1
  erg[["A_2"]] <- A_2
  return(erg)
}

Par = DSGE2VAR.DGP1(
  alpha   = 0.5,
  beta    = 0.99,
  k       = 0.05,
  gamma   = 0.5,
  delta_x = 0.1,
  tau_x   = 0.5,
  tau_pi  = 1.8,
  tau_r   = 0.6,
  rho_x   = 0.5,
  rho_pi  = 0.5,
  rho_r   = 0.5)


# trinity with stochastic volatility
# .DGP1 <- function(Par, Tob, phi, eta_sd, stv_ar, stv_sd, burn = 1000){
#   
#   # DGP Parameters
#   A_1      <- Par$A_1
#   A_2      <- Par$A_2
#   B        <- Par$B
#   
#   # weights of lambda for shock 1 and 2
#   s1 = s2 = 0.5
#   
#   Training <- burn
#   
#     # no volatility break
#     # Initilization
#     Y = epsilon = h = matrix(0, 3, Tob + Training)
#     
#     lambda    = rep(NA, Tob + Training)
#     iv        = rep(NA, Tob + Training)
#     Y[, 1]    = c(1,1,1)
#     Y[, 2]    = c(1.2, 0.4, 0.9)
#     
#     
#     h[3,1:2]  = rnorm(2, sd = stv_sd[1])
#     lambda[1:2] = rnorm(2, sd = stv_sd[2])
# 
#     # h[2,1:2]  = rnorm(2, sd = stv_sd[2])
#     # h[3,1:2]  = rnorm(2, sd = stv_sd[3])
#     
#     # simulation of shocks
#     for(i in 3:(Tob + Training)){
#       
#       # h[,i]   = stv_ar %*% h[,(i-1),drop=F] + diag(stv_sd) %*% rnorm(3)
#       foo_lag   = matrix(c(log(h[3,(i-1)]), lambda[i-1]), nrow = 2, ncol = 1)
#       foo       = stv_ar %*% foo_lag + diag(stv_sd) %*% rnorm(2)
#       lambda[i] = foo[2]
#       h[,i]     = c(s1*exp(lambda[i]), s2*exp(lambda[i]), exp(foo[1]))
#       H_sqrt    = diag(sqrt(h[,i]))
#       e_t       = rnorm(3)
#       e_t       = H_sqrt %*% e_t %>% t()
#       
#       # simulation of variable values and iv
#       e_t            <- matrix(e_t, nrow = 3, ncol = 1) # "reshape" e_t
#       Y[, i]         <- A_1 %*% Y[, (i-1)] + A_2 %*% Y[, (i-2)] + B %*% e_t
#       epsilon[,i]    <- e_t
#       iv[i]         <- phi * e_t[3] + rnorm(1, sd = eta_sd)
#      
#       # i = i+1
#     }
#   
#   
#   # discard the traning data
#   Y          <- Y[, -(1:Training)]        %>% t()
#   epsilon    <- epsilon[, -(1:Training)]  %>% t()
#   h          <- h[, -(1:Training)]        %>% t()
#   iv         <- iv[ -(1:Training)]
#   lambda     <- lambda[ -(1:Training)]
#   
#   # give them names
#   colnames(Y)          <- c("x", "pi", "r")
#   colnames(epsilon)    <- c("demand", "supply", "monetary")
#   
#   # result
#   erg                 <- list()
#   erg[["Y"]]          <- Y
#   erg[["shock"]]      <- epsilon
#   erg[["h"]]          <- h
#   erg[["iv"]]         <- iv
#   erg[["lambda"]]     <- lambda
#   
#   
#   return(erg)
# }


# trinity with stochastic volatility (mp goes first)
.DGP1b <- function(Tob, phi, eta_sd, stv_ar, stv_sd, burn = 1000){
  
  # DGP Parameters
  A_1      <- A_1_s
  A_2      <- A_2_s
  B        <- B_s
  
  # weights of lambda for shock 1 and 2
  s1 = s2 = 0.5
  
  Training <- burn
  
  # no volatility break
  # Initilization
  Y = epsilon = h = matrix(0, 3, Tob + Training)
  
  lambda    = rep(NA, Tob + Training)
  iv        = rep(NA, Tob + Training)
  Y[, 1]    = c(1,1,1)
  Y[, 2]    = c(1.2, 0.4, 0.9)
  
  
  h[1,1:2]  = exp(rnorm(2, sd = stv_sd[1]))
  lambda[1:2] = rnorm(2, sd = stv_sd[2])
  
  # h[2,1:2]  = rnorm(2, sd = stv_sd[2])
  # h[3,1:2]  = rnorm(2, sd = stv_sd[3])
  
  # simulation of shocks
  for(i in 3:(Tob + Training)){
    
    foo_lag   = matrix(c(log(h[1,(i-1)]), lambda[i-1]), nrow = 2, ncol = 1)
    foo       = stv_ar %*% foo_lag + diag(stv_sd) %*% rnorm(2)
    lambda[i] = foo[2]
    h[,i]     = c(exp(foo[1]), s1*exp(lambda[i]), s2*exp(lambda[i]))
    H_sqrt    = diag(sqrt(h[,i]))
    e_t       = rnorm(3)
    e_t       = H_sqrt %*% e_t %>% t()
    
    # simulation of variable values and iv
    e_t            <- matrix(e_t, nrow = 3, ncol = 1) # "reshape" e_t
    Y[, i]         <- A_1 %*% Y[, (i-1)] + A_2 %*% Y[, (i-2)] + B %*% e_t
    epsilon[,i]    <- e_t
    iv[i]         <- phi * e_t[1] + rnorm(1, sd = eta_sd)
    
    # i = i+1
  }
  
  
  # discard the traning data
  Y          <- Y[, -(1:Training)]        %>% t()
  epsilon    <- epsilon[, -(1:Training)]  %>% t()
  h          <- h[, -(1:Training)]        %>% t()
  iv         <- iv[ -(1:Training)]
  lambda     <- lambda[ -(1:Training)]
  
  # give them names
  colnames(Y)          <- c("r", "x", "pi")
  colnames(epsilon)    <- c("monetary", "demand", "supply")
  
  # result
  erg                 <- list()
  erg[["Y"]]          <- Y
  erg[["shock"]]      <- epsilon
  erg[["h"]]          <- h
  erg[["iv"]]         <- iv
  erg[["lambda"]]     <- lambda
  
  
  return(erg)
}

# trinity
# .DGP0 <- function(Par, Tob, phi, eta_sd, burn = 1000){
#   
#   # DGP Parameters
#   A_1      <- Par$A_1
#   A_2      <- Par$A_2
#   B        <- Par$B
#   
#   
#   Training <- burn
#   
#   # no volatility break
#   # Initilization
#   Y = epsilon = matrix(0, 3, Tob + Training)
#   
#   iv        = rep(NA, Tob + Training)
#   Y[, 1]    = c(1,1,1)
#   Y[, 2]    = c(1.2, 0.4, 0.9)
#   
#   # simulation of shocks
#   for(i in 3:(Tob + Training)){
#     
#     e_t       = rnorm(3) %>% t()
#     
#     # simulation of variable values and iv
#     e_t            <- matrix(e_t, nrow = 3, ncol = 1) # "reshape" e_t
#     Y[, i]         <- A_1 %*% Y[, (i-1)] + A_2 %*% Y[, (i-2)] + B %*% e_t
#     epsilon[,i]    <- e_t
#     iv[i]         <- phi * e_t[3] + rnorm(1, sd = eta_sd)
#     
#     # i = i+1
#   }
#   
#   
#   # discard the traning data
#   Y          <- Y[, -(1:Training)]        %>% t()
#   epsilon    <- epsilon[, -(1:Training)]  %>% t()
#   iv         <- iv[ -(1:Training)]
#   
#   # give them names
#   colnames(Y)          <- c("x", "pi", "r")
#   colnames(epsilon)    <- c("demand", "supply", "monetary")
#   
#   # result
#   erg                 <- list()
#   erg[["Y"]]          <- Y
#   erg[["shock"]]      <- epsilon
#   erg[["iv"]]         <- iv
#   
#   return(erg)
# }