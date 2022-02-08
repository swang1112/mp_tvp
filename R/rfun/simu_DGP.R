
DGP1 <- function(Par, Tob, phi, eta_sd, stv_ar, stv_sd, burn = 1000){
  
  # DGP Parameters
  A_1      <- Par$A_1
  A_2      <- Par$A_2
  B        <- Par$B
  
  
  Training <- burn
  
    # no volatility break
    # Initilization
    Y = epsilon = h = iv = matrix(0, 3, Tob + Training)
    Y[, 1]  = c(1,1,1)
    Y[, 2]  = c(1.2, 0.4, 0.9)
    
    
    h[1,1:2]  = rnorm(2, sd = stv_sd[1])
    h[2,1:2]  = rnorm(2, sd = stv_sd[2])
    h[3,1:2]  = rnorm(2, sd = stv_sd[3])
    
    # simulation of shocks
    for(i in 3:(Tob + Training)){
        
      h = stv_ar %*% h[,(i-1),drop=F] + diag(stv_sd) %*% rnorm(3)
      H = diag(exp(h/2))
      e_t     <- rnorm(3)
      e_t     <- H %*% e_t %>% t()
      
      # simulation of variable values and iv
      e_t            <- matrix(e_t, nrow = 3, ncol = 1) # "reshape" e_t
      Y[, i]         <- A_1 %*% Y[, (i-1)] + A_2 %*% Y[, (i-2)] + B %*% e_t
      epsilon[,i]    <- e_t
      iv[,i]         <-  phi * e_t[3] + rnorm(1, sd = eta_sd)
     
    }
  
  
  # discard the traning data
  Y          <- Y[, -(1:Training)]        %>% t()
  epsilon    <- epsilon[, -(1:Training)]  %>% t()
  iv         <- iv[, -(1:Training)]  %>% t()
  # give them names
  colnames(Y)          <- c("x", "pi", "r")
  colnames(epsilon)    <- c("demand", "supply", "monetary")
  

  # result
  erg                 <- list()
  erg[["Y"]]          <- Y
  erg[["shock"]]      <- epsilon
  erg[["iv"]]         <- iv
  
  return(erg)
}