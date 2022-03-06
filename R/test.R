y = dat$Y
u = resid(vars::VAR(y, p = 2, type = 'none'))


u_std = u 
for (i in 1:nrow(u)) {
  H_sqrt_inv = diag(1/sqrt(dat$h[i+2,])) %*% A_0
  
  u_std[i,] = H_sqrt_inv %*% t(u[i,,drop = F])
}
u_std %>% cov

dat$h %>% colMeans()
par(mfrow=c(3,1), mar = c(1,1,0,0))
for (i in 1:3) {
  plot(dat$h[,i], type = "l")
}


par(mfrow=c(3,1), mar = c(1,1,0,0))
for (i in 1:3) {
  plot(dat$Y[,i], type = "l")
}

dat$shock %>% colMeans()
par(mfrow=c(3,1), mar = c(1,1,0,0))
for (i in 1:3) {
  plot(dat$shock[,i], type = "l")
}

par(mfrow=c(1,1), mar = c(1,1,0,0))
plot(dat$shock[,1], type = "l")
lines(dat$iv, col = 2)

par(mfrow=c(1,1), mar = c(2,2,2,0))
plot(dat$lambda, type = "l")
mean(dat$lambda)
mean(exp(dat$lambda)/2)

# within DGP -------------------------------------------------------------


u = resid(vars::VAR(Y, p = 2, type = 'none'))


u_std = u 
for (i in 1:nrow(u)) {
  H_sqrt_inv = diag(1/sqrt(h[i+2,])) %*% A_0
  
  u_std[i,] = H_sqrt_inv %*% t(u[i,,drop = F])
}
u_std %>% cov


h %>% colMeans()
par(mfrow=c(3,1), mar = c(1,1,0,0))
for (i in 1:3) {
  plot(h[,i], type = "l")
}


par(mfrow=c(3,1), mar = c(1,1,0,0))
for (i in 1:3) {
  plot(Y[,i], type = "l")
}

epsilon %>% colMeans()
par(mfrow=c(3,1), mar = c(1,1,0,0))
for (i in 1:3) {
  plot(epsilon[,i], type = "l")
}

par(mfrow=c(1,1), mar = c(1,1,0,0))
plot(epsilon[,1], type = "l")
lines(iv, col = 2)

par(mfrow=c(1,1), mar = c(2,2,2,0))
plot(lambda, type = "l")
mean(lambda)
mean(exp(lambda)/2)
