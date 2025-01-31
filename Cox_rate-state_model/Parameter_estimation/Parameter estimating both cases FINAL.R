setwd("/Users/zhuldyzaybaki/Desktop/Clean code/")
counts = as.matrix(read.csv("counts.csv", header=F))
inW = as.matrix(read.csv("inW.csv", header=F))
mpressW = as.matrix(read.csv("mpressW.csv", header=T))
VprodW = as.matrix(read.csv("VprodW.csv", header=F))


delta = 1
deltaCell = max(inW)
nspatial = dim(inW)[1]/27
N = dim(inW)[1]

################ Functions ###################

likfunc1 = function(theta, counts,inW, VprodW, mpressW, sigma2, delta, nspatial){
  zeta = list(theta1 = theta[1], theta2 = theta[2], alpha = 0, beta= -20)
  res = approxlik(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)
  return(res)
}

likfunc2 = function(theta, counts,inW, VprodW, mpressW, sigma2, delta, nspatial){
  zeta = list(theta1 = theta[1], theta2 = theta[2], alpha = theta[3], beta= theta[4])
  res = approxlik(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)
  return(res)
}



################ Polynomial m(s,t). ###################

sigma2 = (7.17 * 7.17)
mpressW = mpressW/500
sigma2 = sigma2 / (500*500)


# 1. Fit theta1 and theta2
theta = c(-5, 21)
opt1 = optim(par = theta, fn = likfunc1, method = "L-BFGS-B", control = list(fnscale = -1),
             counts = counts, inW = inW, VprodW = VprodW, mpressW = mpressW, 
             sigma2=sigma2, delta = delta, nspatial = nspatial)
opt1$par
zeta = list(theta1 = opt1$par[1], theta2 = opt1$par[2], alpha = 0, beta= -20)
approxlik(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)
approxscore(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)


# 2. Update for alpha and beta
theta = c(opt1$par[1], opt1$par[2], 0.1, -2.3)
opt2 = optim(par = theta, fn = likfunc2, method = "L-BFGS-B", control = list(fnscale = -1),
             counts = counts, inW = inW, VprodW = VprodW, mpressW = mpressW, 
             sigma2=sigma2, delta = delta, nspatial = nspatial)
opt2$par #   -5.288393  13.310838   5.038845 -13.711295

zeta = list(theta1 = opt2$par[1], theta2 = opt2$par[2], alpha = opt2$par[3], beta= opt2$par[4])
approxlik(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)
approxscore(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)

# 3. Improve score
newzeta = newtonstep( zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial )
approxlik(newzeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)  
approxscore(newzeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)
zeta=newzeta

## After one iteration: $theta1 = -5.297343  $theta2 = 13.28657  $alpha = 5.035816  $beta = -15.34813
## Score: $s1 = -0.07103783  $s2 = -0.003998882  $s3 = -0.007763556  $s4 = -3.864789e-05

# 4. Godabme CIs:
G = godambe(zeta, VprodW, mpressW, sigma2, delta, deltaCell, nspatial )  
stds = sqrt(diag(G))
c(zeta$theta1, zeta$theta2, zeta$alpha) - 1.96*stds # -5.484558  9.411129  3.534851
c(zeta$theta1, zeta$theta2, zeta$alpha) + 1.96*stds # -5.110128 17.162011  6.536780


# 5. Bootstrapped CIs:
lam = inW * exp(zeta$theta1 + zeta$theta2*VprodW)/rateGm(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
nn = matrix(NA, nrow = N, ncol = 100)
for (i in 1:N) {
  nn[i,] = rpois(100, lam[i])
}
pars = matrix(NA, nrow = 4, ncol = 100)
for (i in 1:100) {
  print(i)
  parm=c(opt1$par[1], opt1$par[2], 0.1, -2.3) 
  optm = optim(par = parm, fn=likfunc2, method = "L-BFGS-B", control=list(fnscale=-1), hessian=TRUE, 
               lower = c(-20, -20, 0.0, -20), upper = c(20, 20, 20, 20), sigma2 = sigma2,
               counts=nn[,i], inW=inW, VprodW=VprodW, mpressW=mpressW, delta=delta, nspatial=nspatial)
  pars[,i] = optm$par  
}

ci_theta1 <- quantile(pars[1,], c(0.025, 0.975)) # -5.554393 -4.915087 
ci_theta2 <- quantile(pars[2,], c(0.025, 0.975)) #  8.412964 16.880758 
ci_alpha <- quantile(pars[3,], c(0.025, 0.975)) # 0.007717155 0.040000000 
ci_beta <- quantile(pars[4,], c(0.025, 0.975)) # -15.491033  -2.261588 






################ NAM output m(s,t) ###################

mpressW = read.csv("input data/mpressW_nam.csv")
sigma2 = (2.3 * 2.3)/4
mpressW = mpressW/100
sigma2 = sigma2 / (100*100)


# 1. Fit theta1 and theta2
theta = c(-5, 21)
opt1 = optim(par = theta, fn = likfunc1, method = "L-BFGS-B", control = list(fnscale = -1),
             counts = counts, inW = inW, VprodW = VprodW, mpressW = mpressW, 
             sigma2=sigma2, delta = delta, nspatial = nspatial)
opt1$par # -4.892968 15.521052
zeta = list(theta1 = opt1$par[1], theta2 = opt1$par[2], alpha = 0, beta= -20)
approxlik(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)
approxscore(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)


# 2. Update for alpha and beta
theta = c(opt1$par[1], opt1$par[2], 0.1, -2.3)
opt2 = optim(par = theta, fn = likfunc2, method = "L-BFGS-B", control = list(fnscale = -1),
             counts = counts, inW = inW, VprodW = VprodW, mpressW = mpressW, 
             sigma2=sigma2, delta = delta, nspatial = nspatial)
opt2$par #   -5.470231 13.286847  1.299532 -8.865480

zeta = list(theta1 = opt2$par[1], theta2 = opt2$par[2], alpha = opt2$par[3], beta= opt2$par[4])
approxlik(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)
approxscore(zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)

# 3. Improve score
newzeta = newtonstep( zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial )
approxlik(newzeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)  
approxscore(newzeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial)
zeta=newzeta
## After two iteration: $theta1 = -5.470929  $theta2 = 13.3209  $alpha = 1.292514 $beta = -10.71836
## Score: $s1 = -0.006035527 $s2 = -1.564074e-05 $s3 = 0.004394845  $s4 = -0.0002939171


# 4. Godabme CIs:
G = godambe(zeta, VprodW, mpressW, sigma2, delta, deltaCell, nspatial )  
stds = sqrt(diag(G))
c(zeta$theta1, zeta$theta2, zeta$alpha) - 1.96*stds # -5.6705802  9.2035386  0.9653541
c(zeta$theta1, zeta$theta2, zeta$alpha) + 1.96*stds # -5.271277 17.438254  1.619674


# 5. Bootstrapped CIs:
lam = inW * exp(zeta$theta1 + zeta$theta2*VprodW)/rateGm(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
nn = matrix(NA, nrow = N, ncol = 100)
for (i in 1:N) {
  nn[i,] = rpois(100, lam[i])
}
pars = matrix(NA, nrow = 4, ncol = 100)
for (i in 1:100) {
  print(i)
  parm=c(opt1$par[1], opt1$par[2], 0.1, -2.3) 
  optm = optim(par = parm, fn=likfunc2, method = "L-BFGS-B", control=list(fnscale=-1), hessian=TRUE, 
               lower = c(-20, -20, 0.0, -20), upper = c(20, 20, 20, 20), sigma2 = sigma2,
               counts=nn[,i], inW=inW, VprodW=VprodW, mpressW=mpressW, delta=delta, nspatial=nspatial)
  pars[,i] = optm$par  
}

ci_theta1 <- quantile(pars[1,], c(0.025, 0.975)) # -5.909454 -5.222704 
ci_theta2 <- quantile(pars[2,], c(0.025, 0.975)) #  6.24043 17.58266 
ci_alpha <- quantile(pars[3,], c(0.025, 0.975))/100 # 0.009429646 0.041523144


