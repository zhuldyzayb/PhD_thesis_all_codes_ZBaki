############  My functions for parameter estimation  ############

## Estimating values with MC's code and my code match for 32 by 32 grid

## produces exact same result with MC
rateGm <- function( alpha, beta, mpressW, delta, nspatial)
{
  res <- rep(0, length(mpressW))
  for(i in 1:length(mpressW)){
    k = floor((i-1)/nspatial)
    if(k==0) { res[i] = 1 }
    else {
      res[i] = ( res[i - nspatial] + exp(beta)*delta ) * exp( alpha*(mpressW[i] - mpressW[i-nspatial]))
    }
  } 
  return(res)
}

## exact the same function
approxlik <- function( zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial )
{
  gammam <- rateGm(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
  lambdam <- exp( zeta$theta1 + VprodW * zeta$theta2 )
  correction <- exp( - zeta$alpha * zeta$alpha * sigma2 ) 
  lambdaW <- correction * lambdam * inW * delta / gammam
  
  res <- sum( counts[inW>0 & lambdaW > 0] * log(lambdaW[inW>0 & lambdaW > 0]) )
  res <- res - sum( lambdaW[inW>0 & lambdaW > 0] )
  
  return(res)
}


## almost the same 
Elambda <- function( alpha, beta, mpressW, sigma2, delta, nspatial )
{
  res <- 0.0 * mpressW
  pressW <- mpressW
  
  for( n in 1:1000) { 
    pressW <- mpressW + rnorm(n=length(mpressW), mean = 0.0, sd=sqrt(sigma2) )
  
    res <- res + 1.0/rateGm( alpha, beta, pressW, delta, nspatial )
  }
  
  return(res / 1000.0)
}
## similar to MC definition
dalpha <- function( alpha, beta, mpressW, delta, nspatial )
{
  res <- 0*mpressW
  rg <- rateGm( alpha, beta, mpressW, delta, nspatial )
  
  for(i in 1:length(mpressW)){
    k = floor((i-1)/nspatial)
    if(k>0){
      res[i] = res[i-nspatial] * exp( alpha * (mpressW[i] - mpressW[i-nspatial]) ) + (mpressW[i] - mpressW[i-nspatial])*rg[i]
    }
  }
  return(res)
}
## almost the same definition
dbeta <- function( alpha, beta, mpressW, delta, nspatial ) 
{
  res <- 0*mpressW
  
  for(i in 1:length(mpressW)){
    k = floor((i-1)/nspatial)
    if(k>0){
      res[i] = ( res[i-nspatial] + exp(beta)*delta ) * exp( alpha*(mpressW[i] - mpressW[i-nspatial]) )
    }
  }
  
  
  return(res)
}

## exactly the same function
approxscore <- function( zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial )
{
  gammam <- rateGm(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
  Egammam <- Elambda(zeta$alpha, zeta$beta, mpressW, sigma2, delta, nspatial)
  lambdam <- exp( zeta$theta1 + VprodW * zeta$theta2 )
  ElambdaW <- lambdam * inW * delta * Egammam
  
  da <- dalpha(zeta$alpha, zeta$beta, mpressW, delta, nspatial) 
  db <- dbeta(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
  
  b1 <- sum( counts - ElambdaW ) 
  b2 <- sum( (counts - ElambdaW) * VprodW ) 
  b3 <- sum( - (da/gammam) * (counts - ElambdaW) ) 
  b4 <- sum( - (db/gammam) * (counts - ElambdaW) )
  
  return( list( s1=b1, s2=b2, s3=b3, s4=b4 ) )
}

## almost the same
dalphadalpha <- function( alpha, beta, mpressW, delta, nspatial ) 
{
  res <- 0*mpressW
  da <- dalpha( alpha, beta, mpressW, delta, nspatial )
  
  for(i in 1:length(mpressW)){
    k = floor((i-1)/nspatial)
    if(k>0){
      dif = mpressW[i] - mpressW[i-nspatial]
      res[i] = res[i-nspatial] * exp(alpha*dif) + dif * ( da[i-nspatial] * exp(alpha*dif) + da[i] )
    }
  }
  
  return(res)
}

## exactly the same
dbetadbeta <- function( alpha, beta, mpressW, delta, nspatial )
{
  res <- dbeta( alpha, beta, mpressW, delta, nspatial )
  return(res)
}

dalphadbeta <- function( alpha, beta, mpressW, delta, nspatial )
{
  res <- 0*mpressW
  db <- dbeta( alpha, beta, mpressW, delta, nspatial )
  
  for(i in 1:length(mpressW)){
    k = floor((i-1)/nspatial)
    if(k>0){
      res[i] = res[i-nspatial] * exp( alpha*(mpressW[i] - mpressW[i-nspatial]) ) + (mpressW[i] - mpressW[i-nspatial])*db[i]
    }
  }
  
  return(res)
}

## slightly tweeked definition
Edalpha <- function( zeta, VprodW, mpressW, sigma2, delta, nspatial )
{
  res <- 0.0 * mpressW
  pressW <- 0.0 * mpressW
  
  for( n in 1:1000) {
    pressW <- mpressW + rnorm(n=length(mpressW), mean = 0.0, sd=sqrt(sigma2) )
  
    rg <- rateGm( zeta$alpha, zeta$beta, pressW, delta, nspatial )
    res <- res + dalpha( zeta$alpha, zeta$beta, pressW, delta, nspatial ) / (rg*rg)
  }
  
  fac <- exp( zeta$theta1 + zeta$theta2 * VprodW ) 
  return( - fac * res / 1000.0)
}

## slightly tweeked definition
Edbeta <- function( zeta, VprodW, mpressW, sigma2, delta, nspatial )
{
  res <- 0.0 * mpressW
  pressW <- 0.0 * mpressW
  
  for( n in 1:1000) {
    pressW <- mpressW + rnorm(n=length(mpressW), mean = 0.0, sd=sqrt(sigma2) )
    
    rg <- rateGm( zeta$alpha, zeta$beta, pressW, delta, nspatial )
    res <- res + dbeta( zeta$alpha, zeta$beta, pressW, delta, nspatial ) / (rg*rg)
  }
  
  fac <- exp( zeta$theta1 + zeta$theta2 * VprodW )
  return( - fac * res / 1000.0)
}

newtonstep <- function( zeta, counts, inW, VprodW, mpressW, sigma2, delta, nspatial )
{
  gammam <- rateGm(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
  Egammam <- Elambda(zeta$alpha, zeta$beta, mpressW, sigma2, delta, nspatial)
  lambdam <- exp( zeta$theta1 + VprodW * zeta$theta2 )
  ElambdaW <- lambdam * inW * delta * Egammam
  
  EdalambdaW <- Edalpha( zeta, VprodW, mpressW, sigma2, delta, nspatial ) * inW * delta
  EdblambdaW <- Edbeta( zeta, VprodW, mpressW, sigma2, delta, nspatial ) * inW * delta
  
  da <- dalpha(zeta$alpha, zeta$beta, mpressW, delta, nspatial )
  db <- dbeta(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
  daa <- dalphadalpha(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
  dab <- dalphadbeta(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
  dbb <- dbetadbeta(zeta$alpha, zeta$beta, mpressW, delta, nspatial)
  
  Eclv <- counts - ElambdaW 
  b1 <- sum(Eclv)
  b2 <- sum(VprodW * Eclv)
  b3 <- sum( - (da/gammam) * Eclv )
  b4 <- sum( - (db/gammam) * Eclv )
  b <- c( -b1, -b2, -b3, -b4 ) 
  
  A <- matrix( nrow=4, ncol=4 )
  A[1,] <- c( - sum(ElambdaW), - sum( ElambdaW * VprodW ), 
              - sum(EdalambdaW), - sum( EdblambdaW ) )
  
  A[2,1] <- A[1,2]
  A[2,(2:4)] <- c( - sum( ElambdaW * VprodW * VprodW ), 
                   - sum( EdalambdaW * VprodW ), - sum( EdblambdaW * VprodW ) )
  
  A[3,1] <- sum( (da/gammam) * ElambdaW )
  A[3,2] <- sum( (da/gammam) * VprodW * ElambdaW )
  
  A[3,3] <- sum( ( (da/gammam) * (da/gammam) - (daa/gammam) ) * Eclv ) + sum( (da/gammam) * EdalambdaW )
  A[3,4] <- sum( ( (da/gammam) * (db/gammam) - (dab/gammam) ) * Eclv ) + sum( (da/gammam) * EdblambdaW )
  
  A[4,1] <- sum( (db/gammam) * ElambdaW )
  A[4,2] <- sum( (db/gammam) * VprodW * ElambdaW )
  
  A[4,3] <- sum( ( (da/gammam) * (db/gammam) - (dab/gammam) ) * Eclv ) + sum( (db/gammam) * EdalambdaW )
  A[4,4] <- sum( ( (db/gammam) * (db/gammam) - (dbb/gammam) ) * Eclv ) + sum( (db/gammam) * EdblambdaW )
  
  #print(eigen(A)$values)
  print(b)
  print(A)
  nextzeta <- solve(A, b) 
  zeta1 <- list( theta1 = nextzeta[1] + zeta$theta1,
                 theta2 = nextzeta[2] + zeta$theta2,
                 alpha = nextzeta[3] + zeta$alpha, 
                 beta = nextzeta[4] + zeta$beta )
  return( zeta1 )
}


expLam <- function(zeta, VprodW, mpressW, sigma2, nspatial)
{
  val = mpressW * 0.0
  res = mpressW * 0.0
  dt = length(mpressW)/nspatial
  
  for(i in 1:1000){
    press = mpressW + rnorm(n = length(mpressW), mean = 0.0, sd = sqrt(sigma2))
    mpress = rep(press[1:nspatial],dt) - press  # X_0 - X_t
    val = mpress * exp( zeta$theta1 + zeta$theta2*VprodW + zeta$alpha*mpress)
    res = res + val
  }
  return(res/1000)
}



godambe <- function(zeta, VprodW, mpressW, sigma2, delta, deltaCell, nspatial )
{
  U = matrix(nrow=3, ncol=3)
  V = matrix(nrow=3, ncol=3)
  dt = length(mpressW)/nspatial
  
  mpress = rep(mpressW[1:nspatial], dt) - mpressW
  lam = exp(zeta$theta1 + zeta$theta2*VprodW + zeta$alpha*mpress + (zeta$alpha^2)*sigma2) * delta * deltaCell
  dlam = expLam(zeta, VprodW, mpressW, sigma2, nspatial) * delta * deltaCell
  
  U[,1] = c(sum( lam ), sum( VprodW * lam ), sum( mpress * lam ) )
  U[,2] = c(sum( VprodW * lam ), sum( VprodW * VprodW * lam ), sum( mpress * VprodW *lam ) )
  U[,3] = c(sum( dlam ), sum( VprodW * dlam ), sum( mpress * dlam)  )
  
  V[,1] = U[,1]
  V[,2] = U[,2]
  V[,3] = c(sum( mpress * lam ), sum( mpress * VprodW * lam ), sum( mpress * mpress * lam ) )
  
  Ginv = NA
  Ginv = solve(U)%*%V%*%t(solve(U)) 
  
  return(Ginv)
}

