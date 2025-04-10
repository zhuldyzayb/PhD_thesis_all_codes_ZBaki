
### Estimate the log-likelihhood value
## x-observations
## lambda-predicted values for lambda OR: exp(predict(fit))
## This is the same sunction as dpois function built-in in R.
## However it accomodates forthe case when the value of x is not exactly an integer
estimated_likelihood<-function(x, lambda, log=FALSE){
  if(log==FALSE){res=sum(exp(-lambda)*(lambda^x)/gamma(x+1))}
  else {
    res=sum(log((lambda^x)*exp(-lambda)/gamma(x+1)))
  }
  return(res)
}

## For QuasiPoisson case the qAIC=-2*log-Likelihood/chat+2*k
## k=number of parameters estimated 
## L=log-likelihood value
## chat=dispersion parameter; can be obtained from summary(fit)

estimated_AIC<-function(L, chat,k){
  return(2*k-2*L/chat)
}

AIC <- function( pleModel )
{
  fitXw <- ( pleModel$fitted.values ) * ( pleModel$data[,3] )
  logfitXyXw <- ( log( pleModel$fitted.values ) ) * (pleModel$data[,5] ) * ( pleModel$data[,3] )
  res <- -2 * sum( logfitXyXw - fitXw )
  
  res <- res + 2 * length( pleModel$coefficients )
  return( res )
}
