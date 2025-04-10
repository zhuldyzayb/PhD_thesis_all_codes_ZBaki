require(Rcpp)
### use info from folder PCF
dyn.load("/Users/zhuldyzaybaki/Desktop/PhD/Cpp_codes/pcfTime.so") 
PairCorrelationTime <- function( X, minTime, maxTime, lambdaX, epsilon, nh, hmax, smooth )
{
  # shift points and window to fit in [0,T]
  # X is space time point pattern
  # nh the number of lags for the pcf,
  # between 0 and hmax
  # smooth 0: box kernel, otherwise Gaussian kernel
  
  pt <- X - minTime
  npoints <- length(X)
  coef <- lambdaX
  #KernelIntensityTime( X, minTime, maxTime, sigma )
  
  if( smooth ) {
    res <- .C("pcfTimeSmooth", as.integer(nh), as.double(hmax),
              as.integer(npoints), as.double(pt),
              as.double(maxTime-minTime), as.double(coef),
              ans = as.double(numeric(nh+1)),  as.double(epsilon))}
  else {
    res <- .C("pcfTime", as.integer(nh), as.double(hmax),
              as.integer(npoints), as.double(pt),
              as.double(maxTime-minTime), as.double(coef),
              ans = as.double(numeric(nh+1)),  as.double(epsilon))
  }
  res1 <- seq( 0, hmax, by = hmax/nh )
  res2 <- res$ans
  return( list( lag=res1, pcf=res2) )
}
