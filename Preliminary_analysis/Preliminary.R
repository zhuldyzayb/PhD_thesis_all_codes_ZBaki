##############. Chapter 3 model fitting.  ##################

setwd("/Users/zhuldyzaybaki/Desktop/Clean code")

require(lubridate)

## Observed eqs:
eqs = read.csv("filtered_eqs.csv")
eqs = eqs[,-c(1,2,3,7,9)]
eqs$Date = as.Date(eqs$Date, format = "%Y-%m-%d")
eqs$Time  = as.numeric(eqs$Date - as.Date("1995-01-01"))+1
T = as.numeric(as.Date("2021-12-31") - as.Date("1995-01-01") ) + 1


## Production volumes:
gas = read.csv2("Groningen gas production.csv")
gas = gas[,-1]
colnames(gas) = c("well", "date", "days", "prod")
gas$date = as.Date(gas$date, format="%d/%m/%Y")
gas = na.omit(gas)

dates = seq(min(gas$date), max(gas$date), by='month')
monthlygas = rep(0, length(dates))
for(i in 1:length(dates)){
  temp = (gas$date==dates[i])
  monthlygas[i] = sum(gas$prod[temp])
}
monthlygas = monthlygas/(10^9)
monthlygas = data.frame(dates, monthlygas)
colnames(monthlygas) = c('date', 'prod')
monthlygas$days = days_in_month(monthlygas$date)
monthlygas$dprod = monthlygas$prod/monthlygas$days 

#################################
#        Count process          #
#################################

anncount = rep(0,27)
anngas = anncount
cumgas = anngas

# N(t)
for(i in 1995:2021){
  anncount[i-1994] = sum(year(eqs$Date)==i)
}
# V(t-1)
for (i in 1995:2021) {
  anngas[i-1994] = sum(monthlygas$prod[year(monthlygas$date)==i-1]) 
}
# \tilde V(t-1)
cumgas = cumsum(anngas)
# log \tilde V(t-1)
cumgas = log(cumgas)


countmod = glm(anncount~anngas+cumgas, family = poisson(link = 'log'))
summary(countmod)
countmod$coefficients  # param estimates
confint(countmod, level = 0.95)   # CI's
lamfit = countmod$fitted.values   # \hat\lambda fitted intensity
presid = (anncount - lamfit)/sqrt(lamfit)   # Pearson residuals
plot(c(1995:2021), presid, ylim = c(-3,3), ylab = "Pearson residuals", xlab = "", pch=19)

# ACF with CIs
calc_autocorr <- function(n, lambda, h) {
  N <- length(n)
  numerator <- sum(n[(h+1):N] * n[1:(N-h)] / (lambda[(h+1):N] * lambda[1:(N-h)]))
  denominator <- N - h
  return(numerator / denominator)
}

h_range <- 0:20
autocorr <- sapply(h_range, function(h) calc_autocorr(anncount, lamfit, h))   # Observed data
# Simulations
NN <- t(sapply(1:27, function(i) rpois(19, lamfit[i])))
# simulated ACFs
acfmat <- sapply(1:19, function(i) sapply(h_range, function(h) calc_autocorr(NN[,i], lamfit, h)))
# plot ACF
plot(autocorr, ylim = c(0,2), pch = 19, col='red', xlab = "h", ylab = "Empirical autocorrelation function")
points(apply(acfmat, 1, max), pch=19)
points(apply(acfmat, 1, min), pch=19)
# \Sigma = (Cov_\alpha)^(-1);    Cov_\alpha = (-Hessian)

## production in 2021: 6.482033; 2022: 4.559024;  2023:  1.460676
v1 = c(6.482033, 4.559024)
v2 = c( cumsum(anngas)[27]+v1[1], cumsum(anngas)[27]+v1[1]+v1[2])
S = vcov(countmod)

XX = cbind(as.matrix(rep(1, 27)), anngas, cumgas)## first 27 years
Xp = c(1,v1[1], log(v2[1]), 1, v1[2], log(v2[2]))
Xp = matrix(Xp, nrow=3, byrow=TRUE)
XX = rbind(XX,Xp )

sigms = rep(0, 29)
for(i in 1:29){
  sigms[i] = t(XX[i,])%*%S%*%XX[i,]
}

lampred = exp( countmod$coefficients[1] + countmod$coefficients[2]*v1 + countmod$coefficients[3]*log(v2))
lambda = rbind(as.matrix(lamfit), as.matrix(lampred))

upper = lambda*(1 + sqrt(sigms)*1.96)
lower = lambda*(1 - sqrt(sigms)*1.96)

tt = seq(1995, 2023)
par(mar=c(3,3,2,2))
plot(lambda, type = "p", col = "black", pch = 16, xlab = "", ylab = "", ylim=c(0,35), axes=FALSE)
axis(side=1, at = seq(1,29,5), labels = seq(1995, 2023,5))
axis(side=2, at = seq(0,35,5), labels = seq(0, 35,5))
box()
points(length(lambda) - 1, lambda[length(lambda) - 1], col = "red", pch = 16)
points(length(lambda), lambda[length(lambda)], col = "red", pch = 16)
points(lower[1:29], col='grey', pch=16, type='o')
points(upper[1:29], col='grey', pch=16, type='o')


## Cross-validation
CV = 0
for(i in 1:27){
  mod = glm(anncount[-i]~ anngas[-i] + cumgas[-i], family = poisson(link = "log"))
  est = exp( mod$coefficients[1] + mod$coefficients[2] * anngas[i] + mod$coefficients[3]*cumgas[i])
  CV = CV + (anncount[i] - est)^2 
}
CV = CV/27





#################################
#       Point process           #
#################################
win <- owin(c(0, T-1), c(0, 1))  
events <- ppp(x = eqs$Time, y = rep(0, length(eqs$Time)), window = win)
# Covariates
grid_points <- seq(0, T-1, length.out = 3000)

prod_upto_date = function(date, prod_data){
  target_date = as.Date(date, format="%Y-%m-%d")
  start_date = target_date - 365
  
  # in between months
  tot_bwn = sum( prod_data$prod[ (prod_data$date > start_date) & (prod_data$date<= (target_date - months(1)) ) ] )
  # current month                   
  tot_td = prod_data$dprod[prod_data$date==floor_date(target_date, "month")] * day(target_date)   
  # starting  month
  tot_sd = prod_data$dprod[prod_data$date==floor_date(start_date, "month")] * as.numeric(( days_in_month(start_date) - day(start_date)))
  return(tot_bwn+tot_td+tot_sd)
}

cum_prod_upto_date = function(date, prod_data){
  target_date = as.Date(date)
  start_date = as.Date("1995-01-01")
  
  # from Jan 1995
  tot_bf = sum( prod_data$prod[ prod_data$date >=start_date & prod_data$date <= (target_date - months(1)) ] )
  # current month                   
  tot_td = prod_data$dprod[prod_data$date==floor_date(target_date, "month")] * day(target_date)   
  return(tot_bf+tot_td)
}

V1 = grid_points*0
V2 = grid_points*0
for(i in 1:length(grid_points)){
  d = as.Date("1995-01-01") + grid_points[i]
  V1[i] = prod_upto_date(d, monthlygas)
  V2[i] = cum_prod_upto_date(d, monthlygas)
}


x = seq(0,3, 0.001)

y1 = x * exp( - exp( -2 + 1.5 * x)) 
y2 = x * exp( - exp( 0 + 1 * x)) 
y3 = x * exp( - exp( 0.5 + 0.5 * x)) 


par(mar = c(5,5,4,2))
plot(x, y1, xlab = "x", ylab = "f(x)", type = 'l', cex = 0.5, cex.lab = 2)
lines(x, y2, lty = 2)
lines(x, y3, lty = 3)

