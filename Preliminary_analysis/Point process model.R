##### This file is created to work on the section 5 of the Log_gaussian paper:
##### namely the Point processes approach
##### functions are saved in files PseudoLikelihood.R and quasiAIC.R
##### last updated  16/02/2021
require(lubridate)
require(Hmisc)
require(boot)
require(ggplot2)
require(Rcpp)

set.seed(123)
dat=read.csv(file.choose())  ### all the eqs within gas field with mag >=1.5
dat=dat[dat$Year>=1995,] ### trim for starting from Jan 1 1995
dat$Date=as.Date(dat$Date, format='%m/%d/%Y')
dat=as.data.frame(dat)
times=difftime(dat$Date, "1994-12-31")
times=as.numeric(times)
times=as.data.frame(times)

####  Pseudo likelihood functions
mybox=list()  ### create a box: for entire time interval (1995-2002)
mybox=list(zrange=c(0, 9497))  
myfrac=3000   ### divide into smaller pieces
PL=pseudo_likelihood(times, mybox, myfrac)  ### obtain both dummy points and actual incidences



#### Once we have obtained the u_points and the corresponding weights 
### create a daily gas production covariate
gasval=read.csv(file.choose()) ## monthly prod values
colnames(gasval)=c("Time", 'Year', 'Month', 'Gas')  ## originally in m^3
gasval=gasval[gasval$Year>=1993,]   ## trim for time interval including 1994 b/c of definition for C(t,24)
gasval=cbind(gasval, as.matrix(rep(NA, nrow(gasval))), as.matrix(rep(NA, nrow(gasval))))
numd1=c(31,28,31,30,31,30,31,31,30,31,30,31)  ## number of days in each month
numd2=c(31,29,31,30,31,30,31,31,30,31,30,31)  ## for leap year
for(i in 1993:2020){
  if(i%%4==0){gasval[gasval$Year==i, 5]=as.matrix(numd2)}
  else {gasval[gasval$Year==i, 5]=as.matrix(numd1)}
}
gasval[,6]=gasval[,4]/gasval[,5]
colnames(gasval)=c("Time", 'Year', 'Month', 'Gas', "NumDays", "DailyGas")  ## gas still in m^3
dayGas=NA  ## daily gas production values
for(i in 1:336){
  if(i==1){dayGas=as.matrix(rep(gasval[i,6], gasval[i,5]))}
  else {dayGas=rbind(dayGas, as.matrix(rep(gasval[i,6], gasval[i,5])))}
}
daily=seq(as.Date("1993-01-01"), as.Date("2020-12-31"), by="days")
dailydata=as.data.frame(cbind(daily, dayGas))  ### each calendar date and average gas prod value in m^3
dailydata$daily=as.Date(daily)
colnames(dailydata)=c('date','gas')


dates=as.Date("1994-12-31")+PL[[1]]$time  ##all the points considered including dummy pts

### define C(t,12) total over last year for each day
gas1=matrix(NA, nrow=3323, ncol=1)
for(i in 1:3323){
  k=sum(dailydata$date<dates[i])
  temp=dailydata[(k-366):(k-1),]
  gas1[i]=sum(temp$gas)
  }

### define C(t,24)
gas2=matrix(NA, nrow=3323, ncol=1)
for(i in 1:3323){
  k=sum(dailydata$date<dates[i])
  temp=dailydata[(k-731):(k-1),]
  gas2[i]=sum(temp$gas)
}

### we need to convert gas values into bcm
gas1=gas1/(10^9)
gas2=gas2/(10^9)



###########################   Perform the modelling  #########################
u_points=PL[[1]]
t=u_points$time
fit1=glm(formula=yy~1+log(t)+gas1, family = quasipoisson(link='log'), data=u_points, weights = w)
summary(fit1)  ##chat=1.914229 (dispersion parameter)
confint(fit1)
AIC1=AIC(fit1)  ##2752.191


####### Estimated intensity for the whole (0, 9497)
time=c(1:9497)
mydays=dailydata$date[731:10227]
Ct12=matrix(NA, nrow=9497, ncol=1)
for(i in 1:9497){
  k=sum(dailydata$date<mydays[i])
  temp=dailydata[(k-366):(k-1),]
  Ct12[i]=sum(temp$gas)
}
Ct12=Ct12/(10^9)
intens=exp(fit1$coefficients[1]+fit1$coefficients[2]*log(time)+fit1$coefficients[3]*Ct12)
dat1=as.data.frame(cbind(Ct12, intens, c(1:9497)))
colnames(dat1)=c('gas','lambda', 'day')
ggplot(data=dat1, aes(x=day, y=gas))+geom_line()+theme(legend.position = "none")+
  labs(x='Days', y='C(t,12)')

ggplot(data=dat1, aes(x=day, y=lambda))+geom_line()+labs(x='Days', y='Estimated intensity')


#### estimate for Lambda(9498) jan 1 2021
xx=sum(dailydata$gas[year(dailydata$date)==2020])/(10^9) ##C(9498,12)
lam21=exp(fit1$coefficients[1]+fit1$coefficients[2]*log(9498)+fit1$coefficients[3]*xx)

lambda=fit1$fitted.values
H=matrix(NA, nrow=3, ncol=3)
H[1,1]=sum(lambda*u_points$w)
H[1,2]=sum(lambda*log(t)*u_points$w)
H[1,3]=sum(lambda*gas1*u_points$w)
H[2,1]=H[1,2]
H[2,2]=sum(lambda*log(t)*log(t)*u_points$w)
H[2,3]=sum(lambda*gas1*log(t)*u_points$w)
H[3,1]=H[1,3]
H[3,2]=H[2,3]
H[3,3]=sum(lambda*gas1*gas1*u_points$w)

F=solve(H)
Sig=F[c(2:3),c(2:3)]
X21=c(log(9498),xx)
s=t(X21)%*%Sig%*%X21
lam21*(1-s*1.96)  
lam21*(1+s*1.96)  


### Replace the C(t,12) with logC(t,12)
fit2=glm(formula=yy~1+log(t)+log(gas1), family = quasipoisson(link='log'), data=u_points, weights = w)
summary(fit2)  ##chat=1.933066
confint(fit2)
AIC2=AIC(fit2) ### 2753.982


#### For updated model with addition of C(t, 24)
fit3=glm(formula=yy~1+log(t)+gas1+gas2, family = quasipoisson(link='log'), data=u_points, weights = w)
summary(fit3)  ##chat=1.908728
confint(fit3)
AIC3=AIC(fit3)
anova(fit3, test='Chisq')




#######################  Residuals     #########################
## we will only consider residuals for the first fit

y=u_points$yy
w=u_points$w
l=fit1$fitted.values
R=matrix(NA, nrow=3323, ncol=1)
for(i in 1:3323){
  if(y[i]==0){R[i]=-w[i]*l[i]}
  else{ R[i]=1-w[i]*l[i]}
}
R.p=R/sqrt(w*l)

R.d=matrix(NA, nrow=3323, ncol=1)
for(i in 1:3323){
  if(y[i]==0) {R.d[i]=-sqrt(2*w[i]*l[i])}
  else {
    R.d[i]=sign(y[i]-l[i])*sqrt(2*w[i]*(y[i]*log(y[i]/l[i])-(y[i]-l[i])))}
  }

## Smooth the residuals using s(x) function  (raw residuals)

sx=function(x,s,lambda,w, sigma){
  p1=sum(exp(-(rep(x,323)-times)^2/(2*sigma^2)))/(sigma*sqrt(2*pi))
  p2=sum(exp(-(rep(x, 3323)-s)^2/(2*sigma^2))*lambda*w)/(sigma*sqrt(2*pi))
  return(p1-p2)
}

res=matrix(NA, nrow=500, ncol=2)
for(i in 1:500){res[i,1]=9407*i/500
  res[i,2]=sx(res[i,1], u_points$time, fit1$fitted.values, u_points$w, sigma=30)}
rdat=as.data.frame(res)
colnames(rdat)=c('Time', 'Res')
ggplot(data = rdat, aes(x=Time, y=Res))+geom_point()+labs(x='Days', y='Smoothed raw residuals')



#### smoothed deviance residuals 
### use info from folder PCF
dyn.load("/Users/zhuldyzaybaki/Desktop/PhD/Cpp_codes/pcfTime.so") 

SmoothResiduals <- function( pleModel, minTime, maxTime, nbins, sigma, typeDev  )
{
  pt <- pleModel$data[,1] - minTime
  npoints <- length(pt)
  if( typeDev=="deviance" ) { wts <- residuals( pleModel, type="deviance" )  }
  else if( typeDev=="pearson" ) { wts <- residuals( pleModel, type="pearson" ) }
  else { wts <- pleModel$data$w * ( pleModel$data$yy - fitted.values(pleModel) ) }
  
  res <- .C("IntensityBins", as.integer(nbins),
            as.integer(npoints), as.double(pt), as.double(wts),
            as.double(maxTime-minTime),
            ans = as.double(numeric(nbins)),  as.double(sigma))
  
  res1 <- seq( 0.5, nbins - 0.5, by=1) * (maxTime-minTime) / nbins
  res1 <- minTime + res1
  res2 <- res$ans
  return( list( time=res1, values=res2) )
}

R.d=SmoothResiduals(fit1, 0, 9497, 500,30,typeDev = 'deviance')
plot(R.d$time,R.d$values)


tempdat=as.data.frame(cbind(R.d$time, R.d$values))
colnames(tempdat)=c('Time', 'Deviance')
ggplot(data = tempdat, aes(x=Time, y=Deviance))+geom_point()+labs(x='Days', y='Smoothed devance residuals')



#########  PCF #######
### actual PCF
cc=fit1$fitted.values[3001:3323]  ### fitted intensities lambda hat for actual point process
pp=PairCorrelationTime(as.matrix(times), 0, 9497, cc, epsilon=90, nh=700, hmax=7000, smooth=1 )
ggplot(data=as.data.frame(pp), aes(x=lag, y=pcf, col='red'))+geom_line(show.legend = FALSE)+labs(x='Lags', y='Estimated pair correlation function')


####### Simulations for PCF #####
## define a matrix for storing a pcf estimates
sim.pcf=matrix(NA, nrow=701, ncol=19)
## simulate a Poisson process with total intensity lambda in (0:9497)
sim.pois=matrix(NA, nrow=9497, ncol=19)
for(i in 1:9497){
  sim.pois[i,]=rpois(19, intens[i])
}
sim.pois=cbind(as.matrix(c(1:9497)), sim.pois)

for(j in 1:19){
ind=sim.pois[(sim.pois[,j+1]!=0), c(1,j+1)]  ### define eqrthquake times
s.times=NULL
for(i in 1:nrow(ind)){
  if(i==1){s.times=as.matrix(rep(ind[i,1], ind[i,2]))}
  else {s.times=rbind(s.times, as.matrix(rep(ind[i,1], ind[i,2])))}
}
### get weights
s.pl=pseudo_likelihood(as.data.frame(s.times), mybox, myfrac)
s.upt=s.pl[[1]]
n=nrow(s.upt)
## find gas production values
s.dat=as.Date("1994-12-31")+s.upt$time  ##all the points considered including dummy pts
s.gas=matrix(NA, nrow=n, ncol=1)
for(i in 1:n){
  k=sum(dailydata$date<s.dat[i])
  temp=dailydata[(k-366):(k-1),]
  s.gas[i]=sum(temp$gas)
}
s.gas=s.gas/(10^9)
sim.fit=glm(formula=yy~1+log(time)+s.gas, family = quasipoisson(link='log'), data=s.upt, weights = w)
s.t=s.upt$time[3001:n]
s.cc=sim.fit$fitted.values[3001:n]
ppt=PairCorrelationTime(as.matrix(s.t), 0, 9497, s.cc, epsilon=90, nh=700, hmax=7000, smooth=1 )
sim.pcf[,j]=ppt$pcf
}
summary(sim.pcf)
dim(sim.pcf)
up=matrix(NA, nrow=701,ncol=1)
down=matrix(NA, nrow=701,ncol=1)
for(i in 1:701){
  xx=sim.pcf[i,]
  up[i]=max(xx)
  down[i]=min(xx)
}

pcf.dat=as.data.frame(cbind(pp$lag, pp$pcf, up, down))
colnames(pcf.dat)=c('lag', 'pcf', 'max', 'min')
         
gg.pcf=ggplot(data=pcf.dat, aes(x=lag, y=pcf))+geom_line(show.legend = FALSE)+labs(x='Lags', y='Estimated pair correlation function')
gg.pcf+geom_line(aes(y=max, col='red'),show.legend = FALSE)+geom_line(aes(y=min, col='red'),show.legend = FALSE)





######## Cross-Validation using MSE ##### (for each year)
PL=pseudo_likelihood(times, mybox, myfrac)  ### obtain both dummy points and actual incidences



#### Once we have obtained the u_points and the corresponding weights 
### create a daily gas production covariate
PL=pseudo_likelihood(times, mybox, myfrac)  ### obtain both dummy points and actual incidences
dates=as.Date("1994-12-31")+PL[[1]]$time  ##all the points considered including dummy pts

### define C(t,12) total over last year for each day
dim(u_points)  ### first 3000 dummy pts and 323 actual observations
length(gas1)  ### 3323
length(t) ### 3323

id1=(rep(FALSE, 3000))
mse=matrix(NA, nrow=26)
for(i in 1:26){
id2=(year(dat$Date)==1994+i)
id=rbind(as.matrix(id1), as.matrix(id2))
cv.fit=glm(formula=yy[-id]~1+log(t[-id])+gas1[-id], family = quasipoisson(link='log'), data=u_points, weights = w[-id])
mse[i]=mean((u_points$yy[id]-exp(cv.fit$coefficients[1]+cv.fit$coefficients[2]*log(t[id])+cv.fit$coefficients[3]*gas1[id]))^2)
}
MSE=mean(mse)  ### 0.4417








         
         
         
