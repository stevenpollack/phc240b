# Use the exponential, Weibull and gamma distributions to model the time from randomization until death in patients assigned to the experimental treatment group. Report point estimates and approximate 95% confidence intervals for each involved parameter.
require(survival)
require(flexsurv)
require(ggplot2)

testPatients <- subset(veteran,trt==2)
testPatSurvObj <- with(data=testPatients,expr={Surv(time,status,type="right")})

### exponential fit
expFit <- flexsurvreg(testPatSurvObj ~ 1, dist="exp")
expFit$res

### weibull fit
weibullFit <- flexsurvreg(testPatSurvObj ~ 1, dist="weibull")
weibullFit$res

### gamma fit
gammaFit <- flexsurvreg(testPatSurvObj ~ 1, dist="gamma")
gammaFit$res

### generalized gamma fit
genGammaFit <- flexsurvreg(testPatSurvObj ~ 1, dist="gengamma")
genGammaFit$res

plot(expFit,ci=F,col="red")
lines(weibullFit,ci=F,col="green")
lines(gammaFit,ci=F,col="yellow")

# Discuss the adequateness of the exponential, Weibull, gamma and generalized gamma distributions using the likelihood ratio test.

modelLikelihoods <- c(genGammaFit$loglik, gammaFit$loglik, weibullFit$loglik, expFit$loglik)

lrtMat <- outer(modelLikelihoods,modelLikelihoods,FUN=function(x,y){-2*(x-y)})
colnames(lrtMat) <- rownames(lrtMat) <- c("Gen Gamma","Gamma", "Weibull", "Exponential")

# In the Diabetic Retinopathy Study (RDS), patients each had one (random) eye randomized to either photocoagulation therapy or follow-up without treatment. Visual acuity was measured prospectively. A subset of the RDS dataset is available in the timereg package and is referred to as diabetes. Use the Gompertz distribution to analyze the time from administration of photocoagulation therapy until blindness in treated patients. Provide point estimates and approximate 95% confidence intervals for both parameters. Is the exponential distribution an adequate simplification of the Gompertz distribution in this case? Justify your answer formally.


library(flexsurv)

### simulation study
n <- 100
rate <- 12.5

data <- replicate(n=2000,expr=rexp(n=n,rate=rate)) # columns represent experiments, rows, individual death times

sorted.data <- apply(X=data,MARGIN=2,FUN='sort')

# ithOrderCDF <- function(i,rate,n,x) {
#   outer.summand <- sapply(X=i:n,FUN=function(k){
#     inner.summand <- (-1)^(0:k)*choose(k,0:k)*exp(-rate*x*(n-k+0:k))
#     choose(n,k)*sum(inner.summand)
#   })
#   sum(outer.summand)
# }
# 
# ithOrderCDF <- Vectorize(ithOrderCDF,vectorize.args=c("i","x"))

ithOrderCDF2 <- function(i,rate,n,x) {
  sum(choose(n,i:n)*pexp(q=x,rate=rate)^(i:n)*pexp(q=x,rate=rate,lower.tail=F)^(n-(i:n)))
}
  
ithOrderCDF2 <- Vectorize(ithOrderCDF2,vectorize.args=c("i","x"))

library(ggplot2)

# look at 17th order statistic
ggplot.df <- data.frame(X=sorted.data[17,])
ggplot.df2 <- data.frame(X=seq(from=min(ggplot.df),to=max(ggplot.df),length.out=1000))
ggplot.df2 <- transform(ggplot.df2, Y=ithOrderCDF2(i=17,rate=rate,n=n,x=X))

ggplot() + stat_ecdf(data=ggplot.df, aes(x=X)) + geom_line(data=ggplot.df2, aes(x=X,y=Y), colour="red")