library(parallel)
num.of.cores <- detectCores()

library(survival)

#### simulation study
# Lambda_n = Nelson-Aalen estimator
S0 <- function(t) { pexp(q=t,rate=1,lower.tail=F)}
G <- function(x) { pnorm(q=x,mean=0,sd=1)}
# J_0 = (1,2)
grid.res <- 100
J_0.grid <- seq(from=1,to=2,length.out=grid.res)

m <- 1000 # size of empirical average, affects accuracy of H_n(z)
n <- c(200,500,1000) # number of U_i's in \U_n, and O_i affects accuracy of lim \U_n

### generate data
set.seed(5)
latent.data <- rbind(rexp(n=n[1],rate=1),rexp(n=n[1],rate=0.5))
O <- apply(X=latent.data,MARGIN=2,FUN=function(col){
  t <- col[1]; c <- col[2]
  if (t <= c) {
    return(c(t,1))
  } else {
    return(c(c,0))
  }
})
#O <- t(O) # O = (Y,Delta)
O <- O[,order(O[1,])] # make our lives easier by sorting

U <- (replicate(n=m,expr=rnorm(n=n[1],mean=0,sd=1))) # U[,i] = raw data for \U_{i,n}

### process data
Y <- O[1,]; Delta <- O[2,]
n.risk <- sapply(X=Y,FUN=function(y.i){sum(Y >= y.i)})
Pn <- ecdf(Y)
S.hat <- function(t){1-Pn(t)}

### estimate influence curve
IC.Pn <- function(o) {
  y <- o[1]; delta <- o[2]
  function(t) {
    sapply(X=t,FUN=function(t) {
      a <- delta*as.numeric( 0 <= y & y <= t)/S.hat(y)
      m <- which.min(Y <= min(t,y)) - 1
      b <- sum( Delta[1:m] / (S.hat(Y[1:m])*n.risk[1:m]), na.rm=T )  
      return(a-b)
    })
  }
}

IC.mat <- apply(X=O,MARGIN=2,FUN=function(o){IC.Pn(o)(J_0.grid)})
### estimate variance of influence curve
sigma.n <- function(t) {
  sapply(X=t,FUN=function(t){sum( Delta * as.numeric(Y <= t) / (S.hat(Y)*n.risk) )})
}

### estimate \U_{k,n}(t)
mathbb.U <- function(k,n) {
  function(t) {
    a <- apply(X=O,MARGIN=1,FUN=function(o){IC.Pn(o)(t)/sigma.n(t)})
    sum(U[,k]*a)/sqrt(n)
  }
}

### build vector of \sup|\U_k's|
system.time(expr={sup.U.n <- mclapply(X=1:100,mc.cores=num.of.cores,FUN=function(k){
  max(abs(sapply(X=J_0.grid,FUN=function(t){mathbb.U(k,n[1])(t)})))
})})


H.n <- Vectorize(vectorize.args=c("z"),function(z,n) {
  mean( sapply(X=1:m, FUN=function(k){max(abs(mathbb.U(t=J_0.grid,k,n))) <= z }) )
})

             
                 