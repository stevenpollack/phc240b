library(parallel)
num.of.cores <- detectCores()

#library(survival) # don't really need this

#### simulation study
# Lambda_n = Nelson-Aalen estimator
S0 <- function(t) { pexp(q=t,rate=1,lower.tail=F)}
G <- function(x) { pnorm(q=x,mean=0,sd=1)}

m <- 1000 # size of empirical average, affects accuracy of H_n(z)
n <- c(200,500,1000) # number of U_i's in \U_n, and O_i affects accuracy of lim \U_n

### generate data
set.seed(5)
latent.data <- rbind(rexp(n=n[2],rate=1),rexp(n=n[2],rate=0.5))
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

U <- (replicate(n=m,expr=rnorm(n=n[2],mean=0,sd=1))) # U[,i] = raw data for \U_{i,n}

### process data
Y <- O[1,]; Delta <- O[2,]
n.risk <- sapply(X=Y,FUN=function(y.i){sum(Y >= y.i)})
S.hat <- Vectorize(FUN=function(t){mean(Y >= t)},vectorize.args=c("t"))
S.n <- stepfun(x=Y,y=c(1,S.hat(Y)))

### build parsimonious grid
# J_0 = (1,2)
# based on where we actually see failures

# grid.res <- 100 # deprecated
# J_0.grid.uniform <- seq(from=1,to=2,length.out=grid.res) # deprecated

J_0.grid.failure <- Y[(1 <= Y & Y <= 2 & Delta == 1)]
J_0.grid <- unique(c(1,J_0.grid.failure,2))
# check grid: plot(J_0.grid,rep(x=1,times=length(J_0.grid))

### estimate influence curve
### IC.Pn(t;o) = a(t;o) - b(t;o)

b <- cumsum( Delta / (S.n(Y)*n.risk) )
IC.Pn <- function(O,t) {
  ### helper function
  IC.1 <- function(t,y,delta) {
    sapply(X=t,FUN=function(t){
      a <- delta * as.numeric( 0 <= y & y <= t) / S.n(y)
      m <- which.min(Y <= min(t,y)) - 1
      return(a-b[m])
    })
  }
           
  ### logic to handle overloading
  if (!is.null(dim(O))) { # passing matrix of O's
    apply(X=O,MARGIN=2,FUN=function(o,times){IC.1(times,y=o[1],delta=o[2])},times=t)
  } else { # passing just a single observation
    IC.1(t,y=O[1],delta=O[2])
  }
}

system.time(IC.mat <- IC.Pn(O,J_0.grid)) # 1.316 sec.

### estimate variance of influence curve
sigma.n <- Vectorize(FUN=function(t) {
  m <- which.min(Y <= t) - 1
  b[m]
})

### check that this is right:
emp <- apply(X=IC.mat,MARGIN=1,FUN=function(row){mean(row^2)})
theor <- sigma.n(J_0.grid)
max(abs(emp-theor))
system.time(sigma.n(J_0.grid)) # 0.002sec

### estimate \U_{k,n}(t)
mathbb.U <- function(k,n) {
  IC.Pn.k <- IC.mat[,k]
  sapply(X=t,FUN=function(t0){
    a <- apply(X=O,MARGIN=2,FUN=function(o){IC.Pn(o,t0)}/sigma.n(t0)})
    sum(U[,k]*a)/sqrt(n)  
  })
}

### build vector of \sup|\U_k's|
system.time(expr={sup.U.n <- mclapply(X=1:100,mc.cores=num.of.cores,FUN=function(k){
  max(abs(sapply(X=J_0.grid,FUN=function(t){mathbb.U(k,n[2])(t)})))
})})


H.n <- Vectorize(vectorize.args=c("z"),function(z,n) {
  mean( sapply(X=1:m, FUN=function(k){max(abs(mathbb.U(t=J_0.grid,k,n))) <= z }) )
})

             
                 