library(parallel)
num.of.cores <- detectCores()

library(survival)

#### simulation study
# Lambda_n = Nelson-Aalen estimator
S0 <- function(t) { pexp(q=t,rate=1,lower.tail=F)}
G <- function(x) { pnorm(q=x,mean=0,sd=1)}

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
S.hat <- Vectorize(FUN=function(t){mean(Y >= t)},vectorize.args=c("t"))

### build grid
# J_0 = (1,2)
grid.res <- 100
J_0.grid.uniform <- seq(from=1,to=2,length.out=grid.res)
J_0.grid.failure <- Y[(1 <= Y & Y <= 2 & Delta == 1)]
J_0.grid <- sort(unique(c(J_0.grid.uniform,J_0.grid.failure)))


### estimate influence curve
### IC.Pn(t;o) = a(t;o) - b(t;o)

b <- cumsum( Delta / (S.hat(Y)*n.risk) )
IC.Pn <- function(o,t) {
  y <- o[1]; delta <- o[2]
  sapply(X=t,FUN=function(t){
    a <- delta * as.numeric( 0 <= y & y <= t) / S.hat(y)
    m <- which.min(Y <= min(t,y)) - 1
    return(a-b[m])
  })
}

# ---------------------------
# IC.Pn2 <- function(o) {
#   y <- o[1]; delta <- o[2]
#   function(t) {
#     a <- delta * as.numeric(0 <= y & y <= t) / S.hat(y)
#     b <- sum(sapply(X=1:n[1],FUN=function(i){
#       Delta[i] * as.numeric(Y[i] <= min(t,y)) / (n[1] * S.hat(Y[i])^2 )
#     }),na.rm=T)
#     a-b
#   }
# }
# -------------------------

system.time(IC.mat <- apply(X=O,MARGIN=2,FUN=function(o){IC.Pn(o,J_0.grid)})) # 3.7 sec.

### estimate variance of influence curve
sigma.n <- Vectorize(FUN=function(t) {
  m <- which.min(Y <= t) - 1
  b[m]
})

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
  max(abs(sapply(X=J_0.grid,FUN=function(t){mathbb.U(k,n[1])(t)})))
})})


H.n <- Vectorize(vectorize.args=c("z"),function(z,n) {
  mean( sapply(X=1:m, FUN=function(k){max(abs(mathbb.U(t=J_0.grid,k,n))) <= z }) )
})

             
                 