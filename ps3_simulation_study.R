library(parallel)
num.of.cores <- detectCores()

#library(survival) # don't really need this

#### simulation study
# Lambda_n = Nelson-Aalen estimator
# S0 <- function(t) { pexp(q=t,rate=1,lower.tail=F)}
# G <- function(x) { pnorm(q=x,mean=0,sd=1)}

m <- 1000 # size of empirical average, affects accuracy of H_n(z)
n <- c(200,500,1000) # number of U_i's in \U_n, and O_i affects accuracy of lim \U_n

### generate data
### ----------------
#set.seed(5)
genSimData <- function(n,m) {
  latent.data <- rbind(rexp(n=n,rate=1),rexp(n=n,rate=0.5))
  O <- apply(X=latent.data,MARGIN=2,FUN=function(col){
    t <- col[1]; c <- col[2]
    if (t <= c) {
      return(c(t,1))
    } else {
      return(c(c,0))
    }
  })
  O <- O[,order(O[1,])] # make our lives easier by sorting
  U <- (replicate(n=m,expr=rnorm(n=n,mean=0,sd=1))) # U[,i] = raw data for \U_{i,n}
  list(O=O,U=U)
}
### ---------------

### process data
### --------------
sim.data <- genSimData(n[1],m)
Y <-sim.data$O[1,]; Delta <- sim.data$O[2,]
n.risk <- sapply(X=Y,FUN=function(y.i){sum(Y >= y.i)})
S.hat <- Vectorize(FUN=function(t){mean(Y >= t)},vectorize.args=c("t"))
S.n <- stepfun(x=Y,y=c(1,S.hat(Y)))
### -------------

### build parsimonious grid
### --------------
# J_0 = (1,2)
# based on where we actually see failures

J_0.grid.failure <- Y[(1 <= Y & Y <= 2 & Delta == 1)]
J_0.grid <- unique(c(1,J_0.grid.failure,2))
### --------------

### estimate influence curve and
### asymptotic variance.
### IC.Pn(t;o) = a(t;o) - b(t;o)
### ----------------

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

system.time(IC.mat <- IC.Pn(sim.data$O,J_0.grid)) # 1.316 sec., 7.168 sec
### estimate variance of influence curve
sigma.n <- sqrt(rowMeans(IC.mat^2))
# ----
# sigma.n <- Vectorize(FUN=function(t) {
#   m <- which.min(Y <= t) - 1
#   b[m]
# })
# ----
### check that this is right:
# max(rowMeans(IC.mat)) # 8.08e-17
# emp <- rowMeans(IC.mat^2)
# theor <- sigma.n(J_0.grid)
# max(abs(emp-theor)/emp) # 0.508717
# esystem.time(sigma.n(J_0.grid)) # 0.002sec
# ----

### estimate \U_{k,n}(t)
U <- sim.data$U
IC.mat.std <- diag(1/sigma.n) %*% IC.mat
mathbb.U <- n[1]^(-0.5) * IC.mat.std %*% U

### make H.n
W.n <- apply(X=mathbb.U,MARGIN=2,FUN=function(U){max(abs(U))})
H.n <- Vectorize(function(z) {mean( W.n <= z )})
### find 95%ile
c.alpha <- uniroot(f=function(c) H.n(c) - 0.95,interval=c(0,5))$root

### make Nelson-Aalen estimator
b2 <- cumsum( Delta / n.risk ) 
nelsonAalen <- function(t) {
  sapply(X=t,FUN=function(t){b2[which.min( Y <= t )-1]})
}
Lambda.n <- nelsonAalen(J_0.grid)

### make sigma.diamond and confidence band
sigma.diamond <- sigma.n / Lambda.n

confidence.band <- sapply(X=c(1,-1),FUN=function(x){exp(-exp(log(Lambda.n) + x * c.alpha * sigma.diamond / sqrt(n[1])))})

### check coverage of band
S.0 <- sapply(X=J_0.grid,FUN=function(t){pexp(q=t,rate=1,lower.tail=F)})

S0IsCovered <- function(S.0,confidence.band,J_0.grid) {
  l.band <- confidence.band[,1]; u.band <- confidence.band[,2]
  grid.size <- length(J_0.grid)
  checks <- sapply(X=1:(grid.size-1),FUN=function(time.i){
    check.left <- (l.band[time.i] <= S.0[time.i]) & (S.0[time.i] <= u.band[time.i])
    check.right <- (l.band[time.i] <= S.0[time.i+1]) & (S.0[time.i+1] <= u.band[time.i])
    check.left & check.right
  })
  all(checks)
}
print(S0IsCovered(S.0,confidence.band,J_0.grid))

library(ggplot2)
plot.df <- data.frame(X=J_0.grid,S0=S.0,LBand=confidence.band[,1],HBand=confidence.band[,2])
band.plot <- ggplot(data=plot.df, aes(x=X)) + geom_line(aes(y=S0)) + geom_step(aes(y=LBand),color='red') + geom_step(aes(y=HBand),color='red')
show(band.plot)

### make data generation routine
genSimData <- function(n,m) {
  latent.data <- rbind(rexp(n=n,rate=1),rexp(n=n,rate=0.5))
  O <- apply(X=latent.data,MARGIN=2,FUN=function(col){
    t <- col[1]; c <- col[2]
    if (t <= c) {
      return(c(t,1))
    } else {
      return(c(c,0))
    }
  })
  O <- O[,order(O[1,])] # make our lives easier by sorting
  U <- (replicate(n=m,expr=rnorm(n=n,mean=0,sd=1))) # U[,i] = raw data for \U_{i,n}
  list(O=O,U=U)
}

### encapsulate simulation routine
simulationRoutine <-function(n,m,verbose=F) {
  ### make and process data
  ### --------------
  sim.data <- genSimData(n,m)
  Y <-sim.data$O[1,]; Delta <- sim.data$O[2,]
  n.risk <- sapply(X=Y,FUN=function(y.i){sum(Y >= y.i)})
  S.hat <- Vectorize(FUN=function(t){mean(Y >= t)},vectorize.args=c("t"))
  S.n <- stepfun(x=Y,y=c(1,S.hat(Y)))
  ### -------------
  
  ### build parsimonious grid
  # J_0 = (1,2)
  # based on where we actually see failures
  
  J_0.grid.failure <- Y[(1 <= Y & Y <= 2 & Delta == 1)]
  J_0.grid <- unique(c(1,J_0.grid.failure,2))

  
  ### estimate influence curve and
  ### asymptotic variance.
  ### IC.Pn(t;o) = a(t;o) - b(t;o)
  ### ----------------
  
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
  
  t1 <- system.time(IC.mat <- IC.Pn(sim.data$O,J_0.grid))
  if (verbose) {
    cat(paste("built IC.mat in ", t1[3], " secs.\n", sep="" ))
  }
  
  ### estimate variance of influence curve
  sigma.n <- sqrt(rowMeans(IC.mat^2))
  
  ### estimate \U_{k,n}(t)
  U <- sim.data$U
  IC.mat.std <- diag(1/sigma.n) %*% IC.mat
  mathbb.U <- n^(-0.5) * IC.mat.std %*% U
  
  ### make H.n
  W.n <- apply(X=mathbb.U,MARGIN=2,FUN=function(U){max(abs(U))})
  H.n <- Vectorize(function(z) {mean( W.n <= z )})
  ### find 95%ile
  c.alpha <- uniroot(f=function(c) H.n(c) - 0.95,interval=c(0,5))$root
  
  ### make Nelson-Aalen estimator
  b2 <- cumsum( Delta / n.risk ) 
  nelsonAalen <- function(t) {
    sapply(X=t,FUN=function(t){b2[which.min( Y <= t )-1]})
  }
  Lambda.n <- nelsonAalen(J_0.grid)
  
  ### make sigma.diamond and confidence band
  sigma.diamond <- sigma.n / Lambda.n
  
  confidence.band <- sapply(X=c(1,-1),FUN=function(x){exp(-exp(log(Lambda.n) + x * c.alpha * sigma.diamond / sqrt(n[1])))})
  
  ### check coverage of band
  S.0 <- sapply(X=J_0.grid,FUN=function(t){pexp(q=t,rate=1,lower.tail=F)})
  
  S0IsCovered <- function(S.0,confidence.band,J_0.grid) {
    l.band <- confidence.band[,1]; u.band <- confidence.band[,2]
    grid.size <- length(J_0.grid)
    checks <- sapply(X=1:(grid.size-1),FUN=function(time.i){
      check.left <- (l.band[time.i] <= S.0[time.i]) & (S.0[time.i] <= u.band[time.i])
      check.right <- (l.band[time.i] <= S.0[time.i+1]) & (S.0[time.i+1] <= u.band[time.i])
      check.left & check.right
    })
    all(checks)
  }
  
  coverage.status <- S0IsCovered(S.0,confidence.band,J_0.grid)
  if (verbose) {
    cat(paste("S.0 is covered? ", coverage.status, "\n", sep=""))
    require(ggplot2)
    plot.df <- data.frame(X=J_0.grid,S0=S.0,LBand=confidence.band[,1],HBand=confidence.band[,2])
    band.plot <- ggplot(data=plot.df, aes(x=X)) + geom_line(aes(y=S0)) + geom_step(aes(y=LBand),color='red') + geom_step(aes(y=HBand),color='red')
    show(band.plot)
  }
  
  return(coverage.status)
}

sim.study1 <- replicate(n=100,expr={simulationRoutine(n[1],m,F)})