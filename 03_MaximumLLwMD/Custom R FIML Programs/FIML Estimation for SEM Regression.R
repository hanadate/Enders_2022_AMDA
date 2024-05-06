options(scipen = 999)
library(fdir)
library(mvtnorm)

# set working directory to this script's location
fdir::set()

# read data
dat <- read.table("pain.dat")
names(dat) <- c("id", "male", "age", "edugroup", "workhrs", 
                "exercise", "pain", "anxiety", "stress", "control", 
                "interfere", "depress", "disability")

# data and missingness indicator
Y <- cbind(dat$control, dat$depress)
R <- Y[,2] == 999

# evaluate first derivatives at current parameter values
deriv.1 <- function(mu, sigma){
  
  # initialize derivative matrices
  d.mu <- matrix(0, nrow = 2, ncol = 1)
  d.sigma <- rep(0,ncol(duplication.matrix(nrow(sigma))))
  
  # sum individual contributions to derivatives
  for(i in 1:N){
    
    # select elements of data and parameter matrices that correspond to observed data
    if(R[i]){tau <- s <- matrix(c(1,0), nrow = 1)} else {tau <- s <- diag(2)}
    Y.i <- s %*% Y[i,]
    mu.i <- s %*% mu
    sigma.i <- s %*% sigma %*% t(s)
    
    # compute individual contributions to derivatives
    d.mu <- d.mu + t(tau) %*% solve(sigma.i) %*% ((Y.i - mu.i))
    d.sigma <- d.sigma + t(duplication.matrix(nrow(sigma))) %*% vec(-.5*(t(tau) %*% solve(sigma.i) %*% tau - t(tau) %*% solve(sigma.i) %*% matrix((Y.i - mu.i), ncol = 1) %*% t(matrix((Y.i - mu.i), ncol = 1)) %*% solve(sigma.i) %*% tau))

    }
  
  return(list(d.mu, d.sigma))
  
}

# evaluate second derivatives at current parameter values
deriv.2 <- function(mu, sigma){
  
  # initialize derivative matrices
  d.mu <- matrix(0, nrow = length(mu), ncol = length(mu))
  d.sigma <- matrix(0, nrow = ncol(duplication.matrix(nrow(sigma))), ncol = ncol(duplication.matrix(nrow(sigma))))
  d.mu.sigma <- matrix(0, nrow = length(mu), ncol = ncol(duplication.matrix(nrow(sigma))))
  
  # sum individual contributions to derivatives
  for(i in 1:N){
    
    # select elements of data and parameter matrices that correspond to observed data
    if(R[i]){tau <- s <- matrix(c(1,0), nrow = 1)} else {tau <- s <- diag(2)}
    Y.i <- s %*% Y[i,]
    mu.i <- s %*% mu
    sigma.i <- s %*% sigma %*% t(s)
    
    # compute individual contributions to derivatives
    d.mu <- d.mu + t(tau) %*% solve(sigma.i) %*% tau
    d.sigma <- d.sigma + t(duplication.matrix(nrow(sigma))) %*% (t(tau) %*% solve(sigma.i) %*% tau %x% (t(tau) %*% solve(sigma.i) %*% matrix((Y.i - mu.i), ncol = 1) %*% t(matrix((Y.i - mu.i), ncol = 1)) %*% solve(sigma.i) %*% tau  - .5*t(tau) %*% solve(sigma.i) %*% tau)) %*% duplication.matrix(nrow(sigma))
    d.mu.sigma <- d.mu.sigma + (t(tau) %*% solve(sigma.i) %*% tau %x% (t(matrix((Y.i - mu.i), ncol = 1)) %*% solve(sigma.i) %*% tau)) %*% duplication.matrix(nrow(sigma))
  
  }
  
  return(list(-d.mu, -d.sigma, -d.mu.sigma))
  
}

# constants
nvar <- ncol(Y)
nparms <- nvar + (nvar + 1) * nvar / 2
N <- nrow(Y)

# initialize algorithmic features
stopcriterion <- .0000001
stepsize <- 1
maxiterations <- 300

# starting values
# see model-implied moments in Eq. 3.25
g0 <- 5; var.r <- 1; b0 <- 5; b1 <- 0; var.e <- 1
alpha <- c(g0, b0)
beta <- matrix(c(0,b1,0,0), ncol = nvar)
psi <- diag(c(var.r,var.e))

# initialize matrices to hold iterative history
params <- matrix(c(g0,b0,b1,var.r,var.e), nrow = 1) %x% matrix(1, nrow = maxiterations, ncol = 1)
logL <- rep(0,maxiterations)

# begin newton algorithm
iteration <- 1
maxdiff <- 1
while (maxdiff > stopcriterion) {

  # advance iteration index
  iteration <- iteration + 1
  print(paste0("iteration = ", iteration))
  
  # current sem/regression parameters
  alpha <- params[iteration-1,1:nvar]
  beta <- matrix(c(0,params[iteration-1,nvar+1],0,0), ncol = nvar)
  psi <- diag(params[iteration-1,(nparms-1):nparms])
  
  # current model-implied mean vector and covariance matrix
  mu <- solve(diag(nvar) - beta) %*% alpha
  sigma <- solve(diag(nvar) - beta) %*% psi %*% t(solve(diag(nvar) - beta))
  
  # log-likelihood evaluated at the current estimates
  logL[iteration-1] <- sum(dmvnorm(Y[R == F,],mu,sigma, log = T)) + sum(dnorm(Y[R == T,1],mu[1],sqrt(sigma[1,1]), log = T))
  
  # construct derivative matrices
  first.derivs <- deriv.1(mu,sigma)
  second.derivs  <- deriv.2(mu,sigma)
  gradient <- c(first.derivs[[1]], first.derivs[[2]])
  zeroblock <- matrix(0, nrow = nrow(second.derivs[[3]]), ncol = ncol(second.derivs[[3]]))
  hessian <- rbind(cbind(second.derivs[[1]], zeroblock), cbind(t(zeroblock),second.derivs[[2]]))
  
  # derivatives linking model-implied moments to regression parameters (see Table 3.5)
  delta <- matrix(0, nrow = nparms, ncol = nparms)
  delta[1,1] <- delta[2,2] <- delta[3,4] <- delta[5,5] <- 1
  delta[2,1] <- beta[2,1]
  delta[2,3] <- alpha[1]
  delta[4,3] <- psi[1,1]
  delta[4,4] <- beta[2,1]
  delta[5,3] <- 2*psi[1,1]*beta[2,1]
  delta[5,4] <- beta[2,1]^2
  
  # update parameters
  params[iteration,] <- params[iteration-1,] - stepsize * solve(t(delta) %*% hessian %*% (delta)) %*% (t(delta) %*% gradient)

  # calculate max difference between consecutive parameters to assess convergence
  if(iteration < maxiterations){
    maxdiff <- max(abs(params[iteration,] - params[iteration-1,]))
  } else {maxdiff <- stopcriterion}

  # compute final log-likelihood and standard errors after convergence
  if(maxdiff < stopcriterion){
    
    # final sem/regression parameters
    alpha <- params[iteration,1:nvar]
    beta <- matrix(c(0,params[iteration,nvar+1],0,0), ncol = nvar)
    psi <- diag(params[iteration,(nparms-1):nparms])
    
    # final model-implied moments
    mu <- solve(diag(nvar) - beta) %*% alpha
    sigma <- solve(diag(nvar) - beta) %*% psi %*% t(solve(diag(nvar) - beta))
    
    # log likelihood evaluated at the final model-implied estimates
    logL[iteration] <- sum(dmvnorm(Y[R == F,],mu,sigma, log = T)) + sum(dnorm(Y[R == T,1],mu[1],sqrt(sigma[1,1]), log = T))
    
    # final derivatives linking model-implied moments to regression parameters
    delta <- matrix(0, nrow = nparms, ncol = nparms)
    delta[1,1] <- delta[2,2] <- delta[3,4] <- delta[5,5] <- 1
    delta[2,1] <- beta[2,1]
    delta[2,3] <- alpha[1]
    delta[4,3] <- psi[1,1]
    delta[4,4] <- beta[2,1]
    delta[5,3] <- 2*psi[1,1]*beta[2,1]
    delta[5,4] <- beta[2,1]^2
    
    # standard errors with observed information
    second.derivs  <- deriv.2(mu,sigma)
    hessian <- rbind(cbind(second.derivs[[1]], second.derivs[[3]]), cbind(t(second.derivs[[3]]), second.derivs[[2]]))
    information <- t(delta) %*% hessian %*% (delta)
    cov.params <- solve(-information)
    se.params <- sqrt(diag(cov.params))
    
  }
  
}

# summarize
summary <- cbind(params[iteration,c(1,4,2,3,5)], se.params[c(1,4,2,3,5)])
logL <- cbind(seq(1:(iteration))-1, logL[1:iteration], params[1:iteration,c(1,4,2,3,5)])
rownames(summary) <- c("muX","varX","B0","B1","vare")
colnames(summary) <- c("Est.", "Std. Err.")
colnames(logL) <- c("Iteration", "logL", rownames(summary))
print(paste0("Estimation Summary after ", iteration-1, " Newton Iterations:"))
print(round(summary, digits = 3))
print(paste0("Iteration History:"))
print(round(logL, digits = 15))
