# load packages
library(fdir)
library(mvtnorm)
library(MCMCpack)
library(truncnorm)

# set working directory to this script's location
fdir::set()

# read data
dat <- read.table("employeecomplete.dat")
names(dat) = c("employee", "team", "turnover", "male", "empower", "lmx", "jobsat", "climate", "cohesion")
N <- nrow(dat)

# select variables
Y <- dat$jobsat
X <- dat$lmx

# assemble X matrix
X <- cbind(1, X)

# initialize algorithmic features
iterations <- 11000
burnin <- 5000
proposalsd0 <- .035
proposalvar_multiplier <- 1

# initial threshold estimates
highcat <- as.numeric(max(names(table(Y))))
lowcat <- as.numeric(min(names(table(Y))))
proportions <- table(Y)/sum(table(Y))
tau0 <- rep(0,highcat - 1)
for(c in 1:(highcat-1)){tau0[c] <- qnorm(sum(proportions[1:c]),0,1)}
tau0 <- c(-Inf, tau0 - tau0[1], Inf)
  
# store parameter estimates
Ystar <- rep(0, length(Y))
betas <- matrix(0, ncol = 2, nrow = iterations); betas[,1] <- 2
thresholds <- matrix(rep(tau0, each = iterations), nrow = iterations)
threshaccept <- matrix(0, nrow = iterations, ncol = 1)

# gibbs sampler
for(t in 2:iterations){

  if(t%%100 == 0){
    print(paste0("Iteration = " , t))
  }
  
  # tune MH acceptance rates
  tunecheck <- 100
  if(t%%tunecheck == 0 & t <= burnin){
    if(mean(threshaccept[(t-(tunecheck-1)):t]) > .50){proposalvar_multiplier <- proposalvar_multiplier * 1.25}
    if(mean(threshaccept[(t-(tunecheck-1)):t]) < .25){proposalvar_multiplier <- proposalvar_multiplier * .75}
    print(paste0("Iteration = ", t, " | Acceptance rate for Last 50 Adjustment Reps = ", round(mean(threshaccept[(t-(tunecheck-1)):t]),2)))
  }
  
  # draw candidate thresholds from proposal distribution
  thresh_old <- thresholds[t-1,]
  sdMH <- proposalsd0 * proposalvar_multiplier
  thresh_new <- thresh_old
  for(c in 3:highcat){
    thresh_new[c] <- rtruncnorm(1, a = thresh_new[c-1], b = thresh_old[c+1], mean = thresh_old[c], sd = sdMH)
  }
  
  # compute predicted latent scores
  xBi <- X %*% betas[t-1,]

  # compute likeihood part of IR
  ll <- rep(0, length(Y))
  for(i in 1:N){
    ll[i] <- (pnorm(thresh_new[Y[i]+1], xBi[i]) - pnorm(thresh_new[Y[i]], xBi[i])) / (pnorm(thresh_old[Y[i]+1], xBi[i]) - pnorm(thresh_old[Y[i]],xBi[i]))
  } 
  
  # compute normalization offset part of IR
  offset <- rep(1, length(thresh_new))
  for(c in 3:highcat){
    offset[c] <- (pnorm(thresh_old[c+1],thresh_old[c],sdMH) - pnorm(thresh_new[c-1],thresh_old[c],sdMH)) / 
      (pnorm(thresh_new[c+1],thresh_new[c],sdMH) - pnorm(thresh_old[c-1],thresh_new[c],sdMH))
  }
  
  # importance ratio
  ir <- prod(ll) * prod(offset)
  u <- runif(1)
  if(ir > u){
    threshaccept[t] <- 1
    thresholds[t,] <- thresh_new
  } else{
    thresholds[t,] <- thresh_old
  }
  
  # sample latent scores
  for(c in lowcat:highcat){
    Ystar[Y == c] <- rtruncnorm(length(Y[Y == c]), a = thresholds[t,c], b = thresholds[t,c+1], mean = xBi[Y == c], sd = 1)
  }
  
  # sample analysis model coefficients
  betahat <- solve(crossprod(X,X)) %*% crossprod(X,Ystar)
  covbeta <- solve(crossprod(X,X))
  betas[t,] <- rmvnorm(1, mean = betahat, sigma = covbeta)

}


# summarize posterior distributions
regparams <- matrix(0, nrow = 2, ncol = 5)
rownames(regparams) <- c("B0", "B1")
colnames(regparams) <- c("Mean", "StdDev", "2.5%", "50%", "97.5%")

# summarize posterior distribution of the slope
regparams[2,1] <- mean(betas[(burnin+1):iterations,2])
regparams[2,2] <- sd(betas[(burnin+1):iterations,2])
regparams[2,3:5] <- quantile(betas[(burnin+1):iterations,2], c(.025, .50, .975))

# summarize posterior distributions
threshparams <- matrix(0, nrow = 6, ncol = 5)
rownames(threshparams) <- c("tau1","tau2","tau3","tau4","tau5","tau6")
colnames(threshparams) <- c("Mean", "StdDev", "2.5%", "50%", "97.5%")

# summarize posterior distribution of the mean
threshparams[,1] <- colMeans(thresholds[(burnin+1):iterations,2:(highcat)])
threshparams[,2] <- sqrt(diag(cov(thresholds[(burnin+1):iterations,2:highcat])))
threshparams[1,3:5] <- quantile(thresholds[(burnin+1):iterations,2], c(.025, .50, .975))
threshparams[2,3:5] <- quantile(thresholds[(burnin+1):iterations,3], c(.025, .50, .975))
threshparams[3,3:5] <- quantile(thresholds[(burnin+1):iterations,4], c(.025, .50, .975))
threshparams[4,3:5] <- quantile(thresholds[(burnin+1):iterations,5], c(.025, .50, .975))
threshparams[5,3:5] <- quantile(thresholds[(burnin+1):iterations,6], c(.025, .50, .975))
threshparams[6,3:5] <- quantile(thresholds[(burnin+1):iterations,7], c(.025, .50, .975))

print(paste0("Posterior Distribution Summary from ", iterations - burnin, " Iterations:"))
summary <- rbind(regparams,threshparams)
print(round(summary,3))

# kernel density plots of posteriors
params <- cbind(betas,thresholds[,2:7])
for(p in 1:nrow(summary)){
  plot(density(params[(burnin+1):iterations,p]), main = rownames(summary)[p], xlab = "Parameter Value")
}

