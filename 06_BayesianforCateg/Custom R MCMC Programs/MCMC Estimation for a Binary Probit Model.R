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
X <- dat$lmx
Y <- dat$turnover

# assemble X matrix
X <- cbind(1, X)

# initialize algorithmic features
set.seed(90291)
iterations <- 11000
burnin <- 1000

# initialize Ystar
Ystar <- rep(0, length(Y))

# store parameter estimates
betas <- matrix(0, ncol = ncol(X), nrow = iterations)
sigma2e <- 1
R2 <- matrix(1, iterations)
stanbeta <- matrix(1, iterations)

# gibbs sampler
t <- 2
for(t in 2:iterations){
  
  # print iteration history
  if(t %% 500 == 0){print(paste0("Iteration = " , t))}
  
  # sample latent scores
  Ystarhat <- X %*% betas[t-1,]
  Ystar[Y == 0] <- rtruncnorm(length(Y[Y == 0]), a = -Inf, b = 0, mean = Ystarhat[Y == 0], sd = 1)
  Ystar[Y == 1] <- rtruncnorm(length(Y[Y == 1]), a = 0, b = Inf, mean = Ystarhat[Y == 1], sd = 1)
  
  # regression coefficients, conditional on residual variance and latent scores
  betahat <- solve(crossprod(X,X)) %*% crossprod(X,Ystar)
  covbeta <- solve(crossprod(X,X)) * sigma2e
  betas[t,] <- rmvnorm(1, mean = betahat, sigma = covbeta)
  
  # compute R^2 and standardized B1
  R2[t] <- betas[t,2]^2 * var(X[,2]) / (betas[t,2]^2 * var(X[,2]) + 1)
  stanbeta[t] <- betas[t,2] * sd(X[,2]) / sqrt(betas[t,2]^2 * var(X[,2]) + 1)

}

# summarize posterior distributions
summary <- matrix(0, nrow = 5, ncol = 5)
rownames(summary) <- c("B0", "B1", "res. var.", "B1(std.)","R-sq")
colnames(summary) <- c("Mean", "StdDev", "2.5%", "50%", "97.5%")

# summarize posterior distributions
params <- cbind(betas, sigma2e, stanbeta, R2)
for(p in 1:nrow(summary)){
  summary[p,1] <- mean(params[(burnin+1):iterations,p])
  summary[p,2] <- sd(params[(burnin+1):iterations,p])
  summary[p,3:5] <- quantile(params[(burnin+1):iterations,p], c(.025, .50, .975))
  plot(density(params[(burnin+1):iterations,p]), main = rownames(summary)[p], xlab = "Parameter Value")
}

print(paste0("Posterior Distribution Summary from ", iterations - burnin, " Iterations:"))
print(round(summary, 3))


