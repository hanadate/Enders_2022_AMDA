# load packages
library(fdir)
library(mvtnorm)
library(MCMCpack)
library(truncnorm)

# set working directory to this script's location
fdir::set()

# read data
dat <- read.table("employee.dat", na.strings = "999")
names(dat) = c("employee", "team", "turnover", "male", "empower", "lmx", "jobsat", "climate", "cohesion")
N <- nrow(dat)

# select variables
X <- dat$lmx
Y <- dat$turnover

# missing data indicators
Yind <- rep(0, N)
Yind[is.na(Y)] <- 1
Xind <- rep(0, N)
Xind[is.na(X)] <- 1

# initial X imputations at the mean
X[Xind == 1] <- mean(X, na.rm = T)

# initial discrete imputes w random draw from binomial distribution w p = .50
Yinitial <- rbinom(length(Y), 1, .50)
Y[Yind == 1] <- Yinitial[Yind == 1]

# initial latent scores
pY1 <- mean(Y)
Ystarmean <- -1 * qnorm(pY1, lower.tail = F)
Ystar <- rep(0, length(Y))
Ystar[Y == 0] <- rtruncnorm(length(Y[Y == 0]), a = -Inf, b = 0, mean = Ystarmean, sd = 1)
Ystar[Y == 1] <- rtruncnorm(length(Y[Y == 1]), a = 0, b = Inf, mean = Ystarmean, sd = 1)

# inspect Ystar values
plot(density(Ystar))
abline(v = 0)

# assemble X matrix
X <- cbind(1, X)

# initialize algorithmic features
set.seed(90291)
iterations <- 11000
burnin <- 1000

# store parameter estimates
betas <- matrix(0, ncol = ncol(X), nrow = iterations)
sigma2e <- 1
R2 <- matrix(1, iterations)
stanbeta <- matrix(1, iterations)
muX <- matrix(0, iterations)
sigma2X <- matrix(1, iterations)

# gibbs sampler
t <- 2
for(t in 2:iterations){
  
  # print iteration history
  if(t %% 500 == 0){print(paste0("Iteration = " , t))}
  
  # sample latent scores for complete cases
  Ystarhat <- X %*% betas[t-1,]
  Ystar[Y == 0 & Yind == 0] <- rtruncnorm(length(Y[Y == 0 & Yind == 0]), a = -Inf, b = 0, mean = Ystarhat[Y == 0 & Yind == 0], sd = 1)
  Ystar[Y == 1 & Yind == 0] <- rtruncnorm(length(Y[Y == 1 & Yind == 0]), a = 0, b = Inf, mean = Ystarhat[Y == 1 & Yind == 0], sd = 1)
  
  # sample regression coefficients, conditional on current data
  betahat <- solve(crossprod(X,X)) %*% crossprod(X,Ystar)
  covbeta <- solve(crossprod(X,X)) * sigma2e
  betas[t,] <- rmvnorm(1, mean = betahat, sigma = covbeta)
  
  # compute R^2 and standardized B1
  R2[t] <- betas[t,2]^2 * var(X[,2]) / (betas[t,2]^2 * var(X[,2]) + 1)
  stanbeta[t] <- betas[t,2] * sd(X[,2]) / sqrt(betas[t,2]^2 * var(X[,2]) + 1)
  
  # sample X mean, conditional on variance and current data
  Xbar <- mean(X[,2])
  VarXbar <- sigma2X[t-1] / N
  muX[t] <- rnorm(1, mean = Xbar, sd = sqrt(VarXbar))
  
  # sample X variance, conditional on mean and current data
  a <- N / 2
  b <- sum((X[,2] - muX[t])^2) / 2
  sigma2X[t] <- 1 / rgamma(1, a, rate = b)
  
  # impute missing latent Ys, conditional on regression model
  Ystarhat <- X %*% betas[t,]
  varY <- sigma2e
  imps <- rnorm(length(Y), Ystarhat, sqrt(varY))
  Ystar[Yind == 1] <- imps[Yind == 1]
  
  # categorize imputes (if we want to save them)
  Y[Ystar > 0] <- 1
  Y[Ystar <= 0] <- 0

  # impute X conditional on two sets of model parameters
  Xhat <- (sigma2X[t] * betas[t,2] * (Ystar - betas[t,1]) + sigma2e * muX[t]) / (betas[t,2]^2 * sigma2X[t] + sigma2e)
  varX <- (sigma2e * sigma2X[t]) / (betas[t,2]^2 * sigma2X[t] + sigma2e)
  imps <- rnorm(N, mean = Xhat, sd = sqrt(varX))
  X[Xind == 1, 2] <- imps[Xind == 1]
  
}

# summarize posterior distributions
summary <- matrix(0, nrow = 7, ncol = 5)
rownames(summary) <- c("B0", "B1", "res. var.", "B1(std.)","R-sq", "X mean", "X var")
colnames(summary) <- c("Mean", "StdDev", "2.5%", "50%", "97.5%")

# summarize posterior distributiona
params <- cbind(betas, sigma2e, stanbeta, R2, muX, sigma2X)
for(p in 1:nrow(summary)){
  summary[p,1] <- mean(params[(burnin+1):iterations,p])
  summary[p,2] <- sd(params[(burnin+1):iterations,p])
  summary[p,3:5] <- quantile(params[(burnin+1):iterations,p], c(.025, .50, .975))
  plot(density(params[(burnin+1):iterations,p]), main = rownames(summary)[p], xlab = "Parameter Value")
}

print(paste0("Posterior Distribution Summary from ", iterations - burnin, " Iterations:"))
print(round(summary, 3))


