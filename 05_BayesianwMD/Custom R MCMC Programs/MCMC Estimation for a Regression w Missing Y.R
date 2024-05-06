# load packages
library(fdir)
library(mvtnorm)

# set working directory to this script's location
fdir::set()

# read data
datmiss <- read.table("employee.dat", na.strings = "999")
names(datmiss) = c("employee", "team", "turnover", "male", "empower", "lmx", "jobsat", "climate", "cohesion")
N <- nrow(datmiss)

# complete data (to get a complete X)
datcomp <- read.table("employeecomplete.dat")
names(datcomp) = c("employee", "team", "turnover", "male", "empower", "lmx", "jobsat", "climate", "cohesion")

# assemble data
X <- cbind(1, datcomp$lmx)
Y <- datmiss$empower

# missing data indicator
Yind <- rep(0, N)
Yind[is.na(Y)] <- 1

# initial imputations fill in the mean
Y[Yind == 1] <- mean(Y, na.rm = T)

# initialize algorithmic features
set.seed(90291)
iterations <- 11000
burnin <- 1000

# store parameter estimates
betas <- matrix(0, ncol = ncol(X), nrow = iterations)
sigma2e <- matrix(1, iterations)
R2 <- matrix(1, iterations)
stanbeta <- matrix(1, iterations)

# gibbs sampler
for(t in 2:iterations){
  
  # print iteration history
  if(t %% 500 == 0){print(paste0("Iteration = " , t))}

  # regression coefficients, conditional on residual variance and missing values
  betahat <- solve(crossprod(X,X)) %*% crossprod(X,Y)
  covbeta <- solve(crossprod(X,X)) * sigma2e[t-1]
  betas[t,] <- rmvnorm(1, mean = betahat, sigma = covbeta)
  
  # sample variance, conditional on coefficients and missing values
  df <- N / 2
  sumofsq <- sum((Y - X %*% betas[t,])^2) / 2
  sigma2e[t] <- 1 / rgamma(1, df, rate = sumofsq)
  
  # compute R^2 and standardized B1
  R2[t] <- betas[t,2]^2 * var(X[,2]) / var(Y)
  stanbeta[t] <- betas[t,2] * sd(X[,2]) / sd(Y)
  
  # impute Y conditional on analysis model parameters
  Yhat <- X %*% betas[t,]
  varY <- sigma2e[t]
  imps <- rnorm(N, Yhat, sqrt(varY))
  Y[Yind == 1] <- imps[Yind == 1]
  
}

# summarize posterior distributions
summary <- matrix(0, nrow = 5, ncol = 5)
rownames(summary) <- c("B0", "B1", "res. var.", "B1(std.)", "R-sq")
colnames(summary) <- c("Mean", "StdDev", "2.5%", "50%", "97.5%")

# summarize posterior distributiona
params <- cbind(betas, sigma2e, stanbeta, R2)
for(p in 1:nrow(summary)){
  summary[p,1] <- mean(params[(burnin+1):iterations,p])
  summary[p,2] <- sd(params[(burnin+1):iterations,p])
  summary[p,3:5] <- quantile(params[(burnin+1):iterations,p], c(.025, .50, .975))
  plot(density(params[(burnin+1):iterations,p]), main = rownames(summary)[p], xlab = "Parameter Value")
}

print(paste0("Posterior Distribution Summary from ", iterations - burnin, " Iterations:"))
print(round(summary, 3))


