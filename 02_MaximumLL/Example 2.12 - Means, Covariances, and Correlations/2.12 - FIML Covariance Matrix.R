# example 2.12: fiml means, variances, and correlations
# requires fdir and lavaan packages
library(fdir)
library(lavaan)
# set working directory
fdir::set()

# read data from working directory
dat <- read.table("employeecomplete.dat")
names(dat) <- c("employee","team","turnover","male","empower","lmx","worksat","climate","cohesion")
# turnover: Intend to quit job in the next 6 months
# empower: Employee empowerment composite
# lmx: Leader-member exchange (relationship quality w/ sv) composite
# worksat(jobsat): Work satisfaction rating
# climate: Leadership climate composite (team-level)
# cohesion: Team cohesion composite

# specify model
model <- '
  worksat ~~ empower
  worksat ~~ lmx
  empower ~~ lmx
'

# estimate model in lavaan
fit <- lavaan::sem(model, dat, meanstructure = T, fixed.x = F)
summary(fit, standardize = T)
