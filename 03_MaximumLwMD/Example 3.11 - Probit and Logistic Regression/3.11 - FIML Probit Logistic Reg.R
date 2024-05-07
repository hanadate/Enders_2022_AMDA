# example 3.11: fiml logistic and probit regression with factored (sequential) specification for incomplete predictors
# requires fdir and mdmb packages
library(fdir)
library(mdmb)
library(tidyverse)
# set working directory
fdir::set()

# read data from working directory
dat <- read.table("employee.dat", na.strings = "999")
names(dat) <- c("employee","team","turnover","male","empower","lmx","worksat","climate","cohesion")

# summarize incomplete variables to determine ranges for pseudo-imputations
summary(dat[,c("turnover","lmx","empower","male")])

# set ranges (nodes) for pseudo-imputations
nodes.empower <- seq(0, 50, by = 1)
nodes.lmx <- seq(-5, 32, by = 1)
nodes.turnover <- c(0,1)

# model for empower predictor (empower | male)
model.empower <- list( "model"="linreg", "formula" = empower ~ male, nodes = nodes.empower)

# model for lmx predictor (lmx | empower male)
model.lmx <- list( "model"="linreg", "formula" = lmx ~ empower + male, nodes = nodes.lmx)

# logistic model for turnover outcome (turnover | lmx empower male)
logistic.model.turnover <- list( "model"="logistic", "formula" = turnover ~ lmx + empower + male, nodes = nodes.turnover)

# probit model for turnover outcome (turnover | lmx empower male)
probit.model.turnover <- list( "model"="oprobit", "formula" = turnover ~ lmx + empower + male, nodes = nodes.turnover)

# combine predictor models into a list
predictor.models <- list(empower = model.empower, lmx = model.lmx)

# estimate factored logistic regression model w mdmb
fit.logistic <- mdmb::frm_em(dat = dat, dep = logistic.model.turnover, ind = predictor.models) 
summary(fit.logistic)

# estimate factored probit regression model w mdmb
fit.probit <- mdmb::frm_em(dat = dat, dep = probit.model.turnover, ind = predictor.models) 
summary(fit.probit)
