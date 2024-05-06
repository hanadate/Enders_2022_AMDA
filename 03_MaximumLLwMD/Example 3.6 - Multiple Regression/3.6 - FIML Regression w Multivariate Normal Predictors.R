# example 3.6: fiml regression with multivariate normal incomplete predictors
# requires fdir and lavaan packages
library(fdir)
library(lavaan)
library(lavaanPlot)
library(tidyverse)

# set working directory
fdir::set()

# read data from working directory
dat <- read.table("smoking.dat", na.strings = "999")
names(dat) <- c("id","intensity","hvysmoker","age","parsmoke","female","race","income","educ")
# intensity : somoking intensity, cigarettes per day
# hvysmoker : heavy smoking
# parsmoke : parental smoking 
# race : 1white, 2black, 3hispanic, 4other

# specify model with labels
model.labels <- 'intensity ~ b1*parsmoke + b2*age + b3*income'
wald.constraints <- 'b1 == 0; b2 == 0; b3 == 0'

# estimate model in lavaan
# model -> model.lavels
fit <- lavaan::sem(model.labels, dat, fixed.x = F, missing = "fiml")
lavaan::summary(fit, rsquare = T, standardize = T)
lavaanPlot(fit, coefs=TRUE, stand=TRUE, covs=TRUE)
# option ?buildPatshs()

# wald test that all slopes equal 0
lavaan::lavTestWald(fit, constraints = wald.constraints)

# missing data patterns and proportion observed data (coverage)
lavaan::inspect(fit, "patterns")
lavaan::inspect(fit, "coverage")
