# example 3.6: fiml regression with factored (sequential) specification for incomplete predictors
# requires fdir and mdmb packages
library(fdir)
library(mdmb)
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

# summarize incomplete predictors to determine ranges for pseudo-imputations
summary(dat[,c("intensity","parsmoke","income","age")])

# set ranges (nodes) for pseudo-imputations
nodes.income <- seq(1, 20, by = 1)
nodes.parsmoke <- c(0,1)
nodes.intensity <- seq(0, 31, by = 1)

# model for parsmoke predictor (parsmoke | age)
model.parsmoke <- list( "model"="logistic", "formula" = parsmoke ~ age, nodes = nodes.parsmoke)

# model for income predictor (income | parsmoke, age)
model.income <- list( "model"="linreg", "formula" = income ~ parsmoke + age, nodes = nodes.income)

# model for intensity outcome (intensity | parsmoke income age)
# age mean = 21.55
# income mean = 8.43
model.intensity <- list("model" = "linreg", "formula" = intensity ~  parsmoke + I(income - 8.43) + I(age - 21.55), nodes = nodes.intensity)

# combine predictor models into a list
predictor.models <- list(parsmoke = model.parsmoke, income = model.income)

# estimate factored regression model w mdmb
fit <- mdmb::frm_em(dat = dat, dep = model.intensity, ind = predictor.models) 
summary(fit)
