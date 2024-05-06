# example 3.10: regression with auxiliary variables
# requires fdir, lavaan, semtools, and mdmb packages
library(fdir)
library(lavaan)
library(semTools)
library(mdmb)

# set working directory
fdir::set()

# read data from working directory
dat <- read.table("drugtrial.dat", na.strings = "999")
names(dat) <- c("id","male","drug","severity0","severity1","severity3","severity6","dropgrp","edrop","ldrop","dropout","sdrop3","sdrop6")
# drug: medication condition dummy code
# severity0: baseline
# severity1: at 1week
# dropgrp: dropout group 1:completer 2:3-week 3:6-week
# edrop: 3-week dropout
# ldrop: 6-week dropout

# center at the means
dat$severity0.cgm <- dat$severity0 - mean(dat$severity0, na.rm = T)
dat$male.cgm <- dat$male - mean(dat$male, na.rm = T)
summary(dat)

##################################################################
# model without auxiliary variables in lavaan
##################################################################

# specify lavaan model with no auxiliary variables
model <- 'severity6 ~ drug + severity0.cgm + male.cgm'

# estimate model with no auxiliary variables, table3.10
fit.noaux <- lavaan::sem(model, dat, fixed.x = F, missing = "fiml")
lavaan::summary(fit.noaux, rsquare = T, standardize = T)

##################################################################
# saturated correlates model in lavaan
##################################################################

# estimate model with semtools and lavaan, table3.10
fit.satcor <- semTools::sem.auxiliary(model, dat, aux = c("severity1","severity3"))
lavaan::summary(fit.satcor, rsquare = T, standardize = T)

##################################################################
# extra dependent variable model in lavaan
##################################################################

# specify lavaan model with auxiliary variables as extra outcomes
model.extradv <- 'severity6 ~ drug + severity0.cgm + male.cgm
                  severity1 ~ drug + severity0.cgm + male.cgm
                  severity3 ~ drug + severity0.cgm + male.cgm
                  severity6 ~~ severity1
                  severity6 ~~ severity3
                  severity1 ~~ severity3'

# estimate model in lavaan, table3.10
fit.extradv <- lavaan::sem(model.extradv, dat, fixed.x = F, missing = "fiml")
lavaan::summary(fit.extradv, rsquare = T, standardize = T)

##################################################################
# sequential specification in lavaan
##################################################################

# specify lavaan model with factored (sequential) specification for auxiliary variables
model.seqaux <- 'severity6 ~ drug + severity0.cgm + male.cgm
                 severity3 ~ severity6 + drug + severity0.cgm + male.cgm
                 severity1 ~ severity3 + severity6 + drug + severity0.cgm + male.cgm'

# estimate model in lavaan
fit.seqaux <- lavaan::sem(model.seqaux, dat, fixed.x = F, missing = "fiml")
lavaan::summary(fit.seqaux, rsquare = T, standardize = T)

##################################################################
# sequential specification in mdmb
##################################################################

# estimate sample statistics in lavaan for centering
lavaan::inspectSampleCov(fit.satcor, dat, fixed.x = F, missing = "fiml")

# summarize incomplete variables to determine ranges for pseudo-imputations
summary(dat[,c("severity6","severity0","severity3","severity1","male")])
# severity0 mean = 5.367
# male mean = .474
# for model as const
severity0.mean <- round(mean(dat$severity0, na.rm=TRUE),3)
male.mean <- round(mean(dat$male, na.rm=TRUE),3)

# set ranges (nodes) for pseudo-imputations
nodes.severity0 <- seq(1, 7, by = .25)
nodes.severity6 <- seq(1, 7, by = .25)
nodes.severity3 <- seq(1, 7, by = .25)
nodes.severity1 <- seq(1, 7, by = .25)

# model for severity0 predictor (severity1 | drug male)
model.severity0 <- list( "model"="linreg", 
                         "formula" = severity0 ~ drug + male,
                         nodes = nodes.severity0)

# model for severity6 outcome (severity6 | severity0 drug male)
# severity0 mean = 5.367
# male mean = .474
model.severity6 <- list("model" = "linreg", 
                        "formula"=severity6~drug+I(severity0-5.367)+I(male-.474),
                        nodes = nodes.severity6)
model.severity6.2 <- list("model" = "linreg",
                          "formula"=severity6~drug+I(severity0-mean(severity0))+I(male-mean(male)), 
                          nodes = nodes.severity6)
model.severity6.3 <- list("model" = "linreg",
                          "formula"=severity6~drug+I(severity0)+I(male), 
                          nodes = nodes.severity6)

model.severity6
model.severity6.2
model.severity6.3

# model for severity3 auxiliary variable (severity3 | severity6 severity0 drug male)
model.severity3 <- list("model"="linreg",
                        "formula"=severity3~severity6+drug+severity0+male,
                        nodes = nodes.severity3)

# model for severity1 auxiliary variable (severity1 | severity3 severity6 severity0 drug male)
model.severity1 <- list("model"="linreg",
                        "formula"=severity1~severity3+severity6+drug+severity0+male,
                        nodes = nodes.severity1)

# combine predictor models into a list (focal model functions as a predictor model)
predictor.models <- list(severity0 = model.severity0, 
                         severity6 = model.severity6, 
                         severity3 = model.severity3)

# estimate factored regression model w mdmb
fit <- mdmb::frm_em(dat = dat, dep = model.severity1, ind = predictor.models) 
summary(fit)
