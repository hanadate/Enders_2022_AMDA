# example 5.8: bayes regression with factored (sequential) specification for incomplete predictors and auxiliary variables
# requires fdir, lavaan, and mdmb packages

# set working directory
fdir::set()

# read data from working directory
dat <- read.table("drugtrial.dat", na.strings = "999")
names(dat) <- c("id","male","drug","severity0","severity1","severity3","severity6","dropgrp","edrop","ldrop","dropout","sdrop3","sdrop6")

# specify lavaan model for sample statistics
model <- 'severity0 ~~ severity1
          severity0 ~~ severity3
          severity0 ~~ severity6
          severity0 ~~ drug
          severity0 ~~ male
          severity1 ~~ severity3
          severity1 ~~ severity6
          severity1 ~~ drug
          severity1 ~~ male
          severity3 ~~ severity6
          severity3 ~~ drug
          severity3 ~~ male
          drug ~~ male'

# estimate sample statistics in lavaan
lavaan::inspectSampleCov(model, dat, meanstructure = TRUE, missing = "fiml")

# summarize variables to determine ranges for pseudo-imputations
summary(dat[,c("severity6","severity0","severity3","severity1","drug","male")])

# set ranges (nodes) for pseudo-imputations
nodes.severity0 <- seq(1, 7, by = .25)
nodes.severity6 <- seq(1, 7, by = .25)
nodes.severity3 <- seq(1, 7, by = .25)
nodes.severity1 <- seq(1, 7, by = .25)

# model for severity0 predictor (severity1 | drug male)
model.severity0 <- list( "model"="linreg", "formula" = severity0 ~ drug + male, nodes = nodes.severity0)

# model for severity6 outcome (severity6 | severity0 drug male)
# severity0 mean = 5.367
# male mean = .474
model.severity6 <- list("model" = "linreg", "formula" = severity6 ~ drug + I(severity0 - 5.367) + I(male - .474), nodes = nodes.severity6)

# model for severity3 auxiliary variable (severity3 | severity6 severity0 drug male)
model.severity3 <- list( "model"="linreg", "formula" = severity3 ~ severity6 + drug + severity0 + male, nodes = nodes.severity3)

# model for severity1 auxiliary variable (severity1 | severity3 severity6 severity0 drug male)
model.severity1 <- list( "model"="linreg", "formula" = severity1 ~ severity3 + severity6 + drug + severity0 + male, nodes = nodes.severity1)

# combine predictor models into a list (focal model functions as a predictor model)
predictor.models <- list(severity0 = model.severity0, severity6 = model.severity6, severity3 = model.severity3)

# estimate factored regression model w mdmb
fit <- mdmb::frm_fb(dat = dat, dep = model.severity1, ind = predictor.models) 
summary(fit)
