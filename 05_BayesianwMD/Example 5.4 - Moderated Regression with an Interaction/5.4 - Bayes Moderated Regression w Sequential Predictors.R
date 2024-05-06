# example 3.8: moderated regression with factored (sequential) specification for incomplete predictors
# requires fdir, lavaan, and mdmb packages

# set working directory
fdir::set()

# read data from working directory
dat <- read.table("pain.dat", na.strings = "999")
names(dat) <- c("id", "txgrp", "male", "age", "edugroup", "workhrs", "exercise", "paingrps", 
                "pain", "anxiety", "stress", "control", "depress", "interfere", "disability",
                paste0("dep", seq(1:7)), paste0("int", seq(1:6)), paste0("dis", seq(1:6)))

# specify lavaan model for sample statistics
model <- '
  disability ~~ depress
  disability ~~ male
  disability ~~ pain
  depress ~~ male
  depress ~~ pain
  male ~~ pain
'

# estimate sample statistics in lavaan
lavaan::inspectSampleCov(model, dat, meanstructure = TRUE, missing = "fiml")

# summarize incomplete predictors to determine ranges for pseudo-imputations
summary(dat[,c("disability","depress","male","pain")])

# set ranges (nodes) for pseudo-imputations
nodes.pain <- c(0,1)
nodes.depress <- seq(0, 35, by = 1)
nodes.disability <- seq(0, 47, by = 1)

# model for pain predictor (pain | male)
model.pain <- list( "model"="logistic", "formula" = pain ~ male, nodes = nodes.pain)

# model for depression predictor (depress | pain male)
model.depress <- list( "model"="linreg", "formula" = depress ~ pain + male, nodes = nodes.depress)

# model for disability outcome (disability | depress male depress*male pain)
# depression mean = 14.718
model.disability <- list("model" = "linreg", "formula" = disability ~ I(depress - 14.718) + male + I((depress - 14.718)*male) + pain, nodes = nodes.disability)

# combine predictor models into a list
predictor.models <- list(pain = model.pain, depress = model.depress)

# estimate factored regression model w mdmb
fit <- mdmb::frm_fb(dat = dat, dep = model.disability, ind = predictor.models) 
summary(fit)