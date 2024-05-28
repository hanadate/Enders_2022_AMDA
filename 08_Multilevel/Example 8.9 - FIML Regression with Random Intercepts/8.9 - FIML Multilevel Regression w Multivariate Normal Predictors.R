# example 8.9: maximum likelihood estimation for a multilevel regression with random intercepts;
# requires fdir and lavaan packages

# set working directory
fdir::set()

# read imputed data from working directory
dat <- read.table("problemsolving2level.dat", na.strings = "999")
names(dat) <- c("school","student","condition","teachexp","eslpct","ethnic","male","frlunch",
                 "achievegrp","stanmath","efficacy1","efficacy2","probsolve1","probsolve2") 

# specify model
model <- '
    level: 1
    probsolve2 ~ probsolve1  + stanmath + frlunch
    level: 2
    probsolve2 ~ teachexp + condition
'

# estimate model in lavaan
fit <- lavaan::sem(model, dat,cluster = "school", meanstructure = TRUE, fixed.x = FALSE, missing = "fiml")
lavaan::summary(fit, standardize = TRUE)

# missing data patterns and proportion observed data (coverage)
inspect(fit, "patterns")
inspect(fit, "coverage")