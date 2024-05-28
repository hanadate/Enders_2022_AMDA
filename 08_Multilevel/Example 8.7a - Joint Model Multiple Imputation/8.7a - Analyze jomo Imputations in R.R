# example 8.7: analyze joint model multilevel multiple imputations from jomo
# requires fdir and mitml packages

# set working directory
fdir::set()

# load data
load("jomo.imps.Rdata")

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
model <- "probsolve2 ~ probsolve1  + stanmath + frlunch + teachexp + condition + (1 | clus)"
analysis <- with(implist, lme4::lmer(model, REML = T))

# significance tests with barnard & rubin degrees of freedom
df <- 29 - 5 - 1
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = df)
estimates
confint(estimates, level = .95)