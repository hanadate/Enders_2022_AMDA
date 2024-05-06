# example 7.3: analyze joint model imputations from jomo
# requires fdir and mitml packages

# set working directory
fdir::set()

# load data
load("jomo.imps.Rdata")
imps$change <- imps$mathpost - imps$mathpre

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
analysis <- with(implist, lm(change ~ 1))

# significance tests with barnard & rubin degrees of freedom
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = 249)
estimates
confint(estimates, level = .95)
