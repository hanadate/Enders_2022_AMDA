# example 8.8: analyze model-based multiple imputations for a multilevel regression with random intercepts
# requires fdir, rockchalk, lme4, and mitml packages

# set working directory
fdir::set()

# read imputed data from working directory
imps <- read.table("imps.dat")
names(imps) <- c("Imputation","school","student","condition","teachexp","eslpct","ethnic","male","frlunch",
                 "achievegrp","stanmath","efficacy1","efficacy2","probsolve1","probsolve2")

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
model <- "probsolve2 ~ probsolve1  + stanmath + frlunch + teachexp + condition + (1 | school)"
analysis <- with(implist, lme4::lmer(model, REML = T))

# significance tests with barnard & rubin degrees of freedom
df <- 29 - 5 - 1
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = df)
estimates
confint(estimates, level = .95)