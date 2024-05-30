# example 8.2: analyze model-based multiple imputations for a multilevel regression with random intercepts
# requires fdir, rockchalk, lme4, and mitml packages
library(fdir)
library(rockchalk)
library(lme4)
library(mitml)
library(tidyverse)

# set working directory
fdir::set()

# inspect problemsolving2level.dat
edu2lv <- read_table("problemsolving2level.dat",col_names = FALSE,na="999") 
dim(edu2lv) # 982 14
names(edu2lv) <- c("school","student","condition","teachexp","eslpct","ethnic","male","frlunch",
                   "achievegrp","stanmath","efficacy1","efficacy2","probsolve1","probsolve2")
summary(edu2lv)

# read imputed data from working directory
imps <- read.table("imps.dat")
names(imps) <- c("Imputation","school","student","condition","teachexp","eslpct","ethnic","male","frlunch",
                 "achievegrp","stanmath","efficacy1","efficacy2","probsolve1","probsolve2",
                 "ranicept.l2", "frlunch.latent", "frlunch.l2mean","stanmath.l2mean", "probsolv1.l2mean") 
summary(imps)
# distribution of random intercept residuals
plot(density(imps$ranicept.l2))
rockchalk::skewness(imps$ranicept.l2)
rockchalk::kurtosis(imps$ranicept.l2)

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
model <- "probsolve2 ~ probsolve1  + stanmath + frlunch + teachexp + condition + (1 | school)"
analysis <- with(implist, lme4::lmer(model, REML = T))

# significance tests with barnard & rubin degrees of freedom
df <- 29 - 5 - 1
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = df)
estimates
confint(estimates, level = .95)