# example 7.4: analyze fully conditional specification multiple imputations from blimp
# requires fdir and mitml packages
library(fdir)
library(mitml)
library(tidyverse)
# library(rblimp)
library(mice)

# set working directory
fdir::set()

# read math
dat.all <- readr::read_table(file="math.dat",
                             col_names=c("id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost"),
                             na="999")
summary(dat.all)
# for centering mathpre, stanread, (frlunch:binary)
dat.all.means <- apply(dat.all,2,mean,na.rm=TRUE)

# read imputed data from working directory
imps <- read.table("imps.dat", na.strings="999.0000")
names(imps) <- c("Imputation", "id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost")

# compute change score from imputed data
imps$change <- imps$mathpost - imps$mathpre
summary(imps)
glimpse(imps)

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
analysis <- with(implist, lm(change ~ 1))

# significance tests with barnard & rubin degrees of freedom
estimates <- mitml::testEstimates(analysis,
                                  extra.pars = T, 
                                  df.com = 249)
estimates
confint(estimates, level = .95)

# TODO
# table7.2 posterior summary of the FCS
# analyse
fit.mathpost <- with(implist, 
                     exp=lm(mathpost ~ 
                              I(mathpre-dat.all.means["mathpre"])+
                              I(stanread-dat.all.means["stanread"])+
                              I(frlunch-dat.all.means["frlunch"])
                     ))
summary(pool(fit.mathpost))
#------ intercept, mathpre, stanread, frlunch
# mine: 56.74, 0.42, 0.32, -1.80
# true: 56.64, 0.42, 0.31, -1.04

fit.stanread <- with(implist, 
                     exp=lm(stanread ~ 
                              I(mathpre-dat.all.means["mathpre"])+
                              frlunch+
                              I(mathpost-dat.all.means["mathpost"])
                     ))
summary(pool(fit.stanread))
#------ intercept, mathpre, frlunch, mathpost
# mine: 54.50, 0.04, -4.51, 0.48
# true: 56.64, 0.42, 0.31, -1.04

fit.frlunch <- with(implist, exp=glm(frlunch ~ 
                                    I(mathpre-dat.all.means["mathpre"])+
                                    I(mathpost-dat.all.means["mathpost"])+
                                    I(stanread-dat.all.means["stanread"]),
                                  family=binomial))
summary(pool(fit.frlunch))
#------ intercept, mathpre, mathpost, stanread
# mine: -0.47, 0.01, -0.04, -0.06
# true: -0.27, 0.01, -0.02, -0.04 

