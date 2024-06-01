# example 7.4: analyze fully conditional specification multiple imputations from mice
# This script is not included in the orignal source.
library(mice)
library(fdir)
library(tidyverse)
library(doParallel)
library(coda)

fdir::set()

# read math
dat.all <- readr::read_table(file="math.dat",
                             col_names=c("id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost"),
                             na="999")
summary(dat.all)

# centering mathpre, stanread, (frlunch:binary)
dat.all.means <- apply(dat.all,2,mean,na.rm=TRUE)

# select cols, and define frlunch as categorical.
# dat <- dat.all.centerd %>% 
dat <- dat.all %>% 
  dplyr::select(mathpost, mathpre, stanread, frlunch) %>% 
  dplyr::mutate(frlunch=as.factor(frlunch))
glimpse(dat)
summary(dat)

# # detect cores
# detectCores() #16
# # imputes
# t <- now()
# imps <- futuremice(data=dat, parallelseed=1111,
#                    m=100, maxit=1000,
#                    n.core=10,
#                    method=c("norm","norm","norm","logreg"))
# now()-t
# save(imps, file = "mice.imps.Rdata")
load("mice.imps.Rdata")
# 2mins
summary(imps)
imps$formulas

# # conv
# plot(imps)
# conv <- convergence(imps)
# plot(conv)
# # diagnostics 
# bwplot(imps) # box-and-whisker
# stripplot(imps) # for num. red is imp.
# densityplot(imps) # density

# # complete data
# # 100 is numeric. 100L is long-integer.
# dat.comp <- mice::complete(imps, action=100L, include=FALSE)
# summary(dat.comp)

# pool
# TODO: table2
fit.mathpost <- with(imps, 
                     exp=lm(mathpost ~ 
                              I(mathpre-dat.all.means["mathpre"])+
                              I(stanread-dat.all.means["stanread"])+
                              frlunch
                     ))
summary(pool(fit.mathpost))
# mine: 57.46, 0.42, 0.32, -1.70
# true: 56.64, 0.42, 0.31, -1.04

fit.stanread <- with(imps, 
                     exp=lm(stanread ~ 
                              I(mathpre-dat.all.means["mathpre"])+
                              frlunch+
                              I(mathpost-dat.all.means["mathpost"])
                     ))
summary(pool(fit.stanread))
# mine: 54.55, 0.05, -4.73, 0.46
# true: 56.64, 0.42, 0.31, -1.04

fit.frlunch <- with(imps, exp=glm(frlunch ~ 
                                    I(mathpre-dat.all.means["mathpre"])+
                                    I(mathpost-dat.all.means["mathpost"])+
                                    I(stanread-dat.all.means["stanread"]),
                                  family=binomial))
summary(pool(fit.frlunch))
# mine: -0.46, 0.01, -0.03, -0.06
# true: -0.27, 0.01, -0.02, -0.04 
