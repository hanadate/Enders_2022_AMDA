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

# select cols, and define frlunch as categorical.
dat <- dat.all %>% 
  dplyr::select(mathpost, mathpre, stanread, frlunch) %>% 
  # # center
  # dplyr::mutate(mathpre=mathpre-mean(mathpre,na.rm=TRUE),
  #               stanread=stanread-mean(stanread,na.rm=TRUE),
  #               frlunch=frlunch-mean(frlunch,na.rm=TRUE)) %>% 
  dplyr::mutate(frlunch=as.factor(frlunch))
glimpse(dat)
summary(dat)

# detect cores
detectCores() #16
# imputes
t <- now()
imps <- futuremice(data=dat, parallelseed=1111, 
                 m=100, maxit=1000, 
                 n.core=10)
now()-t
save(imps, file = "mice.imps.Rdata") 
load("mice.imps.Rdata")
# 2mins
summary(imps)
imps$formulas

# conv
plot(imps)
conv <- convergence(imps)
plot(conv)
# diagnostics 
bwplot(imps) # box-and-whisker
# stripplot(imps) # for num. red is imp.
densityplot(imps) # density

# complete data
# 100 is numeric. 100L is long-integer.
data.imped <- mice::complete(imps, action=100L, include=FALSE) 
summary(data.imped)

data.imped.long <- mice::complete(imps, action="long", include=FALSE)
data.imped.long %>% glimpse
summary(data.imped.long)

# analyse and pool
fit.mathpost <- with(imps, exp=lm(mathpost ~ mathpre+stanread+frlunch))
summary(fit.mathpost)
pool(fit.mathpost)
summary(pool(fit.mathpost))

fit.stanread <- with(imps, exp=lm(stanread ~ mathpre+frlunch+mathpost))
pool(fit.stanread)
summary(pool(fit.stanread))

fit.frlunch <- with(imps, exp=lm(frlunch ~ mathpre+mathpost+stanread))
summary(pool(fit.frlunch))

