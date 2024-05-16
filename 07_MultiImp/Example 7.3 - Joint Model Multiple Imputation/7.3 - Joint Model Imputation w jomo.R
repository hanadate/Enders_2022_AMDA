# example 7.3: joint model multiple imputation w jomo
# requires fdir and jomo packages
library(fdir)
library(jomo)
library(tidyverse)

# set working directory
fdir::set()

# read data from working directory
dat <- read.table("math.dat", na.strings = "999")
names(dat) <- c("id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost")

# inspect dat
summary(dat)

# define binary variable as factor
# note:  dat$frlunch is correct, not data$frlunch.
# note: to reproduce table7.1, levels=c(1,0) is correct, not c(0,1). just sign direction changes.
dat$frlunch <- factor(dat$frlunch, levels = c(1,0))

# select variables for imputation
vars2impute <- c("frlunch", "stanread","mathpre", "mathpost") 

# joint model imputation with jomo
set.seed(90291)
imps <- jomo::jomo(dat[vars2impute], nburn = 1000, nbetween = 1000, nimp = 100)

# save data
save(imps, file = "jomo.imps.Rdata")

# table7.1
imps %>% glimpse
