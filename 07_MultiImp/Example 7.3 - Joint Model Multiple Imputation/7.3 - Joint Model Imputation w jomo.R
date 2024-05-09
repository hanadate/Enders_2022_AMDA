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
dat$frlunch <- factor(data$frlunch, levels = c(0,1))

# select variables for imputation
vars2impute <- c("frlunch", "stanread","mathpre", "mathpost") 

# joint model imputation with jomo
set.seed(90291)
imps <- jomo::jomo(dat[vars2impute], nburn = 1000, nbetween = 1000, nimp = 100)

# save data
save(imps, file = "jomo.imps.Rdata")