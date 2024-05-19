# example 7.4: plot fully conditional specification imputations
# requires the fdir package
library(fdir)
library(tidyverse)

# set working directory
fdir::set()

# read original data from working directory
original.dat <- read.table("math.dat", na.strings = "999")
names(original.dat) <- c("id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost")
summary(original.dat)
# read imputed data from working directory
imps <- read.table("imps.dat",  na.strings="999.0000")
names(imps) <- c("Imputation", "id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost")
summary(imps)
# add missing data indicator to imputed data
imps$indicator <- rep(is.na(original.dat$mathpost), 100)

# observed and imputed data summaries
summary(imps[imps$indicator == F, "mathpost"])
summary(imps[imps$indicator == T, "mathpost"])

summary(imps[imps$indicator == F, "stanread"])
summary(imps[imps$indicator == T, "stanread"])

# plot observed vs. imputed data
xrange <- c(10,100)
hist(imps[imps$indicator == F, "mathpost"], 
     breaks = seq(from = xrange[1], to = xrange[2], by = 2), 
     col= rgb(0.4,0.4,0.4,0.4), xlim = xrange, 
     xlab = "Math Posttest Scores", ylab = NA, 
     axes = T, main = "Fig7.3 Observed vs. Imputed Data")
hist(imps[imps$indicator == T, "mathpost"], 
     breaks = seq(from = xrange[1], to = xrange[2], by = 2), 
     col = "red", add = T, xlab = NA, ylab = NA, axes = F)

hist(imps[imps$indicator == F, "stanread"], 
     breaks = seq(from = xrange[1], to = xrange[2], by = 2), 
     col= rgb(0.4,0.4,0.4,0.4), xlim = xrange, 
     xlab = "Standardized Reading", ylab = NA, 
     axes = T, main = "Fig7.4 Observed vs. Imputed Data")
hist(imps[imps$indicator == T, "stanread"], 
     breaks = seq(from = xrange[1], to = xrange[2], by = 2), 
     col = "red", add = T, xlab = NA, ylab = NA, axes = F)
