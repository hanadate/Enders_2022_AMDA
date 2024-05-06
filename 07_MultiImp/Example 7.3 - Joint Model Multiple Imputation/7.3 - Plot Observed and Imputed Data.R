# example 7.3: plot joint model imputations from jomo
# requires the fdir package

# set working directory
fdir::set()

# load jomo imputations
load("jomo.imps.RData")

# add missing data indicator to imputed data
original.dat <- imps[imps$Imputation == 0,]
imps <- imps[imps$Imputation > 0,]
imps$indicator <- rep(is.na(original.dat$mathpost), 100)

# observed and imputed data summaries
summary(imps[imps$indicator == F, "mathpost"])
summary(imps[imps$indicator == T, "mathpost"])

# plot observed vs. imputed data
xrange <- c(10,100)
hist(imps[imps$indicator == F, "mathpost"], breaks = seq(from = xrange[1], to = xrange[2], by = 2), col= rgb(0.4,0.4,0.4,0.4), xlim = xrange, xlab = "Math Posttest Scores", ylab = NA, axes = T, main = "Observed vs. Imputed Data")
hist(imps[imps$indicator == T, "mathpost"], breaks = seq(from = xrange[1], to = xrange[2], by = 2), col = "red", add = T, xlab = NA, ylab = NA, axes = F)