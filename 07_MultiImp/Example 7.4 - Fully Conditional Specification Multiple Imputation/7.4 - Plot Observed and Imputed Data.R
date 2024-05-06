# example 7.4: plot fully conditional specification imputations
# requires the fdir package

# set working directory
fdir::set()

# read original data from working directory
original.dat <- read.table("math.dat", na.strings = "999")
names(original.dat) <- c("id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost")

# read imputed data from working directory
imps <- read.table("imps.dat")
names(imps) <- c("Imputation", "id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost")

# add missing data indicator to imputed data
imps$indicator <- rep(is.na(original.dat$mathpost), 100)

# observed and imputed data summaries
summary(imps[imps$indicator == F, "mathpost"])
summary(imps[imps$indicator == T, "mathpost"])

# plot observed vs. imputed data
xrange <- c(10,100)
hist(imps[imps$indicator == F, "mathpost"], breaks = seq(from = xrange[1], to = xrange[2], by = 2), col= rgb(0.4,0.4,0.4,0.4), xlim = xrange, xlab = "Math Posttest Scores", ylab = NA, axes = T, main = "Observed vs. Imputed Data")
hist(imps[imps$indicator == T, "mathpost"], breaks = seq(from = xrange[1], to = xrange[2], by = 2), col = "red", add = T, xlab = NA, ylab = NA, axes = F)