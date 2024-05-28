# example 8.7: joint model mulilevel multiple imputation w jomo
# requires fdir and jomo packages

# set working directory
fdir::set()

# read imputed data from working directory
dat <- read.table("problemsolving2level.dat", na.strings = "999")
names(dat) <- c("school","student","condition","teachexp","eslpct","ethnic","male","frlunch",
                 "achievegrp","stanmath","efficacy1","efficacy2","probsolve1","probsolve2") 

# define binary variable as a factor
dat$frlunch <- factor(dat$frlunch, levels = c(0,1))

# select variables for imputation
dat$icept <- 1
vars2impute.l1 <- c("frlunch","stanmath","probsolve2")
vars2impute.l2 <- c("teachexp")
complete.l1 <- c("icept","probsolve1")
complete.l2 <- c("icept","condition")

# joint model imputation with jomo
set.seed(90291)
imps <- jomo::jomo(dat[vars2impute.l1], Y2 = dat[vars2impute.l2], X = dat[complete.l1], X2 = dat[complete.l2], 
                   clus = dat$school, nburn = 2000, nbetween = 2000, nimp = 100)

# save data
save(imps, file = "jomo.imps.Rdata")


