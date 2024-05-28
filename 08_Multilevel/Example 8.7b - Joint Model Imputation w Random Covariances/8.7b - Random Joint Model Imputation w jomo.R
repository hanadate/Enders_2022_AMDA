# example 8.7: random joint model mulilevel multiple imputation w jomo
# requires fdir and jomo packages

# set working directory
fdir::set()

# read imputed data from working directory
dat <- read.table("diary.dat", na.strings = "999")
names(dat) <- c('person','day','pain','sleep','posaff','negaff','lifegoal','female','educ','diagnoses','activity',
                 'painaccept','catastrophize','stress')

# define binary variable as a factor
dat$female <- factor(dat$female, levels = c(0,1))

# define level-1 and level-2 variables
dat$icept <- 1
vars2impute.l1 <- c("posaff","pain","sleep")
vars2impute.l2 <- c("painaccept")
complete.l2 <- c("icept","female")

# joint model imputation with random covariance matrices in jomo
imps <- jomo::jomo(dat[vars2impute.l1], Y2 = dat[vars2impute.l2], X2 = dat[complete.l2], clus = dat$person, 
             nburn = 2000, nbetween = 2000, nimp = 100, meth = "random")

# save data
save(imps, file = "jomo.imps.Rdata")


