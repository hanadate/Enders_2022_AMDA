# example 8.7: analyze joint model multilevel multiple imputations from jomo
# requires fdir, rockchalk, lme4, and mitml packages

# set working directory
fdir::set()

# load data
load("jomo.imps.Rdata")
imps <- imps[imps$Imputation > 0,]

# create grouping variable that indexes clusters and imputations
imps$clus.by.imp <- imps$Imputation * 1000 + as.numeric(imps$clus)

# add group means to data
imps <- rockchalk::gmc(imps, c("pain"), by = c("clus.by.imp"), FUN = mean, suffix = c(".l2mean", ".cent"), fulldataframe = TRUE)

# within-cluster (group mean) center at the mean of the group means
imps <- rockchalk::gmc(imps, c("pain"), by = c("clus"), FUN = mean, suffix = c(".meanl2mean", ".cwc"), fulldataframe = TRUE)

# center at grand means
imps$sleep.cgm <- imps$sleep - mean(imps$sleep)
imps$pain.l2mean.cgm <- imps$pain.l2mean - mean(imps$pain)
imps$painaccept.cgm <- imps$painaccept - mean(imps$painaccept)

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
model <- "posaff ~ pain.cwc + sleep.cgm + pain.l2mean.cgm + painaccept.cgm + female + (1 + pain.cwc | clus)"
analysis <- with(implist, lme4::lmer(model, REML = T))

# significance tests with barnard & rubin degrees of freedom
df <- 132 - 5 - 1
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = df)
estimates
confint(estimates, level = .95)

# test within- and between-cluster slope difference
slope.diff <- c("pain.l2mean.cgm - pain.cwc")
mitml::testConstraints(analysis, constraints = slope.diff, df.com = df)
