# example 8.3: analyze model-based multiple imputations for a multilevel regression with random slopes
# requires fdir, rockchalk, lme4, and mitml packages

# set working directory
fdir::set()

# read imputed data from working directory
imps <- read.table("imps.dat")
names(imps) <- c('Imputation','person','day','pain','sleep','posaff','negaff','lifegoal','female','educ','diagnoses','activity',
                 'painaccept','catastrophize','stress','ranicept.l2','ranslp.l2','pain.l2mean','sleep.l2mean')

# center at grand means and imputed latent group means
imps$sleep.cgm <- imps$sleep - mean(imps$sleep)
imps$pain.cwc <- imps$pain - imps$pain.l2mean
imps$pain.l2mean.cgm <- imps$pain.l2mean - mean(imps$pain.l2mean)
imps$painaccept.cgm <- imps$painaccept - mean(imps$painaccept)

# distribution of random intercept and slope residuals
plot(density(imps$ranicept.l2))
rockchalk::skewness(imps$ranicept.l2)
rockchalk::kurtosis(imps$ranicept.l2)
plot(density(imps$ranslp.l2))
rockchalk::skewness(imps$ranslp.l2)
rockchalk::kurtosis(imps$ranslp.l2)

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
model <- "posaff ~ pain.cwc + sleep.cgm + pain.l2mean.cgm + painaccept.cgm + female + (1 + pain.cwc | person)"
analysis <- with(implist, lme4::lmer(model, REML = T))

# significance tests with barnard & rubin degrees of freedom
df <- 132 - 5 - 1
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = df)
estimates
confint(estimates, level = .95)

# test within- and between-cluster slope difference
slope.diff <- c("pain.l2mean.cgm - pain.cwc")
mitml::testConstraints(analysis, constraints = slope.diff, df.com = df)
