# example 8.5: analyze model-based multiple imputations for a 3-level regression with a cross-level interaction
# requires fdir, rockchalk, lme4, and mitml packages

# set working directory
fdir::set()

# read imputed data from working directory
imps <- read.table("imps.dat")
names(imps) <- c("Imputation","school","student","wave","condition","teachexp","eslpct","race",
                 "male","frlunch","achievegrp","stanmath","month0", "month7", "probsolve", "efficacy",
                 "ranicept.l2","ranslp.l2","ranicept.l3","ranslp.l3","frlunch.l2mean","frlunch.l3mean",
                 "stanmath.l3mean","residual.l1")

# center at grand means
imps$stanmath.cgm <- imps$stanmath - mean(imps$stanmath)
imps$frlunch.cgm <- imps$frlunch - mean(imps$frlunch)
imps$teachexp.cgm <- imps$teachexp - mean(imps$teachexp)

# compute product of imputed variables
imps$month.by.cond <- imps$month7 * imps$condition

# distribution of random intercept and slope residuals
plot(density(imps$ranicept.l2))
rockchalk::skewness(imps$ranicept.l2)
rockchalk::kurtosis(imps$ranicept.l2)
plot(density(imps$ranslp.l2))
rockchalk::skewness(imps$ranslp.l2)
rockchalk::kurtosis(imps$ranslp.l2)
plot(density(imps$ranicept.l3))
rockchalk::skewness(imps$ranicept.l3)
rockchalk::kurtosis(imps$ranicept.l3)
plot(density(imps$ranslp.l3))
rockchalk::skewness(imps$ranslp.l3)
rockchalk::kurtosis(imps$ranslp.l3)


# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
model <- "probsolve ~ month7  + stanmath.cgm + frlunch.cgm + teachexp.cgm + condition + month.by.cond + (1 + month7 | school/student)"
analysis <- with(implist, lme4::lmer(model, REML = T))

# significance tests with barnard & rubin degrees of freedom
df <- 29 - 6 - 1
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = df)
estimates
confint(estimates, level = .95)

# test conditional effects
growth.control <- c("month7 + month.by.cond*0")
mitml::testConstraints(analysis, constraints = growth.control, df.com = df)
growth.expgrp <- c("month7 + month.by.cond*1")
mitml::testConstraints(analysis, constraints = growth.expgrp, df.com = df)

