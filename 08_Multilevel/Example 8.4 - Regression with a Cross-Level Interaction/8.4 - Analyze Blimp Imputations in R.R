# example 8.4: analyze model-based multiple imputations for a multilevel regression with a cross-level interaction
# requires fdir, rockchalk, lme4, and mitml packages

# set working directory
fdir::set()

# read imputed data from working directory
imps <- read.table("imps.dat")
names(imps) <- c("Imputation", "employee","team","turnover","male","empower","lmx","jobsat","climate","cohesion",
                "ranicept.l2", "ranslp.l2", "lmx.l2mean")

# center at grand means and imputed latent group means
imps$male.cgm <- imps$male - mean(imps$male)
imps$lmx.cwc <- imps$lmx - imps$lmx.l2mean
imps$cohesion.cgm <- imps$cohesion - mean(imps$cohesion)
imps$climate.cgm <- imps$climate - mean(imps$climate)

# compute product of imputed variables
imps$lmx.by.climate <- imps$lmx.cwc * imps$climate.cgm

# distribution of random intercept and slope residuals
plot(density(imps$ranicept.l2))
rockchalk::skewness(imps$ranicept.l2)
rockchalk::kurtosis(imps$ranicept.l2)
plot(density(imps$ranslp.l2))
rockchalk::skewness(imps$ranslp.l2)
rockchalk::kurtosis(imps$ranslp.l2)

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
model <- "empower ~ lmx.cwc  + male.cgm + cohesion.cgm + climate.cgm + lmx.by.climate + (1 + lmx.cwc | team)"
analysis <- with(implist, lme4::lmer(model, REML = T))

# significance tests with barnard & rubin degrees of freedom
df <- 105 - 5 - 1
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = df)
estimates
confint(estimates, level = .95)

# test conditional effects
climate.sd <- sd(imps$climate)
lmx.at.high.climate <- c("lmx.cwc + lmx.by.climate*climate.sd")
mitml::testConstraints(analysis, constraints = lmx.at.high.climate, df.com = df)

lmx.at.low.climate <- c("lmx.cwc + lmx.by.climate*-climate.sd")
mitml::testConstraints(analysis, constraints = lmx.at.low.climate, df.com = df)



