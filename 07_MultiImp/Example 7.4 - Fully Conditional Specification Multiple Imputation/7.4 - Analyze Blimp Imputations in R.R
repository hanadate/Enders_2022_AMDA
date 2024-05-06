# example 7.4: analyze fully conditional specification multiple imputations from blimp
# requires fdir and mitml packages

# set working directory
fdir::set()

# read imputed data from working directory
imps <- read.table("imps.dat")
names(imps) <- c("Imputation", "id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost")

# compute change score from imputed data
imps$change <- imps$mathpost - imps$mathpre

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
analysis <- with(implist, lm(change ~ 1))

# significance tests with barnard & rubin degrees of freedom
estimates <- mitml::testEstimates(analysis, extra.pars = T, df.com = 249)
estimates
confint(estimates, level = .95)