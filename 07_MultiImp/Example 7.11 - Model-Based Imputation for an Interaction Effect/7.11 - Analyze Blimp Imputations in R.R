# example 7.11: analyze model-based multiple imputations from blimp
# requires fdir and mitml packages

# set working directory
fdir::set()

# read imputed data from working directory
imps <- read.table("imps.dat")
names(imps) <- c("Imputation", "id", "txgrp", "male", "age", "edugroup", "workhrs", "exercise", "paingrps", 
                "pain", "anxiety", "stress", "control", "depress", "interfere", "disability",
                paste0("dep", seq(1:7)), paste0("int", seq(1:6)), paste0("dis", seq(1:6)))

# center lower-order variable and compute product of imputed variables
imps$depress.cgm <- imps$depress - mean(imps$depress)
imps$depress.by.male <- imps$depress.cgm * imps$male

# analysis and pooling
implist <- mitml::as.mitml.list(split(imps, imps$Imputation))
analysis <- with(implist, lm(disability ~ depress.cgm + male + depress.by.male + pain))
nullmodel <- with(implist, lm(disability ~ 1))

# significance tests with barnard & rubin degrees of freedom
mitml::testEstimates(analysis, extra.pars = T, df.com = 270)

# wald test that all slopes = 0
mitml::testModels(analysis, nullmodel, df.com = 270, method = "D1")

# likelihood ratio test that all slopes = 0
mitml::testModels(analysis, nullmodel, method = "D3")

# test conditional effects (simple slopes)
slope.male <- c("depress.cgm + depress.by.male*1")
mitml::testConstraints(analysis, constraints = slope.male, df.com = 270)

slope.female <- c("depress.cgm + depress.by.male*0")
mitml::testConstraints(analysis, constraints = slope.female, df.com = 270)
