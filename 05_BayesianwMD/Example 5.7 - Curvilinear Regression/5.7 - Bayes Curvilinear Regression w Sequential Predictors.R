# example 5.7: curvilinear regression with factored (sequential) specification for incomplete predictors
# requires fdir, lavaan, and mdmb packages

# set working directory
fdir::set()

# read data from working directory
dat <- read.table("math.dat", na.strings = "999")
names(dat) <- c("id", "male", "frlunch", "achievegrp", "stanread", "efficacy", "anxiety", "mathpre", "mathpost")

# specify lavaan model for sample statistics
model <- '
  mathpost ~~ anxiety
  mathpost ~~ frlunch
  mathpost ~~ mathpre
  mathpost ~~ male
  anxiety ~~ frlunch
  anxiety ~~ mathpre
  anxiety ~~ male
  frlunch ~~ mathpre
  frlunch ~~ male
  mathpre ~~ male
'

# estimate sample statistics in lavaan
lavaan::inspectSampleCov(model, dat, meanstructure = TRUE, missing = "fiml")

# summarize incomplete predictors to determine ranges for pseudo-imputations
summary(dat[,c("mathpost","anxiety","frlunch","mathpre","male")])

# set ranges (nodes) for pseudo-imputations
nodes.frlunch <- c(0,1)
nodes.anxiety <- seq(-10, 65, by = 1)
nodes.mathpost <- seq(25, 95, by = 1)

# model for frlunch predictor (frlunch | mathpre male)
model.frlunch <- list("model"="logistic", "formula" = frlunch ~ mathpre + male, nodes.frlunch)

# model for anxiety predictor (anxiety | frlunch mathpre male)
model.anxiety <- list( "model"="linreg", "formula" = anxiety ~ frlunch + mathpre + male, nodes = nodes.anxiety)

# model for mathpost outcome (mathpost | anxiety anxiety^2 frlunch mathpre male
# anxiety mean = 18.058
model.mathpost <- list("model" = "linreg", "formula" = mathpost ~ I(anxiety - 18.058) + I((anxiety - 18.058)^2) + frlunch + mathpre + male, nodes = nodes.mathpost)

# combine predictor models into a list
predictor.models <- list(frlunch = model.frlunch, anxiety = model.anxiety)

# estimate factored regression model w mdmb
fit <- mdmb::frm_fb(dat = dat, dep = model.mathpost, ind = predictor.models) 
summary(fit)