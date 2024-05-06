# example 5.4: bayes regression with factored (sequential) specification for incomplete predictors
# requires fdir and mdmb packages

# set working directory
fdir::set()

# read data from working directory
dat <- read.table("employee.dat", na.strings = "999")
names(dat) <- c("employee","team","turnover","male","empower","lmx","worksat","climate","cohesion")

# summarize incomplete variables to determine ranges for pseudo-imputations
summary(dat[,c("empower","lmx","climate","male")])

# set ranges (nodes) for pseudo-imputations
nodes.lmx <- seq(-5, 25, by = 1)
nodes.climate <- seq(5, 40, by = 1)
nodes.empower <- seq(5, 50, by = 1)

# model for climate predictor (climate |male)
model.climate <- list( "model"="linreg", "formula" = climate ~ male, nodes = nodes.climate)

# model for lmx predictor (lmx | climate male)
model.lmx <- list( "model"="linreg", "formula" = lmx ~ climate + male, nodes = nodes.lmx)

# model for empower outcome (empower | male)
model.empower <- list( "model"="linreg", "formula" = empower ~ lmx + climate + male, nodes = nodes.empower)

# combine predictor models into a list
predictor.models <- list(climate = model.climate, lmx = model.lmx)

# estimate factored regression model w mdmb
fit <- mdmb::frm_fb(dat = dat, dep = model.empower, ind = predictor.models) 
summary(fit)

