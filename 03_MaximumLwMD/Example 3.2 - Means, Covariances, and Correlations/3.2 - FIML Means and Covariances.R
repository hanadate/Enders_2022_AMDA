# example 3.2: fiml means, variances, and correlations
# requires fdir and lavaan packages
library(fdir)
library(lavaan)
library(tidyverse)
# set working directory
fdir::set()

# read data from working directory
dat <- read.table("employee.dat", na.strings = "999")
names(dat) <- c("employee","team","turnover","male","empower","lmx","worksat","climate","cohesion")
dat %>% glimpse()
# turnover: Intend to quit job in the next 6 months
# empower: Employee empowerment composite
# lmx: Leader-member exchange (relationship quality w/ sv) composite
# worksat(jobsat): Work satisfaction rating
# climate: Leadership climate composite (team-level)
# cohesion: Team cohesion composite

# specify model
model <- '
  worksat ~~ empower
  worksat ~~ lmx
  empower ~~ lmx
'

# table3.2
# estimate model in lavaan
fit <- lavaan::sem(model, dat, meanstructure = TRUE, fixed.x = FALSE, missing = "fiml")
summary(fit, standardize = TRUE)

# missing data patterns and proportion observed data (coverage)
# lavInspect()
# table3.2
inspect(fit, "patterns")
inspect(fit, "coverage")
inspect(fit, "cor.all")
inspect(fit, "cov.all")

# table3.1
table31 <- dat %>% 
  dplyr::mutate(
    mis.pat=case_when(
      !is.na(worksat) & !is.na(empower) & !is.na(lmx) ~ "p1",
      !is.na(worksat) & !is.na(empower) & is.na(lmx) ~ "p2",
      !is.na(worksat) & is.na(empower) & !is.na(lmx) ~ "p3",
      is.na(worksat) & !is.na(empower) & !is.na(lmx) ~ "p4",
      is.na(worksat) & is.na(empower) & !is.na(lmx) ~ "p5",
      TRUE ~ "else"
    )
  ) %>% 
  dplyr::group_by(mis.pat) %>% 
  dplyr::summarise(n=n()) %>% 
  dplyr::mutate(freq=n/sum(n))
table31

