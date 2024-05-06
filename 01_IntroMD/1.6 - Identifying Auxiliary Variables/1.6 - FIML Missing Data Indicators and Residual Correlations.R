###
# note by S.H. as of Apr 25, 2024
# Edited to be enable to view easily TABLE 1.2 and 1.3 as you see them in the book of AMDA.
###

# identify auxiliary variables with missing data indicators and residual correlations

# load packages
library(fdir)
library(lavaan)
library(semTools)
library(tidyverse)

# set working directory and load data
# set()
dat <- read.table("pain.dat", na.strings = "999")
names(dat) <- c("id", "txgrp", "male", "age", "edugroup", "workhrs", "exercise", "paingrps", 
                "pain", "anxiety", "stress", "control", "depress", "interfere", "disability",
                paste0("dep", seq(1:7)), paste0("int", seq(1:6)), paste0("dis", seq(1:6)))
dat %>% glimpse

#################################################################
# IDENTIFY CORRELATES OF MISSINGNESS USING COHEN'S D EFFECT SIZE
#################################################################

# missing data indicators
dat$mis.depress <- is.na(dat$depress)
dat$mis.interf <- is.na(dat$interfere)
dat$mis.pain <- is.na(dat$pain)

# auxiliary variables regressed on missing data indicators
model.depress <- '
     age ~ mis.depress
     exercise ~ mis.depress
     anxiety ~ mis.depress
     stress ~ mis.depress
     control ~ mis.depress
     disability ~ mis.depress
'

model.interfere <- '
     age ~ mis.interf
     exercise ~ mis.interf
     anxiety ~ mis.interf
     stress ~ mis.interf
     control ~ mis.interf
     disability ~ mis.interf
'

model.pain <- '
     age ~ mis.pain
     exercise ~ mis.pain
     anxiety ~ mis.pain
     stress ~ mis.pain
     control ~ mis.pain
     disability ~ mis.pain
'

# estimate model in lavaan with standardized auxiliary variables (outcomes)
fit.depress <- sem(model.depress, dat, fixed.x = T, missing = "fiml")
est.depress <- standardizedSolution(fit.depress, type = "std.nox") %>% 
  dplyr::filter(rhs=="mis.depress"&lhs!="mis.depress")


# estimate model in lavaan with standardized auxiliary variables (outcomes)
fit.interfere <- sem(model.interfere, dat, fixed.x = T, missing = "fiml")
est.interfere <- standardizedSolution(fit.interfere, type = "std.nox") %>% 
  dplyr::filter(rhs=="mis.interf"&lhs!="mis.interf")

# estimate model in lavaan with standardized auxiliary variables (outcomes)
fit.pain <- sem(model.pain, dat, fixed.x = T, missing = "fiml")
est.pain <- standardizedSolution(fit.pain, type = "std.nox") %>% 
  dplyr::filter(rhs=="mis.pain"&lhs!="mis.pain")

# bind res: TABLE 1.2 is this..
res.binded_table1_2 <- dplyr::bind_rows(est.depress, est.interfere, est.pain) %>% 
  dplyr::select(lhs, rhs, est.std) %>% 
  dplyr::mutate(est.std=round(est.std,2)) %>% 
  tidyr::pivot_wider(names_from=rhs, values_from=est.std) %>% 
  dplyr::select(lhs, mis.pain, mis.interf, mis.depress)
res.binded_table1_2

#################################################################
# IDENTIFY AUXILIARY VARIABLES BASED ON RESIDUAL CORRELATIONS
#################################################################

# incomplete variables regressed on all other analysis variables
model.depress <- 'depress ~ interfere + pain'
model.interf <- 'interfere ~ pain + depress'
model.pain <- 'pain ~ depress + interfere'

# potential auxiliary variables
auxvars <- c("age","exercise","anxiety","stress","control","disability")

# estimate model with lavaan and semTools
fit.depress <- sem.auxiliary(model.depress, dat, fixed.x = F, missing = "fiml", aux = auxvars)
# summary(fit.depress, rsquare = T, standardize = T)
est.depress <- standardizedSolution(fit.depress, type="std.nox") %>% 
  dplyr::filter(rhs=="depress"&lhs!="depress")

# estimate model with lavaan and semTools
fit.interf <- sem.auxiliary(model.interf, dat, fixed.x = F, missing = "fiml", aux = auxvars)
# summary(fit.interf, rsquare = T, standardize = T)
est.interf <- standardizedsolution(fit.interf, type="std.nox") %>% 
  dplyr::filter(rhs=="interfere"&lhs!="interfere")

# estimate model with lavaan and semTools
fit.pain <- sem.auxiliary(model.pain, dat, fixed.x = F, missing = "fiml", aux = auxvars)
# summary(fit.pain, rsquare = T, standardize = T)
est.pain <- standardizedsolution(fit.pain, type="std.nox") %>% 
  dplyr::filter(rhs=="pain"&lhs!="pain")

# bind res: TABLE 1.3 is this.
res.binded_table1_3 <- dplyr::bind_rows(est.depress, est.interf, est.pain) %>% 
  dplyr::select(lhs, rhs, est.std, se, z, pvalue) %>% 
  dplyr::mutate(across(where(is.double), ~round(.x, 2)))
res.binded_table1_3
