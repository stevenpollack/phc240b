# Lab 2-3. 
# * Installing flexsurv package
# * Reading datasets in R
# * Working with flexsurv package

#============================================================
# Installation (2 ways)
#============================================================
# (I) Easier. 
# Update R to version 3.0.1 by going to http://cran.cnr.berkeley.edu/
install.packages("flexsurv")	
library(flexsurv)

# (II) Harder
# Install required packages, download the source for flexsurv package and install it manually

install.packages("muhaz")
install.packages("mvtnorm")
install.packages("eha")
install.packages("survival")
install.packages("splines")

# REPLACE WITH THE PATH TO DIRECTORY OF YOUR DOWNLOADED FLEXSURV PACKAGE
PATH_flexsurv <- "/Users/monikaizano/Dropbox/F_13/GSI_PH240B_Survival"

pack_path <- paste(PATH_flexsurv, "/flexsurv_0.2.tar.gz", sep="")
install.packages(pack_path, repos = NULL, type="source")
library(flexsurv)


#============================================================
# Reading/Working w Datasets (http://data.princeton.edu/wws509/datasets)
#============================================================
# First install package "foreign" to be able to read STATA dataset
install.packages("foreign")
library(foreign)


#------------------------------------------------------
# Time to Ph.D.
#------------------------------------------------------
# year: coded 1 to 14, representing years of graduate school.
# university: coded 1 for Berkeley, 2 for Columbia, 3 for Princeton.
# residence: coded 1 for permanent residents, 2 for temporary residents.
# events: number of students graduating in this category.
# exposure: number of person-years of exposure to graduation in this category.


phd_dat <- read.dta("http://data.princeton.edu/wws509/datasets/phd.dta")
phd_dat
class(phd_dat)
# [1] "data.frame"

#------------------------------------------------------
# Gehan-Freirich Survival Data
#------------------------------------------------------
# data show the length of remission in weeks for two groups of leukemia patients, treated and control, and were analyzed by Cox in his original proportional hazards paper.
# Treatment: coded Treated (drug) or Control (placebo),
# Time: weeks of remission,
# Failure: coded 1 if a failure (relapse), 0 if censored
leuk_dat <- read.dta("http://data.princeton.edu/wws509/datasets/gehan.dta")
leuk_dat
class(leuk_dat)
# [1] "data.frame"

#------------------------------------------------------
# Ovarian cancer dataset that comes with flexsurv package
#------------------------------------------------------
data(ovarian)
ovarian

#------------------------------------------------------
# Fit parametric survival distributions to above dataset
#------------------------------------------------------

# Fit the Exponential model
fite <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="exp")
# Look at the fit
fite
# Plot the fit vs. KM
plot(fite)

# Fit the Weibull model
fitw <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="weibull")
# Look at the fit
fitw
# Plot the fit vs. KM
plot(fitw)
# Plot the fit vs. KM and exponential (blue)
plot(fitw)
lines(fite, col="blue", lwd.ci=1, lty.ci=1)
legend(900,1,c("Wei", "Exp"), lty=c(1,1,1),
	lwd=c(2,2,2),
	col=c("red","blue"))

# Fit the Generalized Gamma model
fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="gengamma")
# Look at the fit
fitg
plot(fitg)
# Plot the fit vs. KM, Exponential (blue), Weibull (green)
plot(fitg)
lines(fite, col="blue", lwd.ci=1, lty.ci=1)
lines(fitw, col="green", lwd.ci=1, lty.ci=1)
legend(900,1,c("Gen Gamma", "Exp","Wei"), lty=c(1,1,1),
	lwd=c(2,2,2),
	col=c("red","blue","green"))

## Identical AIC, probably not enough data in this simple example for a
## very flexible model to be worthwhile.
