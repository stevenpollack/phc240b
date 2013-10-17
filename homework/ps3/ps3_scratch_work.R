# Malignant lesions are often excised in patients with melanoma. The survival and characteristics of 205 patients having undergone this surgical procedure at Odense University Hospital in Denmark is provided in the dataset melanoma. This dataset is available as part of the timereg package. We are interested in studying the distribution of the time from surgery until death. Consider patients alive at the end of follow-up or dead from other causes to be right-censored.

library(timereg)
data(melanoma)
attach(melanoma)

# (a) Plot the Kaplan-Meier estimator of the survival function along with (pointwise) 95% confidence intervals. Provide the associated estimate and 95% confidence interval of the 5-year survival probability. Also compute the pointwise probability of surviving at least 5 years given survival until 2 years.

library(survival)

# convert status to a binary situation, indicating
# whether a patient died or not.
is.dead <- rep(F,times=length(melanoma$status))
is.dead[which(melanoma$status %in% c(1,3))] <- T

# make survival object and fit K-M curve
mel.surv.obj <- Surv(time=melanoma$days,event=is.dead,type="right")
km.curve <- survfit(formula=mel.surv.obj~1,conf.int=0.95,conf.type='plain')

# plot curve
source("homework/ps3/ggsurv.R")
plot1a <- ggsurv(km.curve) + labs(x="Days", y=expression(paste(hat(S),"(",x,")",sep="")),title="Kaplan-Meier curve for melanoma data")
show(plot1a)

# get estimate of S(5 years), assuming right continous
five.yrs <- 365*5
indx.of.nearest.time <- max(which(km.curve$time <= five.yrs))
survival.estimate <- km.curve$surv[indx.of.nearest.time] # 73.2%
survival.estimate.ci <- c(km.curve$lower[indx.of.nearest.time],km.curve$upper[indx.of.nearest.time]) # (67.1%, 79.4%)

# if left continuous, use min(which(km.curve$time >= five.yrs))

# estimate P( S >= 5 | S >= 2):
two.yrs <- 365*2
indx.of.nearest.time2 <- max(which(km.curve$time <= two.yrs))
survival.estimate2 <- km.curve$surv[indx.of.nearest.time2]
conditional.estimate <- survival.estimate / survival.estimate2 # 81.6%

# (b) On a single figure, plot the Kaplan-Meier estimator of the survival function corresponding to the subgroups of patients with and without ulcerated lesions. Perform a logrank test to determine whether survivorship differs by ulceration status, report results and interpret your findings.

# stratify data and plot curves
km.curve2 <- survfit(formula=mel.surv.obj~ulc,conf.type='plain')
plot1b <- ggsurv(km.curve2,cens.col="black",CI=F) + scale_color_discrete(name="Ulcer status",breaks=c(0,1),labels=c("Absent","Present")) + scale_linetype_manual(guide=F,values=c(1,1)) + labs(x="Days",y=expression(paste(hat(S),"(",x,")",sep="")), title="K-M curves for patients with varying ulcer stati")
show(plot1b)

# do log-rank test
logrankTest.results <- survdiff(formula=mel.surv.obj~ulc,rho=0) # p-value < 1e-6
### there's a good reason to believe ulcerated people have different hazards.

# (c) Fit a Cox proportional hazards model for the survival time, including ulceration status, lesion thickness and gender as main-effect terms. Report a point estimate, standard error, 95% confidence interval and p-value for each parameter in your model. Carefully interpret each parameter estimate.

cox.model <- coxph(formula=mel.surv.obj ~ ulc + sex + thick)
summary(cox.model)

library(xtable)

# scrap data
beta <- cox.model$coefficients
se <- sqrt(diag(cox.model$var))
z <- qnorm(0.975)
pvalues <- 2*pnorm(beta/se,lower.tail=F) # two-sided alternative

# place in matrix for output
prob1c.out <- cbind(beta,se,beta - z*se, beta + z*se,pvalues)
colnames(prob1c.out) <- c("coefficient", "standard error", "lower 95%", "upper 95%", "p-value")
rownames(prob1c.out) <- c("Ulceration status", "Sex", "Lesion thickness")
prob1c.out

xtable(prob1c.out, digits=4)

# just changing a person's ulceration status from absent to present increases their hazard 2.6x; men have 1.59x higher hazard than women; as a tumor becomes .01mm thicker, a persons hazard increases by a factor 1.001.

detach(melanoma)