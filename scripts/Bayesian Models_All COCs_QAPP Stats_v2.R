# This script is a modified version of the "Bayesian Models_All COCs.R" script, specifically 
#  designed to run additional statistics required for the QAPP.  These include:

#  1. If the COC data includes potentially influential outliers, determine whether the model conclusions
#     are sensitive to these outliers.
#  2. Does changing the variance structure result in a quantitative change to the parameters for spatial
#     landscape predictors?
#  3. If our selected random effects structure is incorrect, does a simpler structure result in qualitative
#     differences to the model outcome?
#  4. Assess model precision by comparing the absolute distance between model estimates and the S8 data
#     used to fit the model (root mean square error -- RMSE)

# This script uses v10 lme model results (v12 for zinc; all are nested random effects: agency/location) to generate 
#  the Bayesian outputs.

# Author: Eva Dusek Jennings
# Revised: Oct 7, 2026
#---------------------------------------------------------------------------------------

#options(mc.cores = 2)

#if having trouble with installing packages, install them from binary, like this:
#  install.packages("igraph", type="binary")

#devtools::install_github("paul-buerkner/brms")
library(brms)
library(nlme)
library(ggplot2)
library(loo)
library(dplyr)
#devtools::install_github("rmcelreath/rethinking")  #this may not work b/c dependency "cmdstanr" isn't available

#install.packages("performance")
#library(performance)
#install.packages("HLMdiag")
library(HLMdiag)


#methods(class="brmsfit")  #complete list of methods available for brmsfit models

#----------------#
#  Total Copper  #
#----------------#

load(file="results/Frequentist_Copper Models.RData")
Cu.Form4  #lme model equation
Cu.r1X  #random effect in lme model
Cu.vf1X  #variance structure for lme model
Cu.lme <- lme(Cu.Form4, data=Cu.coc2, method="REML", random = Cu.r1X, weights=Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))


#------
# 1. Outliers
#------

par(mfrow=c(1,1))
plot(Cu.coc2$result)
#two possible high outliers, one possible low outlier

Cu.lme_inf <- hlm_influence(Cu.lme, level=1)
dotplot_diag(Cu.lme_inf$cooksd, name = "cooks.distance", cutoff = "internal")
Cu.coc2[c(198, 401, 46, 180, 483),]  #top 5 observations in terms of Cook's distance

Cu.coc3 <- Cu.coc2[-c(198),]  #remove the top one, which is particularly influential
Cu.lme.rmOutliers <- lme(Cu.Form4, data=Cu.coc3, method="REML", random = Cu.r1X, weights=Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

coefficients(Cu.lme)
coefficients(Cu.lme.rmOutliers)

fixed.effects(Cu.lme)
fixed.effects(Cu.lme.rmOutliers)

plot(Cu.lme)
plot(Cu.lme.rmOutliers)

round((abs(fixed.effects(Cu.lme) - fixed.effects(Cu.lme.rmOutliers) ) / fixed.effects(Cu.lme))*100, 1)

#while there is one possible outlier, it has a negligible effect on the landscape predictors and intercept (0.1 to 0.5% change)


#------
# End of 1
#------

#Copper lme model summary
Cu.lme <- lme(Cu.Form4, data=Cu.coc2, method="REML", random = Cu.r1X, weights=Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(Cu.lme)

par(mfrow=c(3,1))
Cu.lme0 <- lme(Cu.Form4, data=Cu.coc2, method="REML", random = Cu.r1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E.lme0 <- resid(object = Cu.lme0, type = "normalized")
plot(fitted(Cu.lme0), E.lme0, main="no var struct", xlab="fitted", ylab="std residuals", col="gray", pch=16)

E.lme <- resid(object = Cu.lme, type = "normalized")
plot(fitted(Cu.lme), E.lme, main="var cov = location", xlab="fitted", ylab="std residuals", col="gray", pch=16)

Cu.lme1 <- lme(Cu.Form4, data=Cu.coc2, method="REML", random = Cu.r1X, weights=varIdent(form= ~1|agency), control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E.lme1 <- resid(object = Cu.lme1, type = "normalized")
plot(fitted(Cu.lme1), E.lme1, main="var cov = agency", xlab="fitted", ylab="std residuals", col="gray", pch=16)
#var cov = agency or location are both improvements on no variance structure; var cov=location is the best.
#  explanation for variance covariate = agency is that some agencies have more diversity in their types of sites, so having the variance
#  covariate set to agency allows us to compensate for this (residual error E(ijk) would be expected to be higher for some 
#  agencies and lower for others)


#Bayesian Mixed Model - check various variance structures to see if the one selected by lme is the best
#   agency/location as nested random effect; no variance structure
# fit0 <- brm(formula= result ~ summer + rain + sqrt_traffic + devAge2 + (1|agency/location),
fit0 <- brm(formula= update(Cu.Form4, .~. + (1|agency/location)),
            data=Cu.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency)
fit1 <- brm(bf(formula= update(Cu.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
            data=Cu.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|location) -- best model from LME 
fit2 <- brm(bf(formula= update(Cu.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=Cu.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 2),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)


#model validation using approximate leave-one-out cross-validation
#loo package (developed by Vehtari, Gelman and Gabry (2017a, 2017b)) allows calculation of LOOIC, similar to AIC.
#  We are looking to make sure that the Pareto shape k parameter for each data point (used to test reliability and convergence
#  rate of the PSIS-based estimates) is below 0.7
fit0 <- add_criterion(fit0, criterion=c("loo"))
fit1 <- add_criterion(fit1, criterion=c("loo"))
fit2 <- add_criterion(fit2, criterion=c("loo"))
fit2 <- add_criterion(fit2, criterion=c("loo"), moment_match=TRUE)
loo_compare(fit0, fit1, fit2, criterion="loo")  #top one in the output gives the best model
loo_compare(fit0, fit1, criterion="loo")  #top one in the output gives the best model

# loo(fit0)  #highest looic is the best model; lowest elpd (expected log predictive density) is the best model
# loo(fit1)
# loo(fit2)

# PSIS diagnostic tool - look for points (potential influential outliers) above 0.5, and especially above 0.7
plot(loo(fit0, cores=getOption("mc.cores", 1)), main="no var struct")
plot(loo(fit1, cores=getOption("mc.cores", 1)), main="varcov=agency")
plot(loo(fit2, cores=getOption("mc.cores", 1)), main="varcov=location")

pareto_k_ids(loo(fit1), threshold=0.5)  #22, 120
pareto_k_ids(loo(fit2), threshold=0.5)  #109, 120, 186, 256, 274



### FOR SOME ABOVE MODELS, THERE WERE SOME PROBLEMS WITH PARETO-K VALUES BEING TOO HIGH --
#      TRY A STUDENT-T DISTRIBUTION RATHER THAN A GAUSSIAN (NORMAL) DISTRIBUTION
#by using the student-t distribution instead of the gaussian, we are allowing fatter tails, and the influential
#  observation (the highest SNO-HDR value) becomes less of an outlier and its pareto_k value is not longer > 0.7
fit2.t <- brm(bf(formula= update(Cu.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=Cu.coc2,
            family=student,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 2),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)



#------
# 2. Changing Variance Structure
#------

Cu.fit0 <- summary(fit0)$fixed
Cu.fit1 <- summary(fit1)$fixed
Cu.fit2 <- summary(fit2)$fixed

abs(Cu.fit0$Estimate - Cu.fit2$Estimate[-2]) / Cu.fit2$Estimate[-2]
abs(Cu.fit1$Estimate[-2] - Cu.fit2$Estimate[-2]) / Cu.fit2$Estimate[-2]
#comparison between [no variance structure] or [agency as variance covariate] and [location as variance covariate] shows 
#  1% or less difference in values for global intercept or landscape predictors (items 1, 4, 5)


#agency/location as nested random effect; no variance structure
fit0.t <- brm(formula= update(Cu.Form4, .~. + (1|agency/location)),
            data=Cu.coc2,
            family=student,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency)
fit1.t <- brm(bf(formula= update(Cu.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
            data=Cu.coc2,
            family=student,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20),
            cores = getOption("mc.cores", 2)
)

fit0.t <- add_criterion(fit0.t, criterion=c("loo"))
fit1.t <- add_criterion(fit1.t, criterion=c("loo"))
fit2.t <- add_criterion(fit2.t, criterion=c("loo"))
fit2.t <- add_criterion(fit2.t, criterion=c("loo"), moment_match=TRUE)

# PSIS diagnostic tool - look for points (potential influential outliers) above 0.5, and especially above 0.7
par(mfrow=c(3,1), oma=c(3,3,3,1))
plot(loo(fit0.t, cores=getOption("mc.cores", 1)), main="no var struct")
plot(loo(fit1.t, cores=getOption("mc.cores", 1)), main="varcov=agency")
plot(loo(fit2.t, cores=getOption("mc.cores", 1)), main="varcov=location")  #this one has one data point just over 0.7.  I don't love that...


#I explored removing the agency nesting piece, but I'd really have to do it for all COCs, and I don't think that will
#   hold for every single one.  The sd_agency_intercept is close to 1 (not 0, which would indicate removal is ok), 
#   while the elpd_diff is small (indicating I could consider removing it).  I'm not sure what effect agency has here, but should keep

loo_compare(fit0.t, fit1.t, fit2.t, criterion="loo")  #top one in the output gives the best model


Cu.fit0.t <- summary(fit0.t)$fixed
Cu.fit1.t <- summary(fit1.t)$fixed
Cu.fit2.t <- summary(fit2.t)$fixed
round((abs(Cu.fit0.t$Estimate - Cu.fit2.t$Estimate[-2]) / Cu.fit2.t$Estimate[-2])*100, 1)
round((abs(Cu.fit1.t$Estimate[-2] - Cu.fit2.t$Estimate[-2]) / Cu.fit2.t$Estimate[-2])*100, 1)
#compare models that ALL use a student-t distribution for residual variance

# abs(Cu.fit2.t$Estimate[-2] - Cu.fit2$Estimate[-2]) / Cu.fit2.t$Estimate[-2]
#comparison between the predictors assuming data are from a Gaussian vs Student-t distribution 
#  2.9% or less difference in values for global intercept or landscape predictors (items 1, 4, 5)

#------
# End of 2
#------


# fit2.t <- add_criterion(fit2.t, criterion=c("loo"))
# loo_compare(fit0, fit1, fit2.t, criterion="loo")  #top one in the output gives the best model
# #best fit is for fit2.t, the model with sigma = (~1|location), and student-t distribution for data

plot(loo(fit2.t, cores=getOption("mc.cores", 1)), main="varcov=location")
pareto_k_ids(loo(fit2.t), threshold=0.5)  #109, 186, 256
pareto_k_ids(loo(fit2.t), threshold=0.7)  #186


#look at the relative amount of variability within each location, and relative amount of variability within each agency
par(mfrow=c(2,1))
boxplot(result ~ location, data=Cu.coc2)
boxplot(result ~ agency, data=Cu.coc2)
#I don't think I'd use the bottom plot to argue for using agency as a variance covariate (as there are other things that
#  could contribute to the spread, such as amount of sqrt_traffic for an agencies 1 to 3 sites), but the top plot indicates
#  that variability within a location is similar (except for KIC_HDR) for all locations

summary(fit1.t)$fixed  
summary(fit2.t)$fixed


#look at residuals vs fitted values for each model - are any markedly better than others?
par(mfrow=c(4,1), mar=c(4,4,4,2))
resid.0t <- residuals(fit0.t, type="ordinary")
fitted.0t <- fitted(fit0.t, scale="response")
plot(resid.0t[,1] ~ fitted.0t[,1], ylab="residuals", xlab="fitted values", main="Model 0")

resid.1t <- residuals(fit1.t, type="ordinary")
fitted.1t <- fitted(fit1.t, scale="response")
plot(resid.1t[,1] ~ fitted.1t[,1], ylab="residuals", xlab="fitted values", main="Model 1")

resid.2t <- residuals(fit2.t, type="ordinary")
fitted.2t <- fitted(fit2.t, scale="response")
plot(resid.2t[,1] ~ fitted.2t[,1], ylab="residuals", xlab="fitted values", main="Model 2.t")
#these all look pretty similar...


#-------------------------------------------
#https://tem11010.github.io/regression_brms/

#graphical posterior predictive checking. Compare observed data to simulated data from the posterior predictive distribution. 
#  This is a density plot, where the observed y values are plotted with expected values from the posterior distribution
pp_check(fit1, ndraws=20)
pp_check(fit2.t, ndraws=20)
#both of these seem to be similar, just with different scales

#Look at the fit based on the grouping variable. Here are scatter-plots with the observed chemical concentrations (log scale) 
#  on the y-axis and the average model predictions (across all posterior samples) on the x-axis.
#  Red line is the 1:1 line, indicating perfect fit of model predictions to data.  Any locations where model doesn't fit?
pp_check(fit2.t, type = "scatter_avg_grouped", group = "location") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
#-------------------------------------------

Cu.brm <- fit2.t  #best Bayesian model for copper, so far


#------
# 3. No agency as a random effect
#------

fit2.t.RE_loc <- brm(bf(formula= update(Cu.Form4, .~. + (1|location)), 
                 sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
              data=Cu.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              cores = getOption("mc.cores", 2),
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

Cu.fit2.t.RE_loc <- summary(fit2.t.RE_loc)$fixed

Cu.fit2.t.RE_loc
Cu.fit2.t

round( (abs(Cu.fit2.t$Estimate[-2] - Cu.fit2.t.RE_loc$Estimate[-2]) / Cu.fit2.t$Estimate[-2]) * 100, 1)
#comparison between the predictors assuming random effects is LOCATION only (not nested inside agency)
#  0.1% difference in global intercept parameter, 8.2% difference in sqrt_traffic parameter, 2.8% difference in devAge2 parameter

#------
# End of 3
#------


#------
# 4. RMSE
#------

#load(file="../results/Bayesian_Copper.Rdata")

# 
# library(performance)
# performance_rmse(Cu.brm)
# 

Cu_resid2 <- (residuals(Cu.brm)[, 1])^2  #1. take the residuals & square them
Cu_RMSE <- sqrt(mean(Cu_resid2))  #2 & 3. take the mean of the squared residuals, then take the sqrt of the result
Cu_RMSE
# 0.5195612

#RSR = RMSE/SD
Cu_RSR <- Cu_RMSE / sd(Cu.coc2$result)
Cu_RSR
# 0.4630891

#------
# End of 4
#------


save(Cu.brm, file="results/Bayesian_Copper.Rdata")



#--------------------------#
#  Total Suspended Solids  #
#--------------------------#


load(file="results/Frequentist_TSS Models.RData")
TSS.Form4
TSS.r1X
TSS.vf1X
TSS.lme <- lme(TSS.Form4, data=TSS.coc2, method="REML", random = TSS.r1X, weights=TSS.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
#summary(TSS.lme)

#------
# 1. Outliers
#------

plot(TSS.coc2$result)
#two possible low outliers

TSS.lme_inf <- hlm_influence(TSS.lme, level=1)
dotplot_diag(TSS.lme_inf$cooksd, name = "cooks.distance", cutoff = "internal")
TSS.coc2[c(182, 270, 69, 360, 171),]  #top 5 observations in terms of Cook's distance

TSS.coc3 <- TSS.coc2[-c(182, 270),]  #remove the top two, which are particularly influential
TSS.lme.rmOutliers <- lme(TSS.Form4, data=TSS.coc3, method="REML", random = TSS.r1X, weights=TSS.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

# coefficients(TSS.lme)
# coefficients(TSS.lme.rmOutliers)

fixed.effects(TSS.lme)
fixed.effects(TSS.lme.rmOutliers)

# plot(TSS.lme)
# plot(TSS.lme.rmOutliers)

round(abs(fixed.effects(TSS.lme) - fixed.effects(TSS.lme.rmOutliers) ) / fixed.effects(TSS.lme) *100, 1)

#three potential influential outliers based on Cook's distance.  Effect on landscape predictors are low:
#  0% for global intercept, 2.3% for oldDev predictor, and 2.4% for sqrt_CO2_road

#------
# End of 1
#------



par(mfrow=c(3,1))
TSS.lme0 <- lme(TSS.Form4, data=TSS.coc2, method="REML", random = TSS.r1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E.lme0 <- resid(object = TSS.lme0, type = "normalized")
plot(fitted(TSS.lme0), E.lme0, main="no var struct", xlab="fitted", ylab="std residuals", col="gray", pch=16)

E.lme <- resid(object = TSS.lme, type = "normalized")
plot(fitted(TSS.lme), E.lme, main="var cov = location", xlab="fitted", ylab="std residuals", col="gray", pch=16)

TSS.lme1 <- lme(TSS.Form4, data=TSS.coc2, method="REML", random = TSS.r1X, weights=varIdent(form= ~1|agency), control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E.lme1 <- resid(object = TSS.lme1, type = "normalized")
plot(fitted(TSS.lme1), E.lme1, main="var cov = agency", xlab="fitted", ylab="std residuals", col="gray", pch=16)

AIC(TSS.lme0, TSS.lme, TSS.lme1) #AIC also best for var cov=location
#var cov = location looks best, and also has lowest AIC
#  explanation for variance covariate = agency is that some agencies have more diversity in their types of sites, so having the variance
#  covariate set to agency allows us to compensate for this (residual error E(ijk) would be expected to be higher for some 
#  agencies and lower for others)


#fit1.t <- brm(bf(formula= update(Cu.Form4, .~. + (1|agency/location)), 
                 

#Bayesian Mixed Model - check various variance structures to see if the one selected by lme is the best
#   agency/location as nested random effect; no variance structure
fit0 <- brm(bf(formula= update(TSS.Form4, .~. + (1|agency/location))), 
            data=TSS.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency)
fit1 <- brm(bf(formula= update(TSS.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
            data=TSS.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 2),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|location) -- best model from LME 
fit2 <- brm(bf(formula= update(TSS.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=TSS.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 2),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)


#model validation using approximate leave-one-out cross-validation
#loo package (developed by Vehtari, Gelman and Gabry (2017a, 2017b)) allows calculation of LOOIC, similar to AIC.
#  We are looking to make sure that the Pareto shape k parameter for each data point (used to test reliability and convergence
#  rate of the PSIS-based estimates) is below 0.7
fit0 <- add_criterion(fit0, criterion=c("loo"))
fit1 <- add_criterion(fit1, criterion=c("loo"))
fit2 <- add_criterion(fit2, criterion=c("loo"))
fit2 <- add_criterion(fit2, criterion=c("loo"), moment_match=TRUE)
loo_compare(fit0, fit1, fit2, criterion="loo")  #top one in the output gives the best model

# PSIS diagnostic tool - look for points (potential influential outliers) above 0.5, and especially above 0.7
plot(loo(fit0, cores=getOption("mc.cores", 1)), main="no var struct")
plot(loo(fit1, cores=getOption("mc.cores", 1)), main="varcov=agency")
plot(loo(fit2, cores=getOption("mc.cores", 1)), main="varcov=location")

pareto_k_ids(loo(fit1), threshold=0.5)  #none!
pareto_k_ids(loo(fit2), threshold=0.5)  #5 different samples
pareto_k_ids(loo(fit2), threshold=0.7)  #none



### FOR SOME ABOVE MODELS, THERE WERE SOME PROBLEMS WITH PARETO-K VALUES BEING TOO HIGH --
#      TRY A STUDENT-T DISTRIBUTION RATHER THAN A GAUSSIAN (NORMAL) DISTRIBUTION
#by using the student-t distribution instead of the gaussian, we are allowing fatter tails, and the influential
#  observation (the highest SNO-HDR value) becomes less of an outlier and its pareto_k value is not longer > 0.7
fit2.t <- brm(bf(formula= update(TSS.Form4, .~. + (1|agency/location)), 
                 sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
              data=TSS.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              cores = getOption("mc.cores", 2),
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)


fit2.t <- add_criterion(fit2.t, criterion=c("loo"))
pareto_k_ids(loo(fit2.t), threshold=0.5)  #one sample
loo_compare(fit2, fit2.t, criterion="loo")  #top one in the output gives the best model


#------
# 2. Changing Variance Structure
#------

#agency/location as nested random effect; no variance structure
fit0.t <- brm(formula= update(TSS.Form4, .~. + (1|agency/location)),
              data=TSS.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              cores = getOption("mc.cores", 2),
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency)
fit1.t <- brm(bf(formula= update(TSS.Form4, .~. + (1|agency/location)), 
                 sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
              data=TSS.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              cores = getOption("mc.cores", 2),
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

fit0.t <- add_criterion(fit0.t, criterion=c("loo"))
fit1.t <- add_criterion(fit1.t, criterion=c("loo"))
fit2.t <- add_criterion(fit2.t, criterion=c("loo"))

# PSIS diagnostic tool - look for points (potential influential outliers) above 0.5, and especially above 0.7
plot(loo(fit0.t, cores=getOption("mc.cores", 1)), main="no var struct")
plot(loo(fit1.t, cores=getOption("mc.cores", 1)), main="varcov=agency")
plot(loo(fit2.t, cores=getOption("mc.cores", 1)), main="varcov=location")

TSS.fit0.t <- summary(fit0.t)$fixed
TSS.fit1.t <- summary(fit1.t)$fixed
TSS.fit2.t <- summary(fit2.t)$fixed
round((abs(TSS.fit0.t$Estimate - TSS.fit2.t$Estimate[-2]) / TSS.fit2.t$Estimate[-2])*100, 1)
round((abs(TSS.fit1.t$Estimate[-2] - TSS.fit2.t$Estimate[-2]) / TSS.fit2.t$Estimate[-2])*100, 1)
#compare models that ALL use a student-t distribution for residual variance

# abs(TSS.fit2.t$Estimate[-2] - TSS.fit2$Estimate[-2]) / TSS.fit2.t$Estimate[-2]
# #comparison between the predictors assuming data are from a Gaussian vs Student-t distribution 
# #  2.9% or less difference in values for global intercept or landscape predictors (items 1, 4, 5)

#------
# End of 2
#------


fit2.t <- add_criterion(fit2.t, criterion=c("loo"))
loo_compare(fit0.t, fit1.t, fit2.t, criterion="loo")  #top one in the output gives the best model
#best fit is for fit2.t, the model with sigma = (~1|location), and student-t distribution for data

plot(loo(fit2.t, cores=getOption("mc.cores", 1)), main="varcov=location")
pareto_k_ids(loo(fit2.t), threshold=0.5)  #353
pareto_k_ids(loo(fit2.t), threshold=0.7)  #no values above 0.7

#look at residuals vs fitted values for each model - are any markedly better than others?
par(mfrow=c(4,1), mar=c(4,4,4,2))
resid.0 <- residuals(fit0, type="ordinary")
fitted.0 <- fitted(fit0, scale="response")
plot(resid.0[,1] ~ fitted.0[,1], ylab="residuals", xlab="fitted values", main="Model 0")

resid.1 <- residuals(fit1, type="ordinary")
fitted.1 <- fitted(fit1, scale="response")
plot(resid.1[,1] ~ fitted.1[,1], ylab="residuals", xlab="fitted values", main="Model 1")

resid.2 <- residuals(fit2, type="ordinary")
fitted.2 <- fitted(fit2, scale="response")
plot(resid.2[,1] ~ fitted.2[,1], ylab="residuals", xlab="fitted values", main="Model 2")

resid.2t <- residuals(fit2.t, type="ordinary")
fitted.2t <- fitted(fit2.t, scale="response")
plot(resid.2t[,1] ~ fitted.2t[,1], ylab="residuals", xlab="fitted values", main="Model 2")
#these plots all look remarkably similar - all with a bit of tapering at high TSS values


#Bayesian Mixed Model; agency/location as random effect; variance structure = varIdent(form = ~1|location)
TSS.brm <- fit2.t

#graphical posterior predictive checking. Compare observed data to simulated data from the posterior predictive distribution. 
#  This is a density plot, where the observed y values are plotted with expected values from the posterior distribution
pp_check(TSS.brm, ndraws=200)

#Look at the fit based on the grouping variable. Here are scatter-plots with the observed chemical concentrations (log scale) 
#  on the y-axis and the average model predictions (across all posterior samples) on the x-axis.
#  Red line is the 1:1 line, indicating perfect fit of model predictions to data.  Any locations where model doesn't fit?
pp_check(TSS.brm, type = "scatter_avg_grouped", group = "location") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
#-------------------------------------------


#------
# 3. No agency as a random effect
#------

fit2.t.RE_loc <- brm(bf(formula= update(TSS.Form4, .~. + (1|location)), 
               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=TSS.coc2,
            family=student,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 2),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

TSS.fit2.t.RE_loc <- summary(fit2.t.RE_loc)$fixed

TSS.fit2.t.RE_loc
TSS.fit2.t


round((abs(TSS.fit2.t$Estimate[-2] - TSS.fit2.t.RE_loc$Estimate[-2]) / TSS.fit2.t$Estimate[-2])*100, 1)
#comparison between the predictors assuming random effects is LOCATION only (not nested inside agency)
#  0.2% difference in global intercept parameter, 10.5% difference in sqrt_traffic parameter, 1.1% difference in devAge2 parameter

#------
# End of 3
#------


#------
# 4. RMSE
#------

#load(file="../results/Bayesian_TSS.Rdata")

TSS_resid2 <- (residuals(TSS.brm)[, 1])^2  #1. take the residuals & square them
TSS_RMSE <- sqrt(mean(TSS_resid2))  #2 & 3. take the mean of the squared residuals, then take the sqrt of the result
TSS_RMSE
# 0.9600283

#RSR = RMSE/SD
TSS_RSR <- TSS_RMSE / sd(TSS.coc2$result)
TSS_RSR
# 0.7770575

#------
# End of 4
#------

save(TSS.brm, file="results/Bayesian_TSS.Rdata")




#--------------------#
#  Total Phosphorus  #
#--------------------#

load(file="results/Frequentist_Total Phosphorus Models.RData")
P.Form4
P.r1X
P.vf1X
#  Phosphorus lme model doesn't have too strong of a landscape predictor; try a model with only rain + summer
P.lme <- lme(P.Form4, data=P.coc2, method="REML", random = P.r1X, weights=P.vf1X)
#P.lme <- lme(P.Form4, data=P.coc2, method="ML", random = P.r1X, weights=P.vf1X)

#------
# 1. Outliers
#------

plot(P.coc2$result)
#hmmmm.... maybe a low outlier??

# #identify the location of the lowest point/ potential outlier:
# which(P.coc2$result==min(P.coc2$result))  #its row 213...
# 
# #Total Phosphorus lme model summary -- one outlier removed (row 213)
# P.lme.outlierRemoved <- lme(P.Form4, data=P.coc2[-213,], method="REML", random = P.r1X, weights=P.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(P.lme.outlierRemoved)
# 
# #Phosphorus lme model summary -- all data points
# P.lme <- lme(P.Form4, data=P.coc2, method="REML", random = P.r1X, weights=P.vf1X)
# summary(P.lme)
# 
# #fixed coefficients for the lme model with all data points, vs the one with the outlier removed
# P.lme$coefficients$fixed
# P.lme.outlierRemoved$coefficients$fixed
# abs(P.lme$coefficients$fixed - P.lme.outlierRemoved$coefficients$fixed)/P.lme$coefficients$fixed  #ratio of difference btwn two models/original model
# #fixed effects for the model with all data points differs from the one with the outlier removed by 1.4% or less.
# #  the greatest difference is in rain fixed effect (4%) followed by summer (1.4%), sqrt_CO2_road (1.3%) and global intercept (0.1%)


P.lme_inf <- hlm_influence(P.lme, level=1)
dotplot_diag(P.lme_inf$cooksd, name = "cooks.distance", cutoff = "internal")
P.coc2[c(342, 102, 70, 120, 31),]  #top 5 observations in terms of Cook's distance
P.coc3 <- P.coc2[-c(342, 102, 70),]  #remove the top three, which looks particularly influential
P.lme.rmOutliers <- lme(P.Form4, data=P.coc3, method="REML", random = P.r1X, weights=P.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

coefficients(P.lme)
coefficients(P.lme.rmOutliers)

fixed.effects(P.lme)
fixed.effects(P.lme.rmOutliers)

plot(P.lme)
plot(P.lme.rmOutliers)

round(abs(fixed.effects(P.lme) - fixed.effects(P.lme.rmOutliers) ) / fixed.effects(P.lme) * 100, 1)


#------
# End of 1
#------


par(mfrow=c(3,1))
P.lme0 <- lme(P.Form4, data=P.coc2, method="REML", random = P.r1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E.lme0 <- resid(object = P.lme0, type = "normalized")
plot(fitted(P.lme0), E.lme0, main="no var struct", xlab="fitted", ylab="std residuals", col="gray", pch=16)

E.lme <- resid(object = P.lme, type = "normalized")
plot(fitted(P.lme), E.lme, main="var cov = location", xlab="fitted", ylab="std residuals", col="gray", pch=16)

P.lme1 <- lme(P.Form4, data=P.coc2, method="REML", random = P.r1X, weights=varIdent(form= ~1|agency), control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E.lme1 <- resid(object = P.lme1, type = "normalized")
plot(fitted(P.lme1), E.lme1, main="var cov = agency", xlab="fitted", ylab="std residuals", col="gray", pch=16)

AIC(P.lme0, P.lme, P.lme1) #AIC best for var cov=location
#var cov = location looks best, and also has lowest AIC
#  explanation for variance covariate = agency is that some agencies have more diversity in their types of sites, so having the variance
#  covariate set to agency allows us to compensate for this (residual error E(ijk) would be expected to be higher for some 
#  agencies and lower for others)


#Bayesian Mixed Model - check various variance structures to see if the one selected by lme is the best

fit0 <- brm(bf(formula= update(P.Form4, .~. + (1|agency/location))), 
#               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=P.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 2),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency)
fit1 <- brm(bf(formula= update(P.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
            data=P.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 2),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|location) -- best model from LME 
fit2 <- brm(bf(formula= update(P.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=P.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 2),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)


#------
# 2. Changing Variance Structure
#------

P.fit0 <- summary(fit0)$fixed
P.fit1 <- summary(fit1)$fixed
P.fit2 <- summary(fit2)$fixed

round(abs(P.fit0$Estimate - P.fit2$Estimate[-2]) / P.fit2$Estimate[-2] * 100, 1)
round(abs(P.fit1$Estimate[-2] - P.fit2$Estimate[-2]) / P.fit2$Estimate[-2] * 100, 1)
#comparison between [no variance structure] or [agency as variance covariate] and [location as variance covariate] shows 
#  0.2 to 0.3% difference in global intercept and 0.3 to 1.8% difference in sqrt_CO2 road predictor (items 1, 4)

#------
# End of 2
#------


#compare the Leave One Out (loo) criterion for these three models; LOO is a cross-validation technique to validate the model
fit2 <- add_criterion(fit2, criterion=c("loo")) #, moment_match=TRUE)
fit2 <- add_criterion(fit2, criterion=c("loo"), moment_match=TRUE)
plot(loo(fit2, cores=getOption("mc.cores", 1)))  #two points above 0.7
pareto_k_ids(loo(fit2), threshold=0.5)  #4 samples
pareto_k_ids(loo(fit2), threshold=0.7)  #1 sample 

#single predictor (sqrt_CO2_road) with student-t distribution
fit2.t <- brm(bf(formula= update(P.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=P.coc2,
            family=student,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 1),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#------
# 2. Changing Variance Structure
#------

#agency/location as nested random effect; no variance structure
fit0.t <- brm(bf(formula= update(P.Form4, .~. + (1|agency/location))), 
              data=P.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency)
fit1.t <- brm(bf(formula= update(P.Form4, .~. + (1|agency/location)), 
                 sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
              data=P.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

fit0.t <- add_criterion(fit0.t, criterion=c("loo"))
fit1.t <- add_criterion(fit1.t, criterion=c("loo"))
fit2.t <- add_criterion(fit2.t, criterion=c("loo"))

# PSIS diagnostic tool - look for points (potential influential outliers) above 0.5, and especially above 0.7
par(mfrow=c(3,1))
plot(loo(fit0.t, cores=getOption("mc.cores", 1)), main="no var struct")
plot(loo(fit1.t, cores=getOption("mc.cores", 1)), main="varcov=agency")
plot(loo(fit2.t, cores=getOption("mc.cores", 1)), main="varcov=location")

P.fit0.t <- summary(fit0.t)$fixed
P.fit1.t <- summary(fit1.t)$fixed
P.fit2.t <- summary(fit2.t)$fixed
round((abs(P.fit0.t$Estimate - P.fit2.t$Estimate[-2]) / P.fit2.t$Estimate[-2])*100, 1)
round((abs(P.fit1.t$Estimate[-2] - P.fit2.t$Estimate[-2]) / P.fit2.t$Estimate[-2])*100, 1)
#compare models that ALL use a student-t distribution for residual variance

pareto_k_ids(loo(fit2.t), threshold=0.5)  #1 samples
pareto_k_ids(loo(fit2.t), threshold=0.7)  #0 samples 
#student-t distribution better than Gaussian distribution

# abs(P.fit2$Estimate[-2] - P.fit2.t$Estimate[-2]) / P.fit2.t$Estimate[-2]
#comparison between the predictors assuming data are from a Gaussian vs Student-t distribution 
#  0.4% difference in global intercept and 2.1% difference in sqrt_CO2 road predictor (items 1, 4)

#------
# End of 2
#------


fit2.t <- add_criterion(fit2.t, criterion=c("loo"))
plot(loo(fit2.t, cores=getOption("mc.cores", 1)))  #one point above 0.5 (below 0.7)
pareto_k_ids(loo(fit2.t), threshold=0.5)  #213
#t-distribution is a definite improvement on the fit

loo_compare(fit2, fit2.t, criterion="loo")




#look at residuals vs fitted values for candidate model - are any markedly better than others?
par(mfrow=c(2,1), mar=c(4,4,4,2))
resid.2t <- residuals(fit2.t, type="ordinary")
fitted.2t <- fitted(fit2.t, scale="response")
plot(resid.2t[,1] ~ fitted.2t[,1], ylab="residuals", xlab="fitted values", main="Model 2 - student's t distribution")

resid.2 <- residuals(fit2, type="ordinary")
fitted.2 <- fitted(fit2, scale="response")
plot(resid.2[,1] ~ fitted.2[,1], ylab="residuals", xlab="fitted values", main="Model 2 - Gaussian distribution")
#these plots look about the same...

#Bayesian Mixed Model; oldDev + agency/location as random effect; variance structure = varIdent(form = ~1|location)
P.brm <- fit2.t

#Look at the fit based on the grouping variable. Here are scatter-plots with the observed chemical concentrations (log scale) 
#  on the y-axis and the average model predictions (across all posterior samples) on the x-axis.
#  Red line is the 1:1 line, indicating perfect fit of model predictions to data.  Any locations where model doesn't fit?
pp_check(P.brm, type = "scatter_avg_grouped", group = "location") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)

# pp_check(fit2.noPred.t, type = "scatter_avg_grouped", group = "location") + 
#   geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
# #either way we go, the fit at TAC_IND isn't good.  High values are under-predicted, low values are over-predicted


#------
# 3. No agency as a random effect
#------

#what is the AIC for comparing best model with random effects = location only?
P.r2 <- formula( ~ 1 | location)  #non-nested random effects -- location only
P.lme.nonNestedRE <- lme(P.Form4, data=P.coc2, method="REML", random = P.r2, weights=varIdent(form= ~1|location), control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
AIC(P.lme, P.lme.nonNestedRE) #AIC best for var cov=location

fit2.t.RE_loc <- brm(bf(formula= update(P.Form4, .~. + (1|location)), 
                      sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
                   data=P.coc2,
                   prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                   cores = getOption("mc.cores", 4),
                   control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

summary(fit2.t.RE_loc)
# plot(fit2.t.RE_loc, variable=c("b_Intercept", "b_rain", "b_summer1", "b_sqrt_CO2_road"))  #recognizable covariates)
#plot(P.brm, 

P.fit2.t.RE_loc <- summary(fit2.t.RE_loc)$fixed

P.fit2.t.RE_loc
P.fit2.t

round(abs(P.fit2.t$Estimate[-2] - P.fit2.t.RE_loc$Estimate[-2]) / P.fit2.t$Estimate[-2] *100, 1)
#comparison between the predictors assuming random effects is LOCATION only (not nested inside agency)

plot(fit2.t)
fit2.t.RE_loc <- add_criterion(fit2.t.RE_loc, criterion=c("loo"))
fit2.t.RE_loc <- add_criterion(fit2.t.RE_loc, criterion=c("loo"), moment_match=TRUE)
loo_compare(fit2.t, fit2.t.RE_loc, criterion="loo")

# to decide whether agency can be removed as the nesting random effect, look at several things:
#  1.) plot of sd_agency__intercept: is it close to 0 or does it comfortably overlap 0?  (here, sd_agency__intercept is btwn 0 and 1)
#  2.) is elpd_diff small (within 1-2 std errors)?  (here, it is 5)
# looks like we should keep agency in there.

#------
# End of 3
#------


#------
# 4. RMSE
#------

#load(file="../results/Bayesian_Phosphorus.Rdata")

P_resid2 <- (residuals(P.brm)[, 1])^2  #1. take the residuals & square them
P_RMSE <- sqrt(mean(P_resid2))  #2 & 3. take the mean of the squared residuals, then take the sqrt of the result
P_RMSE
# 0.6213033

#RSR = RMSE/SD
P_RSR <- P_RMSE / sd(P.coc2$result)
P_RSR
# 0.6564339

#------
# End of 4
#------

save(P.brm, file="results/Bayesian_Phosphorus.Rdata")

#load(file="../results/Bayesian_Phosphorus.RData")


#--------------#
#  Total Zinc  #
#--------------#

load(file="results/Frequentist_Total Zinc Models.RData")
totZn.Form4
totZn.r1X
totZn.vf1X

# #remove any zinc values over 800 ug/L (these should have been removed in the frequentist model stage also)
# totZn.coc2 <- totZn.coc2 %>%
#   filter(!result > log(800))
# max(exp(totZn.coc2$result))

totZn.lme <- lme(totZn.Form4, data=totZn.coc2, method="REML", random = totZn.r1X, weights=totZn.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

#------
# 1. Outliers
#------

plot(totZn.coc2$result)

totZn.lme_inf <- hlm_influence(totZn.lme, level=1)
dotplot_diag(totZn.lme_inf$cooksd, name = "cooks.distance", cutoff = "internal")
totZn.coc2[c(253, 437, 249, 192, 146),]  #top 5 observations in terms of Cook's distance

totZn.coc3 <- totZn.coc2[-c(253, 437, 249, 192, 146),]  #remove the top five, all of which are particularly influential
totZn.lme.rmOutliers <- lme(totZn.Form4, data=totZn.coc3, method="REML", random = totZn.r1X, weights=totZn.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

coefficients(totZn.lme)
coefficients(totZn.lme.rmOutliers)

fixed.effects(totZn.lme)
fixed.effects(totZn.lme.rmOutliers)

plot(totZn.lme)
plot(totZn.lme.rmOutliers)

round(abs(fixed.effects(totZn.lme) - fixed.effects(totZn.lme.rmOutliers) ) / fixed.effects(totZn.lme) * 100, 1)
#these five outliers are not very influential.  Effect on landscape predictor is low:


#------
# End of 1
#------


#compare models with simpler variance structures (no variance covariate, agency=var cov, location=var cov)
par(mfrow=c(3,1))
totZn.lme0 <- lme(totZn.Form4, data=totZn.coc2, method="REML", random = totZn.r1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E.lme0 <- resid(object = totZn.lme0, type = "normalized")
plot(fitted(totZn.lme0), E.lme0, main="no var struct", xlab="fitted", ylab="std residuals", col="gray", pch=16)

E.lme <- resid(object = totZn.lme, type = "normalized")
plot(fitted(totZn.lme), E.lme, main="var cov = location", xlab="fitted", ylab="std residuals", col="gray", pch=16)

totZn.lme1 <- lme(totZn.Form4, data=totZn.coc2, method="REML", random = totZn.r1X, weights=varIdent(form= ~1|agency), control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E.lme1 <- resid(object = totZn.lme1, type = "normalized")
plot(fitted(totZn.lme1), E.lme1, main="var cov = agency", xlab="fitted", ylab="std residuals", col="gray", pch=16)

AIC(totZn.lme0, totZn.lme, totZn.lme1) #AIC best for var cov=location
#var cov = location looks best, and also has lowest AIC
#  explanation for variance covariate = agency is that some agencies have more diversity in their types of sites, so having the variance
#  covariate set to agency allows us to compensate for this (residual error E(ijk) would be expected to be higher for some 
#  agencies and lower for others)


#Bayesian Mixed Model - check various variance structures to see if the one selected by lme is the best
fit0 <- brm(bf(formula= update(totZn.Form4, .~. + (1|agency/location))), 
            #               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=totZn.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            cores = getOption("mc.cores", 4),
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency)
fit1 <- brm(bf(formula= update(totZn.Form4, .~. + (1|agency/location)), 
               sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
            data=totZn.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20),
            cores = getOption("mc.cores", 4),
            save_pars=save_pars(all=TRUE)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|location) -- best model from LME 
fit2 <- brm(bf(formula= update(totZn.Form4, .~. + (1|agency/location)),                
               sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
            data=totZn.coc2,
            prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
            control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20),
            cores = getOption("mc.cores", 4),
            save_pars=save_pars(all=TRUE)
)
#there was 1 divergent transition after warmup.  Try t-distribution to check for improvement

fit0 <- add_criterion(fit0, criterion=c("loo"))
fit1 <- add_criterion(fit1, criterion=c("loo"))
fit2 <- add_criterion(fit2, criterion=c("loo"))
fit2 <- add_criterion(fit2, criterion=c("loo"), moment_match=TRUE)
loo_compare(fit0, fit1, fit2, criterion="loo")

# PSIS diagnostic tool - look for points (potential influential outliers) above 0.5, and especially above 0.7
plot(loo(fit0, cores=getOption("mc.cores", 1)), main="no var struct")
plot(loo(fit1, cores=getOption("mc.cores", 1)), main="varcov=agency")
plot(loo(fit2, cores=getOption("mc.cores", 1)), main="varcov=location")

pareto_k_ids(loo(fit2), threshold=0.5)  #6 samples
pareto_k_ids(loo(fit2), threshold=0.7)  #3 samples
#try t-distribution for improved fit


# #------
# # 2. Changing Variance Structure   ### NOTE THAT THE BEST FIT ZINC MODEL USES A STUDENT T-DISTRIBUTION RATHER THAN NORMAL!
# #------
# 
# totZn.fit0 <- summary(fit0)$fixed
# totZn.fit1 <- summary(fit1)$fixed
# totZn.fit2 <- summary(fit2)$fixed
# 
# abs(totZn.fit0$Estimate - totZn.fit2$Estimate[-2]) / totZn.fit2$Estimate[-2]
# abs(totZn.fit1$Estimate[-2] - totZn.fit2$Estimate[-2]) / totZn.fit2$Estimate[-2]
# #comparison between [no variance structure] or [agency as variance covariate] and [location as variance covariate] shows 
# #  0.07% to 0.3% difference in global intercept (item 1), 0.9 to 3.2% difference in not_greenBE (item 4), and 1.0 to 0.1% difference in sqrt_CO2_transport (item 5)
# 
# #------
# # End of 2
# #------


#agency/location as nested random effect; variance structure = varIdent(form = ~1|location) -- student t distribution
fit2.t <- brm(bf(formula= update(totZn.Form4, .~. + (1|agency/location)), 
                 sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
              data=totZn.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              cores = getOption("mc.cores", 4),
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)


#------
# 2. Changing Variance Structure
#------

#agency/location as nested random effect; no variance structure
fit0.t <- brm(bf(formula= update(totZn.Form4, .~. + (1|agency/location))), 
              data=totZn.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              cores = getOption("mc.cores", 4),
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency)
fit1.t <- brm(bf(formula= update(totZn.Form4, .~. + (1|agency/location)), 
                 sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
              data=totZn.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              cores = getOption("mc.cores", 4),
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

fit0.t <- add_criterion(fit0.t, criterion=c("loo"))
fit1.t <- add_criterion(fit1.t, criterion=c("loo"))
fit2.t <- add_criterion(fit2.t, criterion=c("loo"))

# PSIS diagnostic tool - look for points (potential influential outliers) above 0.5, and especially above 0.7
plot(loo(fit0.t, cores=getOption("mc.cores", 1)), main="no var struct")
plot(loo(fit1.t, cores=getOption("mc.cores", 1)), main="varcov=agency")
plot(loo(fit2.t, cores=getOption("mc.cores", 1)), main="varcov=location")

pareto_k_ids(loo(fit2.t), threshold=0.5)  #2 samples
pareto_k_ids(loo(fit2.t), threshold=0.7)  #none!
#use t-distribution

totZn.fit0.t <- summary(fit0.t)$fixed
totZn.fit1.t <- summary(fit1.t)$fixed
totZn.fit2.t <- summary(fit2.t)$fixed
round((abs(totZn.fit0.t$Estimate - totZn.fit2.t$Estimate[-2]) / totZn.fit2.t$Estimate[-2])*100, 1)
round((abs(totZn.fit1.t$Estimate[-2] - totZn.fit2.t$Estimate[-2]) / totZn.fit2.t$Estimate[-2])*100, 1)
#compare models that ALL use a student-t distribution for residual variance

# abs(totZn.fit2$Estimate[-2] - totZn.fit2.t$Estimate[-2]) / totZn.fit2.t$Estimate[-2]
#comparison between the predictors assuming data are from a Gaussian vs Student-t distribution 
#  0.0 to 0.4% difference in global intercept,  0.6% to 1.6% difference in not_greenBE, 4.9 to 4.8% difference in sqrt_CO2_transport (items 1, 4, 5)

#------
# End of 2
#------


# #compare the Leave One Out (loo) criterion for these three models; LOO is a cross-validation technique to validate the model
# fit1 <- add_criterion(fit1, criterion=c("loo"), moment_match=TRUE)
# fit2 <- add_criterion(fit2, criterion=c("loo"), moment_match=TRUE)
# fit1.t <- add_criterion(fit1, criterion=c("loo"), moment_match=TRUE)
# fit2.t <- add_criterion(fit2.t, criterion=c("loo"))

loo_compare(fit1, fit1.t, fit2, fit2.t, criterion="loo")  #top one in the output gives the best model

# # PSIS diagnostic tool - look for points above 0.5, and especially above 0.7
# par(mfrow=c(3,1), mar=c(2,4,4,2))
# plot(loo(fit1, cores=getOption("mc.cores", 1)))  #one point above 1.0
# plot(loo(fit2, cores=getOption("mc.cores", 1)))  #3 points above 0.5; 1 point above 0.7
# plot(loo(fit2.t, cores=getOption("mc.cores", 1)))  #1 point above 0.5; no points above 0.7
#using the student's t-distribution gives an improvement

#posterior predictive check of our two models
# pp_check(fit2, ndraws=100)
pp_check(fit2.t, ndraws=100)
# use the student t-distribution; it deals with the 1 high pareto K value.

#look at the relative amount of variability within each location, and relative amount of variability within each agency
par(mfrow=c(2,1))
boxplot(result ~ location, data=totZn.coc2)
boxplot(result ~ agency, data=totZn.coc2)
#There is variability at both the location AND the agency scale.  I think I'd prefer to use location, as agency sometimes has
#  only 1 location, sometimes 3

summary(fit2.t)$fixed

#look at residuals vs fitted values for candidate model - are any markedly better than others?
par(mfrow=c(2,1), mar=c(4,4,4,2))
resid.2t <- residuals(fit2.t, type="ordinary")
fitted.2t <- fitted(fit2.t, scale="response")
plot(resid.2t[,1] ~ fitted.2t[,1], ylab="residuals", xlab="fitted values", main="Model 2t - student-T distribution")

totZn.brm <- fit2.t

#-------------------------------------------
#https://tem11010.github.io/regression_brms/

#graphical posterior predictive checking. Compare observed data to simulated data from the posterior predictive distribution. 
#  This is a density plot, where the observed y values are plotted with expected values from the posterior distribution
pp_check(totZn.brm, ndraws=200)

#Look at the fit based on the grouping variable. Here are scatter-plots with the observed chemical concentrations (log scale) 
#  on the y-axis and the average model predictions (across all posterior samples) on the x-axis.
#  Red line is the 1:1 line, indicating perfect fit of model predictions to data.  Any locations where model doesn't fit?
pp_check(totZn.brm, type = "scatter_avg_grouped", group = "location") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
#-------------------------------------------


#------
# 3. No agency as a random effect
#------

fit2.t.RE_loc <- brm(bf(formula= update(totZn.Form4, .~. + (1|location)), 
                        sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
                     data=totZn.coc2,
                     prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                     cores = getOption("mc.cores", 4),
                     control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

totZn.fit2.t.RE_loc <- summary(fit2.t.RE_loc)$fixed

totZn.fit2.t.RE_loc
totZn.fit2.t

round(abs(totZn.fit2.t$Estimate[-2] - totZn.fit2.t.RE_loc$Estimate[-2]) / totZn.fit2.t$Estimate[-2] * 100, 1)
#comparison between the predictors assuming random effects is LOCATION only (not nested inside agency)
#  0.3% difference in global intercept parameter, 0.3% difference in not_greenBE, 0.8% difference in sqrt_CO2_transport

#------
# End of 3
#------


#------
# 4. RMSE
#------

#load(file="../results/Bayesian_TotalZinc.Rdata")

totZn_resid2 <- (residuals(totZn.brm)[, 1])^2  #1. take the residuals & square them
totZn_RMSE <- sqrt(mean(totZn_resid2))  #2 & 3. take the mean of the squared residuals, then take the sqrt of the result
totZn_RMSE
# 0.5531018

#RSR = RMSE/SD
totZn_RSR <- totZn_RMSE / sd(totZn.coc2$result)
totZn_RSR
# 0.4890489

#------
# End of 4
#------


save(totZn.brm, file="results/Bayesian_TotalZinc.Rdata")
#load(file="../results/Bayesian_TotalZinc_not_greenBE.Rdata")


#---------------------------#
#  Total Kjeldahl Nitrogen  #
#---------------------------#

load(file="results/Frequentist_Total Kjeldahl Nitrogen Models_censtat.RData")
TKN.Form4
TKN.r1X
TKN.vf1X
TKN.lme <- lme(TKN.Form4, data=TKN.coc2, method="REML", random = TKN.r1X, weights=TKN.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))


#below, we will use the brms censored models rather than running brms models on data generated by regression on order statistics (ROS).
#  From a Google search: "In brms, censored models treat censored observations as latent/missing values constrained to a known threshold 
#  or interval, evaluated directly via the likelihood's cumulative distribution function (CDF) or survival function. This is better than 
#  regression on order statistics (ROS) because ROS uses ad-hoc imputations or approximations for censored points, whereas likelihood-based 
#  censoring in Stan/brms provides exact, unbiassed marginalization of the censored range."
#  Furthermore, brms allows censoring combined with random effects.  ROS breaks down in hierarchical settings.

#------
# 1. Outliers
#------

#check for outliers
plot(TKN.coc2$result)   #hmmm.... no real outliers here, but lets see if removing the top point makes a difference
TKN.coc2$result[order(TKN.coc2$result)]
hist(TKN.coc2$result)

#identify the location of the highest point/ potential outlier:
which(TKN.coc2$result==max(TKN.coc2$result))  #its row 287...

#Total Kjeldahl Nitrogen lme model summary -- one outlier removed (row 265)
TKN.lme.outlierRemoved <- lme(TKN.Form4, data=TKN.coc2[-which(TKN.coc2$result==max(TKN.coc2$result)),], 
                              method="REML", random = TKN.r1X, weights=TKN.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(TKN.lme.outlierRemoved)

#Total Kjeldahl Nitrogen lme model summary -- all data points
TKN.lme <- lme(TKN.Form4, data=TKN.coc2, method="REML", random = TKN.r1X, weights=TKN.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(TKN.lme)

#fixed coefficients for the lme model with all data points, vs the one with the outlier removed
TKN.lme$coefficients$fixed
TKN.lme.outlierRemoved$coefficients$fixed
abs(TKN.lme$coefficients$fixed - TKN.lme.outlierRemoved$coefficients$fixed)/TKN.lme$coefficients$fixed  #ratio of difference btwn two models/original model
#fixed effects for the model with all data points differs from the one with the outlier removed by 2% or less.
#  the greatest difference is in sqrt_traffic fixed effect (2%) followed by devAge2 (0.5%) and global intercept (0.05%)



plot(TKN.coc2$result)
#two possible outliers?

TKN.lme_inf <- hlm_influence(TKN.lme, level=1)
dotplot_diag(TKN.lme_inf$cooksd, name = "cooks.distance", cutoff = "internal")
TKN.coc2[c(31, 390, 118, 60, 126),]  #top 5 observations in terms of Cook's distance

min(TKN.coc2$result)
max(TKN.coc2$result)

TKN.coc3 <- TKN.coc2[-c(31, 390),]  #remove the top two, which are particularly influential
TKN.lme.rmOutliers <- lme(TKN.Form4, data=TKN.coc3, method="REML", random = TKN.r1X, weights=TKN.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

coefficients(TKN.lme)
coefficients(TKN.lme.rmOutliers)

fixed.effects(TKN.lme)
fixed.effects(TKN.lme.rmOutliers)

plot(TKN.lme)
plot(TKN.lme.rmOutliers)

round(abs(fixed.effects(TKN.lme) - fixed.effects(TKN.lme.rmOutliers) ) / fixed.effects(TKN.lme) * 100, 1)

#two potential influential outliers based on Cook's distance.  
# NOTE: neither of the two potentially influential observations are the min or max!  One is close to the max, other is somewhere in the middle and 
#       was generated as a replacement for censored data.  I don't think there is compelling reason to eliminate them.

#------
# End of 1
#------


#look at the relative amount of variability within each location, and relative amount of variability within each agency
par(mfrow=c(2,1))
boxplot(result ~ location, data=TKN.coc2)
boxplot(result ~ agency, data=TKN.coc2)
#There is variability at both the location AND the agency scale.  I'd prefer to use location, as agency sometimes has
#  only 1 location, sometimes 3

#Try Bayesian censored model methods
#add a column to coc2, which indicates whether there is no censoring "none", or left-censoring "left" of data points
TKN.coc2 <- TKN.coc2 %>% 
  mutate(cen1 = if_else(cen==TRUE, "left", "none"))



TKN.Form4.cen <- as.formula(bquote(log(oconc) | cens(cen1) ~ .(TKN.Form4[[3]]) ))  #Form4, but starting with "log(oconc) | cens(cen1)"


#Bayesian Mixed Model using brms built-in censored data function; variance structure = varIdent(form = ~1|location)
fit2.cen <- brm(bf(update(TKN.Form4.cen, .~. + (1|agency/location)), 
                   sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
                data=TKN.coc2,
                prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                cores = getOption("mc.cores", 4),
                control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#Bayesian Mixed Model using brms built-in censored data function; variance structure = varIdent(form = ~1|agency)
fit1.cen <- brm(bf(update(TKN.Form4.cen, .~. + (1|agency/location)), 
                   sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
                data=TKN.coc2,
                prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                cores = getOption("mc.cores", 4),
                control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#Bayesian Mixed Model using brms built-in censored data function; variance structure = varIdent(form = ~1|agency)
fit0.cen <- brm(bf(update(TKN.Form4.cen, .~. + (1|agency/location)) ), 
                   #sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
                data=TKN.coc2,
                prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                cores = getOption("mc.cores", 4),
                control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#### ROS models should give (almost) the same results as the CEN models; just double check below.

#Bayesian Mixed Model using the ROS data that were also used for frequentist lme models; variance structure = varIdent(form = ~1|agency)
fit1.ros <- brm(bf(update(TKN.Form4, .~. + (1|agency/location)), 
                     sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
                  data=TKN.coc2,
                  prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                  cores = getOption("mc.cores", 4),
                  control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#Bayesian Mixed Model using the ROS data that were also used for frequentist lme models; variance structure = varIdent(form = ~1|location)
fit2.ros <- brm(bf(update(TKN.Form4, .~. + (1|agency/location)), 
                     sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
                  data=TKN.coc2,
                  prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                  cores = getOption("mc.cores", 4),
                  control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#check the pareto k values for the brms censored model
fit2.cen <- add_criterion(fit2.cen, criterion=c("loo"))#, moment_match=TRUE)
fit2.cen <- add_criterion(fit2.cen, criterion=c("loo"), moment_match=TRUE)
fit1.cen <- add_criterion(fit1.cen, criterion=c("loo"))#, moment_match=TRUE)
fit2.ros <- add_criterion(fit2.ros, criterion=c("loo"))#, moment_match=TRUE)
fit1.ros <- add_criterion(fit1.ros, criterion=c("loo"))#, moment_match=TRUE)
loo_compare(fit2.cen, fit2.ros, fit1.cen, fit1.ros, criterion="loo")
#WOAH!  cen vs ros give very different answers here.  Stick with cen, which integrates the censored data into the CDF.
#  Also note that the model with variance covariate of agency fits better than with var cov = location


par(mfrow=c(2,2), mar=c(4,4,4,2))
plot(loo(fit2.cen, cores=getOption("mc.cores", 1)), main="fit2.cen")  #one point above 0.7, 3 points above 1.0
plot(loo(fit1.cen, cores=getOption("mc.cores", 1)), main="fit1.cen")  #one point above 0.7, 3 points above 1.0
#several pareto k-values for the two .cen models are really high - try a student-t distribution model.
#  note that the ros models don't have any high pareto k-values.  I'd prefer to stick with the Bayesian censored
#  methods that are built-in, though.

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency) -- student t distribution
fit1.t.cen <- brm(bf(update(TKN.Form4.cen, .~. + (1|agency/location)), 
                     sigma ~ (1|agency)),  #equivalent to varIdent(form= ~1|agency)
                  data=TKN.coc2,
                  family=student,
                  prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                  cores = getOption("mc.cores", 4),
                  control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; variance structure = varIdent(form = ~1|location) -- student t distribution
fit2.t.cen <- brm(bf(update(TKN.Form4.cen, .~. + (1|agency/location)), 
                     sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
                  data=TKN.coc2,
                  family=student,
                  prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                  cores = getOption("mc.cores", 4),
                  control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

#agency/location as nested random effect; no variance structure
fit0.t.cen <- brm(bf(update(TKN.Form4.cen, .~. + (1|agency/location))), 
                  data=TKN.coc2,
                  family=student,
                  prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                  cores = getOption("mc.cores", 4),
                  control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

fit0.t.cen <- add_criterion(fit0.t.cen, criterion=c("loo"))
fit1.t.cen <- add_criterion(fit1.t.cen, criterion=c("loo"))
fit2.t.cen <- add_criterion(fit2.t.cen, criterion=c("loo"))
loo_compare(fit2.t.cen, fit1.t.cen, fit0.t.cen, criterion="loo")
#fit1.cen is best fit!

loo_compare(fit1.cen, fit1.t.cen, criterion="loo")

# PSIS diagnostic tool - look for points (potential influential outliers) above 0.5, and especially above 0.7
par(mfrow=c(3,1))
plot(loo(fit0.t.cen, cores=getOption("mc.cores", 1)), main="no var struct")
plot(loo(fit1.t.cen, cores=getOption("mc.cores", 1)), main="varcov=agency")
plot(loo(fit2.t.cen, cores=getOption("mc.cores", 1)), main="varcov=location")




#------
# 2. Changing Variance Structure
#------

TKN.fit0.t.cen <- summary(fit0.t.cen)$fixed
TKN.fit1.t.cen <- summary(fit1.t.cen)$fixed   #best model is fit1.t!
TKN.fit2.t.cen <- summary(fit2.t.cen)$fixed

round(abs(TKN.fit0.t.cen$Estimate - TKN.fit1.t.cen$Estimate[-2]) / TKN.fit1.t.cen$Estimate[-2] * 100, 1)
round(abs(TKN.fit2.t.cen$Estimate[-2] - TKN.fit1.t.cen$Estimate[-2]) / TKN.fit1.t.cen$Estimate[-2] * 100, 1)
#comparison between [no variance structure] or [location as variance covariate] and [agency as variance covariate] shows 
#  0.01 to 0.4% difference in global intercept, 4 to 8% difference in sqrt_traffic, & 0.6 to 1.5% difference in devAge2 (items 1, 4,5)

#------
# End of 2
#------


#this will be important in #3 below, since using 1/agency as a variance covariate wouldn't make sense if it wasn't a nested RE 
TKN.fit2.t.cen <- summary(fit2.t.cen)$fixed



#posterior predictive check of our two models
pp_check(fit1.t.cen, ndraws=100)
pp_check(fit1.cen, ndraws=100)
pp_check(fit1.ros, ndraws=100)
#while both models would be fine, choose the brms cen model (fit1.t.cen); this uses censored
#  methods within the Bayesian context, making it perhaps more suitable than the ROS method
#  that was used outside of the Bayesian context.

#look at residuals vs fitted values for candidate model - are any markedly better than others?
par(mfrow=c(2,1), mar=c(4,4,4,2))
resid.1t.cen <- residuals(fit1.t.cen, type="ordinary")
fitted.1t.cen <- fitted(fit1.t.cen, scale="response")
plot(resid.1t.cen[,1] ~ fitted.1t.cen[,1], ylab="residuals", xlab="fitted values", main="Model 1 - t-distr censored model")

TKN.brm <- fit1.t.cen   #this is the model we are choosing -- using censored methods within the Bayesian context

#Look at the fit based on the grouping variable. Here are scatter-plots with the observed chemical concentrations (log scale) 
#  on the y-axis and the average model predictions (across all posterior samples) on the x-axis.
#  Red line is the 1:1 line, indicating perfect fit of model predictions to data.  Any locations where model doesn't fit?
pp_check(TKN.brm, type = "scatter_avg_grouped", group = "location") + 
  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
#-------------------------------------------


#------
# 3. No agency as a random effect
#------

#agency/location as nested random effect; variance structure = varIdent(form = ~1|agency) -- student t distribution
fit2.t.cen.RE_loc <- brm(bf(update(TKN.Form4.cen, .~. + (1|location)), 
                     sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
                  data=TKN.coc2,
                  family=student,
                  prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                  cores = getOption("mc.cores", 4),
                  control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

TKN.fit2.t.cen.RE_loc <- summary(fit2.t.cen.RE_loc)$fixed

TKN.fit2.t.cen.RE_loc
TKN.fit2.t.cen

round(abs(TKN.fit2.t.cen$Estimate[-2] - TKN.fit2.t.cen.RE_loc$Estimate[-2]) / TKN.fit2.t.cen$Estimate[-2] *100, 1)
#comparison between the predictors assuming random effects is LOCATION only (not nested inside agency)

### NOTE: this comparison was run for the SECOND-best model, which uses location as a variance covariate.  
#         The model that uses agency as a variance covariate is not appropriate for comparison to a model
#         where agency was dropped as the random effect.


#------
# End of 3
#------


#------
# 4. RMSE
#------

#load(file="../results/Bayesian_TotalKjeldahlNitrogen.Rdata")

TKN_resid2 <- (residuals(TKN.brm)[, 1])^2  #1. take the residuals & square them
TKN_RMSE <- sqrt(mean(TKN_resid2))  #2 & 3. take the mean of the squared residuals, then take the sqrt of the resulting mean
TKN_RMSE
# 0.6239021

#RSR = RMSE/SD
TKN_RSR <- TKN_RMSE / sd(TKN.coc2$result)
TKN_RSR
# 0.6821426


#------
# End of 4
#------

save(TKN.brm, file="results/Bayesian_TotalKjeldahlNitrogen.Rdata")


# #---------------------------------------------------#
# #  Save Model Results Used in Sensitivity Analysis  #   #####  oops!  I must have cleared the memory after running half the script.
# #---------------------------------------------------#
# 
# save(Cu.fit0, Cu.fit1, Cu.fit2, Cu.fit0.t, Cu.fit1.t, Cu.fit2.t, Cu.fit2.t.RE_loc,
#      TSS.fit0, TSS.fit1, TSS.fit2, TSS.fit2.RE_loc,
#      P.fit0, P.fit1, P.fit2, P.fit0.t, P.fit1.t, P.fit2.t, P.fit2.t.RE_loc,
#      totZn.fit0, totZn.fit1, totZn.fit2, totZn.fit0.t, totZn.fit1.t, totZn.fit2.t, totZn.fit2.t.RE_loc,
#      TKN.fit0.cen, TKN.fit1.cen, TKN.fit2.cen, TKN.fit0.t.cen, TKN.fit1.t.cen, TKN.fit2.t.cen, TKN.fit1.t.cen.RE_loc, TKN.fit2.t.cen.RE_loc,
#      file="../results/Bayesian_Model_Results_Sensitivity_Analysis.Rdata")
# 
# #load(file="../results/Bayesian_Model_Results_Sensitivity_Analysis.Rdata")

