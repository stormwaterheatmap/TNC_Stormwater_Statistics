
#run Frequentist Copper Model_v9.R first...
myForm <- my.formulas.m4[[5]]

#these lines of code assess fit of this particular model in terms of COC vs. individual predictors, and predictor correlation
myModel.traffic <- lme(data=coc2, my.formulas.m4[[5]], method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(myModel.traffic)

myModel.sqrt_traffic <- lme(data=coc2, my.formulas.m4[[16]], method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(myModel.sqrt_traffic)

library(brms)
library(loo)

Cu.coc2 <- coc2
#fit2.t.sqrt_traffic <- fit2.t.traffic
fit2.t.sqrt_traffic <- brm(bf(result ~ summer + rain + sqrt_traffic + devAge2 + (1|agency/location), 
                 sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
              data=Cu.coc2,
              family=student,
              prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
              cores = getOption("mc.cores", 1),
              control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)


fit2.t.traffic <- brm(bf(result ~ summer + rain + traffic + devAge2 + (1|agency/location), 
                         sigma ~ (1|location)),  #equivalent to varIdent(form= ~1|location)
                      data=Cu.coc2,
                      family=student,
                      prior = c(set_prior("normal(0,10)", class="b")),   #non-informative priors on all predictors
                      cores = getOption("mc.cores", 1),
                      control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 20)
)

source("Bayesian Model Plot Functions.R")  #script with plotting functions for Bayesian figures

#read in standardized values; this will be used in plots of prediction intervals vs raw data
std.vals <- read.csv(file="../processed_data/spatial_predictor_standardization_values.csv")
colnames(std.vals) <- c("predictor", "mean", "sd")

#read in backtransform powers for spatial predictors, used for plotting raw spatial predictor values
xform.pwr <- read.csv("../processed_data/spatial_predictor_backtransform_power.csv", header=TRUE)


predInt.vs.rawPreds.2preds(Cu.coc2, fit2.t.sqrt_traffic, "Total Copper vs sqrt_traffic", preds=c("sqrt_traffic", "devAge2"))

predInt.vs.rawPreds.2preds(Cu.coc2, fit2.t.traffic, "Total Copper vs traffic", preds=c("traffic", "devAge2"))
