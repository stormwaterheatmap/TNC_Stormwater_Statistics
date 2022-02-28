# This generalized script assesses linear relationships between total phosphorous and landscape predictors
# using mixed effects models

# Note: This version takes the updated set of predictors (from July 7, 2021) and examines each set
#       with the COC of interest, to see which set of predictors may be most correlated with
#       the COC.  Based on linear fits, a select set of predictors (best_predictors) is tested in
#       various combinations within LME models

#Four models are generated for the COC of interest (note the model numbers are 1,3,4,5 (no 2)):
#   1.) median value for all locations, irrespective of date
#   3.) COC value based only on land use (categorical) + rain, precip, summer (as applicable)
#   4.) COC value based on up to 3 landscape parameters NOT including land use + rain


# Eva Dusek Jennings
# Feb 18, 2022  -- with sqrt_traffic
#----------------------------------------------

rm(list = ls(all = T))

# Load libraries ----------------------------------------------
library(effects)
library(EnvStats)  #note the objects that are masked from other packages...
library(tidyverse)
library(ggplot2)
library(graphics)
library(stats)
library(nlme)
library(grid)
library(gridExtra)
library(lme4)  #for lmer models
library(lattice)
library(gstat)  #for spatial correlation exploration
library(sp)  #for spatial correlation exploration
library(sjPlot)
library(sjlabelled)
library(TMB)
library(glmmTMB)
library(here)
library(lattice)
library(sp)
library(NADA)
library(ggplot2)
library(knitr)
library(texreg)
library(huxtable)
#run script with functions specific to COC analysis
source(here("..", "functions", "COC analysis functions.R"))
source(here("..", "functions", "HighstatLibV10.R")) #Zuur library, incl panel.smooth2, VIF

#toggle for whether all exploratory parts of code should be run (TRUE) or just essential parts (FALSE) 
run_exploratory_code <- FALSE

#read in stormwater data with spatial predictors-------------------------------
s8 <- read.csv(here("..", "processed_data", "s8data_with_spatial_predictors.csv"))
s8$start_date <- as.Date(s8$start_date)
s8$conc <- s8$result
s8$month <- as.factor(s8$month)
s8$location_id <- as.factor(s8$location_id)
s8$loc <- as.factor(s8$loc)
s8$agency <- as.factor(s8$agency)
s8$season <- as.factor(s8$season)
s8$land_use <- as.factor(s8$land_use)

#list of landscape predictors-------------------------------
aa <- which(str_detect(names(s8), "daymet"))  #after the last predictor column is a series of "daymet" columns
bb <- aa[which(aa > 34)][1] - 1  #identify the last predictor column
predictors <- names(s8)[34:bb]  #first predictor is always column 34; get all predictors
n.preds <- length(predictors)

#select the chemical of interest
this_param <- 'Total Phosphorus - Water - Total'
this_param_short <- "Phosphorus"
coc <- s8[which(s8$parameter==this_param),]  #create the "coc" dataframe for this chemical of interest only

#Non Detect Handing -------------------------------

#look at censored points
ggplot(coc)+geom_jitter(aes(x=agency,y=conc,color=nondetect_flag))+scale_y_log10()

#which values are ND?  min/max detection limit?  where and when were these samples recorded?
coc.nd <- coc[which(coc$nondetect_flag==TRUE),]
pct.cen <- 100*nrow(coc.nd)/nrow(coc)  #for phosphorus, 8 samples are ND (1.86%); 3 of these (TAC) were "estimated"
nrow(coc.nd)
min(coc.nd$conc)
max(coc.nd$conc)
coc.nd[,c("location_id", "start_date", "conc")]  #location, date, and detection limit of ND samples

#NOTE: Two TAC samples (TAC003S8D 2010-02-15 and TFWFD1 2011-11-16) had results of 118.0 and 107.0, but were classified
#      as non-detect.  Additionally, a sample at TFWFD1 from 2010-02-15 was 16.6 ug/L,but classified as ND.  Christian
#      looked into these samples, and saw that they were actually flagged as "estimated".  For this study, assume that 
#      these are the correct values, and remove the ND flag.  All other ND's were at POT (4 total) and SNO_COM (1 total)
coc[which(coc$location_id=="TFWFD1" & coc$nondetect_flag==TRUE),]$nondetect_flag <- FALSE
coc[which(coc$location_id=="TAC003S8D_OF245" & coc$nondetect_flag==TRUE),]$nondetect_flag <- FALSE

#Since there are only five true ND's for Phosphorus, substitute the ND's with half the detection limit
coc <- coc %>% 
  mutate(conc = ifelse(nondetect_flag, conc * 0.5, conc))

coc.nd <- coc[which(coc$nondetect_flag==TRUE),]
pct.cen <- 100*nrow(coc.nd)/nrow(coc)  #following sample adjustment, 5 samples are ND (1.17%)


# distribution exploration ----------------------------------------------
if(run_exploratory_code==TRUE) {
  #explore the data; look for underlying distribution
  par(mfrow=c(2,2))
  hist(coc$conc, breaks=50)
  qqnorm((coc$conc), main=paste("QQ-Normal plot of", this_param_short))
  qqline((coc$conc))
  hist(log(coc$conc), breaks=50)
  qqnorm(log(coc$conc), main=paste("QQ-Normal plot of log(", this_param_short, ")", sep=""))
  qqline(log(coc$conc))
  
  #consider other distribution (log-normal and gamma).  Compare QQ plots for these distributions
  par(mfrow=c(2,2))
  qqPlot(log(coc$conc), dist="norm", estimate.params=TRUE, add.line=TRUE, main=paste("QQ log-normal plot of", this_param_short))
  qqPlot(log(coc$conc), dist="norm", param.list=list(mean=-3.3, sd=2.09), add.line=TRUE, main=paste("QQ log-normal plot of", this_param_short))
  qqPlot(coc$conc, dist="gamma", param.list=list(shape=.7, scale=1), points.col="blue",
         add.line=TRUE, main=paste("QQ-gamma plot of", this_param_short))
  qqPlot(sqrt(coc$conc), dist="norm", estimate.params=TRUE, add.line=TRUE, main=paste("QQ plot of sqrt-transformed", this_param_short))
  #qqPlot(1/(coc$conc), dist="norm", estimate.params=TRUE, add.line=TRUE)
}

# gamma distribution & log-normal are both decent.  Use log-normal.  Ln-transform concentration data in column "result"
coc$result <- log(coc$conc)


#-------------------------------------------------------------------------------------------#
#  Explore Cleveland Dotplots and Boxplots of COC data or predictors conditional on Agency  #
#-------------------------------------------------------------------------------------------#

colors_agency <- c("red", "orange", "yellow", "green", "blue", "purple")
colors_location <- colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)]

if(run_exploratory_code==TRUE) {
  #Cleveland Dotplot - use to detect violations of homogeneity: We are looking for whether the spread of data values
  #   differs between sampling locations (or agencies).  If so, it indicates heterogeneity, and that there may be problems with violation of
  #   homogeneity in a linear regression model applied on these data.  We're also looking for outliers.
  par(mfrow=c(1,1))
  dotchart(coc$result, groups=coc$location_id, pch=19, col=colors_agency[as.numeric(coc$agency)],
           xlab="concentration", main=paste("Cleveland Dotplot:", this_param_short))
  
  #pairs plots allow visualization of interactions between possible predictors.  Look for relationships between "result"
  #   and all of the predictor variables
  pairs(coc %>% select(result, predictors[1:ceiling(n.preds/3)]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  pairs(coc %>% select(result, predictors[ (ceiling(n.preds/3)+1) : (ceiling(n.preds/3)*2) ]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  pairs(coc %>% select(result, predictors[ (ceiling(n.preds/3)*2+1) : n.preds]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  
  
  #look also at interactions between result and rainfall
  pairs(coc %>% select(result, antecedant_dry_days_std, daymet_precip_std, daymet_3day_std,
                       daymet_7day_std, daymet_14day_std, daymet_21day_std, daymet_28daySR_std), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  #for phosphorus, important predictors are likely: 
  #  traffic, sqrt_CO2_road
  #  daymet 28day (-0.2), daymet 21day (-0.2), daymet 14day (-0.2)
}


#-----------------------------------------------------------------------------------------------------#
#  Plot COC vs various predictors: look for potential predictors to model & sources of heterogeneity  #
#-----------------------------------------------------------------------------------------------------#

lp_plots()
pr_plots()
tlp_plots()
mp_plots()


#best predictors for this COC; make sure to only have one version (transformed or not) of each predictor!
best_predictors <- c("urbRES", "roofs", "sqrt_traffic", "sqrt_popn", "pm25_na", "sqrt_slope", "sqrt_CO2_res", "sqrt_CO2_tot", 
                     "sqrt_CO2_road", "devAge2", "grass", "paved")
#NOTE: grass & paved are only added b/c they were part of the best fit model for old version of predictors!

pred_i <- which(predictors %in% best_predictors)
# grid.arrange(
#   for(i in pred_i) {
#     get(paste("lp", i, sep=""))
#   },
#   nrow=4, ncol=4)

lp6 <- my.ggplot(6)
lp8 <- my.ggplot(8)
lp13 <- my.ggplot(13)
lp14 <- my.ggplot(14)
lp18 <- my.ggplot(18)
lp19 <- my.ggplot(19)
lp20 <- my.ggplot(20)
lp22 <- my.ggplot(22)
lp24 <- my.ggplot(24)
grid.arrange(lp6, lp8, lp13, lp14, lp18, lp19, lp20, lp22, lp24, nrow=4, ncol=4)

bp_coefs <- bp_signs <- rep(NA, length(best_predictors))
names(bp_coefs) <- names(bp_signs) <- best_predictors
for(i in 1:length(best_predictors)) {
  aa <- lm(as.formula(paste("result ~ ", best_predictors[i])), data=coc)
  bp_coefs[i] <- coefficients(aa)[2]
}
bp_signs <- sign(bp_coefs)

#examine correlations between predictors; generate a second vector of only predictors that aren't highly correlated
#  this second vector will be used to generate FormX
pairs(coc[best_predictors], lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)

elements_to_remove <- c("sqrt_popn", "sqrt_CO2_res", "sqrt_CO2_road")  #these predictors are highly correlated with others
best_predictors2 <- best_predictors[!(best_predictors %in% elements_to_remove)]


#--------------------------------------------------#
#  Any evidence of changes in COC conc over time?  #
#--------------------------------------------------#

if(run_exploratory_code==TRUE) {
  #for each agency, is there a trend in COC concentration over time?
  xyplot(result~start_date|agency, data=coc,
         xlab="time", ylab="log(conc)",
         strip=function(bg="white",...)
         strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
  #   For phosphorus, no evidence of patterns in concentration changing over time
}


#---------------------------------------------------#
#  How does rainfall change COC conc at each site?  #
#---------------------------------------------------#

if(run_exploratory_code==TRUE) {
  
  #daymet precip
  xyplot(result~daymet_precip_std|location_id, data=coc,
         xlab="daymet precip", ylab="log(conc)",
         strip=function(bg="white",...)
           strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
  
  #3-day daymet precip
  xyplot(result~daymet_3day_std|location_id, data=coc,
         xlab="daymet 3-day precip", ylab="log(conc)",
         strip=function(bg="white",...)
           strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
  
  #7-day daymet precip
  xyplot(result~daymet_7day_std|location_id, data=coc,
         xlab="daymet 7-day precip", ylab="log(conc)",
         strip=function(bg="white",...)
           strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
  
  #14-day daymet precip
  xyplot(result~daymet_14day_std|location_id, data=coc,
         xlab="daymet 14-day precip", ylab="log(conc)",
         strip=function(bg="white",...)
           strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
  
  #21-day daymet precip
  xyplot(result~daymet_21day_std|location_id, data=coc,
         xlab="daymet 21-day precip", ylab="log(conc)",
         strip=function(bg="white",...)
           strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
  
  #28-day daymet precip
  xyplot(result~daymet_28daySR_std|location_id, data=coc,
         xlab="daymet 28-day sqrt xform precip", ylab="log(conc)",
         strip=function(bg="white",...)
           strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
  
  #antecedant dry days
  xyplot(result~antecedant_dry_days_std|location_id, data=coc,
         xlab="antecedant dry days", ylab="log(conc)",
         strip=function(bg="white",...)
           strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
#for total phosphorus, start to see a decrease in total phosphorus by location
#   with 7-day daymet precip.  Signal is really strong for 21- and 28-day precip
}


#--------------------------------------------#
#  Create coc2: only predictors we will use  #
#--------------------------------------------#

#see lines 126-128 for selection of predictors that are important for this COC!
coc2 <- coc %>%
  select(result,
         location=loc, 
         latitude, longitude,
         agency,
         land_use,
         year, 
         month, season, start_date,
         all_of(best_predictors),
         dry=antecedant_dry_days_std,
         rain=daymet_21day_std,  #choose whichever precip measure makes the most sense for this coc!
         mPrecip=mPrecip  #note: if one of the other mPrecip measures looked better, use it here instead!
  ) %>%
  mutate(year=as.factor(year))  %>%
  mutate(landuse = case_when(
    land_use == "LDR" ~ "1.LDR",
    land_use == "HDR" ~ "2.HDR",
    land_use == "COM" ~ "3.COM",
    land_use == "IND" ~ "4.IND")) %>%
  mutate(landuse = as.factor(landuse)) %>%
  mutate(summer = case_when(
    season=="1" ~ "0",
    season=="2" ~ "0",
    season=="3" ~ "1",
    season=="4" ~ "0")) %>%
  mutate(summer=as.factor(summer))

monthlyAvgPrecip <- data.frame(month=c(1:12), 
                               precip=c(5.24, 4.09, 3.92, 2.75, 2.03, 1.55, 0.93, 1.16, 1.61, 3.24, 5.67, 6.06))
coc2$monthlyPrecip <- monthlyAvgPrecip[coc2$month,2]

rain <- "21-day precip, standardized"  #precip measure that was used for this COC


#------------------------------------------------#
#  Formulas for Possible Predictor Combinations  #
#------------------------------------------------#

weather <- c("rain", "summer")

FormX <- as.formula( paste("result ~ landuse + ", paste(weather, collapse=" + "), " + ", paste(best_predictors2, collapse=" + ")) )
#NOTE: it is not practical at this point to add interactions -- that will come once we find best model using just main effects
#      Also, don't include pairs of highly correlated predictors here (e.g.: nodev:devAge, paved:impervious, impervious:grass)


#-------------------------------------------------------#
#  Zuur et al Method for Mixed Effects Model Selection  #
#-------------------------------------------------------#

#------------------#
#  All Predictors  #
#------------------#

# Steps 1-2: fit LM for ALL PREDICTORS to data using gls; look for heterogeneity
#-----------

if(run_exploratory_code==TRUE) {
  M.gls1X <- gls(FormX, data=coc2, method="REML")  #use REML for comparing diff't random variances
  E.1X <- resid(object=M.gls1X, type="normalized")  
  boxplots.resids2(M.gls1X, E.1X, "X")  #check for heterogeneity in residuals vs fitted values; look for source in other plots
  
  #is there heteroscedasticity in any of these potential covariates?
  plot(coc2$result ~ coc2$location)
  plot(coc2$result ~ coc2$agency)
  plot(coc2$result ~ coc2$month)  #unlikely to be considered for use as a covariate
  plot(coc2$result ~ coc2$season)
  plot(coc2$result ~ coc2$year)
  plot(coc2$result ~ coc2$dry)
  abline(lm(coc2$result ~ coc2$dry))
  plot(coc2$result ~ coc2$rain)
  abline(lm(coc2$result ~ coc2$rain))
  plot(coc2$result ~ coc2$monthlyPrecip)
  plot(coc2$result ~ coc2$mPrecip)
  #
  
  
  #Step 3:  Choose a variance structure (if there was heterogeneity)
  #-------  for selecting random structure, use REML to compare - its better at capturing random structure
  #         REML estimates the random effects by considering linear combinations of the data that remove the 
  #         fixed effects. If these fixed effects are changed, the likelihoods of the two models will not be directly comparable.
  #         Compare AIC for the various beyond-optimal models with different variance structures.
  #         Alternately, can use likelihood ratio test (anova) to compare nested models
  
  
  M.gls2X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency))
  # M.gls3X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|location))
  M.gls4X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|season))
  # M.gls5X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency * season))
  # M.gls6X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency * year))
  M.gls7X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varExp(form= ~dry)))
  M.gls8X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varConstPower(form= ~dry)))
  # M.gls14X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varExp(form= ~monthlyPrecip)))  #agency & monthly average precip (same for all years & locations)
  # M.gls15X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varConstPower(form= ~monthlyPrecip)))
  # M.gls16X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varExp(form= ~mPrecip)))  #agency & monthly precip (diff't by year & location)
  # M.gls17X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varConstPower(form= ~(mPrecip/10))))
  M.gls18X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varExp(form= ~rain)))
  M.gls19X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varConstPower(form= ~rain)))
  M.gls20X <- gls(FormX, data=coc2, method="REML", weights=varConstPower(form= ~rain))
  M.gls21X <- gls(FormX, data=coc2, method="REML", weights=varExp(form= ~rain))
  M.gls22X <- gls(FormX, data=coc2, method="REML", weights=varFixed(~rain))
  
  anova(M.gls1X, M.gls2X,
        M.gls4X, M.gls7X, M.gls8X,
        M.gls18X, M.gls19X, M.gls20X, M.gls21X, M.gls22X)
  #for phosphorus: AIC is best for agency & dry (7, AIC=1048.6), followed by agency (2).  If the predictor "dry" doesn't
  #     fit well into the model, use only agency as the variance covariate.

  #likelihood ratio tests for nested versions of the best fit model
  anova(M.gls1X, M.gls2X, M.gls7X)  #note: for phosphorus, the fit of model 2 is not much worse than model 7 (delta AIC = 5.2)
  
  #residual plots for best fit models -- look for homogeneity of residuals
  E.7X <- resid(object=M.gls7X, type="normalized")  
  plot.resids(M.gls7X, E.7X, "X")
}

#variance functions for best fit model so far
#vf1X <- varComb(varIdent(form= ~1|agency), varExp(form= ~dry))
vf1X <- varIdent(form=~1|agency)

# Steps 4-6: Find the proper random effects structure; look for temporal and spatial autocorrelation
#-----------  note: look for spatial correlation AFTER setting random effects!

if(run_exploratory_code==TRUE) {
  M.vf1X <- gls(FormX, data=coc2, method="REML", weights=vf1X)
  E.vf1X <- resid(object=M.vf1X, type="normalized")  
  
  #plots of residuals: no variance structure vs. variance structure
  par(mfrow=c(2,2), mar=c(4,4,4,1))
  plot(fitted(M.gls1X), E.1X, main="no variance structure", xlab="fitted values", xaxt="s", ylab="normalized residuals", col="gray", pch=16)
  plot(fitted(M.vf1X), E.vf1X, main="variance covariate: agency", xlab="fitted values", xaxt="s", ylab="normalized residuals", col="gray", pch=16)
  
  
    
  #Random intercept model; this is nested in the best fit model from Step 2, so can compare with a likelhiood ratio test
  M1.lme1X <- lme(data=coc2, FormX, random = ~1|agency, method="REML", weights=vf1X,
                  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  anova(M.vf1X, M1.lme1X)  #the likelihood ratio test p-value compares the two nested models, and whether they are significantly different
  #For phosphorus, the model with random intercept is significantly better than without!
  ### NOTE: with all best_predictors in place, the best fit model doesn't have a random intercept.
  #         however, once only 3 or less predictors are used, the lme model fits much better than gls
}

#random structure for best fit model so far
r1X <- formula(~1|agency)   ##### NOTE: try also GLS models with no random intercept -- they may fit better!!

#### see note on line 668 re: gls vs lme model with random intercept.  the random intercept model ends up
#       fitting a LOT better than the gls model when number of predictors is limited to 3

if(run_exploratory_code==TRUE) {
  #best random effects structure & residual plots to test it
  M1.rX <- lme(data=coc2, FormX, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  summary(M1.rX)
  E1.rX <- residuals(object=M1.rX, type="normalized")
  boxplots.resids2(M1.rX, E1.rX, "X")
  
  #Check for temporal correlation; look for patterns in Auto-correlation plot for residuals
  par(mfrow=c(2,2))
  plot(E1.rX~coc2$start_date, pch=16, col="cadet blue")
  acf(E1.rX, na.action=na.pass, main="Auto-correlation plot for residuals")  #look for lines extending past blue dashed range
  #for phosphorus, no indication of temporal autocorrelation
  
  #Check for spatial correlation
  mydata2 <- data.frame(E1.rX, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
  coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
  bubble(mydata2,"E1.rX",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
  
  #variogram to test for spatial correlation
  Vario1 = variogram(E1.rX ~ 1, mydata2)
  plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
  Vario2 <- variogram(E1.rX ~ 1, mydata2, alpha = c(0, 45, 90,135) )
  plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?
  #for phosphorus, variograms look good!
}

# Steps 7-8: Find the proper fixed effects structure; use ML for model comparisons
#-----------  

#-------------------------------------------------------------#
#  Exhaustive Search for Model with Best 3 Predictors + Rain  #
#-------------------------------------------------------------#

if(run_exploratory_code==TRUE) {
  #Construct all possible models with 1, 2 and 3 landscape predictors, both including and not including rain;
  #   use AIC to determine which models provide the best fit
  preds.1 <- best_predictors
  
  #one predictor
  dd.1 <- data.frame(Var1=preds.1, Var2=NA, Var3=NA)
  
  #two predictors
  dd.2 <- data.frame(NULL)
  for (i in 1:length(preds.1)) {
    preds.2 <- preds.1[-c(1:i)]
    dd.2 <- rbind(dd.2, expand.grid(preds.1[i], preds.2))
  }  
  dd.2$Var3 <- NA
  
  #three predictors
  dd.3 <- data.frame(NULL)
  for (j in 1:(length(preds.1)-2)) {
    aa <- preds.1[j]
    preds.1a <- preds.1[-c(1:j)]  #remove current predictor "j", and all predictors before it
    for (i in 1:(length(preds.1a)-1)) {
      preds.2a <- preds.1a[-c(1:i)]
      bb <- expand.grid(preds.1a[i], preds.2a)
      dd.3 <- rbind(dd.3, data.frame(bb, Var3=aa))
    }
  }
  
  #combine one, two and three-predictor options into one dataframe
  dd <- rbind(dd.1, dd.2, dd.3)  #, dd.1, dd.2, dd.3)
  
  dd.str <- paste(dd$Var1, dd$Var2, dd$Var3)  #vector with strings of the predictor sets
  
  #identify and eliminate situations with predictors that don't go together, or are highly correlated
  rem1 <- which(str_detect(dd.str, "sqrt_CO2_tot")==TRUE & str_detect(dd.str, "sqrt_CO2_road")==TRUE)
  rem2 <- which(str_detect(dd.str, "sqrt_popn")==TRUE & str_detect(dd.str, "sqrt_CO2_res")==TRUE)
  rem.all <- unique(c(rem1, rem2))  #keep only non-repetitive instances of indexes that should be removed
  
  dd <- dd[-rem.all, ]  #remove rows with predictors that don't go together!
  
  #generate an exhaustive set of formulas with one, two and three predictors, including and excluding rain
  ee.1 <- list(NULL)  #list of formulas based on the sets of predictors in "dd"; including rain
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i])) ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1]))
    else if(is.na(dd$Var3[i]))  ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1], "+", dd[i,2]))
    else  ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  # ee.2 <- list(NULL) #list of formulas based on the sets of predictors in "dd"; NOT including rain
  # for(i in 1:nrow(dd)) {
  #   if(is.na(dd$Var2[i]))  ee.2[[i]] <- as.formula(paste("result~ summer + ", dd[i,1]))
  #   else if(is.na(dd$Var3[i]))  ee.2[[i]] <- as.formula(paste("result~ summer + ", dd[i,1], "+", dd[i,2]))
  #   else ee.2[[i]] <- as.formula(paste("result~ summer + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  # }
  # ee.3 <- list(NULL) #list of formulas based on the sets of predictors in "dd"; NOT including rain
  # for(i in 1:nrow(dd)) {
  #   if(is.na(dd$Var2[i]))  ee.3[[i]] <- as.formula(paste("result~dry + ", dd[i,1]))
  #   else if(is.na(dd$Var3[i]))  ee.3[[i]] <- as.formula(paste("result~dry + ", dd[i,1], "+", dd[i,2]))
  #   else ee.3[[i]] <- as.formula(paste("result~dry + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  # }
  # ee.4 <- list(NULL) #list of formulas based on the sets of predictors in "dd"; NOT including rain
  # for(i in 1:nrow(dd)) {
  #   if(is.na(dd$Var2[i]))  ee.4[[i]] <- as.formula(paste("result~dry + rain + ", dd[i,1]))
  #   else if(is.na(dd$Var3[i]))  ee.4[[i]] <- as.formula(paste("result~dry + rain + ", dd[i,1], "+", dd[i,2]))
  #   else ee.4[[i]] <- as.formula(paste("result~dry + rain + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  # }
  
  ee <- ee.1  #c(ee.1, ee.2, ee.3, ee.4) 
  
  ###NOTE: Phosphorus is tricky - the linear relationship between individual predictors and phosphorus
  #        is fairly flat for many predictors, so comparing the slight positive or negative relationship
  #        with what the model predicts could be incorrect.
  
#  MyEqns <- list(formulas=ee, aics=rep(0, length(ee)), do_signs_match=rep(FALSE, length(ee)))
  
  #run the lme for each formula, saving the AIC values; takes about 7 minutes to run
  ptm <- proc.time()  #time the code!
  my.aics <- rep(0, length(ee))
  for (i in 1:length(ee)) {
    bb <- lme(data=coc2, ee[[i]], random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
    my.aics[i] <- AIC(bb)
    # MyEqns$aics[i] <- AIC(bb)
    # #extract landscape predictor coefficients 
    # aa <- summary(bb)$coefficients$fixed
    # aa <- aa[names(aa) %in% preds.1]
    # aa <- sign(aa)
    # aa <- aa[sort(names(aa))]
    # #extract coefficient signs for linear models of result~landscape predictors
    # cc <- bp_signs[names(bp_coefs) %in% names(aa)]
    # cc <- cc[sort(names(cc))]
    # #compare coefficient signs from this model, to those of single-predictor linear models
    # MyEqns$do_signs_match[i] <- identical(aa, cc)  #if coefficient signs match, TRUE; (default=FALSE)
    
  }
  proc.time() - ptm  #stop the clock;  9 minutes for 386 models
  
  # xx <- which(MyEqns$do_signs_match==TRUE)
  # 
  # ee.all <- ee
  # ee <- ee[xx]
  # my.aics.all <- my.aics
  # my.aics <- my.aics[xx]
  
  landuseAIC <- AIC(lme(data=coc2, as.formula(paste("result~landuse+", paste(weather, collapse="+"))), random = r1X, 
                        method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8)))
  
  par(mfrow=c(4,4), mar=c(1,2,2,0), oma=c(1,1,1,1))
  # gg <- rep("black", length(ee))
  # gg[c(which(is.na(dd$Var3)), which(is.na(dd$Var3))+nrow(dd), which(is.na(dd$Var3))+2*nrow(dd), which(is.na(dd$Var3))+3*nrow(dd))] <- "turquoise"   #make formulas with only 2 landscape predictors red
  # gg[c(which(is.na(dd$Var2)), which(is.na(dd$Var2))+nrow(dd), which(is.na(dd$Var2))+2*nrow(dd), which(is.na(dd$Var2))+3*nrow(dd) )] <- "pink"  #make formulas with only 1 landscape predictor yellow
  # plot(my.aics, col=gg, pch=16, main="1 pred=pink, 2 preds=turquoise", xaxt="n", ylab="AIC")
  gg <- rep("black", length(ee))
  gg[which(str_detect(as.character(ee), "rain"))] <- "light blue"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(ee), "dry"))] <- "goldenrod"   #make formulas with "dry" goldenrod
  gg[which(str_detect(as.character(ee), "dry") & str_detect(as.character(ee), "rain"))] <- "green"   #make formulas with "dry" goldenrod
  plot(my.aics, col=gg, xaxt="n", pch=16, main="rain=blue, dry=gold, both=green")
  abline(h=landuseAIC, col="gray")
  #  plot(my.aics, col=ifelse(str_detect(as.character(ee), "rain"), "light blue", "black"), xaxt="n", yaxt="n", pch=16, main="rain=light blue")
  
  colAIC <- c("red", "orange", "yellow", "light green", "green3", "cadet blue", "purple1", "tan3", "salmon", "maroon1", "goldenrod", "dark green")
  for (i in 1:length(best_predictors)) {
    plot(my.aics, col=ifelse(str_detect(as.character(ee), best_predictors[i]), colAIC[i], "black"), xaxt="n", yaxt="n", pch=16, main=paste(best_predictors[i], "=", colAIC[i]) )
    abline(h=landuseAIC, col="gray")
  }
  
  plot.new()
  text(x=0.5, y=0.5, labels=paste("Landuse AIC", round(landuseAIC, 1)), cex=2)
  
  #arrange formulae and AICs according to order of AIC value, from smallest to largest
  my.formulas.m4 <- ee[order(my.aics)]
  my.formulas.m4[1:10]
  my.aics.m4 <- my.aics[order(my.aics)]
  my.aics.m4[1:10]  #what are the top 10 AIC's?
  
  #save gls formulas -- no random intercept -- check if these have the best fit once only a few predictors are added!
#  save(my.formulas.m4, my.aics.m4, file=(here("scripts", "formulas_and_aics", "phosphorus_m4_gls.RData")))
  save(my.formulas.m4, my.aics.m4, file=(here("scripts", "formulas_and_aics", "phosphorus_m4_lme.RData")))
  
  
  #obtain formulas and AIC's where "rain" is the sole precip predictor, and only 1 or 2 predictors are used!
  rain.aics <- my.aics[1:72]
  rain.formulas <- ee[1:72]
  rain.formulas.ord <- rain.formulas[order(rain.aics)]
  rain.aics.ord <- rain.aics[order(rain.aics)]
  rain.formulas.ord[1:6]
  rain.aics.ord[1:6]
  
}

#NOTE: for phosphorus, adding "dry" improves fit slightly; check later to see if we should keep "dry" or remove (from variance fxn too)

#NOTE: The early test for whether to use a random intercept model or NO random intercept model indicated
#      best fit for NO random intercept.  HOWEVER, upon running the section above using both gls and lme,
#      it is clear that:  1.) best fit model is the same for both gls and lme
#                         2.) lme model fits MUCH better than gls model
#      stick with the lme model

if (run_exploratory_code==FALSE) {
  load(file=here("scripts", "formulas_and_aics", "phosphorus_m4_lme.RData"))
}




#--------------------------------------------#
#  Quantiles List for Plotting Interactions  #
#--------------------------------------------#

# #list of quantiles for all landscape predictors
# qList <- list(paved=quantile(coc2$paved, probs=c(0,.25,.50,.75,1)),
#               impervious=quantile(coc2$impervious, probs=c(0,.25,.50,.75,1)),
#               roofs=quantile(coc2$roofs, probs=c(0,.25,.50,.75,1)),
#               roof_RES=quantile(coc2$roof_RES, probs=c(0,.25,.50,.75,1)),
#               roof_nonRES=quantile(coc2$roof_nonRES, probs=c(0,.25,.50,.75,1)),
#               grass=quantile(coc2$grass, probs=c(0,.25,.50,.75,1)),
#               trees=quantile(coc2$trees, probs=c(0,.25,.50,.75,1)),
#               traffic=quantile(coc2$traffic, probs=c(0,.25,.50,.75,1)),
#               pm25=quantile(coc2$pm25, probs=c(0,.25,.50,.75,1)),
#               slope=quantile(coc2$slope, probs=c(0,.25,.50,.75,1)),
#               nodev=quantile(coc2$nodev, probs=c(0,.25,.50,.75,1)),
#               devAge=quantile(coc2$devAge, probs=c(0,.25,.50,.75,1)),
#               dry=quantile(coc2$dry, probs=c(0,.25,.50,.75,1)),
#               rain=quantile(coc2$rain, probs=c(0,.25,.50,.75,1)),
#               RES=quantile(coc2$RES, probs=c(0,.25,.50,.75,1)),
#               COM_IND=quantile(coc2$COM_IND, probs=c(0,.25,.50,.75,1)),
#               AG=quantile(coc2$AG, probs=c(0,.25,.50,.75,1)),
#               COM=quantile(coc2$COM, probs=c(0,.25,.50,.75,1)),
#               IND=quantile(coc2$IND, probs=c(0,.25,.50,.75,1)) )


#------------------------#
#  Models 1, 3, 4 and 5  #
#------------------------#

vf1X <- varIdent(form= ~1|agency)  #varComb(varIdent(form= ~1|agency), varExp(form= ~dry))
#I tested Model 4a with and without dry (as covariate and as variance covariate); difference in AIC is -4 pts when dry is included,
#  but BIC is +4 pts when dry is included.  Remove "dry", as it seems like a toss-up whether to keep it or not!  Change 
#  variance function to compensate!

#------ Model 1: median value for all locations --------#

Model1 <- gls(data=coc2, result~1, method="ML")  #here, we use the actual concentration, rather than the transformed one
E1 <- residuals(object=Model1, type="normalized")


#------ Model 3: land use + rain with variance structure & random effects ---------#

Form3 <- as.formula(paste("result ~ landuse + ", paste(weather, collapse=" + ")))
Model3 <- lme(data=coc2, Form3, random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
#Model3 <- gls(data=coc2, Form3, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E3 <- residuals(object=Model3, type="normalized")
boxplots.resids2(Model3, E3, "X")

#------ Model 4: no land use; up to 3 predictors + rain ---------#

#run through the top formulas when only best predictors are on the table;
#  keep only formulas that make sense
myForm <- my.formulas.m4[[2]]
#myForm <- Form4a

#these lines of code assess fit of this particular model in terms of COC vs. individual predictors, and predictor correlation
myModel <- lme(data=coc2, myForm, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
#myModel <- gls(data=coc2, myForm, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(myModel)
AIC(myModel)
plot.single.preds(myModel)
check.cor(myModel)
boxplots.resids2(myModel, residuals(myModel, type="normalized"), "X")

# #compare mixed effects model to linear model for CO2_road
# myLM <- gls(data=coc2, myForm, method="ML", weights=vf1X)
# boxplots.resids2(myLM, residuals(myLM, type="normalized"), "X")
# AIC(myModel, myLM)

#formulas that are worth considering (single plots of predictors make sense)
Form4a <- formula(result ~ rain + summer + grass + paved + sqrt_CO2_road) #my.formulas.m4[[1]]  #AIC=830.3; lme with grass + paved + sqrt_CO2_road
Form4b <- formula(result ~ rain + summer + grass + paved) #my.formulas.m4[[3]]  #AIC=839.8; lme with grass + paved
Form4c <- formula(result ~ rain + summer + sqrt_CO2_road)
# #for GLS model:
# Form4a <- my.formulas.m4[[1]]  #AIC=962.8; gls with grass + paved + sqrt_CO2_road
# Form4b <- my.formulas.m4[[2]]  #AIC=964.9; gls with grass + devAge + sqrt_CO2_road
# Form4c <- my.formulas.m4[[5]]  #AIC=980.8; gls with grass + sqrt_CO2_road

#####  Best fit Model4  ####

Form4 <- Form4a
Model4 <- lme(data=coc2, Form4, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4 <- residuals(object = Model4, type = "normalized")

#try adding interactions and see if they are significant
M4.sub1 <- update(Model4, . ~ . + rain:grass) #not significant (p=0.6869)
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, . ~ . + rain:paved)  #not significant (p=0.9962)
anova(Model4, M4.sub2)
M4.sub3 <- update(Model4, . ~ . + grass:paved)  #not significant (p=0.334)
anova(Model4, M4.sub3)


Form4.alt <- Form4a
Model4.alt <- lme(data=coc2, Form4.alt, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4.alt <- residuals(object = Model4.alt, type = "normalized")


#------ Compare models and plot ---------#
AIC(Model1, Model3, Model4, Model4.alt)

AIC(Model1, Model3, Model4, Model4.alt)[,2] - AIC(Model4)


if(run_exploratory_code==TRUE) {
  #plot model predictions for each model (above)
  M1.preds <- predict(Model1)
  plot.preds.vs.results(M1.preds)
  
  M3.preds <- predict(Model3)
  plot.preds.vs.results(M3.preds)
  
  M4.preds <- predict(Model4)
  plot.preds.vs.results(M4.preds)
  
  M5.preds <- predict(Model5)
  plot.preds.vs.results(M5.preds)
}


#----------------------------------#
#  Look at the Fits of the Models  #
#----------------------------------#

#check residuals for each model, to see if they tell us anything important
par(mfrow=c(2,2), mar=c(2,4,4,1), oma=c(0,0,0,0))
plot(coc2$location, E1, main="Null Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,5,7,9,12,15), labels=c("King", "Pierce", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
plot(coc2$location, E3, main="Categorical Landuse Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,5,7,9,12,15), labels=c("King", "Pierce", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
plot(coc2$location, E4, main="Landscape Predictor Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,5,7,9,12,15), labels=c("King", "Pierce", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
# plot(coc2$location, E4.alt, main="Model 4a (grass + paved + CO2_road)", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,2,2,3,4,4,4,5,5,5,6,6,6)])
# axis(side=1, at=c(2,5,7,9,12,15), labels=c("King", "Pierce", "POT", "Sea", "Sno", "Tac"))
# abline(h=0, col="gray")


# Step 9: Refit models with REML and apply graphical model validation; check for homogeneity, 
#--------    normality and independence

#------ Model 3: landuse only, with variance structure & random structure ---------#

M3.final <- lme(data=coc2, Form3, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E3.final <- residuals(object=M3.final, type="normalized")

#plot residuals against agency, year, month, date and other unused predictors; check if model is missing something important!
boxplots.resids2(M3.final, E3.final, "X")
#is temporal auto-correlation present?  Look for patterns in Auto-correlation plot for residuals
par(mfrow=c(2,2))
plot(E3.final~coc$start_date, pch=16, col="cadet blue")
acf(E3.final, na.action=na.pass, main="Auto-correlation plot for residuals")  #look for lines extending past blue dashed range
#Check for spatial correlation
mydata2 <- data.frame(E3.final, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
bubble(mydata2,"E3.final",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
#variogram to test for spatial correlation
Vario1 = variogram(E3.final ~ 1, mydata2)
plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
Vario2 <- variogram(E3.final ~ 1, mydata2, alpha = c(0, 45, 90,135) )
plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?
#there maybe evidence of spatial auto-correlation with this model!  Evidence that it is lacking some predictors


#------ Model 4: no land use; up to 3 predictors + rain ---------#

M4.final <- lme(data=coc2, Form4, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4.final <- residuals(object=M4.final, type="normalized")

#plot residuals against agency, year, month, date and other unused predictors; check if model is missing something important!
boxplots.resids2(M4.final, E4.final, "X")
#is temporal auto-correlation present?  Look for patterns in Auto-correlation plot for residuals
par(mfrow=c(2,2))
plot(E4.final~coc$start_date, pch=16, col="cadet blue")
acf(E4.final, na.action=na.pass, main="Auto-correlation plot for residuals")  #look for lines extending past blue dashed range
#Check for spatial correlation
mydata2 <- data.frame(E4.final, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
bubble(mydata2,"E4.final",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
#variogram to test for spatial correlation
Vario1 = variogram(E4.final ~ 1, mydata2)
plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
Vario2 <- variogram(E4.final ~ 1, mydata2, alpha = c(0, 45, 90,135) )
plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?
#for Phosphorus, this all looks fine; bubble plot isn't perfect, but not too bad


#------ Model 4, alt: up to 3 predictors + rain ---------#

M4alt.final <- lme(data=coc2, Form4.alt, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4alt.final <- residuals(object=M4alt.final, type="normalized")

#plot residuals against agency, year, month, date and other unused predictors; check if model is missing something important!
boxplots.resids2(M4alt.final, E4alt.final, "X")
#is temporal auto-correlation present?  Look for patterns in Auto-correlation plot for residuals
par(mfrow=c(2,2))
plot(E4alt.final~coc$start_date, pch=16, col="cadet blue")
acf(E4alt.final, na.action=na.pass, main="Auto-correlation plot for residuals")  #look for lines extending past blue dashed range
#Check for spatial correlation
mydata2 <- data.frame(E4alt.final, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
bubble(mydata2,"E4alt.final",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
#variogram to test for spatial correlation
Vario1 = variogram(E4alt.final ~ 1, mydata2)
plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
Vario2 <- variogram(E4alt.final ~ 1, mydata2, alpha = c(0, 45, 90,135) )
plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?


#-----------------#
#  Save Formulas  #
#-----------------#

#save important items with phosphorus-specific names
P.coc2 <- coc2
P.r1X <- r1X
P.vf1X <- vf1X
P.Form3 <- Form3
P.Form4 <- Form4
P.rain <- rain

save(P.coc2, P.r1X, P.vf1X, P.Form3, P.Form4, rain, file=here("results", "Total Phosphorus Models.RData"))


#-------------------------------#
#  Plot Predictor Coefficients  #
#-------------------------------#

#generate plots that are shown in Rmarkdown script
P.null <- gls(data = coc2, result ~ 1, method = "REML") 
P.M3 <- lme(data = coc2, Form3, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
P.M4 <- lme(data = coc2, Form4, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
#P.M4.alt <- lme(data = coc2, Form4.alt, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

P_models <- list(
  Null_Model = P.null,
  Categorical_Landuse_Model = P.M3,
  Landscape_Predictor_Model = P.M4
)

huxtablereg(P_models,
            single.row = TRUE, custom.model.names = names(P_models)) %>%
  set_bottom_border(1, -1, 0.4) %>%
  set_bold(1, -1, TRUE) 

theme_set(theme_sjplot())
plotreg(P_models, custom.title = "Regression Results, Total Phosphorus", custom.model.names = names(P_models))
plot_models(P_models, m.labels = names(P_models),legend.title = "Models", show.values = TRUE,show.intercept = TRUE)


#----------------------------------#
#  Plots to support Eva's summary  #
#----------------------------------#

# NONE :)





