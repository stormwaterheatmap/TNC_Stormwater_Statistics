# This generalized script assesses linear relationships between copper and landscape predictors
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
# July 16, 2021
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
source(here("functions", "COC analysis functions.R"))
source(here("functions", "HighstatLibV10.R")) #Zuur library, incl panel.smooth2, VIF


#toggle for whether all exploratory parts of code should be run (TRUE) or just essential parts (FALSE) 
run_exploratory_code <- FALSE


#read in stormwater data with spatial predictors-------------------------------
s8 <- read.csv(here("processed_data", "s8data_with_spatial_predictors.csv"))
s8$start_date <- as.Date(s8$start_date)
s8$conc <- s8$result
s8$month <- as.factor(s8$month)
s8$location_id <- as.factor(s8$location_id)
s8$loc <- as.factor(s8$loc)
s8$agency <- as.factor(s8$agency)
s8$season <- as.factor(s8$season)
s8$land_use <- as.factor(s8$land_use)

#list of landscape predictors-------------------------------
predictors <- names(s8)[34:61]

#select the chemical of interest
this_param <- "Copper - Water - Total"
this_param_short <- "Copper"
coc <- s8[which(s8$parameter == this_param), ]  #create the "coc" dataframe for this chemical of interest only

#Non Detect Handing -------------------------------
#look at censored points
ggplot(coc)+geom_jitter(aes(x=agency,y=conc,color=nondetect_flag))+scale_y_log10()

#look into which locations/ dates had ND data
coc.nd <- coc[which(coc$nondetect_flag == TRUE), ]  #after removing PIE_LDR, PIE_HDR: 4 ND results, all 0.1
pct.cen <- 100*nrow(coc.nd)/nrow(coc)  #for copper, 1.75% of samples are ND's
nrow(coc.nd)
min(coc.nd$conc)
max(coc.nd$conc)
coc.nd[, c("location_id", "start_date", "conc")]  #location, date, and detection limit of ND samples

#NOTE: For copper, all April-Sept 2009 samples at SNO_HDR are ND's at 0.1ppm.  Samples collected
#      at this location in Oct, Nov and Dec are 3.9, 4.5, 4.6 and 27.9.  All samples for 2010-2013
#      collected at SNO_HDR are > 2.  Furthermore, other metals tested at SNO_HDR between April-Sept
#      2009 were also ND, when no later samples were ND.  Remove these samples -- they are likely 
#      field blanks that got misclassified as samples

#remove likely field blanks
coc <- coc[-(which(coc$location=="SNO_HDR" & coc$nondetect_flag==TRUE & coc$year==2009)),]

#Since only 5 ND's remain (1% of all samples samples), substitute with half of the detection limit
## AFTER REMOVING PIE_LDR and PIE_HDR, no ND's left!
# coc <- coc %>% 
#   mutate(conc = ifelse(nondetect_flag, conc * 0.5, conc))

# distribution exploration ----------------------------------------------
if (run_exploratory_code==TRUE) {
  #explore the data; look for underlying distribution
  par(mfrow = c(2, 2))
  hist(coc$conc, breaks = 50)
  qqnorm((coc$conc), main = paste("QQ-Normal plot of", this_param_short))
  qqline((coc$conc))
  hist(log(coc$conc), breaks = 50)
  qqnorm(log(coc$conc), main = paste("QQ-Normal plot of log(", this_param_short, ")", sep = ""))
  qqline(log(coc$conc))
  
  #consider other distribution (log-normal and gamma).  Compare QQ plots for these distributions
  par(mfrow = c(2, 2))
  qqPlot(log(coc$conc), dist = "norm", estimate.params = TRUE, add.line = TRUE, 
         main = paste("QQ log-normal plot of", this_param_short))
  qqPlot(log(coc$conc), dist = "norm", param.list = list(mean = -3.3, sd = 2.09), add.line = TRUE,
         main = paste("QQ log-normal plot of", this_param_short))
  qqPlot(coc$conc, dist = "gamma", param.list = list(shape = .7, scale = 1), points.col = "blue", add.line = TRUE, 
         main = paste("QQ-gamma plot of", this_param_short))
  qqPlot(sqrt(coc$conc), dist = "norm", estimate.params = TRUE, add.line = TRUE, 
         main = paste("QQ plot of sqrt-transformed", this_param_short))
  #qqPlot(1/(coc$conc), dist="norm", estimate.params=TRUE, add.line=TRUE)
  
  # for both the log-normal and gamma distributions, the high tail is skewed.  lets find out where the
  # high data came from...  ###note: this skew decreased after removing PIE_LDR, PIE_HDR
  coc %>%
    dplyr::filter(parameter == this_param, result > quantile(coc$conc, 0.95)) %>%
    dplyr::select(c(location_id, result, month, year)) %>%
    dplyr::arrange(location_id, result)
  # for copper, 18 of 24 come from SEAC1S8D_OUT; 5 come from TAC001s8D_OF235.  No trend in month & year
}

# best distribution is log-normal.  Ln-transform concentration data in column "result"
coc$result <- log(coc$conc)


#-------------------------------------------------------------------------------------------#
#  Explore Cleveland Dotplots and Boxplots of COC data or predictors conditional on Agency  #
#-------------------------------------------------------------------------------------------#

colors_agency <- c("red", "orange", "yellow", "green", "blue", "purple")

if (run_exploratory_code==TRUE) {
  #Cleveland Dotplot - use to detect violations of homogeneity: We are looking for whether the spread of data values
  #   differs between sampling locations (or agencies).  If so, it indicates heterogeneity, and that there may be problems 
  #   with violation ofhomogeneity in a linear regression model applied on these data.  We're also looking for outliers.
  par(mfrow = c(1, 1))
  dotchart(coc$result, groups = coc$loc, pch = 19, col = colors_agency[as.numeric(coc$agency)],
           xlab = "concentration", main = paste("Cleveland Dotplot:", this_param_short))
  
  #pairs plots allow visualization of interactions between possible predictors.  Look for relationships between "result"
  #   and all of the predictor variables
  pairs(coc %>% select(result, predictors[1:12]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  pairs(coc %>% select(result, predictors[13:24]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  pairs(coc %>% select(result, predictors[25:35]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)

  #look also at interactions between result and rainfall
  pairs(coc %>% select(result, antecedant_dry_days_std, daymet_precip_std, daymet_3day_std,
                       daymet_7day_std, daymet_14day_std, daymet_21day_std, daymet_28daySR_std), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  #for copper, AFTER removal of PIE_HDR & PIE_COM, most important predictors are:
  #  intURB, totRES, nodev, traffic, pm25_na, CO2_tot, CO2_com, CO2_road, devAge, roof_intURB, roof_intURB_IND
  #  daymet14, 21 or 28day (all -0.2)
}

#-----------------------------------------------------------------------------------------------------#
#  Plot COC vs various predictors: look for potential predictors to model & sources of heterogeneity  #
#-----------------------------------------------------------------------------------------------------#

#look for relationships between coc and landscape predictors
lp1 <- my.ggplot(1)
lp2 <- my.ggplot(2)
lp3 <- my.ggplot(3)
lp4 <- my.ggplot(4)
lp5 <- my.ggplot(5)
lp6 <- my.ggplot(6)
lp7 <- my.ggplot(7)
lp8 <- my.ggplot(8)
lp9 <- my.ggplot(9)
lp10 <- my.ggplot(10)
lp11 <- my.ggplot(11)
lp12 <- my.ggplot(12)
lp13 <- my.ggplot(13)
lp14 <- my.ggplot(14)
lp15 <- my.ggplot(15)
lp16 <- my.ggplot(16)
lp17 <- my.ggplot(17)
lp18 <- my.ggplot(18)
lp19 <- my.ggplot(19)
lp20 <- my.ggplot(20)
lp21 <- my.ggplot(21)
lp22 <- my.ggplot(22)
lp23 <- my.ggplot(23)
lp24 <- my.ggplot(24)
lp25 <- my.ggplot(25)
lp26 <- my.ggplot(26)
lp27 <- my.ggplot(27)
lp28 <- my.ggplot(28)

#relationship between COC and precipitation
pr1 <- ggplot(coc, aes(agency, result)) + geom_boxplot()
pr2 <- ggplot(coc, aes(daymet_precip_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_precip_std))
pr3 <- ggplot(coc, aes(daymet_3day_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_3day_std))
pr4 <- ggplot(coc, aes(daymet_7day_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_7day_std))
pr5 <- ggplot(coc, aes(daymet_14day_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_14day_std))
pr6 <- ggplot(coc, aes(daymet_21day_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_21day_std))
pr7 <- ggplot(coc, aes(daymet_28daySR_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_28daySR_std))
pr8 <- ggplot(coc, aes(antecedant_dry_days_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$antecedant_dry_days_std))
# NOTE: I have left out precip & inches_rain_per_hour b/c they are missing ~68 data points

#relationship between COC and temporal/ location predictors
tlp1 <- ggplot(coc, aes(as.factor(year), result)) + geom_boxplot()  #starts in Feb, 2009, ends in April, 2013; trends could be due to months with data?
tlp2 <- ggplot(coc, aes(month, result)) + geom_boxplot(position=position_dodge(width=0.9))
a <- coc %>% count(month)
tlp2 <- tlp2 + annotate(geom="text", x=c(1:12), y=rep(-1.4,12), label=paste("n=", a[,2], sep=""),
                      color="red", size=3.5, angle=90)  #note the small sample size in June-Sept
tlp3 <- ggplot(coc, aes(as.factor(season), result)) + geom_boxplot()
tlp4 <- ggplot(coc, aes(land_use, result)) + geom_boxplot()
tlp5 <- ggplot(coc, aes(location_id, result)) + geom_boxplot()
  
if(run_exploratory_code==TRUE) {
    #Look at various transformations for monthly precip by location; which is most evenly spread out?
    par(mfrow=c(2,2), mar=c(4,4,4,2))
    plot(coc$mPrecip, coc$result)
    plot(coc$mPrecipSR, coc$result)
    plot(coc$mPrecipCR, coc$result)
    plot(coc$mPrecipLog, coc$result)
    #for copper, mPrecip is better than the rest!
    
    #look for relationships between coc's and some predictors; 
    # look especially for sources of heterogeneity
    grid.arrange(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, nrow=3)  #14-day & 21-day look the best -- most evenly spaced data...
    grid.arrange(pr1, tlp5, tlp1, tlp2, tlp3, tlp4, nrow=2)
    #for copper, 28-day prior rainfall does not show heteroskedasticity.
    #   agency, landuse, do show heteroskedasticity
    
    grid.arrange(lp1, lp2, lp3, lp4, lp5, lp6, lp7, lp8, lp9, lp10, lp11, lp12, nrow=3, ncol=4)
    grid.arrange(lp13, lp14, lp15, lp16, lp17, lp18, lp19, lp20, lp21, lp22, lp23, lp24, nrow=3, ncol=4)
    grid.arrange(lp26, lp27, lp28, nrow=3, ncol=4)
}
#best predictors for this COC; make sure to only have one version (transformed or not) of each predictor!
best_predictors <- c("intURB", "intURB_IND", "totRES", "grass", "greenery", "impervious", "nodev",
                     "traffic", "sqrt_popn", "pm25_na", "sqrt_CO2_tot", "sqrt_CO2_com", "sqrt_CO2_road", 
                     "sqrt_CO2_nonroad", "devAge2", "roof_intURB_IND")

pred_i <- which(predictors %in% best_predictors)
# grid.arrange(
#   for(i in pred_i) {
#     get(paste("lp", i, sep=""))
#   },
#   nrow=4, ncol=4)

grid.arrange(lp3, lp4, lp5, lp6, lp7, lp10, lp11, lp13, lp14, lp16, lp20, lp21, 
             lp22, lp23, lp24, lp27, nrow=4, ncol=4)

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

elements_to_remove <- c("intURB", "sqrt_roof_intURB_IND", "greenery", "impervious", 
                        "sqrt_CO2_road", "sqrt_CO2_com", "sqrt_CO2_nonroad", "devAge2")  #these predictors are highly correlated with others
best_predictors2 <- best_predictors[!(best_predictors %in% elements_to_remove)]


#--------------------------------------------------#
#  Any evidence of changes in COC conc over time?  #
#--------------------------------------------------#

if (run_exploratory_code==TRUE) {
  #for each agency, is there a trend in COC concentration over time?
  library(lattice)
  xyplot(result~start_date|agency, data=coc,
         xlab="time", ylab="log(conc)",
         strip=function(bg="white",...)
           strip.default(bg="white",...),
         panel=function(x,y) {
           panel.grid(h=-1, v=2)
           I1 <- order(x)
           llines(x[I1], y[I1], col=1)
         })
  #   For copper, no long-term trend in concentration over time
  #   We also see that the log(conc) is generally higher in some agencies than others
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
  #for copper...
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

rain <- "daymet 21-day precip, standardized"  #precip measure that was used for this COC


#------------------------------------------------#
#  Formulas for Possible Predictor Combinations  #
#------------------------------------------------#

#precip / season predictors that should be added to the model
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
if (run_exploratory_code==TRUE) {
  M.gls1X <- gls(FormX, data=coc2, method="REML")  #use REML for comparing diff't random variances
  E.1X <- resid(object=M.gls1X, type="normalized")  
  boxplots.resids2(M.gls1X, E.1X, "X")  #check for heterogeneity in residuals vs fitted values; look for source in other plots

  par(mfrow=c(1,1), mar=c(4,4,4,1))
  plot(fitted(M.gls1X), E.1X, main="", xlab="fitted values", xaxt="s", ylab="normalized residuals", col="gray", pch=16)
    
  #is there heteroscedasticity in any of these potential covariates?
  plot(coc2$result ~ coc2$location)
  plot(coc2$result ~ coc2$agency)
  plot(coc2$result ~ coc2$month)  #unlikely to be considered for use as a covariate
  plot(coc2$result ~ coc2$season)
  plot(coc2$result ~ coc2$year)
  plot(coc2$result ~ coc2$dry)
  abline(lm(coc2$result ~ coc2$dry))
  plot(coc2$result ~ coc2$rain7)
  plot(coc2$result ~ coc2$rain28)
  plot(coc2$result ~ coc2$rain)
  abline(lm(coc2$result ~ coc2$rain))
  plot(coc2$result ~ coc2$monthlyPrecip)
  plot(coc2$result ~ coc2$mPrecip)
  plot(coc2$result ~ coc2$paved)
  plot(coc2$result ~ coc2$grass)
  #for copper, rain is the only covariate that might be a source of heteroskedasticity
  
  
  #Step 3:  Choose a variance structure (if there was heterogeneity)
  #-------  for selecting random structure, use REML to compare - its better at capturing random structure
  #         REML estimates the random effects by considering linear combinations of the data that remove the 
  #         fixed effects. If these fixed effects are changed, the likelihoods of the two models will not be directly comparable.
  #         Compare AIC for the various beyond-optimal models with different variance structures.
  #         Alternately, can use likelihood ratio test (anova) to compare nested models
  
  M.gls2X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency))
#  M.gls3X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|location))
  M.gls4X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|season))
#  M.gls5X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency * season))
  M.gls18X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varExp(form= ~rain)))  #agency & previous 7-day rainfall (SR-xformed), unstandardized
  M.gls19X <- gls(FormX, data=coc2, method="REML", weights=varComb(varIdent(form= ~1|agency), varConstPower(form= ~rain)))
  M.gls20X <- gls(FormX, data=coc2, method="REML", weights=varConstPower(form= ~rain))
  M.gls21X <- gls(FormX, data=coc2, method="REML", weights=varExp(form= ~rain))
  # M.gls22X <- gls(FormX, data=coc2, method="REML", weights=varExp(form= ~paved))
  # M.gls23X <- gls(FormX, data=coc2, method="REML", weights=varConstPower(form= ~paved))
  M.gls24X <- gls(FormX, data=coc2, method="REML", weights=varFixed(~rain))

  anova(M.gls1X, M.gls2X,
        M.gls4X, 
        M.gls18X, M.gls19X, M.gls20X, M.gls21X, M.gls24X)
  #for copper: AIC and BIC are best for agency & rain (#19), or for just agency (#2)
  
  #likelihood ratio tests for nested versions of the best fit models
  anova(M.gls1X, M.gls2X, M.gls19X)  ### NOTE: since we're not using dry in the model, don't include it in the variance structure!
  
  #residual plots for best fit models -- look for homogeneity of residuals
  E.19X <- resid(object = M.gls19X, type = "normalized")
  boxplots.resids2(M.gls19X, E.19X, "X")
}

#variance function for best fit models so far
vf1X <- varComb(varIdent(form= ~1|agency), varConstPower(form= ~rain))
#stick with variance structure that only uses variance covariates -- since not using dry in model, don't use in V.S.

# Steps 4-6: Find the proper random effects structure; look for temporal and spatial autocorrelation
#-----------  note: look for spatial correlation AFTER setting random effects!

if (run_exploratory_code==TRUE) {
  #----Model with Variance Function 1----#
  M.vf1X <- gls(FormX, data = coc2, method = "REML",  weights = vf1X)
  E.vf1X <- resid(object=M.vf1X, type="normalized")  
  
  #plots of residuals: no variance structure vs. variance structure
  par(mfrow=c(2,2), mar=c(4,4,4,1))
  plot(fitted(M.gls1X), E.1X, main="no variance structure", xlab="fitted values", xaxt="s", ylab="normalized residuals", col="gray", pch=16)
  plot(fitted(M.vf1X), E.vf1X, main="variance covariates: agency & rain", xlab="fitted values", xaxt="s", ylab="normalized residuals", col="gray", pch=16)
  
  #Random intercept model; this is nested in the best fit model from Step 2, so can compare with a likelhiood ratio test
  M1.lme1X <- lme(data=coc2, FormX, random = ~1|agency, method="REML", weights=vf1X,
                  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  anova(M.vf1X, M1.lme1X)  #the likelihood ratio test p-value compares the two nested models, and whether they are significantly different
  #For copper, the model with random intercept is significantly better than the model with no random effects
}

#random structure for best fit model so far
r1X <- formula( ~ 1 | agency)

if (run_exploratory_code==TRUE) {
  #best random effects structure & residual plots to test it
  M1.rX <- lme(data=coc2, FormX, method="REML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  summary(M1.rX)
  E1.rX <- residuals(object=M1.rX, type="normalized")
  boxplots.resids2(M1.rX, E1.rX, "X")
  
  #Check for temporal correlation; look for patterns in Auto-correlation plot for residuals
  par(mfrow=c(2,2))
  plot(E1.rX~coc2$start_date, pch=16, col="cadet blue")
  acf(E1.rX, na.action=na.pass, main="Auto-correlation plot for residuals")  #look for lines extending past blue dashed range
  #for copper, no indication of temporal autocorrelation
  
  #Check for spatial correlation
  mydata2 <- data.frame(E1.rX, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
  coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
  bubble(mydata2,"E1.rX",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
  
  #variogram to test for spatial correlation
  Vario1 = variogram(E1.rX ~ 1, mydata2)
  plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
  Vario2 <- variogram(E1.rX ~ 1, mydata2, alpha = c(0, 45, 90,135) )
  plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?
  #for copper, no indication of temporal autocorrelation
}

# Steps 7-8: Find the proper fixed effects structure; use ML for model comparisons
#-----------

#-------------------------------------------------------------#
#  Exhaustive Search for Model with Best 3 Predictors + Rain  #
#-------------------------------------------------------------#

if (run_exploratory_code==TRUE) {
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
  
  # dd.str <- paste(dd$Var1, dd$Var2, dd$Var3)  #vector with strings of the predictor sets
  
  for (i in 1:length(best_predictors)) {
    aa <- which(cor(coc2[, best_predictors]) [, i] >= 0.85)
    aa <- names(aa)
    for (j in 1:length(aa)) {
      bb <- which(dd[, 1] == best_predictors[i] & dd[, 2] == aa[j] | 
                  dd[, 1] == best_predictors[i] & dd[, 3] == aa[j])
      if (length(bb) > 0) {
        dd <- dd[-bb, ]
      }
    }
  }
  

  # #identify and eliminate situations with predictors that don't go together, or are highly correlated
  # rem1 <- which(str_detect(dd.str, "sqrt_CO2_tot")==TRUE & str_detect(dd.str, "sqrt_CO2_road")==TRUE)
  # rem2 <- which(str_detect(dd.str, "sqrt_popn")==TRUE & str_detect(dd.str, "sqrt_CO2_res")==TRUE)
  # rem3 <- which(str_detect(dd.str, "grass")==TRUE & str_detect(dd.str, "impervious")==TRUE)
  # rem4 <- which(str_detect(dd.str, "paved")==TRUE & str_detect(dd.str, "impervious")==TRUE)
  # rem5 <- which(str_detect(dd.str, "nodev")==TRUE & str_detect(dd.str, "devAge")==TRUE)
  # rem6 <- which(str_detect(dd.str, "sqrt_nodev")==TRUE & str_detect(dd.str, "devAge")==TRUE)
  # rem7 <- which(str_detect(dd.str, "sqrt_CO2_tot")==TRUE & str_detect(dd.str, "sqrt_CO2_com")==TRUE)
  # rem8 <- which(str_detect(dd.str, "grass")==TRUE & str_detect(dd.str, "greenery")==TRUE)
  # rem9 <- which(str_detect(dd.str, "greenery")==TRUE & str_detect(dd.str, "impervious")==TRUE)
  # rem10 <- which(str_detect(dd.str, "greenery")==TRUE & str_detect(dd.str, "totRES")==TRUE)
  # rem11 <- which(str_detect(dd.str, "intURB")==TRUE & str_detect(dd.str, "sqrt_intURB")==TRUE)
  # rem.all <- unique(c(rem1, rem2, rem3, rem4, rem5, rem6, rem7, rem8, rem9, rem10))  #keep only non-repetitive instances of indexes that should be removed
  # 
  # dd <- dd[-rem.all, ]  #remove rows with predictors that don't go together!

  #generate an exhaustive set of formulas with one, two and three predictors, including and excluding rain
  ee.1 <- list(NULL)  #list of formulas based on the sets of predictors in "dd"; including rain
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i])) ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1]))
    else if(is.na(dd$Var3[i]))  ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1], "+", dd[i,2]))
    else  ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  # ee.2 <- list(NULL) #list of formulas based on the sets of predictors in "dd"; NOT including rain
  # for(i in 1:nrow(dd)) {
  #   if(is.na(dd$Var2[i]))  ee.2[[i]] <- as.formula(paste("result~ ", dd[i,1]))
  #   else if(is.na(dd$Var3[i]))  ee.2[[i]] <- as.formula(paste("result~ ", dd[i,1], "+", dd[i,2]))
  #   else ee.2[[i]] <- as.formula(paste("result~ ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
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
  
  ee <- ee.1 #c(ee.1, ee.2)  #c(ee.1, ee.2, ee.3, ee.4) 
  
  MyEqns <- list(formulas=ee, aics=rep(0, length(ee)), do_signs_match=rep(FALSE, length(ee)))

  #run the lme for each formula, saving the AIC values; takes about 7 minutes to run
  ptm <- proc.time()  #time the code!
  my.aics <- rep(0, length(ee))
  for (i in 1:length(ee)) {
    bb <- lme(data=coc2, ee[[i]], random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
    my.aics[i] <- AIC(bb)
    MyEqns$aics[i] <- AIC(bb)
    #extract landscape predictor coefficients 
    aa <- summary(bb)$coefficients$fixed
    aa <- aa[names(aa) %in% preds.1]
    aa <- sign(aa)
    aa <- aa[sort(names(aa))]
    #extract coefficient signs for linear models of result~landscape predictors
    cc <- bp_signs[names(bp_coefs) %in% names(aa)]
    cc <- cc[sort(names(cc))]
    #compare coefficient signs from this model, to those of single-predictor linear models
    MyEqns$do_signs_match[i] <- identical(aa, cc)  #if coefficient signs match, TRUE; (default=FALSE)

  }
  proc.time() - ptm  #stop the clock;  9 minutes for 386 models

 
  xx <- which(MyEqns$do_signs_match==TRUE)

  ee.all <- ee
  ee <- ee[xx]
  my.aics.all <- my.aics
  my.aics <- my.aics[xx]
  
  # ee <- ee.all
  # my.aics <- my.aics.all
      
  par(mfrow=c(5,5), mar=c(1,2,2,0), oma=c(1,1,1,1))
  gg <- rep("black", length(ee))
  gg[c(which(is.na(dd$Var3)), which(is.na(dd$Var3))+nrow(dd), which(is.na(dd$Var3))+2*nrow(dd), which(is.na(dd$Var3))+3*nrow(dd))] <- "turquoise"   #make formulas with only 2 landscape predictors red
  gg[c(which(is.na(dd$Var2)), which(is.na(dd$Var2))+nrow(dd), which(is.na(dd$Var2))+2*nrow(dd), which(is.na(dd$Var2))+3*nrow(dd) )] <- "pink"  #make formulas with only 1 landscape predictor yellow
  plot(my.aics, col=gg, pch=16, main="1 pred=pink, 2 preds=turquoise", xaxt="n", ylab="AIC")
  gg <- rep("black", length(ee))
  gg[which(str_detect(as.character(ee), "impervious"))] <- "red"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(ee), "paved"))] <- "orange"   #make formulas with "rain" light blue
  if(sum(str_detect("paved", preds.1) | str_detect("impervious", preds.1))>0)  plot(my.aics, col=gg, xaxt="n", yaxt="n", pch=16, main="impervious=red, paved=orange")
  gg <- rep("black", length(ee))
  gg[which(str_detect(as.character(ee), "rain"))] <- "light blue"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(ee), "dry"))] <- "goldenrod"   #make formulas with "dry" goldenrod
  gg[which(str_detect(as.character(ee), "dry") & str_detect(as.character(ee), "rain"))] <- "green"   #make formulas with "dry" goldenrod
  plot(my.aics, col=gg, xaxt="n", yaxt="n", pch=16, main="rain=blue, dry=gold, both=green")
  #  plot(my.aics, col=ifelse(str_detect(as.character(ee), "rain"), "light blue", "black"), xaxt="n", yaxt="n", pch=16, main="rain=light blue")

  colAIC <- c("red", "orange", "yellow", "light green", "green3", "cadet blue", "purple1", "tan3", "salmon", "maroon1", "goldenrod", "light blue", "pink",
              "red", "orange", "yellow", "light green", "green3", "cadet blue", "purple1", "tan3", "salmon", "maroon1", "goldenrod", "light blue", "pink")
  for (i in 1:length(best_predictors)) {
    plot(my.aics, col=ifelse(str_detect(as.character(ee), best_predictors[i]), colAIC[i], "black"), xaxt="n", yaxt="n", pch=16, main=paste(best_predictors[i], "=", colAIC[i]) )
  }

  plot.new()
  text(x=0.5, y=0.5, labels=paste("Landuse AIC", 
                                  round(AIC(lme(data=coc2, as.formula(paste("result~landuse+", paste(weather, collapse="+"))), random = r1X, 
                                                method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))), 1)), cex=2)
  
  #arrange formulae and AICs according to order of AIC value, from smallest to largest
  my.formulas.m4 <- ee[order(my.aics)]
  my.formulas.m4[1:10]
  my.aics.m4 <- my.aics[order(my.aics)]
  my.aics.m4[1:10]  #what are the top 10 AIC's?

  save(my.formulas.m4, my.aics.m4, file=(here("scripts", "formulas_and_aics", "total_copper_m4.RData")))
}

if (run_exploratory_code==FALSE) {
  load(file=here("scripts", "formulas_and_aics", "total_copper_m4.RData"))
}

#identify which (ranked) formulas have "rain" in them
aa <- grep("rain", my.formulas.m4)  #ordered list (by AIC value) of all formulae that include rain
my.aics.m4[aa]
my.formulas.m4[aa]

#identify which (ranked) formulas have 2 or less predictors
bb <- as.character(my.formulas.m4[aa])
formulas.m4.two <- my.formulas.m4[aa][which(lengths(gregexpr("\\W+", bb))< 4)]
aics.m4.two <- my.aics.m4[aa][which(lengths(gregexpr("\\W+", bb))< 4)]
FA.two <- cbind(as.character(formulas.m4.two), round(aics.m4.two, 1))



#--------------------------------------------#
#  Quantiles List for Plotting Interactions  #
#--------------------------------------------#

#list of quantiles for all landscape predictors
qList <- list(
  totRES=quantile(coc2$totRES, probs=c(0,.25,.50,.75,1)),
  traffic=quantile(coc2$traffic, probs=c(0,.25,.50,.75,1)),
  pm25_na=quantile(coc2$pm25_na, probs=c(0,.25,.50,.75,1)),
  devAge2=quantile(coc2$devAge2, probs=c(0,.25,.50,.75,1)),
  rain=quantile(coc2$rain, probs=c(0,.25,.50,.75,1))
)
  
  
  # paved=quantile(coc2$paved, probs=c(0,.25,.50,.75,1)),
  #             impervious=quantile(coc2$impervious, probs=c(0,.25,.50,.75,1)),
  #             roofs=quantile(coc2$roofs, probs=c(0,.25,.50,.75,1)),
  #             roof_RES=quantile(coc2$roof_RES, probs=c(0,.25,.50,.75,1)),
  #             roof_nonRES=quantile(coc2$roof_nonRES, probs=c(0,.25,.50,.75,1)),
  #             grass=quantile(coc2$grass, probs=c(0,.25,.50,.75,1)),
  #             trees=quantile(coc2$trees, probs=c(0,.25,.50,.75,1)),
  #             traffic=quantile(coc2$traffic, probs=c(0,.25,.50,.75,1)),
  #             pm25=quantile(coc2$pm25, probs=c(0,.25,.50,.75,1)),
  #             slope=quantile(coc2$slope, probs=c(0,.25,.50,.75,1)),
  #             nodev=quantile(coc2$nodev, probs=c(0,.25,.50,.75,1)),
  #             devAge=quantile(coc2$devAge, probs=c(0,.25,.50,.75,1)),
  #             dry=quantile(coc2$dry, probs=c(0,.25,.50,.75,1)),
  #             rain=quantile(coc2$rain, probs=c(0,.25,.50,.75,1)),
  #             RES=quantile(coc2$RES, probs=c(0,.25,.50,.75,1)),
  #             COM_IND=quantile(coc2$COM_IND, probs=c(0,.25,.50,.75,1)),
  #             AG=quantile(coc2$AG, probs=c(0,.25,.50,.75,1)),
  #             COM=quantile(coc2$COM, probs=c(0,.25,.50,.75,1)),
  #             IND=quantile(coc2$IND, probs=c(0,.25,.50,.75,1)) )
  #                                       

#--------------------#
#  Models 1, 3, & 4  #
#--------------------#

#------ Model 1: median value for all locations --------#

Model1 <- gls(data=coc2, result~1, method="ML")  #here, we use the actual concentration, rather than the transformed one
E1 <- residuals(object=Model1, type="normalized")


#------ Model 3: land use + rain with variance structure & random effects ---------#

Form3 <- as.formula(paste("result ~ landuse + ", paste(weather, collapse=" + ")))
Model3 <- lme(data=coc2, Form3, random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E3 <- residuals(object=Model3, type="normalized")


#------ Model 4: up to 3 predictors + precip/ season components ---------#

#run through the top formulas when only best predictors are on the table;
#  keep only formulas that make sense
myForm <- my.formulas.m4[[7]]
myForm <- update(Form4b, .~. +rain:pm25_na)

#these lines of code assess fit of this particular model in terms of COC vs. individual predictors, and predictor correlation
myModel <- lme(data=coc2, myForm, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(myModel)
AIC(myModel)
plot.single.preds(myModel)
check.cor(myModel)
boxplots.resids2(myModel, residuals(myModel, type="normalized"), "X")

#formulas that are worth considering (single plots of predictors make sense)
Form4a <- formula(result ~ rain + summer + totRES + traffic) #AIC=743.4; model #3; first two had a 3rd predictor only slightly significant
Form4b <- formula(result ~ rain + summer + traffic + devAge2 + pm25_na)  #AIC=758.0; model #7
Form4c <- formula(result ~ rain + summer + traffic + devAge2 + sqrt_CO2_road)  #AIC=768.8; model #15

#####  Best fit Model4  ####

Form4 <- Form4a
Model4 <- lme(data=coc2, Form4, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
Model4a <- Model4
E4a <- residuals(object = Model4a, type = "normalized")

#try adding interactions and see if they are significant
M4.sub1 <- update(Model4, . ~ . + rain:traffic) #not significant (p=0.4437)
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, . ~ . + rain:totRES)  #not significant (p=0.4697)
anova(Model4, M4.sub2)
M4.sub3 <- update(Model4, . ~ . + traffic:totRES)  #highly significant (p<0.0001)
anova(Model4, M4.sub3)

#Plot formula 4 interaction
Plot.Quantile("totRES", "traffic", M4.sub3, yEqn=5)
#This relationship seems plausible, both on basis of data theoretically.  It also makes a big difference to AIC
#  places that have the most traffic have a more dramatic reduction in copper as residential percentage land use increases
#  places that have little traffic are less affected by whether the properties are zoned as residential or not

Form4a.int <- update(Form4, .~. + totRES:traffic)
M4a.int <- update(Model4, .~. + totRES:traffic)
E4a.int <- residuals(object = M4a.int, type = "normalized")
plot.single.preds(M4a.int)
boxplots.resids2(M4a.int, E4a.int, "X")


#####  Alternate Model4  ####

Form4 <- Form4b
Model4 <- lme(data=coc2, Form4, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
Model4b <- Model4
E4b <- residuals(object = Model4b, type = "normalized")

#try adding interactions and see if they are significant
M4.sub1 <- update(Model4, . ~ . + rain:traffic) #not significant (p=0.1576)
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, . ~ . + rain:devAge2)  #not significant (p=0.7504)
anova(Model4, M4.sub2)
M4.sub3 <- update(Model4, . ~ . + rain:pm25_na)  #highly significant (p=0.0021)
anova(Model4, M4.sub3)

M4.sub4 <- update(Model4, . ~ . + traffic:devAge2) #highly significant (p<0.0001)
anova(Model4, M4.sub4)
M4.sub5 <- update(Model4, . ~ . + pm25_na:devAge2)  #highly significant (p=0.0002)
anova(Model4, M4.sub5)
M4.sub6 <- update(Model4, . ~ . + traffic:pm25_na)  #highly significant (p<0.0001)
anova(Model4, M4.sub6)


#Plot formula 4 interaction
Plot.Quantile("rain", "pm25_na", M4.sub3, yEqn=5)
#This relationship makes sense.  In watersheds where atmospheric pollution is high (red dots), 
#  increasing amounts of rain flush the pm25 out of the air, so you get a steeper decline in copper
#  into stormwater than in watersheds where there is little pollution.

Plot.Quantile("traffic", "devAge2",  M4.sub4, yEqn=5)  
#not buying it!  too little low devAge data

Plot.Quantile("devAge2", "pm25_na",  M4.sub5, yEqn=5)  
#This relationship makes theoretical sense, whereby more/older developed areas that are more heavily
#  air-polluted will have higher copper in stormwater.  However, the data don't compellingly
#  appear to have two different trends based on devAge or pm25 (maybe try viewing with pm25 on x-axis?)

Plot.Quantile("traffic", "pm25_na",  M4.sub6, yEqn=5)  
#not buying it!  traffic and pm25 seem too linked together, and the data don't show compelling evidence
#  for a different slope with traffic, when pm25 is low vs high (maybe try viewing with pm25 on x-axis?)


Form4b.int <- update(Form4b, .~. +rain:pm25_na) #+ devAge2:pm25_na)
M4b.int <- update(Model4, .~. +rain:pm25_na) # + devAge2:pm25_na)
E4b.int <- residuals(object = M4b.int, type = "normalized")
plot.single.preds(M4b.int)
boxplots.resids2(M4b.int, E4b.int, "X")

#----------------------#


Plot.Quantile("devAge2", "pm25_na",  M4.sub5, yEqn=5)  
bb <- "devAge2"
aa <- "pm25_na"
myModel <- M4.sub5
xEqn <- 0.15
yEqn <- 5

#create a list for quantiles of the two landscape predictors aa and bb; this is used to generate IQ
cc <- list(unlist(qList[[aa]]), unlist(qList[[bb]]))
names(cc) <- c(aa, bb)

#generate dataframe of interaction fits for landscape predictors aa and bb
IQ <- effect(paste(bb, "*", aa, sep=""), myModel, xlevels= cc, se=TRUE, confidence.level=0.95, typical=mean)
IQ <- as.data.frame(IQ)
IQ[bb] <- factor(as.numeric(unlist(IQ[bb])), levels=qList[[bb]], labels=c("0%", "25%", "50%", "75%", "100%"))

#linear equation for the interaction, to be printed in the plot
eqnText <- paste("ln(", this_param_short, ") = ", 
                 round(myModel$coefficients[[1]]["(Intercept)"], 2), " + ",
                 round(myModel$coefficients[[1]][aa], 2), "*", aa, " + ",
                 round(myModel$coefficients[[1]][bb], 2), "*", bb, " + ", 
                 round(myModel$coefficients[[1]][paste(bb,":",aa, sep="")], 2), "*", aa, "*", bb, sep="")


ggplot() + 
  geom_line(data=IQ, size=1.5, aes(x=IQ[,aa], y=fit, group=IQ[,bb], color=IQ[,bb]))+
  ylab("COC")+ xlab(names(IQ)[2])+
  ggtitle(paste("Interaction between", names(IQ)[1], "and", names(IQ)[2]))+
  scale_color_manual(values=c("blue", "purple", "yellow", "orange", "red"))+  #colors for lines (color)
  geom_point(data=coc2, aes(x=coc2[,aa], y=result, fill=coc2[,bb]), size=3, shape=21, stroke=0) + 
  scale_fill_gradient2(low = "blue", mid="yellow", high = "red")+  #colors for points (fill)
  geom_text(aes(x=xEqn, y=yEqn, label=eqnText), cex=4.5, color="black")+
  labs(color=names(IQ)[2], fill=names(IQ)[2])








#------ Compare models and plot ---------#

#compare all of the models
AIC(Model1, Model3, Model4a, M4a.int, Model4b, M4b.int)
#obtain the delta AIC value
AIC(Model1, Model3, Model4a, M4a.int,  Model4b, M4b.int)[, 2] - AIC(Model3)


#plot model predictions for each model (above)
if (F) {
  #change from F to T if you want to plot these
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
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
plot(coc2$location, E3, main="Categorical Landuse Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
# plot(coc2$location, E4a, main="Model 4a", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
# axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
# abline(h=0, col="gray")
# plot(coc2$location, E4a.int, main="Model 4a.int", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
# axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
# abline(h=0, col="gray")
# plot(coc2$location, E4b, main="Model 4b", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
# axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
# abline(h=0, col="gray")
plot(coc2$location, E4b.int, main="Landscape Predictor Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")



#------ Model 3: landuse only, with variance structure & random structure ---------#

M3.final <- lme(data=coc2, result~landuse+rain, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
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


#------ Model 4.int ---------#

M4int.final <- lme(data=coc2, Form4.int, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4int.final <- residuals(object=M4int.final, type="normalized")

#plot residuals against agency, year, month, date and other unused predictors; check if model is missing something important!
boxplots.resids2(M4int.final, E4int.final, "X")
#is temporal auto-correlation present?  Look for patterns in Auto-correlation plot for residuals
par(mfrow=c(2,2))
plot(E4int.final~coc$start_date, pch=16, col="cadet blue")
acf(E4int.final, na.action=na.pass, main="Auto-correlation plot for residuals")  #look for lines extending past blue dashed range
#Check for spatial correlation
mydata2 <- data.frame(E4int.final, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
bubble(mydata2,"E4int.final",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
#variogram to test for spatial correlation
Vario1 = variogram(E4int.final ~ 1, mydata2)
plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
Vario2 <- variogram(E4int.final ~ 1, mydata2, alpha = c(0, 45, 90,135) )
plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?


#-----------------#
#  Save Formulas  #
#-----------------#

#save important items with copper-specific names
Cu.coc2 <- coc2
Cu.r1X <- r1X
Cu.vf1X <- vf1X
Cu.Form3 <- Form3
Cu.Form4 <- Form4b
Cu.Form4.int <- Form4b.int
Cu.rain <- rain

save(Cu.coc2, Cu.r1X, Cu.vf1X, Cu.Form3, Cu.Form4, Cu.Form4.int, file=here("results", "Copper Models.RData"))


#-------------------------------#
#  Plot Predictor Coefficients  #
#-------------------------------#

#generate plots that are shown in Rmarkdown script
Cu.null <- gls(data = Cu.coc2, result ~ 1, method = "REML") 
Cu.M3 <- lme(data = Cu.coc2, result ~ landuse + rain + summer, random = Cu.r1X, method = "REML", weights = Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
Cu.M4 <- lme(data = Cu.coc2, Cu.Form4, random = Cu.r1X, method = "REML", weights = Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
Cu.M4.int <- lme(data = Cu.coc2, Cu.Form4.int, random = Cu.r1X, method = "REML", weights = Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

Cu_models <- list(
  Null_Model = Cu.null,
  Categorical_Landuse_Model = Cu.M3,
#  Cu.M4 = Cu.M4,
  Landscape_Predictor_Model = Cu.M4.int)

huxtablereg(Cu_models,
            single.row = TRUE, custom.model.names = names(Cu_models)) %>%
  set_bottom_border(1, -1, 0.4) %>%
  set_bold(1, -1, TRUE) 

#plot parameter estimates for each model
theme_set(theme_sjplot())
plotreg(Cu_models, custom.title = "Regression Results, total Copper", custom.model.names = names(Cu_models))
plot_models(Cu_models,m.labels = names(Cu_models),legend.title = "Models", show.values = TRUE,show.intercept = TRUE)
#plot_models(Cu_models[-c(1)],m.labels = names(Cu_models[-c(1)]),legend.title = "Models", show.values = TRUE,show.intercept = TRUE)


#----------------------------------#
#  Plots to support Eva's summary  #
#----------------------------------#

# Model 3 has the lowest AIC; however, there are only 2 IND sites and 2 LDR sites; This model may not hold up when applied to other watersheds in Puget Sound area  
par(mfrow=c(2,1), mar=c(4,4,2,1))
landuseCol <- c("pink", "yellow", "gray", "green")
boxplot(coc2$result ~ coc2$land_use, ylab="ln(Copper)", xlab="", col=landuseCol)
boxplot(coc2$result ~ coc2$location, xlab="", ylab="ln(Copper)", col=landuseCol[c(1,2,4,1,1,1,2,3,1,2,4,1,2,3)], cex.axis=0.6)

#try using CO2_tot instead of sqrt_CO2_tot -- does it make an improvement that is worth considering?
coc3 <- cbind(coc2, CO2_tot=coc$CO2_tot)
Form4.alt <- formula(result ~ rain + summer + traffic + CO2_tot + totRES)
Model4.alt <- lme(data=coc3, Form4.alt, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4.alt <- residuals(object = Model4, type = "normalized")
#trying out CO2_tot instead of sqrt_CO2_tot does NO favors for the model and continues to not fit well
#  to CO2_tot.