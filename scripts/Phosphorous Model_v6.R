# This generalized script assesses linear relationships between total phosphorous and landscape predictors
# using mixed effects models

# Note: This version takes the final set of predictors (from Jan 29, 2021) and examines each set
#       with the COC of interest, to see which set of predictors may be most correlated with
#       the COC.
#       Also note that this version uses landuse percentages for landuse in Model 6, rather than 
#       landuse as a factor

#Four models are generated for the COC of interest (note the model numbers are 1,3,4,6 (no 2 or 5)):
#   1.) median value for all locations, irrespective of date
#   3.) COC value based only on self-reported land use (as a factor) + rain
#   4.) COC value based on up to 3 landscape parameters + rain; does NOT include land use
#   6.) COC value based on up to 3 landscape parameters + land use (as a continuous variable) + rain


# Eva Dusek Jennings
# Jun 18, 2021
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
predictors <- c("impervious_std", "paved_std", "roofs_std", "no2_std", "grass_std", "logTrees_std", 
                "sqrtNodev_std", "pm25_std", "sqrtSlope_std", "logTraffic_std",
                "roof_RES_std", "roof_nonRES_std", "logCOM_IND_std", "logCOM_std", "sqrtRES_std", "IND_std", "sqrtAG_std")

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

if(run_exploratory_code==TRUE) {
  #Cleveland Dotplot - use to detect violations of homogeneity: We are looking for whether the spread of data values
  #   differs between sampling locations (or agencies).  If so, it indicates heterogeneity, and that there may be problems with violation of
  #   homogeneity in a linear regression model applied on these data.  We're also looking for outliers.
  par(mfrow=c(1,1))
  dotchart(coc$result, groups=coc$location_id, pch=19, col=colors_agency[as.numeric(coc$agency)],
           xlab="concentration", main=paste("Cleveland Dotplot:", this_param_short))
  
  #pairs plots allow visualization of interactions between possible predictors.  Look for relationships between "result"
  #   and all of the predictor variables
  pairs(coc %>% select(result, predictors), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  #look also at interactions between result and rainfall
  pairs(coc %>% select(result, antecedant_dry_days_std, daymet_precip_std, daymet_3day_std,
                       daymet_7day_std, daymet_14day_std, daymet_21day_std, daymet_28daySR_std), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  #for phosphorus, important predictors are likely: 
  #  traffic (0.4), roof_nonRES (0.2), COM, slope
  #  daymet 28day (-0.2), daymet 21day (-0.2), daymet 14day (-0.2)
  
  #likely most important predictors
  pairs(coc[, c("result", "agency", "logTraffic_std", "roof_nonRES_std", "logCOM_IND_std", "sqrtSlope_std",
                "daymet_28daySR_std")], 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
}


#-----------------------------------------------------------------------------------------------------#
#  Plot COC vs various predictors: look for potential predictors to model & sources of heterogeneity  #
#-----------------------------------------------------------------------------------------------------#

#look for relationships between coc and some predictor variables
p1 <- ggplot(coc, aes(agency, result)) + geom_boxplot()
p2 <- ggplot(coc, aes(daymet_precip_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_precip_std))
p3 <- ggplot(coc, aes(daymet_3day_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_3day_std))
p4 <- ggplot(coc, aes(daymet_7day_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_7day_std))
p5 <- ggplot(coc, aes(daymet_14day_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_14day_std))
p6 <- ggplot(coc, aes(daymet_21day_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_21day_std))
p7 <- ggplot(coc, aes(daymet_28daySR_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$daymet_28daySR_std))
p8 <- ggplot(coc, aes(antecedant_dry_days_std, result)) + geom_point() + geom_smooth(method = "lm") + annotation_custom(slope_grob(coc$antecedant_dry_days_std))
# NOTE: I have left out precip & inches_rain_per_hour b/c they are missing ~68 data points

p9 <- my.ggplot(1)
p10 <- my.ggplot(2)
p11 <- my.ggplot(3)
p12 <- my.ggplot(4)
p13 <- my.ggplot(5)
p14 <- my.ggplot(6)
p15 <- my.ggplot(7)
p16 <- my.ggplot(8)
p17 <- my.ggplot(9)
p18 <- my.ggplot(10)
p24 <- my.ggplot(11)
p25 <- my.ggplot(12)
p26 <- my.ggplot(13)
p27 <- my.ggplot(14)
p28 <- my.ggplot(15)

p19 <- ggplot(coc, aes(as.factor(year), result)) + geom_boxplot()  #starts in Feb, 2009, ends in April, 2013; trends could be due to months with data?
p20 <- ggplot(coc, aes(month, result)) + geom_boxplot(position=position_dodge(width=0.9))
a <- coc %>% count(month)
p20 <- p20 + annotate(geom="text", x=c(1:12), y=rep(-1.4,12), label=paste("n=", a[,2], sep=""),
               color="red", size=3.5, angle=90)  #note the small sample size in June-Sept
p21 <- ggplot(coc, aes(as.factor(season), result)) + geom_boxplot()
p22 <- ggplot(coc, aes(land_use, result)) + geom_boxplot()
p23 <- ggplot(coc, aes(location_id, result)) + geom_boxplot()

if(run_exploratory_code==TRUE) {
  #Look at various transformations for monthly precip by location; which is most evenly spread out?
  par(mfrow=c(2,2), mar=c(4,4,4,2))
  plot(coc$mPrecip, coc$result)
  plot(coc$mPrecipSR, coc$result)
  plot(coc$mPrecipCR, coc$result)
  plot(coc$mPrecipLog, coc$result)
  #for phosphorus, mPrecip is better than the rest!
  
  #look for relationships between coc's and some predictors; 
  # look especially for sources of heterogeneity
  grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow=3)  #14-day, 21-day and 28-day precip all look good
  grid.arrange(p9, p10, p11, p12, p13, p14, p15, p16, p17, nrow=3)
  grid.arrange(p18, p24, p25, p26, p27, p28, nrow=2)
  grid.arrange(p1, p23, p19, p20, p21, p22, nrow=2)
  #for phosphorus, 14, 21 and 28-day precip don't show evidence of heteroskedasticity; 21-day is most evenly spread out
  #   location, agency, landuse, and month show heteroskedasticity
}


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
         paved=paved_std,
         impervious=impervious_std,
         roofs=roofs_std,
         grass=grass_std,
         trees=logTrees_std,
         traffic=logTraffic_std,
         pm25=pm25_std,
         slope=sqrtSlope_std,
         nodev=sqrtNodev_std,
         RES=sqrtRES_std,
         COM_IND=logCOM_IND_std,
         IND=IND_std,
         COM=logCOM_std,
         AG=sqrtAG_std,
         roof_RES=roof_RES_std,
         roof_nonRES=roof_nonRES_std,
         devAge=devAge2_std,
         dry=antecedant_dry_days_std,
         rain=daymet_21day_std,  #choose 21-day b/c its spread over time is very nice
         mPrecip=mPrecip  #note: if one of the other mPrecip measures looked better, use it here instead!
  ) %>%
  mutate(year=as.factor(year)) %>%
  mutate(landuse = case_when(
    land_use == "LDR" ~ "1.LDR",
    land_use == "HDR" ~ "2.HDR",
    land_use == "COM" ~ "3.COM",
    land_use == "IND" ~ "4.IND")) %>%
  mutate(landuse=as.factor(landuse)) %>%
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

FormX <- formula(result ~ paved + roof_RES + roof_nonRES + grass + trees + nodev + pm25 + traffic + slope + rain)
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
  M.gls3X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|location))
  M.gls4X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|season))
  M.gls5X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency * season))
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
  
  anova(M.gls1X, M.gls2X, M.gls3X,
        M.gls4X, M.gls5X, M.gls7X, M.gls8X,
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
vf1X <- varComb(varIdent(form= ~1|agency), varExp(form= ~dry))


# Steps 4-6: Find the proper random effects structure; look for temporal and spatial autocorrelation
#-----------  note: look for spatial correlation AFTER setting random effects!

if(run_exploratory_code==TRUE) {
  #----Model with Variance Function 1----#
  M.vf1X <- gls(FormX, data=coc2, method="REML", weights=vf1X)
  
  #Random intercept model; this is nested in the best fit model from Step 2, so can compare with a likelhiood ratio test
  M1.lme1X <- lme(data=coc2, FormX, random = ~1|agency, method="REML", weights=vf1X,
                  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  anova(M.vf1X, M1.lme1X)  #the likelihood ratio test p-value compares the two nested models, and whether they are significantly different
  #For phosphorus, the model with random intercept is significantly better than without!
  
  #Random slope models; with our dataset, these only makes sense for rain or dry.  Other predictors don't have enough data points along x-axis to generate a random slope reliably.
  M1.lme2X <- lme(data=coc2, FormX, random= ~1 + rain|agency, method="REML", weights=vf1X,
                  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  anova(M.vf1X, M1.lme1X, M1.lme2X)  #likelihood ratio test
  #for phosphorus, the random slope model has a singular convergence error
}
  
#random structure for best fit model so far
r1X <- formula(~1|agency)

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
  preds.1 <- c("impervious", "paved", "roofs", "roof_RES", "roof_nonRES", "trees", "grass", "nodev", "devAge", "pm25", "traffic", "slope")
  
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
  rem1 <- which(str_detect(dd.str, "roof_nonRES")==TRUE & str_detect(dd.str, "roofs")==TRUE)
  rem2 <- which(str_detect(dd.str, "roof_RES")==TRUE & str_detect(dd.str, "roofs")==TRUE)
  rem3 <- which(str_detect(dd.str, "roof_RES")==TRUE & str_detect(dd.str, "roof_nonRES")==TRUE)
  rem4 <- which(str_detect(dd.str, "paved")==TRUE & str_detect(dd.str, "impervious")==TRUE)
  rem5 <- which(str_detect(dd.str, "nodev")==TRUE & str_detect(dd.str, "devAge")==TRUE)
  rem6 <- which(str_detect(dd.str, "grass")==TRUE & str_detect(dd.str, "impervious")==TRUE)
  rem.all <- unique(c(rem1, rem2, rem3, rem4, rem5, rem6))  #keep only non-repetitive instances of indexes that should be removed
  
  dd <- dd[-rem.all, ]  #remove rows with predictors that don't go together!
  
  #generate an exhaustive set of formulas with one, two and three predictors, including and excluding rain
  ee.1 <- list(NULL)  #list of formulas based on the sets of predictors in "dd"; including rain
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i])) ee.1[[i]] <- as.formula(paste("result~rain + ", dd[i,1]))
    else if(is.na(dd$Var3[i]))  ee.1[[i]] <- as.formula(paste("result~rain + ", dd[i,1], "+", dd[i,2]))
    else  ee.1[[i]] <- as.formula(paste("result~rain + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  ee.2 <- list(NULL) #list of formulas based on the sets of predictors in "dd"; NOT including rain
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i]))  ee.2[[i]] <- as.formula(paste("result~ ", dd[i,1]))
    else if(is.na(dd$Var3[i]))  ee.2[[i]] <- as.formula(paste("result~ ", dd[i,1], "+", dd[i,2]))
    else ee.2[[i]] <- as.formula(paste("result~ ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  ee.3 <- list(NULL) #list of formulas based on the sets of predictors in "dd"; NOT including rain
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i]))  ee.3[[i]] <- as.formula(paste("result~dry + ", dd[i,1]))
    else if(is.na(dd$Var3[i]))  ee.3[[i]] <- as.formula(paste("result~dry + ", dd[i,1], "+", dd[i,2]))
    else ee.3[[i]] <- as.formula(paste("result~dry + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  ee.4 <- list(NULL) #list of formulas based on the sets of predictors in "dd"; NOT including rain
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i]))  ee.4[[i]] <- as.formula(paste("result~dry + rain + ", dd[i,1]))
    else if(is.na(dd$Var3[i]))  ee.4[[i]] <- as.formula(paste("result~dry + rain + ", dd[i,1], "+", dd[i,2]))
    else ee.4[[i]] <- as.formula(paste("result~dry + rain + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  
  ee <- c(ee.1, ee.2, ee.3, ee.4) 
  
  
  #run the lme for each formula, saving the AIC values; takes about 7 minutes to run
  ptm <- proc.time()  #time the code!
  my.aics <- rep(0, length(ee))
  for (i in 1:length(ee)) {
    bb <- lme(data=coc2, ee[[i]], random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
    my.aics[i] <- AIC(bb)
  }
  proc.time() - ptm  #stop the clock;  13 minutes for 940 models
  
  par(mfrow=c(3,5), mar=c(1,2,2,0), oma=c(1,1,1,1))
  gg <- rep("black", length(ee))
  gg[c(which(is.na(dd$Var3)), which(is.na(dd$Var3))+nrow(dd), which(is.na(dd$Var3))+2*nrow(dd), which(is.na(dd$Var3))+3*nrow(dd))] <- "turquoise"   #make formulas with only 2 landscape predictors red
  gg[c(which(is.na(dd$Var2)), which(is.na(dd$Var2))+nrow(dd), which(is.na(dd$Var2))+2*nrow(dd), which(is.na(dd$Var2))+3*nrow(dd) )] <- "pink"  #make formulas with only 1 landscape predictor yellow
  plot(my.aics, col=gg, pch=16, main="1 pred=pink, 2 preds=turquoise", xaxt="n", ylab="AIC")
  gg <- rep("black", length(ee))
  gg[which(str_detect(as.character(ee), "impervious"))] <- "red"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(ee), "paved"))] <- "orange"   #make formulas with "rain" light blue
  plot(my.aics, col=gg, xaxt="n", yaxt="n", pch=16, main="impervious=red, paved=orange")
  gg <- rep("black", length(ee))
  gg[which(str_detect(as.character(ee), "rain"))] <- "light blue"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(ee), "dry"))] <- "goldenrod"   #make formulas with "dry" goldenrod
  gg[which(str_detect(as.character(ee), "dry") & str_detect(as.character(ee), "rain"))] <- "green"   #make formulas with "dry" goldenrod
  plot(my.aics, col=gg, xaxt="n", yaxt="n", pch=16, main="rain=blue, dry=gold, both=green")
  #  plot(my.aics, col=ifelse(str_detect(as.character(ee), "rain"), "light blue", "black"), xaxt="n", yaxt="n", pch=16, main="rain=light blue")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "slope"), "tan3", "black"), xaxt="n", yaxt="n", pch=16, main="slope=light brown")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "trees"), "green3", "black"), xaxt="n", yaxt="n", pch=16, main="trees=green")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "grass"), "light green", "black"), xaxt="n", pch=16, main="grass=light green")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "nodev"), "purple1", "black"), xaxt="n", yaxt="n", pch=16, main="nodev=light purple")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "roofs"), "cadet blue", "black"), xaxt="n", yaxt="n", pch=16, main="roofs=blue")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "pm25"), "yellow", "black"), xaxt="n", yaxt="n", pch=16, main="pm25=yellow")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "traffic"), "gray", "black"), xaxt="n", yaxt="n", pch=16, main="traffic=gray")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "roof_nonRES"), "salmon", "black"), xaxt="n", pch=16, main="roof_nonRES=salmon")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "roof_RES"), "maroon1", "black"), xaxt="n", yaxt="n", pch=16, main="roof_RES=maroon")
  plot(my.aics, col=ifelse(str_detect(as.character(ee), "devAge"), "turquoise", "black"), xaxt="n", yaxt="n", pch=16, main="devAge=turquoise")
  #for phosphorus, best fitting predictors are: roof_RES, slope
  
  #arrange formulae and AICs according to order of AIC value, from smallest to largest
  my.formulas.m4 <- ee[order(my.aics)]
  my.formulas.m4[1:15]
  my.aics.m4 <- my.aics[order(my.aics)]
  my.aics.m4[1:15]  #what are the top 10 AIC's?
  
  save(my.formulas.m4, my.aics.m4, file=(here("scripts", "formulas_and_aics", "phosphorus_m4.RData")))

  #obtain formulas and AIC's where "rain" is the sole precip predictor, and only 1 or 2 predictors are used!
  rain.aics <- my.aics[1:72]
  rain.formulas <- ee[1:72]
  rain.formulas.ord <- rain.formulas[order(rain.aics)]
  rain.aics.ord <- rain.aics[order(rain.aics)]
  rain.formulas.ord[1:6]
  rain.aics.ord[1:6]
  
}

#NOTE: for phosphorus, adding "dry" improves fit slightly; check later to see if we should keep "dry" or remove (from variance fxn too)
Form4a <- formula(result~rain + dry + roof_RES + slope + pm25)  #AIC = 846.3
Form4b <- formula(result~rain + dry + roof_RES + slope + impervious)  #AIC = 852.1
Form4c <- formula(result~rain + dry + roof_RES + slope + grass)  #AIC = 853.2
Form4d <- formula(result~rain + dry + roof_RES + devAge + pm25)  #AIC = 854.5
Form4e <- formula(result~rain + dry + grass + pm25 + paved)  #AIC = 859.2
Form4f <- formula(result~rain + dry + grass + roofs + paved)  #AIC = 860.7
# NONE OF THESE FORMULAS ARE SUITABLE!  Perhaps too many predictors.  Try formulas with only 1-2 predictors

Form4i <- formula(result~rain + paved + grass)  #AIC=867.2   NOTE: grass:paved are highly correlated... I am skeptical...
Form4j <- formula(result~rain + roof_RES + slope)  #AIC=873.8
Form4k <- formula(result~rain + roofs + trees)  #AIC=874.1
Form4m <- formula(result~rain + trees + nodev)  #AIC=891.1

Form4x <- formula(result ~ rain + roof_RES)
Form4y <- formula(result ~ rain + trees)

#-----------------------------------------------------------------------#
#  Exhaustive Search for Model with Landuse + Best 3 Predictors + Rain  #
#-----------------------------------------------------------------------#

if(run_exploratory_code==TRUE) {
  
  #generate an exhaustive list of formulas based on the sets of predictors in "dd"; including landuse & rain
  ff.1 <- list(NULL)
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i]))  ff.1[[i]] <- as.formula(paste("result~rain + COM + ", dd[i,1]))
    else if(is.na(dd$Var3[i]))  ff.1[[i]] <- as.formula(paste("result~rain + COM + ", dd[i,1], "+", dd[i,2]))
    else  ff.1[[i]] <- as.formula(paste("result~rain + COM + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  ff.2 <- list(NULL)
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i])) ff.2[[i]] <- as.formula(paste("result~rain + RES + ", dd[i,1]))
    else if(is.na(dd$Var3[i])) ff.2[[i]] <- as.formula(paste("result~rain + RES + ", dd[i,1], "+", dd[i,2]))
    else ff.2[[i]] <- as.formula(paste("result~rain + RES + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  ff.3 <- list(NULL)
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i])) ff.3[[i]] <- as.formula(paste("result~rain + IND + ", dd[i,1]))
    else if(is.na(dd$Var3[i])) ff.3[[i]] <- as.formula(paste("result~rain + IND + ", dd[i,1], "+", dd[i,2]))
    else ff.3[[i]] <- as.formula(paste("result~rain + IND + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  ff.4 <- list(NULL)
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i])) ff.4[[i]] <- as.formula(paste("result~rain + AG + ", dd[i,1]))
    else if(is.na(dd$Var3[i])) ff.4[[i]] <- as.formula(paste("result~rain + AG + ", dd[i,1], "+", dd[i,2]))
    else ff.4[[i]] <- as.formula(paste("result~rain + AG + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  ff <- c(ff.1, ff.2, ff.3, ff.4)
  
  #run the lme for each formula, saving the AIC values; takes about 7 minutes to run
  ptm <- proc.time()  #time the code!
  my.aics.ff <- rep(0, length(ff))
  for (i in 1:length(ff)) {
    bb <- lme(data=coc2, ff[[i]], method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
    my.aics.ff[i] <- AIC(bb)
  }
  proc.time() - ptm  #stop the clock; takes about 13 minutes
  
  #combine the formulas and AIC values for both landscape and non-landscape formulas
  jj <- ff #c(ee, ff)
  my.aics.jj <- my.aics.ff  #c(my.aics, my.aics.ff)
  
  par(mfrow=c(3,5), mar=c(1,2,2,0), oma=c(1,1,1,1))
  gg <- rep("black", length(jj))
  gg[c(which(is.na(dd$Var3)), which(is.na(dd$Var3))+nrow(dd), which(is.na(dd$Var3)) + 2*nrow(dd), which(is.na(dd$Var3)) + 3*nrow(dd) )] <- "turquoise"   #make formulas with only 2 landscape predictors turquoise
  gg[c(which(is.na(dd$Var2)), which(is.na(dd$Var2))+nrow(dd), which(is.na(dd$Var2)) + 2*nrow(dd), which(is.na(dd$Var2)) + 3*nrow(dd) )] <- "pink"  #make formulas with only 1 landscape predictor pink
  plot(my.aics.jj, col=gg, pch=16, main="1 pred=pink, 2 preds=turquoise", xaxt="n", ylab="AIC")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "rain"), "light blue", "black"), xaxt="n", pch=16, main="rain=light blue")
  gg <- rep("black", length(jj))
  gg[which(str_detect(as.character(jj), " IND "))] <- "gray"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(jj), " COM "))] <- "pink"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(jj), " RES "))] <- "green"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(jj), " AG "))] <- "goldenrod"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(jj), " RES ") & str_detect(as.character(jj), " COM "))] <- "magenta"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(jj), " RES ") & str_detect(as.character(jj), " IND "))] <- "cadet blue"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(jj), " COM ") & str_detect(as.character(jj), " IND "))] <- "tan3"   #make formulas with "rain" light blue
  plot(my.aics.jj, col=gg, xaxt="n", yaxt="n", pch=16, main="COM, RES, IND, AG")
  # abline(h=c(AIC(M1.LUr), AIC(M1.LU)), col="magenta", lty=2, lwd=2)
  # text(x=300, y=AIC(M1.LUr)+10, label="AIC for landuse w/random structure", col="magenta")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "slope"), "tan3", "black"), xaxt="n", yaxt="n", pch=16, main="slope=light brown")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "trees"), "green3", "black"), xaxt="n", yaxt="n", pch=16, main="trees=green")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "grass"), "light green", "black"), xaxt="n", pch=16, main="grass=light green")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "nodev"), "purple1", "black"), xaxt="n", yaxt="n", pch=16, main="nodev=light purple")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "roofs"), "cadet blue", "black"), xaxt="n", yaxt="n", pch=16, main="roofs=blue")
  gg <- rep("black", length(jj))
  gg[which(str_detect(as.character(jj), "impervious"))] <- "red"   #make formulas with "rain" light blue
  gg[which(str_detect(as.character(jj), "paved"))] <- "orange"   #make formulas with "rain" light blue
  plot(my.aics.jj, col=gg, xaxt="n", yaxt="n", pch=16, main="impervious=red, paved=orange")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "pm25"), "yellow", "black"), xaxt="n", yaxt="n", pch=16, main="pm25=yellow")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "traffic"), "gray", "black"), xaxt="n", pch=16, main="traffic=gray")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "roof_nonRES"), "salmon", "black"), xaxt="n", yaxt="n", pch=16, main="roof_nonRES=salmon")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "roof_RES"), "maroon1", "black"), xaxt="n", yaxt="n", pch=16, main="roof_RES=maroon")
  plot(my.aics.jj, col=ifelse(str_detect(as.character(jj), "devAge"), "turquoise", "black"), xaxt="n", yaxt="n", pch=16, main="devAge=turquoise")
  #for phosphorus, looks like AG is a slightly better predictor, with COM, RES and IND tied
  
  #arrange formulae and AICs according to order of AIC value, from smallest to largest
  my.formulas.m5 <- jj[order(my.aics.jj)]
  my.formulas.m5[1:10]
  my.aics.m5 <- my.aics.jj[order(my.aics.jj)]
  my.aics.m5[1:10]
  
  save(my.formulas.m5, my.aics.m5, file=(here("scripts", "formulas_and_aics", "phosphorus_m5.RData")))
  
  #obtain formulas with only one and two predictors + landuse percentage
  two.formulas.m5 <- jj[c(1:72, 235+(1:72), 235*2+(1:72), 235*3+(1:72))]
  two.aics.m5 <- my.aics.jj[c(1:72, 235+(1:72), 235*2+(1:72), 235*3+(1:72))]
  two.formulas.m5.ord <- two.formulas.m5[order(two.aics.m5)]
  two.aics.m5.ord <- two.aics.m5[order(two.aics.m5)]
  two.formulas.m5.ord[1:10]
  two.aics.m5.ord[1:10]
  
  one.aics.m5 <- my.aics.jj[c(1:12, 235+(1:12), 235*2+(1:12), 235*3+(1:12))]
  one.formulas.m5 <- jj[c(1:12, 235+(1:12), 235*2+(1:12), 235*3+(1:12))]
  one.formulas.m5.ord <- one.formulas.m5[order(one.aics.m5)]
  one.aics.m5.ord <- one.aics.m5[order(one.aics.m5)]
  one.formulas.m5.ord[1:10]
  one.aics.m5.ord[1:10]
  
}

#NOTE: these formulas don't include dry!  include it below, to test if it is worth including (also as a variance covariate)
# Form5a <- formula(result ~ rain + AG + roofs + grass + paved)  #AIC = 846.8
# Form5b <- formula(result ~ rain + AG + pm25 + slope + roof_RES)  #AIC = 850.8
# Form5c <- formula(result ~ rain + RES + grass + nodev + paved)  #AIC = 854.5  #this is #7
# Form5d <- formula(result ~ rain + AG + trees + grass + paved)  #AIC = 855.7   #this is #9
### based on my experience with Model4 and with plotting Model5a, there are WAY too many predictors in these models.
#   lets try again with just 2 or 1 predictors...

Form5a <- formula(result ~ rain + AG + paved + grass)  #AIC=857.3; can swap out AG for IND, COM or RES and still get good AIC
Form5b <- formula(result ~ rain + RES + roof_RES + slope)  #AIC=861; can swap out RES for any other land use percentage
Form5c <- formula(result ~ rain + AG + roof_RES + pm25)  #AIC=872.5
Form5d <- formula(result ~ rain + COM + roofs + trees)  #AIC=873.4

Form5i <- formula(result ~ rain + AG + roofs)  #AIC=884.9
Form5j <- formula(result ~ rain + AG + roof_RES)  #AIC=891.3
Form5k <- formula(result ~ rain + AG + nodev)  #AIC=897.9
Form5l <- formula(result ~ rain + RES + trees)  #AIC=898.9
Form5m <- formula(result ~ rain + RES + paved)  #AIC=900.9

#--------------------------------------------#
#  Quantiles List for Plotting Interactions  #
#--------------------------------------------#

#list of quantiles for all landscape predictors
qList <- list(paved=quantile(coc2$paved, probs=c(0,.25,.50,.75,1)),
              impervious=quantile(coc2$impervious, probs=c(0,.25,.50,.75,1)),
              roofs=quantile(coc2$roofs, probs=c(0,.25,.50,.75,1)),
              roof_RES=quantile(coc2$roof_RES, probs=c(0,.25,.50,.75,1)),
              roof_nonRES=quantile(coc2$roof_nonRES, probs=c(0,.25,.50,.75,1)),
              grass=quantile(coc2$grass, probs=c(0,.25,.50,.75,1)),
              trees=quantile(coc2$trees, probs=c(0,.25,.50,.75,1)),
              traffic=quantile(coc2$traffic, probs=c(0,.25,.50,.75,1)),
              pm25=quantile(coc2$pm25, probs=c(0,.25,.50,.75,1)),
              slope=quantile(coc2$slope, probs=c(0,.25,.50,.75,1)),
              nodev=quantile(coc2$nodev, probs=c(0,.25,.50,.75,1)),
              devAge=quantile(coc2$devAge, probs=c(0,.25,.50,.75,1)),
              dry=quantile(coc2$dry, probs=c(0,.25,.50,.75,1)),
              rain=quantile(coc2$rain, probs=c(0,.25,.50,.75,1)),
              RES=quantile(coc2$RES, probs=c(0,.25,.50,.75,1)),
              COM_IND=quantile(coc2$COM_IND, probs=c(0,.25,.50,.75,1)),
              AG=quantile(coc2$AG, probs=c(0,.25,.50,.75,1)),
              COM=quantile(coc2$COM, probs=c(0,.25,.50,.75,1)),
              IND=quantile(coc2$IND, probs=c(0,.25,.50,.75,1)) )


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

Form3 <- formula(result ~ landuse + rain + summer)
Model3 <- lme(data=coc2, Form3, random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E3 <- residuals(object=Model3, type="normalized")
boxplots.resids2(Model3, E3, "X")

#------ Model 4: no land use; up to 3 predictors + rain ---------#

#NOTE: for phosphorus, adding "dry" improves fit, but only slightly (~4 AIC pts);
#   fit was made worse based on BIC (~4 BIC pts).  As a result, I have removed "dry" from all of these equations.
#   AIC values listed here are for the equations when the INCLUDE "dry"
Form4a <- formula(result~rain + roof_RES + slope + pm25)  #AIC = 846.3
Form4b <- formula(result~rain + roof_RES + slope + impervious)  #AIC = 852.1
Form4c <- formula(result~rain + roof_RES + slope + grass)  #AIC = 853.2
Form4d <- formula(result~rain + roof_RES + devAge + pm25)  #AIC = 854.5
Form4e <- formula(result~rain + grass + pm25 + paved)  #AIC = 859.2
Form4f <- formula(result~rain + grass + roofs + paved)  #AIC = 860.7

# #####  Model4a  #####
# Model4a <- lme(data=coc2, Form4a, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model4a)
# 
# # #Test whether "dry" is a worthwhile covariate to include
# # Form4a.1 <- update(Form4a, .~. -dry)
# # Model4a.1 <- lme(data=coc2, Form4a.1, method="ML", random=r1X,
# #                weights=varIdent(form= ~1|agency),
# #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# # summary(Model4a.1)
# # 
# # AIC(Model4a, Model4a.1)  #AIC without "dry" is 4 points higher than WITH "dry"
# # BIC(Model4a, Model4a.1)  #BIC without "dry" is 4 points lower than WITH "dry" -- its a toss-up.  Exclude "dry" to simplify! 
# 
# #check to see whether these predictors are highly correlated or not
# z4a <- coc2 %>%
#   dplyr::select(c(slope, pm25, roof_RES))
# pairs(z4a, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #highest correlation is between pm25:roof_RES (0.4)
# 
# ppp1 <- Plot.One.Predictor.wAgency("slope", "agency", Model4a, yEqn = 7.5, xEqn=0.5)
# ppp2 <- Plot.One.Predictor.wAgency("pm25", "agency", Model4a, yEqn = 7.5, xEqn=-0.5)
# ppp3 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model4a, yEqn = 7.5, xEqn=0.2)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2)
# #these plots look terrible!!
# 
# E4a <- residuals(Model4a, type="normalized")
# boxplots.resids2(Model4a, E4a, "X")
# #fits look good, though!  Don't use this model with such terrible plots!

# #####  Model4b  #####
# Model4b <- lme(data=coc2, Form4b, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model4b)
# 
# #check to see whether these predictors are highly correlated or not
# z4b <- coc2 %>%
#   dplyr::select(c(slope, impervious, roof_RES))
# pairs(z4b, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #highest correlation is impervious:roof_RES (0.6)
# 
# ppp1 <- Plot.One.Predictor.wAgency("slope", "agency", Model4b, yEqn = 7.5, xEqn=-0.5)  
# ppp2 <- Plot.One.Predictor.wAgency("impervious", "agency", Model4b, yEqn = 7.5, xEqn=-0.5)  
# ppp3 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model4b, yEqn = 7.5, xEqn=0.5)  
# grid.arrange(ppp1, ppp2, ppp3, nrow=2)  
# #fit for slope is ridiculous (negative relationship)
# 
# E4b <- residuals(Model4b, type="normalized")
# boxplots.resids2(Model4b, E4b, "X")


# #####  Model4c  #####
# Model4c <- lme(data=coc2, Form4c, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# 
# #check to see whether these predictors are highly correlated or not
# z4c <- coc2 %>%
#   dplyr::select(c(slope, grass, roof_RES))
# pairs(z4c, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #highest correlation is grass:roof_RES (0.5)
# 
# ppp1 <- Plot.One.Predictor.wAgency("slope", "agency", Model4c, yEqn = 7, xEqn=-0.5)
# ppp2 <- Plot.One.Predictor.wAgency("grass", "agency", Model4c, yEqn = 7, xEqn=0.5)
# ppp3 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model4c, yEqn = 7, xEqn=0.5)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2)
# #fit to slope data is really questionable, in my opinion.  Don't use this model!
# 
# E4c <- residuals(Model4c, type="normalized")
# boxplots.resids2(Model4c, E4c, "X")

# #####  Model4d  #####
# Model4d <- lme(data=coc2, Form4d, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# 
# #check to see whether these predictors are highly correlated or not
# z4d <- coc2 %>%
#   dplyr::select(c(devAge, pm25, roof_RES))
# pairs(z4d, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #correlation between devAge:pm25 = 0.5
# 
# ppp1 <- Plot.One.Predictor.wAgency("devAge", "agency", Model4d, yEqn = 7.5, xEqn=-0.5)
# ppp2 <- Plot.One.Predictor.wAgency("pm25", "agency", Model4d, yEqn = 7.5, xEqn=0.5)
# ppp3 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model4d, yEqn = 7.5, xEqn=0.5)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2)
# # pm25 fit is terrible!  Don't use this model!
# 
# E4d <- residuals(Model4d, type="normalized")
# boxplots.resids2(Model4d, E4d, "X")


# #####  Model4e  #####
# Model4e <- lme(data=coc2, Form4e, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model4e)
# 
# #check to see whether these predictors are highly correlated or not
# z4e <- coc2 %>%
#   dplyr::select(c(grass, pm25, paved))
# pairs(z4e, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #correlation coefficient is highest for grass:paved (-0.8)
# 
# ppp1 <- Plot.One.Predictor.wAgency("grass", "agency", Model4e, yEqn = 7.5, xEqn=-0.5)
# ppp2 <- Plot.One.Predictor.wAgency("pm25", "agency", Model4e, yEqn = 7.5, xEqn=0.5)
# ppp3 <- Plot.One.Predictor.wAgency("paved", "agency", Model4e, yEqn = 7.5, xEqn=0.5)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2, ncol=2)
# #fits are generally awful - don't use this model
# 
# E4e <- residuals(Model4e, type="normalized")
# boxplots.resids2(Model4e, E4e, "X")


# #####  Model4f  #####
# Model4f <- lme(data=coc2, Form4f, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model4f)
# 
# #check to see whether these predictors are highly correlated or not
# z4f <- coc2 %>%
#   dplyr::select(c(grass, roofs, paved))
# pairs(z4f, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #correlation coefficient is highest for grass:paved (-0.8)
# 
# ppp1 <- Plot.One.Predictor.wAgency("grass", "agency", Model4f, yEqn = 7.5, xEqn=-0.5)
# ppp2 <- Plot.One.Predictor.wAgency("roofs", "agency", Model4f, yEqn = 7.5, xEqn=0.5)
# ppp3 <- Plot.One.Predictor.wAgency("paved", "agency", Model4f, yEqn = 7.5, xEqn=0.5)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2, ncol=2)
# #fits are generally awful - don't use this model




Form4i <- formula(result~rain + paved + grass)  #AIC=867.2   NOTE: grass:paved are highly correlated... I am skeptical...
Form4j <- formula(result~rain + roof_RES + slope)  #AIC=873.8
Form4k <- formula(result~rain + roofs + trees)  #AIC=874.1
#Form4l <- formula(result~rain + trees + pm25)  #AIC=881.6  TERRIBLE FIT FOR pm25!!
Form4m <- formula(result~rain + roof_RES + devAge)  #AIC=891.1  this model is OK - can't quite explain it, but it is ok.
Form4n <- formula(result~rain + trees + nodev)  #AIC=891.1  this model is OK - can't quite explain it, but it is ok.

Form4x <- formula(result ~ rain + roof_RES)
Form4y <- formula(result ~ rain + trees)


#####  Model4i  #####
Model4i <- lme(data=coc2, Form4i, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(Model4i)

#check to see whether these predictors are highly correlated or not
z4i <- coc2 %>%
  dplyr::select(c(grass, paved))
pairs(z4i, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
#correlation coefficient is highest for grass:paved (-0.8)

ppp1 <- Plot.One.Predictor.wAgency("grass", "agency", Model4i, yEqn = 7.5, xEqn=-0.5)
ppp2 <- Plot.One.Predictor.wAgency("paved", "agency", Model4i, yEqn = 7.5, xEqn=0.5)
grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
#fits are pretty good, except for city of Seattle...

E4i <- residuals(Model4i, type="normalized")
boxplots.resids2(Model4i, E4i, "X")

#try adding SUMMER as a predictor (categorical, 0 for all seasons except summer, which is 1)
Model4i.1 <- update(Model4i, .~. +summer)
summary(Model4i.1)
ppp1 <- Plot.One.Predictor.wAgency("grass", "agency", Model4i.1, yEqn = 7.5, xEqn=-0.5)
ppp2 <- Plot.One.Predictor.wAgency("paved", "agency", Model4i.1, yEqn = 7.5, xEqn=0.5)
grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
E4i.1 <- residuals(Model4i.1, type="normalized")
boxplots.resids2(Model4i.1, E4i.1, "X")

AIC(Model4i, Model4i.1)  #ADDING SUMMER AS A PREDICTOR MAKES A HUGE DIFFERENCE!  (AIC goes down by 32 pts)


# #try swapping out "traffic" for "paved"
# Model4i.2 <- update(Model4i.1, .~. -paved +traffic)
# summary(Model4i.2)
# ppp1 <- Plot.One.Predictor.wAgency("grass", "agency", Model4i.2, yEqn = 7.5, xEqn=-0.5)
# ppp2 <- Plot.One.Predictor.wAgency("traffic", "agency", Model4i.2, yEqn = 7.5, xEqn=0.5)
# grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
# E4i.2 <- residuals(Model4i.2, type="normalized")
# boxplots.resids2(Model4i.2, E4i.2, "X")
# #this model is terrible!  The fit to grass is almost completely flat!


# #####  Model4j  #####
# Model4j <- lme(data=coc2, 
#                update(Form4j, .~. +summer), 
#                method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model4j)
# 
# #check to see whether these predictors are highly correlated or not
# z4j <- coc2 %>%
#   dplyr::select(c(roof_RES, slope))
# pairs(z4j, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #correlation coefficient is 0.1
# 
# ppp1 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model4j, yEqn = 7.5, xEqn=-0.5)
# ppp2 <- Plot.One.Predictor.wAgency("slope", "agency", Model4j, yEqn = 7.5, xEqn=0.5)
# grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
# #I don't like the fit for slope; don't use this model!
# 
# E4j <- residuals(Model4j, type="normalized")
# boxplots.resids2(Model4j, E4j, "X")


#####  Model4k  #####
Model4k <- lme(data=coc2, 
               update(Form4k, .~. +summer), 
               method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(Model4k)

#check to see whether these predictors are highly correlated or not
z4k <- coc2 %>%
  dplyr::select(c(roofs, trees))
pairs(z4k, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
#correlation coefficient is -0.2

ppp1 <- Plot.One.Predictor.wAgency("roofs", "agency", Model4k, yEqn = 7.5, xEqn=-0.5)
ppp2 <- Plot.One.Predictor.wAgency("trees", "agency", Model4k, yEqn = 7.5, xEqn=0.5)
grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
#I like the fit for trees; not totally sure about roofs, but its better than many others.  Considering different intercept for
#  each agency, the roofs fit makes more sense.

E4k <- residuals(Model4k, type="normalized")
boxplots.resids2(Model4k, E4k, "X")


# #####  Model4m  #####    #NOTE: I tried swapping out devAge for nodev, but fitting to nodev was tough
# Model4m <- lme(data=coc2, 
#                update(Form4m, .~. +summer), 
#                method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model4m)
# 
# #check to see whether these predictors are highly correlated or not
# z4m <- coc2 %>%
#   dplyr::select(c(roof_RES, devAge))
# pairs(z4m, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #correlation coefficient is -0.3
# 
# ppp1 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model4m, yEqn = 7.5, xEqn=-0.5)
# ppp2 <- Plot.One.Predictor.wAgency("devAge", "agency", Model4m, yEqn = 7.5, xEqn=0.5)
# grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
# #fit to roof_RES is good; fit to devAge is okay.  Not entirely sure how to explain this all... its a little
# #  convoluted with roof_RES AND devAge as the predictors here...
# 
# E4m <- residuals(Model4m, type="normalized")
# boxplots.resids2(Model4m, E4m, "X")


# #####  Model4n  #####
# Model4n <- lme(data=coc2, 
#                update(Form4n, .~. +summer), 
#                method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model4n)
# 
# #check to see whether these predictors are highly correlated or not
# z4n <- coc2 %>%
#   dplyr::select(c(trees, nodev))
# pairs(z4n, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #correlation coefficient is 0.5
# 
# ppp1 <- Plot.One.Predictor.wAgency("trees", "agency", Model4n, yEqn = 7.5, xEqn=-0.5)
# ppp2 <- Plot.One.Predictor.wAgency("nodev", "agency", Model4n, yEqn = 7.5, xEqn=0.5)
# grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
# #I like the fit for trees; doesn't entirely make sense that nodev increases with phosphorus...
# 
# E4n <- residuals(Model4n, type="normalized")
# boxplots.resids2(Model4n, E4n, "X")


#####  Best fit Model4  ####

Form4 <- update(Form4i, .~. +summer)
Model4 <- Model4i.1

#try adding interactions and see if they are significant
M4.sub1 <- update(Model4, . ~ . + rain:paved) #not significant (p=0.9962)
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, . ~ . + rain:grass)  #not significant (p=0.6869)
anova(Model4, M4.sub2)
M4.sub4 <- update(Model4, . ~ . + grass:paved)  #significant (p=0.334)
anova(Model4, M4.sub4)
#no significant interactions to include

E4 <- residuals(Model4, type="normalized")
boxplots.resids2(Model4, E4, "X")


#------ Model 5: up to 3 predictors + landuse + rain ---------#

#NOTE: these formulas have way too many predictors, and hte "best fit" lines 
# Form5a <- formula(result ~ rain + AG + roofs + grass + paved)  #AIC = 846.8 -- WOW -- HORRID!  Fit lines make NO sense!
# Form5b <- formula(result ~ rain + AG + pm25 + slope + roof_RES)  #AIC = 850.8
# Form5c <- formula(result ~ rain + RES + grass + nodev + paved)  #AIC = 854.5  #this is #7
# Form5d <- formula(result ~ rain + AG + trees + grass + paved)  #AIC = 855.7   #this is #9

#NOTE: AIC's listed here are for formulas NOT including summer.  I've added it here to simplify...
Form5a <- formula(result ~ rain + AG + paved + grass + summer)  #AIC=857.3; can swap out AG for IND, COM or RES and still get good AIC
Form5b <- formula(result ~ rain + RES + roof_RES + slope + summer)  #AIC=861; can swap out RES for any other land use percentage
Form5c <- formula(result ~ rain + AG + roof_RES + pm25 + summer)  #AIC=872.5
Form5d <- formula(result ~ rain + COM + roofs + trees + summer)  #AIC=873.4

Form5i <- formula(result ~ rain + AG + roofs + summer)  #AIC=884.9
Form5j <- formula(result ~ rain + AG + roof_RES + summer)  #AIC=891.3
Form5k <- formula(result ~ rain + AG + nodev + summer)  #AIC=897.9
Form5l <- formula(result ~ rain + RES + trees + summer)  #AIC=898.9
Form5m <- formula(result ~ rain + RES + paved + summer)  #AIC=900.9


# #####  Model5a  #####
# Model5a <- lme(data=coc2, Form5a, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model5a)
# 
# #check to see whether these predictors are highly correlated or not
# z5a <- coc2 %>%
#   dplyr::select(c(AG, grass, paved))
# pairs(z5a, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #highest correlation coefficient AG:paved AND grass:paved (both -0.8) -- too high!
# 
# ppp1 <- Plot.One.Predictor.wAgency("AG", "agency", Model5a, yEqn = 7, xEqn=0.2)
# ppp2 <- Plot.One.Predictor.wAgency("grass", "agency", Model5a, yEqn = 7, xEqn=-0.5)
# ppp3 <- Plot.One.Predictor.wAgency("paved", "agency", Model5a, yEqn = 7, xEqn=0.2)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2)
# #relationship with AG is really questionable in my mind.  Don't use this relationship
# 
# E5a <- residuals(Model5a, type="normalized")
# boxplots.resids2(Model5a, E5a, "X")


# #### Model5b ####
# Model5b <- lme(data=coc2, Form5b, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model5b)
# 
# #check to see whether these predictors are highly correlated or not
# z5b <- coc2 %>%
#   dplyr::select(c(RES, slope, roof_RES))
# pairs(z5b, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #max correlation is 0.4 (RES:slope)
# 
# ppp1 <- Plot.One.Predictor.wAgency("RES", "agency", Model5b, yEqn = 7.5, xEqn=0.2)
# ppp2 <- Plot.One.Predictor.wAgency("slope", "agency", Model5b, yEqn = 7.5, xEqn=-0.2)
# ppp3 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model5b, yEqn = 7.5, xEqn=0.2)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2)
# #fit to slope is dreadful.  Don't use this model
# 
# E5b <- residuals(Model5b, type="normalized")
# boxplots.resids2(Model5b, E5b, "X")


# #### Model5c ####
# Model5c <- lme(data=coc2, Form5c, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model5c)
# 
# #check to see whether these predictors are highly correlated or not
# z5c <- coc2 %>%
#   dplyr::select(c(AG, pm25, roof_RES))
# pairs(z5c, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #max correlation is 0.7 (trees:roof_RES)
# 
# ppp1 <- Plot.One.Predictor.wAgency("AG", "agency", Model5c, yEqn = 7.5, xEqn=0.2)
# ppp2 <- Plot.One.Predictor.wAgency("pm25", "agency", Model5c, yEqn = 7.5, xEqn=0.5)
# ppp3 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model5c, yEqn = 7.5, xEqn=0.2)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2)
# #fit to pm25 is dreadful.  Eliminate this model!
# 
# E5c <- residuals(Model5c, type="normalized")
# boxplots.resids2(Model5c, E5c, "X")


# #### Model5d ####
# Model5d <- lme(data=coc2, Form5d, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model5d)
# 
# #check to see whether these predictors are highly correlated or not
# z5d <- coc2 %>%
#   dplyr::select(c(COM, roofs, trees))
# pairs(z5d, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
# #max correlation is -0.3 (COM:trees)
# 
# ppp1 <- Plot.One.Predictor.wAgency("COM", "agency", Model5d, yEqn = 7.5, xEqn=0.2)
# ppp2 <- Plot.One.Predictor.wAgency("roofs", "agency", Model5d, yEqn = 7.5, xEqn=-0.2)
# ppp3 <- Plot.One.Predictor.wAgency("trees", "agency", Model5d, yEqn = 7.5, xEqn=0.2)
# grid.arrange(ppp1, ppp2, ppp3, nrow=2)
# #fit to COM seems wrong; I also have some hesitance about the fit to "roofs".  This is not my favorite model...
# 
# E5d <- residuals(Model5d, type="normalized")
# boxplots.resids2(Model5d, E5d, "X")


#### Model5i ####
Model5i <- lme(data=coc2, Form5i, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(Model5i)

#check to see whether these predictors are highly correlated or not
z5i <- coc2 %>%
  dplyr::select(c(AG, roofs))
pairs(z5i, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
#highest correlation is -0.1

ppp1 <- Plot.One.Predictor.wAgency("AG", "agency", Model5i, yEqn = 7.5, xEqn=0.2)
ppp2 <- Plot.One.Predictor.wAgency("roofs", "agency", Model5i, yEqn = 7.5, xEqn=0.5)
grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
#fit looks good, but the mechanism doesn't make sense.  Why should AG, which uses fertilizer, have low phosphorus??

E5i <- residuals(Model5i, type="normalized")
boxplots.resids2(Model5i, E5i, "X")
#this fit is about as good as it gets for phosphorus!  


#### Model5j ####
Model5j <- lme(data=coc2, Form5j, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(Model5j)

#check to see whether these predictors are highly correlated or not
z5j <- coc2 %>%
  dplyr::select(c(AG, roof_RES))
pairs(z5j, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
#highest correlation is 0.5

ppp1 <- Plot.One.Predictor.wAgency("AG", "agency", Model5j, yEqn = 7.5, xEqn=0.2)
ppp2 <- Plot.One.Predictor.wAgency("roof_RES", "agency", Model5j, yEqn = 7.5, xEqn=0.5)
grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
#fit looks good, but the mechanism doesn't make sense.  Why should AG, which uses fertilizer, have low phosphorus??

E5j <- residuals(Model5j, type="normalized")
boxplots.resids2(Model5j, E5j, "X")


#### Model5l ####
Model5l <- lme(data=coc2, Form5l, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(Model5l)

#check to see whether these predictors are highly correlated or not
z5l <- coc2 %>%
  dplyr::select(c(RES, trees))
pairs(z5l, lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
#highest correlation is 0.4

ppp1 <- Plot.One.Predictor.wAgency("RES", "agency", Model5l, yEqn = 7.5, xEqn=0.2)
ppp2 <- Plot.One.Predictor.wAgency("trees", "agency", Model5l, yEqn = 7.5, xEqn=0.5)
grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
#fit to RES is only meh...

E5l <- residuals(Model5l, type="normalized")
boxplots.resids2(Model5l, E5l, "X")


###### models based on literature #####
# Model5x <- lme(data=coc2, 
#                result~rain + summer + grass*RES + paved*COM_IND, 
#                method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# summary(Model5x)
# 
# ppp1 <- Plot.One.Predictor.wAgency("RES", "agency", Model5x, yEqn = 7.5, xEqn=0.2)
# ppp2 <- Plot.One.Predictor.wAgency("grass", "agency", Model5x, yEqn = 7.5, xEqn=0.5)
# ppp3 <- Plot.One.Predictor.wAgency("COM_IND", "agency", Model5x, yEqn = 7.5, xEqn=0.2)
# ppp4 <- Plot.One.Predictor.wAgency("paved", "agency", Model5x, yEqn = 7.5, xEqn=0.5)
# grid.arrange(ppp1, ppp2, ppp3, ppp4, nrow=2)
# 
# Plot.Quantile("grass", "RES", Model5x, yEqn=7.5)
# Plot.Quantile("paved", "COM_IND", Model5x, yEqn=7.5)
# #hmmmmm... These interaction plots are not convincing in any sort of way.  Scratch this model.


#### BEST FIT MODEL 5
Form5 <- Form5i
Model5 <- Model5i

#try adding interactions and see if they are significant
M5.sub1 <- update(Model5, . ~ . + AG:rain)   #not significant (p=0.6573)
anova(Model5, M5.sub1)
M5.sub2 <- update(Model5, . ~ . + roofs:rain)  #not significant (p=0.8641)
anova(Model5, M5.sub2)
M5.sub5 <- update(Model5, . ~ . + AG:roofs)  #slightly significant (p<0.0001)
anova(Model5, M5.sub5)

#Plot formula 5 interaction
Plot.Quantile("AG", "roofs", M5.sub5, yEqn=7.5)
#I have NO idea how to explain this interaction... The data is not obviously compelling, but not uncompelling either.

E5 <- residuals(Model5, type="normalized")
boxplots.resids2(Model5, E5, "X")



#------ Compare models and plot ---------#
AIC(Model1, Model3, Model4, Model4k, Model5)

AIC(Model1, Model3, Model4, Model4k, Model5)[,2] - AIC(Model4)

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
plot(coc2$location, E1, main="Model 1", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
plot(coc2$location, E3, main="Model 3", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
plot(coc2$location, E4, main="Model 4i", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
# plot(coc2$location, E4e, main="Model 4e", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
# axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
# abline(h=0, col="gray")
plot(coc2$location, E4k, main="Model 4k", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")


# Step 9: Refit models with REML and apply graphical model validation; check for homogeneity, 
#--------    normality and independence

#------ Model 3: landuse only, with variance structure & random structure ---------#

M3.final <- lme(data=coc2, result~landuse+rain+summer, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E3.final <- residuals(object=M3.final, type="normalized")

#plot best model residuals - check for goodness of model fit
plot.resids(M3.final, E3.final, "X")
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

Form4 <- update(Form4i, .~. +summer)

M4.final <- lme(data=coc2, Form4, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4.final <- residuals(object=M4.final, type="normalized")

#plot best model residuals - check for goodness of model fit
plot.resids(M4.final, E4.final, "X")
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
#for Phosphorus, this all looks great


#------ Model 5: up to 3 predictors + landuse + rain ---------#

M5.final <- lme(data=coc2, Form5, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E5.final <- residuals(object=M5.final, type="normalized")

#plot best model residuals - check for goodness of model fit
plot.resids(M5.final, E5.final, "X")
#plot residuals against agency, year, month, date and other unused predictors; check if model is missing something important!
boxplots.resids2(M5.final, E5.final, "X")
#is temporal auto-correlation present?  Look for patterns in Auto-correlation plot for residuals
par(mfrow=c(2,2))
plot(E5.final~coc$start_date, pch=16, col="cadet blue")
acf(E5.final, na.action=na.pass, main="Auto-correlation plot for residuals")  #look for lines extending past blue dashed range
#Check for spatial correlation
mydata2 <- data.frame(E5.final, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
bubble(mydata2,"E5.final",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
#variogram to test for spatial correlation
Vario1 = variogram(E5.final ~ 1, mydata2)
plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
Vario2 <- variogram(E5.final ~ 1, mydata2, alpha = c(0, 45, 90,135) )
plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?
#for Phosphorus, the bubble plot looks suspicious, but the variogram looks okay.  Hmmm...


#plot parameter estimates for each model
theme_set(theme_sjplot())

plot_model(M3.final, sort.est=TRUE, show.values = TRUE, value.offset =.3, title="Phosphorus: Parameter Estimates, Model 3")
plot_model(M4.final, sort.est=TRUE, show.values = TRUE, value.offset =.3, title="Phosphorus: Parameter Estimates, Model 4")
plot_model(M5.final, sort.est=TRUE, show.values = TRUE, value.offset =.3, title="Phosphorus: Parameter Estimates, Model 5")
#plot_model(bestM, transform="exp")  #transform all estimates by applying exp() function

# plot_model(M4.final, type="slope")
# plot_model(M4.final, type="pred")
# plot_model(M4.final, type = "diag")



#generate plots that are shown in Rmarkdown script
M4k.final <- lme(data=coc2, P.Form4k, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4k.final <- residuals(object=M4k.final, type="normalized")

P.null <- gls(data = coc2, result ~ 1, method = "REML") 
P.M3 <- lme(data = coc2, Form3, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
P.M4i <- lme(data = coc2, Form4, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
P.M4k <- lme(data = coc2, update(Form4k, .~. +summer), random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
P.M5 <- lme(data = coc2, Form5, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

P_models <- list(
  P.null = P.null,
  P.M3 = P.M3,
  P.M4i = P.M4i,
  P.M4k = P.M4k,   ###maybe eliminate this one or the other M4?
  P.M5 = P.M5)

huxtablereg(P_models,
            single.row = TRUE, custom.model.names = names(P_models)) %>%
  set_bottom_border(1, -1, 0.4) %>%
  set_bold(1, -1, TRUE) 

plotreg(P_models, custom.title = "Regression Results, Total Phosphorus", custom.model.names = names(P_models))
plot_models(P_models[-c(1)],m.labels = names(P_models[-c(1)]),legend.title = "Models", show.values = TRUE,show.intercept = TRUE)
            

#save important items with phosphorus-specific names
P.coc2 <- coc2
P.r1X <- r1X
P.vf1X <- vf1X
P.FormLU <- Form3
P.FormBest <- Form4
P.rain <- rain
  
save(P.coc2, P.r1X, P.vf1X, P.FormLU, P.FormBest, rain, file=here("results", "Total Phosphorus Models.RData"))



