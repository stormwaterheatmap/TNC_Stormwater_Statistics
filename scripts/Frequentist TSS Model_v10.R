# This generalized script assesses linear relationships between TSS and landscape predictors
# using mixed effects models

# Note: This version takes the updated set of predictors with updated polygons (from Sep 19, 2026) and examines each set
#       with the COC of interest, to see which set of predictors may be most correlated with the COC.  Based on linear fits, 
#       a select set of predictors (best_predictors) is tested in various combinations within LME models.

#Four models are generated for the COC of interest (note the model numbers are 1,3,4,5 (no 2)):
#   1.) median value for all locations, irrespective of date
#   3.) COC value based only on land use (categorical) + rain, precip, summer (as applicable)
#   4.) COC value based on up to 3 landscape parameters NOT including land use + rain


# Eva Dusek Jennings
# Oct 5, 2026
#----------------------------------------------

rm(list = ls(all = T))

# Load libraries ----------------------------------------------
library(effects)
library(EnvStats)  #note the objects that are masked from other packages...
#library(tidyverse)
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
library(stringr)
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
aa <- which(str_detect(names(s8), "daymet"))  #after the last predictor column is a series of "daymet" columns
bb <- aa[which(aa > 34)][1] - 1  #identify the last predictor column
predictors <- names(s8)[34:bb]  #first predictor is always column 34; get all predictors
n.preds <- length(predictors)

#select the chemical of interest
this_param <- "Total Suspended Solids - Water - Total"
this_param_short <- "TSS"
this_param_abbrev <- "TSS"
coc <- s8[which(s8$parameter==this_param),]  #create the "coc" dataframe for this chemical of interest only

#Non Detect Handing -------------------------------

#look at censored points
ggplot(coc)+geom_jitter(aes(x=agency,y=conc,color=nondetect_flag))+scale_y_log10()

#which values are ND?  min/max detection limit?  where and when were these samples recorded?
coc.nd <- coc[which(coc$nondetect_flag==TRUE),]
pct.cen <- 100*nrow(coc.nd)/nrow(coc)  #for TSS, 0.83% of samples are ND's
nrow(coc.nd)
min(coc.nd$conc)
max(coc.nd$conc)
coc.nd[,c("location_id", "start_date", "conc")]  #location, date, and detection limit of ND samples
# for TSS, no pattern in ND samples; 3 total: SNO_COM on 4/14/10, KICLDR on 1/26/11 and 3/13/13

#Since there are only 3 ND's for TSS, substitute the ND's with half the detection limit
coc <- coc %>% 
  mutate(conc = ifelse(nondetect_flag, conc * 0.5, conc))

ggplot(coc)+geom_jitter(aes(x=agency,y=conc,color=nondetect_flag))+scale_y_log10()


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

# best distribution is log-normal.  Ln-transform concentration data in column "result"
coc$result <- log(coc$conc)


#-------------------------------------------------------------------------------------------#
#  Explore Cleveland Dotplots and Boxplots of COC data or predictors conditional on Agency  #
#-------------------------------------------------------------------------------------------#

colors_agency <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
colors_location <- colors_agency[c(1,1,1,2,2,2,3,4,4,4,5,5,5,6,6,6)]

if(run_exploratory_code==TRUE) {
  #Cleveland Dotplot - use to detect violations of homogeneity: We are looking for whether the spread of data values
  #   differs between sampling locations (or agencies).  If so, it indicates heterogeneity, and that there may be problems with violation of
  #   homogeneity in a linear regression model applied on these data.  We're also looking for outliers.
  par(mfrow=c(1,1))
  dotchart(coc$result, groups=coc$location_id, pch=19, col=colors_agency[as.numeric(coc$agency)],
           xlab="concentration", main=paste("Cleveland Dotplot:", this_param_short))
  
  #pairs plots allow visualization of interactions between possible predictors.  Look for relationships between "result"
  #   and all of the predictor variables
  pairs(coc %>% dplyr::select(result, predictors[1:ceiling(n.preds/3)]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  pairs(coc %>% dplyr::select(result, predictors[ (ceiling(n.preds/3)+1) : (ceiling(n.preds/3)*2) ]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  pairs(coc %>% dplyr::select(result, predictors[ (ceiling(n.preds/3)*2+1) : n.preds]), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  
  #look also at interactions between result and rainfall
  pairs(coc %>% dplyr::select(result, antecedant_dry_days_std, daymet_precip_std, daymet_3day_std,
                       daymet_7day_std, daymet_14day_std, daymet_21day_std, daymet_28daySR_std), 
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)
  #for TSS, important predictors are likely: 
  #  sqrt_traffic, sqrt_CO2_road, devAge
  #  daymet_precip (0.2)
}

#-----------------------------------------------------------------------------------------------------#
#  Plot COC vs various predictors: look for potential predictors to model & sources of heterogeneity  #
#-----------------------------------------------------------------------------------------------------#

lp_plots()
pr_plots()
tlp_plots()
mp_plots()
#for TSS, daymet_precip shows evidence of heteroskedasticity; also has best fit  
#   location, agency, landuse, and year also show heteroskedasticity

#best predictors for this COC; make sure to only have one version (transformed or not) of each predictor!
best_predictors <- c("totRES", "grass", "trees", "roofs", "impervious", "nodev", "ln_traffic_k100", #"sqrt_traffic_k100", "ln_traffic_k50", 
                     "sqrt_popn", "sqrt_CO2_res", "sqrt_CO2_com", "sqrt_CO2_road", "sqrt_CO2_transport", "sqrt_CO2_total",
                     "devAge", "oldDev")

pred_i <- which(predictors %in% best_predictors)
lp_plots2(pred_i)


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

elements_to_remove <- c("totRES", "devAge", "sqrt_CO2_total", "sqrt_CO2_res", "nodev", "ln_traffic_k100")  #these predictors are highly correlated with others
best_predictors2 <- best_predictors[!(best_predictors %in% elements_to_remove)]

pairs(coc[best_predictors2], lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)


#--------------------------------------------------#
#  Any evidence of changes in COC conc over time?  #
#--------------------------------------------------#

if(run_exploratory_code==TRUE) {
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
  #   For TSS, maybe some weird changes over time for Tacoma and Port of Tacoma, 
  #       but time-series is not long enough to tell for sure
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
  #for TSS, the only consistent pattern is increase in TSS with increased daymet_precip on sampling day 
  #  (more intense rain = higher TSS)
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
         rain1=daymet_precip_std,  #choose whichever precip measure makes the most sense for this coc!
         rain14=daymet_14day_std,
         rain21=daymet_21day_std,
         mPrecip=mPrecip  #note: if one of the other mPrecip measures looked better, use it here instead!
  ) %>%
  mutate(year=as.factor(year))  %>%
  mutate(landuse = case_when(
    land_use == "LDR" ~ "1.LDR",
    land_use == "HDR" ~ "2.HDR",
    land_use == "COM" ~ "3.COM",
    land_use == "IND" ~ "4.IND")) %>%
  mutate(landuse = as.factor(landuse)) %>%
  mutate(spring = case_when(
    season=="2" ~ "1",
    TRUE ~ "0")) %>%
  mutate(summer = case_when(
    season=="3" ~ "1",
    TRUE ~ "0")) %>%
  mutate(spring=as.factor(spring),
         summer = as.factor(summer))

monthlyAvgPrecip <- data.frame(month=c(1:12), 
                               precip=c(5.24, 4.09, 3.92, 2.75, 2.03, 1.55, 0.93, 1.16, 1.61, 3.24, 5.67, 6.06))
coc2$monthlyPrecip <- monthlyAvgPrecip[coc2$month,2]

#rain <- "daymet precip, standardized"  #precip measure that was used for this COC


# plot(fitted(Model4), E4, main="residuals", col="gray", pch=16)
# qqplot(E4)
# hist(E4)


#------------------------------------------------#
#  Formulas for Possible Predictor Combinations  #
#------------------------------------------------#

weather <- c("rain1", "spring", "rain21")

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
  
  # #is there heteroscedasticity in any of these potential covariates?
  # plot(coc2$result ~ coc2$location)
  # plot(coc2$result ~ coc2$agency)
  # plot(coc2$result ~ coc2$month)  #unlikely to be considered for use as a covariate
  # plot(coc2$result ~ coc2$season)
  # plot(coc2$result ~ coc2$year)
  # plot(coc2$result ~ coc2$dry)
  # abline(lm(coc2$result ~ coc2$dry))
  # plot(coc2$result ~ coc2$rain)
  # abline(lm(coc2$result ~ coc2$rain))
  # plot(coc2$result ~ coc2$monthlyPrecip)
  # plot(coc2$result ~ coc2$mPrecip)
  # #for TSS, rain is the most likely source of heteroskedasticity
  # 
  
  #Step 3:  Choose a variance structure (if there was heterogeneity)
  #-------  for selecting random structure, use REML to compare - its better at capturing random structure
  #         REML estimates the random effects by considering linear combinations of the data that remove the 
  #         fixed effects. If these fixed effects are changed, the likelihoods of the two models will not be directly comparable.
  #         Compare AIC for the various beyond-optimal models with different variance structures.
  #         Alternately, can use likelihood ratio test (anova) to compare nested models
  
  M.gls2X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency))
  M.gls3X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|location))
  
  anova(M.gls1X, M.gls2X, M.gls3X)
  #for TSS: AIC is best for location, and BIC is best for agency; go with location for now
  
  #residual plots for best fit models -- look for homogeneity of residuals
  E.3X <- resid(object=M.gls3X, type="normalized")  
  boxplots.resids2(M.gls3X, E.3X, "X")

  par(mfrow=c(2,2), mar=c(4.5, 4.5, 4, 1), xaxt="s")
  plot(fitted(M.gls1X), resid(object=M.gls1X, type="normalized"), main="no variance structure", xlab="fitted values", ylab="normalized residuals", col="gray", pch=16)
  plot(fitted(M.gls2X), resid(object=M.gls2X, type="normalized"), main="variance covariate: agency", xlab="fitted values", ylab="normalized residuals", col="gray", pch=16)
  plot(fitted(M.gls3X), resid(object=M.gls3X, type="normalized"), main="variance covariate: location", xlab="fitted values", ylab="normalized residuals", col="gray", pch=16)
  
  #look at the relative amount of variability within each location, and relative amount of variability within each agency
  par(mfrow=c(2,1), mar=c(4, 4, 1, 1))
  boxplot(result ~ location, data=coc2, col=colors_location)
  boxplot(result ~ agency, data=coc2, col=colors_agency)
  
}

#variance functions for best fit model so far
vf1X <- varIdent(form= ~1|location) 


# Steps 4-6: Find the proper random effects structure; look for temporal and spatial autocorrelation
#-----------  note: look for spatial correlation AFTER setting random effects!
#because of the data structure (repeat sampling at each location over time), location MUST be an RE.  Question is, should it be nested w/in agency?
#  check for whether the agency piece accounts for >1% of the total variance once best model is chosen.  If agency contribution is minimal, remove it then.

if (run_exploratory_code==TRUE) {
  #----Model with Variance Function 1----#
  M.vf1X <- gls(FormX, data = coc2, method = "REML",  weights = vf1X)
  E.vf1X <- resid(object=M.vf1X, type="normalized")  
  
  #plots of residuals: no variance structure vs. variance structure
  par(mfrow=c(2,2), mar=c(4,4,4,1))
  plot(fitted(M.gls1X), E.1X, main="no variance structure", xlab="fitted values", xaxt="s", ylab="normalized residuals", col="gray", pch=16)
  plot(fitted(M.vf1X), E.vf1X, main="variance covariates: location", xlab="fitted values", xaxt="s", ylab="normalized residuals", col="gray", pch=16)
  
  #Random intercept model; this is nested in the best fit model from Step 2, so can compare with a likelihood ratio test
  M1.lme1X <- lme(data=coc2, FormX, random = ~1|agency/location, method="REML", weights=vf1X,
                  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  anova(M.vf1X, M1.lme1X)  #the likelihood ratio test p-value compares the two nested models, and whether they are significantly different
  #For TSS, the model with no random effects has a lower AIC, BUT I think this is b/c I would have to have location as a fixed effect in the gls model.
  #   regardless, location has to be a random effect and for now, it should be nested w/in agency.
}
  
#random structure for best fit model so far        ##### NOTE: the "best fit" is with no random structure, but the data dictates location as random effect
r1X <- formula(~1|agency/location)

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
  #for TSS, no indication of temporal autocorrelation
  
  #Check for spatial correlation
  mydata2 <- data.frame(E1.rX, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
  coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
  bubble(mydata2,"E1.rX",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
  #variogram to test for spatial correlation
  Vario1 = variogram(E1.rX ~ 1, mydata2)
  plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
  Vario2 <- variogram(E1.rX ~ 1, mydata2, alpha = c(0, 45, 90,135) )
  plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?
  #for TSS, no indication of spatial autocorrelation
}

# Steps 7-8: Find the proper fixed effects structure; use ML for model comparisons
#-----------  


#-------------------------------#
#  Save Data & Model Structure  #
#-------------------------------#

TSS.coc2 <- coc2
TSS.r1X <- r1X
TSS.vf1X <- vf1X
#TSS.rain <- rain
TSS.preds <- best_predictors

save(TSS.coc2, TSS.r1X, TSS.vf1X, #TSS.rain, 
     TSS.preds, file=here("results", "S8 Data & Model Structure_TSS.RData"))




#-------------------------------------------------------------#
#  Exhaustive Search for Model with Best 3 Predictors + Rain  #
#-------------------------------------------------------------#

preds.1 <- best_predictors


if(run_exploratory_code==TRUE) {
  #Construct all possible models with 1, 2 and 3 landscape predictors, both including and not including rain;
  #   use AIC to determine which models provide the best fit
  
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
  
  #combine one and two-predictor options into one dataframe
  dd <- rbind(dd.1, dd.2) # , dd.3)  #, dd.1, dd.2, dd.3)
  
  #identify which predictors have strong correlation with others; eliminate models with high inter-predictor correlation
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
  
  #generate an exhaustive set of formulas with one, two and three predictors, including weather-related predictors (season, rain)
  ee.1 <- list(NULL)  #list of formulas based on the sets of predictors in "dd"; including rain
  for(i in 1:nrow(dd)) {
    if(is.na(dd$Var2[i])) ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1]))
    else if(is.na(dd$Var3[i]))  ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1], "+", dd[i,2]))
    else  ee.1[[i]] <- as.formula(paste("result~ ", paste(weather, collapse=" + "), " + ", dd[i,1], "+", dd[i,2], "+", dd[i,3]))
  }
  ee <- ee.1
  
  #list of equations and whether the signs match what we expect from linear models
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
  proc.time() - ptm  #stop the clock;  ~9 minutes for 386 models
  
  #identify which of the above models have predictors that are in the expected (+ or -) direction
  xx <- which(MyEqns$do_signs_match==TRUE)
  
  ee.all <- ee
  ee <- ee[xx]
  my.aics.all <- my.aics
  my.aics <- my.aics[xx]
  
  par(mfrow=c(5,5), mar=c(1,2,2,0), oma=c(1,1,1,1))
  gg <- rep("black", length(ee)) ### NOTE: the plot with 1, 2 or 3 predictors doesn't work!  it refers to dd, and ee no longer matches it!
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
  
  save(my.formulas.m4, my.aics.m4, file=(here("scripts", "formulas_and_aics", "TSS_m4.RData")))
}

if (run_exploratory_code==FALSE) {
  load(file=here("scripts", "formulas_and_aics", "TSS_m4.RData"))
}

if (F) {
  load(file=here("results", "S8 Data & Model Structure_TSS.RData"))
  coc2 <- TSS.coc2
  r1X <- TSS.r1X
  vf1X <- TSS.vf1X
#  rain <- TSS.rain
  best_predictors <- TSS.preds
}


#--------------------------------------------#
#  Quantiles List for Plotting Interactions  #
#--------------------------------------------#

#list of quantiles for all landscape predictors
qList <- list(oldDev=quantile(coc2$oldDev, probs=c(0,.25,.50,.75,1)),
              sqrt_CO2_road=quantile(coc2$sqrt_CO2_road, probs=c(0, .25, .5, .75, 1))) #,
              # rain = quantile(coc2$rain, probs=c(0, .25, .5, .75, 1)))



#------------------------#
#  Models 1, 3, 4 and 5  #
#------------------------#

#------ Model 1: median value for all locations --------#

Model1 <- gls(data=coc2, result~1, method="ML")
E1 <- residuals(object=Model1, type="normalized")


#------ Model 3: land use + rain with variance structure & random effects ---------#

Form3 <- as.formula(paste("result ~ landuse + ", paste(weather, collapse=" + ")))
Model3 <- lme(data=coc2, Form3, random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E3 <- residuals(object=Model3, type="normalized")

#------ Model 4: up to 3 predictors + precip/ season components ---------#

#run through the top formulas when only best predictors are on the table;
#  keep only formulas that make sense
myForm <- my.formulas.m4[[1]]   #1 and 2 make sense

#these lines of code assess fit of this particular model in terms of COC vs. individual predictors, and predictor correlation
myModel <- lme(data=coc2, myForm, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(myModel)
AIC(myModel)
plot.single.preds(myModel)
check.cor(myModel)
boxplots.resids2(myModel, residuals(myModel, type="normalized"), "X")

#myModel <- update(myModel, .~. -rain21 +rain14)  #while rain14 is still significant, it has a higher AIC by 2 pts

#formulas that are worth considering (single plots of predictors make sense)
Form4a <- formula(result ~ rain1 + spring + rain21 + oldDev + sqrt_CO2_road)  #my.formulas.m4[[1]]  
Form4b <- formula(result ~ rain1 + spring + rain21 + oldDev)  #my.formulas.m4[[2]]  

#just out of curiousity, try summer as a predictor, rather than spring
Form4a.1 <- update(Form4a, .~. + summer - spring)
Model4a.1 <- lme(data=coc2, Form4a.1, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(Model4a.1)

AIC(myModel, Model4a.1)  #spring is the better predictor for TSS


#####  Model4a  ####

Form4 <- Form4a  #the model with sqrt_CO2_road is preferable, as it more directly explains relationship with TSS 
Model4 <- lme(data=coc2, Form4, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4 <- residuals(object = Model4, type = "normalized")

#try adding interactions and see if they are significant
M4.sub1 <- update(Model4, . ~ . + rain1:oldDev) #just barely significant (p=0.0443
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, .~. +rain1:sqrt_CO2_road)  #not significant (p=0.537)
anova(Model4, M4.sub2)

#none of these are significant enough to include
summary(M4.sub1)

M4.sub1 <- update(Model4, . ~ . + rain21:oldDev) #not significant (p=0.495
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, .~. +rain21:sqrt_CO2_road)  #not significant (p=0.3974)
anova(Model4, M4.sub2)


#update Model4a with interaction
Form4.int <- update(Form4, .~. +rain1:oldDev)
Model4.int <- update(Model4, .~. +rain1:oldDev)
E4.int <- residuals(object = Model4.int, type = "normalized")
summary(Model4.int)
boxplots.resids2(Model4.int, E4.int, "X")

Plot.Quantile("rain1", "oldDev", Model4.int, yEqn=13.1, xEqn=1)
plot.single.preds(Model4.int, xEqn=0)

#This relationship seems plausible for TSS.  More daily rain means more TSS in stormwater, but for highly developed (old dev) areas
#  the TSS starts higher but doesn't increase as fast as for undeveloped areas, where TSS starts lower (during sprinkling rain) but
#  really increases as the rain drives hard (think of dirt running off during big rains).


#model 4b is as follows:
Model4b.int <- update(Model4.int, .~. -sqrt_CO2_road)

AIC(Model4.int, Model4b.int)

#-------  Best fit Model  --------#

Form4 <- Form4.int
Model4 <- Model4.int
E4 <- E4.int



#------ Compare models and plot ---------#
AIC(Model1, Model3, Model4)

AIC(Model1, Model3, Model4)[,2] - AIC(Model4)


if(run_exploratory_code==TRUE) {
  #plot model predictions for each model (above)
  M1.preds <- predict(Model1)
  plot.preds.vs.results(M1.preds)
  
  M3.preds <- predict(Model3)
  plot.preds.vs.results(M3.preds)
  
  M4.preds <- predict(Model4)
  plot.preds.vs.results(M4.preds)
}


#----------------------------------#
#  Look at the Fits of the Models  #
#----------------------------------#

label_loc <- c(2,5,7,9,12,15)
label_text <- c("King", "Pie", "POT", "Sea", "Sno", "Tac")

#check residuals for each model, to see if they tell us anything important
par(mfrow=c(2,2), mar=c(2,4,4,1), oma=c(0,0,0,0))
plot(coc2$location, E1, main="Null Model", ylab="Residuals", xaxt="n", col=colors_location)
axis(side=1, at=label_loc, labels=label_text)
abline(h=0, col="gray")
plot(coc2$location, E3, main="Categorical Land Use Model", ylab="Residuals", xaxt="n", col=colors_location)
axis(side=1, at=label_loc, labels=label_text)
abline(h=0, col="gray")
plot(coc2$location, E4, main="Spatial Predictor Model", ylab="Residuals", xaxt="n", col=colors_location)
axis(side=1, at=label_loc, labels=label_text)
abline(h=0, col="gray")
# plot(coc2$location, E4.alt, main="Model 4, alternate", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,2,2,3,4,4,4,5,5,5,6,6,6)])
# axis(side=1, at=c(2,5,7,9,12,15), labels=c("King", "Pierce", "POT", "Sea", "Sno", "Tac"))
# abline(h=0, col="gray")


# Step 9: Refit models with REML and apply graphical model validation; check for homogeneity, 
#--------    normality and independence

#------ Model 3: land use + rain ---------#

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
#variogram to test for spatial correlation; for Nitrite-Nitrate, there may be spatial autocorrelation with Model3
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
#everything looks good!


#-----------------#
#  Save Formulas  #
#-----------------#

#save important items with TSS-specific names
TSS.coc2 <- coc2
TSS.r1X <- r1X
TSS.vf1X <- vf1X
TSS.Form3 <- Form3
TSS.Form4 <- Form4
#TSS.rain <- rain
TSS.preds <- best_predictors

save(TSS.coc2, TSS.r1X, TSS.vf1X, TSS.Form3, TSS.Form4, #TSS.rain, 
     TSS.preds, file=here("results", "Frequentist_TSS Models.RData"))



#-------------------------------#
#  Plot Predictor Coefficients  #
#-------------------------------#

#generate plots that are shown in Rmarkdown script
TSS.null <- gls(data = coc2, result ~ 1, method = "REML") 
TSS.M3 <- lme(data = coc2, Form3, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
TSS.M4 <- lme(data = coc2, Form4, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

TSS_models <- list(
  Null_Model = TSS.null,
  Categorical_Land_Use_Model = TSS.M3,
  Spatial_Predictor_Model = TSS.M4)

huxtablereg(TSS_models,
            single.row = TRUE, custom.model.names = names(TSS_models)) %>%
  set_bottom_border(1, -1, 0.4) %>%
  set_bold(1, -1, TRUE) 

#plot parameter estimates for each model
theme_set(theme_sjplot())
plotreg(TSS_models, custom.title = "Regression Results, TSS", custom.model.names = names(TSS_models))
plot_models(TSS_models,m.labels = names(TSS_models),
            # colors = c("#482173", "#25858e", "#86d549"),
            colors = c("#0072B2", "#D55E00", "goldenrod2"), 
            legend.title = "Models", show.values = TRUE,show.intercept = TRUE)



