# This generalized script assesses linear relationships between Total Kjeldahl Nitrogen and landscape predictors
#  Because TKN has a relatively high percentage of ND data (12%), this script does not perform mixed effects
#  models, which would be skewed by the ND data.

# Eva Dusek Jennings
# Aug 10, 2022
#----------------------------------------------

rm(list = ls(all = T))


# Load libraries ----------------------------------------------
library(effects)
library(EnvStats)  #note the objects that are masked from other packages...
#library(tidyverse)
library(stringr)
library(ggplot2)
library(graphics)
library(stats)
library(nlme)
library(grid)
library(gridExtra)
#library(lme4)  #for lmer models
library(lattice)
library(gstat)  #for spatial correlation exploration
library(sp)  #for spatial correlation exploration
library(sjPlot)
library(sjlabelled)
library(TMB)
library(glmmTMB)
library(here)
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
this_param <- "Total Kjeldahl Nitrogen - Water - Total"
this_param_short <- "TKN"
coc <- s8[which(s8$parameter==this_param),]  #create the "coc" dataframe for this chemical of interest only

#Non Detect Handing -------------------------------
#look at censored points
ggplot(coc)+geom_jitter(aes(x=agency,y=conc,color=nondetect_flag))+scale_y_log10()
ggplot(coc)+geom_jitter(aes(x=loc,y=conc,color=nondetect_flag))+scale_y_log10()

#which values are ND?  min/max detection limit?  where and when were these samples recorded?
coc.nd <- coc[which(coc$nondetect_flag==TRUE),]
pct.cen <- 100*nrow(coc.nd)/nrow(coc)  #for TKN, 41 samples (9.76%) are ND.  ND level varies greatly (30 to 500)
nrow(coc.nd)
min(coc.nd$conc)
max(coc.nd$conc)
coc.nd[,c("location_id", "start_date", "conc")]  #location, date, and detection limit of ND samples

coc[which(coc$loc=="TAC_HDR" & coc$conc==255),]  #what is going on with this sample that has a strange
#  detection limit??  I looked it up in the raw data, and it is flagged with a "result data qualifier = U"...
#  I think this is an "estimated" value -- lets keep this value of 255 for the censored statistics.
coc[which(coc$loc=="TAC_HDR" & coc$conc==255), "nondetect_flag"] <- FALSE 

#save a copy of the original concentrations
coc$oconc <- coc$result

#For comparison's sake, save a copy where ND's are substituted with half the detection limit
coc <- coc %>%
  mutate(result_halfND = ifelse(nondetect_flag, conc * 0.5, conc))

#use Regression on Order Statistics (ros) to generate missing values (NADA package).  Go location-by-location for this part
coc1 <- NULL  #this will be the new DF that is generated from ROS statistics
for (i in 1:length(unique(coc$loc))) {
  this.loc <- coc[which(coc$loc==levels(coc$loc)[i]), ]  #identify all rows for the current location
  ros.this.loc <- ros(obs=this.loc$result, censored=this.loc$nondetect_flag)  #run the ROS statistics; get modeled values
  df.this.loc <- as.data.frame(ros.this.loc)  #this gives all modeled data points, arranged lowest to highest
  loc.sorted <- this.loc[order(this.loc$result), ]  #sort this.loc by result value, to match ros df
  loc.ros <- cbind(loc.sorted, df.this.loc)  #bind modeled ros values to all columns for the current location
  coc1 <- rbind(coc1, loc.ros)  #build the coc1 data frame, location by location
}

ggplot(coc1)+geom_jitter(aes(x=loc,y=modeled,color=nondetect_flag))+scale_y_log10() #plot new values, by location

#rename
coc <- coc1  #for ease of use down the line, change the name of coc1 to coc 
coc$result <- coc$modeled  #for ease of use down the line, put modeled values into the "result" column

# distribution exploration ----------------------------------------------
if (run_exploratory_code ==TRUE) {
  #explore the data; look for underlying distribution
  par(mfrow=c(2,2))
  hist(coc$result, breaks=50)
  qqnorm((coc$result), main=paste("QQ-Normal plot of", this_param_short))
  qqline((coc$result))
  hist(log(coc$result), breaks=50)
  qqnorm(log(coc$result), main=paste("QQ-Normal plot of log(", this_param_short, ")", sep=""))
  qqline(log(coc$result))
  
  #consider other distribution (log-normal and gamma).  Compare QQ plots for these distributions
  par(mfrow=c(2,2))
  qqPlot(log(coc$result), dist="norm", estimate.params=TRUE, add.line=TRUE, main=paste("QQ log-normal plot of", this_param_short))
  qqPlot(log(coc$result), dist="norm", param.list=list(mean=-3.3, sd=2.09), add.line=TRUE, main=paste("QQ log-normal plot of", this_param_short))
  qqPlot(coc$result, dist="gamma", param.list=list(shape=.7, scale=1), points.col="blue",
         add.line=TRUE, main=paste("QQ-gamma plot of", this_param_short))
  qqPlot(sqrt(coc$result), dist="norm", estimate.params=TRUE, add.line=TRUE, main=paste("QQ plot of sqrt-transformed", this_param_short))
  #qqPlot(1/(coc$conc), dist="norm", estimate.params=TRUE, add.line=TRUE)
}

# log-normal distribution works.  Let's ln-transform concentration data in column "result"
coc$result <- log(coc$result)


#-------------------------------------------------------------------------------------------#
#  Explore Cleveland Dotplots and Boxplots of COC data or predictors conditional on Agency  #
#-------------------------------------------------------------------------------------------#

colors_agency <- c("red", "orange", "yellow", "green", "blue", "purple")
colors_location <- colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)]

if (run_exploratory_code ==TRUE) {
  #Cleveland Dotplot - use to detect violations of homogeneity: We are looking for whether the spread of data values
  #   differs between sampling locations (or agencies).  If so, it indicates heterogeneity, and that there may be problems with violation of
  #   homogeneity in a linear regression model applied on these data.  We're also looking for outliers.
  par(mfrow=c(1,1))
  #show dotchart for (ln-transformed) original concentrations, incl reporting limits
  # dotchart(coc$ln.oconc, groups=coc$location_id, pch=19, col=colors_agency[as.numeric(coc$agency)],
  #          xlab="concentration", main=paste("Cleveland Dotplot:", this_param_short))
  #show dotchart for "result", where reporting limits were halved for non-detects
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
}

#-----------------------------------------------------------------------------------------------------#
#  Plot COC vs various predictors: look for potential predictors to model & sources of heterogeneity  #
#-----------------------------------------------------------------------------------------------------#

lp_plots()
pr_plots()
tlp_plots()
mp_plots()

#for TKN:
#14-day & 21-day precip look the best -- most evenly spaced data...
#mPrecip is better than the rest!
#   location, agency(!), landuse(!!), month(!!), season(!) and year do show heteroskedasticity

#best predictors for this COC; make sure to only have one version (transformed or not) of each predictor!
best_predictors <- c("roofs", "nodev", "sqrt_traffic", "sqrt_popn", "sqrt_CO2_res", "sqrt_CO2_tot", 
                     "sqrt_CO2_road", "devAge2", "roof_intURB", "roof_intURB_IND")

pred_i <- which(predictors %in% best_predictors)
lp_plots(pred_i)

#examine correlations between predictors; generate a second vector of only predictors that aren't highly correlated
#  this second vector will be used to generate FormX
pairs(coc[best_predictors], lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist)

elements_to_remove <- c("sqrt_CO2_road", "nodev", "sqrt_nodev", "sqrt_popn", "roof_intURB_IND")  #these predictors are highly correlated with others
best_predictors2 <- best_predictors[!(best_predictors %in% elements_to_remove)]

bp_coefs <- bp_signs <- rep(NA, length(best_predictors))
names(bp_coefs) <- names(bp_signs) <- best_predictors
for(i in 1:length(best_predictors)) {
  aa <- lm(as.formula(paste("result ~ ", best_predictors[i])), data=coc)
  bp_coefs[i] <- coefficients(aa)[2]
}
bp_signs <- sign(bp_coefs)



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
  #   For TKN, there might be a seasonal trend in Snohomish?
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
}


#--------------------------------------------#
#  Create coc2: only predictors we will use  #
#--------------------------------------------#

#see lines 126-128 for selection of predictors that are important for this COC!
coc2 <- coc %>%
  dplyr::select(result,
         oconc,
         cen=nondetect_flag,
         location=loc, 
         latitude, longitude,
         agency,
         land_use,
         year, 
         month, season, start_date,
         all_of(best_predictors),
         dry=antecedant_dry_days_std,
         rain=daymet_14day_std,  #choose whichever precip measure makes the most sense for this coc!
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
  mutate(summer=as.factor(summer)) %>%
  mutate(season3 = case_when(
    season=="1" ~ "winter_spring",
    season=="2" ~ "winter_spring",
    season=="3" ~ "summer",
    season=="4" ~ "fall")) %>%
  mutate(season3 = as.factor(season3))

monthlyAvgPrecip <- data.frame(month=c(1:12), 
                               precip=c(5.24, 4.09, 3.92, 2.75, 2.03, 1.55, 0.93, 1.16, 1.61, 3.24, 5.67, 6.06))
coc2$monthlyPrecip <- monthlyAvgPrecip[coc2$month,2]

rain <- "daymet 14-day, standardized"  #precip measure that was used for this COC


#------------------------------------------------#
#  Formulas for Possible Predictor Combinations  #
#------------------------------------------------#

if (run_exploratory_code==TRUE) {
  #looking to see if there are relationships between TKN and various seasonal items.  Could try 10-day avg daymet temperature??
  par(mfrow=c(2,2))
  boxplot(coc2$result~coc2$season)  #winter & spring are pretty indistinguishable.  Combine into one, using season3
  boxplot(coc2$result~coc2$month)
  boxplot(coc2$result~coc2$monthlyPrecip)
  boxplot(coc2$result~coc2$dry)
}

weather <- c("rain", "summer", "dry")

#FormX <- as.formula( paste("result ~ landuse + ", paste(weather, collapse=" + "), " + ", paste(best_predictors2, collapse=" + ")) )
FormX <- as.formula( paste("result ~ ", paste(weather, collapse=" + "), " + ", paste(best_predictors2, collapse=" + ")) )
#NOTE: it is not practical at this point to add interactions -- that will come once we find best model using just main effects


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
  
  
  #Step 3:  Choose a variance structure (if there was heterogeneity)
  #-------  for selecting random structure, use REML to compare - its better at capturing random structure
  #         REML estimates the random effects by considering linear combinations of the data that remove the 
  #         fixed effects. If these fixed effects are changed, the likelihoods of the two models will not be directly comparable.
  #         Compare AIC for the various beyond-optimal models with different variance structures.
  #         Alternately, can use likelihood ratio test (anova) to compare nested models
  
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
  plot(coc2$result ~ coc2$mPrecip, col=colors_location[as.numeric(coc2$location)])
  #for 
  
  
  #Step 3:  Choose a variance structure (if there was heterogeneity)
  #-------  for selecting random structure, use REML to compare - its better at capturing random structure
  #         REML estimates the random effects by considering linear combinations of the data that remove the 
  #         fixed effects. If these fixed effects are changed, the likelihoods of the two models will not be directly comparable.
  #         Compare AIC for the various beyond-optimal models with different variance structures.
  #         Alternately, can use likelihood ratio test (anova) to compare nested models
  
  M.gls2X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency))
  M.gls3X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|location))
  
  anova(M.gls1X, M.gls2X, M.gls3X)
  #for TKN: AIC is slightly better (dAIC = 1.2) for location (3); BIC is much better for agency (2)
  
  #look at residuals for the three variance structure options
  par(mfrow=c(2,2), mar=c(4.5, 4.5, 4, 1), xaxt="s")
  plot(fitted(M.gls1X), resid(object=M.gls1X, type="normalized"), main="no variance structure", xlab="fitted values", ylab="normalized residuals", col="gray", pch=16)
  plot(fitted(M.gls2X), resid(object=M.gls2X, type="normalized"), main="variance covariate: agency", xlab="fitted values", ylab="normalized residuals", col="gray", pch=16)
  plot(fitted(M.gls3X), resid(object=M.gls3X, type="normalized"), main="variance covariate: location", xlab="fitted values", ylab="normalized residuals", col="gray", pch=16)
}

#variance functions for best fit models so far
vf1X <- varIdent(form= ~1|agency)

# Steps 4-6: Find the proper random effects structure; look for temporal and spatial autocorrelation
#-----------  note: look for spatial correlation AFTER setting random effects!

if (run_exploratory_code==TRUE) {

  #----Model with Variance Function 1----#
  M.vf1X <- gls(FormX, data=coc2, method="REML", weights=vf1X)
  
  #Random intercept model; this is nested in the best fit model from Step 2, so can compare with a likelhiood ratio test
  M1.lme1X <- lme(data=coc2, FormX, random = ~1|agency/location, method="REML", weights=vf1X,
                  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  anova(M.vf1X, M1.lme1X)  #the likelihood ratio test p-value compares the two nested models, and whether they are significantly different
  #For TKN, with the set of best_predictors2 (best_predictors with highly correlated predictors removed), 
  #  the model with random intercept is significantly better than without
}

#random structure for best fit model so far
r1X <- formula(~1|agency/location)

if (run_exploratory_code==TRUE) {
  #best random effects structure & residual plots to test it
  M1.rX <- lme(data=coc2, FormX, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  summary(M1.rX)
  E1.rX <- residuals(object=M1.rX, type="normalized")
  boxplots.resids2(M1.rX, E1.rX, "X")
  
  #Check for temporal correlation; look for patterns in Auto-correlation plot for residuals
  #  NOTE:  currently, data are sorted by magnitude, so there may be high auto-correlation artifact!
  par(mfrow=c(2,2))
  plot(E1.rX~coc2$start_date, pch=16, col="cadet blue")
  acf(E1.rX, na.action=na.pass, main="Auto-correlation plot for residuals")  #look for lines extending past blue dashed range
  #NOTE: we should fix this, as acf assumes the DF is sorted by date (NOT by value!!)
  
  #Check for spatial correlation
  mydata2 <- data.frame(E1.rX, coc2.longitude=jitter(coc2$longitude, amount=0.05), coc2.latitude=jitter(coc2$latitude, amount=0.05))  #jitter the X-Y coordinates, so that data points aren't on top of each other
  coordinates(mydata2)<-c("coc2.longitude","coc2.latitude")
  bubble(mydata2,"E1.rX",col=c("black","grey"), main="Residuals",xlab="X-coordinates", ylab="Y-coordinates")
  
  #variogram to test for spatial correlation
  Vario1 = variogram(E1.rX ~ 1, mydata2)
  plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
  Vario2 <- variogram(E1.rX ~ 1, mydata2, alpha = c(0, 45, 90,135) )
  plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?
  #for TKN, no indication of spatial autocorrelation
}


# Steps 7-8: Find the proper fixed effects structure; use ML for model comparisons
#-----------  

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
  
  # 
  # 
  # 
  #   
  # #run the lme for each formula, saving the AIC values; takes about 7 minutes to run
  # ptm <- proc.time()  #time the code!
  # my.aics <- rep(0, length(ee))
  # for (i in 1:length(ee)) {
  #   bb <- lme(data=coc2, ee[[i]], random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  #   my.aics[i] <- AIC(bb)
  # }
  # proc.time() - ptm  #stop the clock;  9 minutes for 386 models
  
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
  
  save(my.formulas.m4, my.aics.m4, file=(here("formulas_and_aics", "total_kjeldahl_nitrogen_m4_censtat.RData")))
}

if (run_exploratory_code==FALSE) {
  load(file=here("formulas_and_aics", "total_kjeldahl_nitrogen_m4_censtat.RData"))
}



#--------------------------------------------#
#  Quantiles List for Plotting Interactions  #
#--------------------------------------------#

# #list of quantiles for all landscape predictors
qList <- list(devAge2=quantile(coc2$devAge2, probs=c(0,.25,.50,.75,1)),
              sqrt_traffic=quantile(coc2$sqrt_traffic, probs=c(0,.25,.50,.75,1)),
              sqrt_CO2_road=quantile(coc2$sqrt_CO2_road, probs=c(0,.25,.50,.75,1))
)


#------------------------#
#  Models 1, 3, 4 and 5  #
#------------------------#

#NOTE: when I ran these including "dry", it was NOT significant.  Remove from fixed effects AND variance structure!

#------ Model 1: median value for all locations --------#

Model1 <- gls(data=coc2, result~1, method="ML")  #here, we use the actual concentration, rather than the transformed one
E1 <- residuals(object=Model1, type="normalized")

#------ Model 3: land use + rain with variance structure & random effects ---------#

Form3 <- update(as.formula(paste("result ~ landuse + ", paste(weather, collapse=" + "))), .~. -dry)
#Form3 <- as.formula(paste("result ~ landuse + ", paste(weather, collapse=" + ")))
Model3 <- lme(data=coc2, Form3, random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E3 <- residuals(object=Model3, type="normalized")

#------ Model 4: no land use; up to 3 predictors + rain ---------#

#run through the top formulas when only best predictors are on the table;
#  keep only formulas that make sense
#  NOTE:  removing "dry" from these equations makes sense, as dry is just barely (or not at all) significant
#  also note that with the censored data statistics, we generated values (based on a distribution) for
#  censored data points, without consideration for wet/dry/season dynamics.  This makes me less confident
#  about keeping "dry" in the equation.
#myForm <- my.formulas.m4[[2]] 
myForm <- update(my.formulas.m4[[2]], .~. -dry)

#these lines of code assess fit of this particular model in terms of COC vs. individual predictors, and predictor correlation
myModel <- lme(data=coc2, myForm, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
#myModel <- gls(data=coc2, myForm, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(myModel)
AIC(myModel)
plot.single.preds(myModel)
check.cor(myModel)
boxplots.resids2(myModel, residuals(myModel, type="normalized"), "X")

#formulas that are worth considering (single plots of predictors make sense)
Form4a <- formula(result ~ rain + summer + sqrt_traffic + devAge2) #my.formulas.m4[[2]]; AIC=796.6
Form4b <- formula(result ~ rain + summer + sqrt_CO2_road + devAge2) #my.formulas.m4[[5]]; AIC=799.7
Form4c <- formula(result ~ rain + summer + sqrt_traffic + nodev) #my.formulas.m4[[8]]; AIC=801.4


#####  Best fit Model4  ####

Form4 <- Form4a
Model4 <- lme(data=coc2, Form4, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4 <- residuals(object = Model4, type = "normalized")

#try adding interactions and see if they are significant
M4.sub1 <- update(Model4, . ~ . + rain:sqrt_traffic) #not significant (p=0.1844)
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, . ~ . + rain:devAge2)  #not significant (p=0.5384)
anova(Model4, M4.sub2)


Form4 <- Form4b
Model4 <- lme(data=coc2, Form4, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4 <- residuals(object = Model4, type = "normalized")

#try adding interactions and see if they are significant
M4.sub1 <- update(Model4, . ~ . + rain:sqrt_CO2_road) #not significant (p=0.7386)
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, . ~ . + rain:devAge2)  #not significant (p=0.5145)
anova(Model4, M4.sub2)
M4.sub3 <- update(Model4, . ~ . + sqrt_CO2_road:devAge2)  #not significant (p=0.2135)
anova(Model4, M4.sub3)


Form4 <- Form4c
Model4 <- lme(data=coc2, Form4, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4 <- residuals(object = Model4, type = "normalized")

#try adding interactions and see if they are significant
M4.sub1 <- update(Model4, . ~ . + rain:sqrt_traffic) #not significant (p=0.2144)
anova(Model4, M4.sub1)
M4.sub2 <- update(Model4, . ~ . + rain:nodev)  #not significant (p=0.3398)
anova(Model4, M4.sub2)


#best model
Form4 <- Form4a
Model4 <- lme(data=coc2, Form4, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))


#------ Compare models and plot ---------#
AIC(Model1, Model3, Model4)

AIC(Model1, Model3, Model4)[,2] - AIC(Model4)


#plot model predictions for each model (above)
if (run_exploratory_code==TRUE) {
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

Model4a <- lme(data=coc2, Form4a, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4a <- residuals(object = Model4a, type = "normalized")

# Model4b <- lme(data=coc2, Form4b, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# E4b <- residuals(object = Model4b, type = "normalized")

#check residuals for each model, to see if they tell us anything important
par(mfrow=c(2,2), mar=c(2,4,4,1), oma=c(0,0,0,0))
plot(coc2$location, E1, main="Null Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
plot(coc2$location, E3, main="Categorical Landuse Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
plot(coc2$location, E4a, main="Landscape Predictor Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
# plot(coc2$location, E4b, main="Model 4b: sqrt_CO2_road + devAge2", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
# axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
# abline(h=0, col="gray")



# Step 9: Refit models with REML and apply graphical model validation; check for homogeneity, 
#--------    normality and independence

#------ Model 3: land use + rain ---------#    ##### NOTE: for some reason, the acf is giving what I think are completely
#                                                          bogus plots of auto-correlation.  Looking at lattice plots from
#                                                          earlier in code, doesn't look like that level of autocorrelation
#                                                          could possibly be present.  Maybe something is off with the work
#                                                          we did to employ censtat methods?  Will look into this with non-censtat TKN


M3.final <- lme(data=coc2, result~landuse+rain+summer, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
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
#variogram to test for spatial correlation; 
Vario1 = variogram(E3.final ~ 1, mydata2)
plot(Vario1)  #if there is no spatial correlation, will see a horizontal bar of points at the top of the plot
Vario2 <- variogram(E3.final ~ 1, mydata2, alpha = c(0, 45, 90,135) )
plot(Vario2)  #do we see any different patterns in the different directions, or roughly the same pattern?
# for Nitrite-Nitrate, there may be spatial autocorrelation with Model3 (south end/ SW site seem biased)


#------ Model 4: no land use; up to 3 predictors + rain ---------#  ##### NOTE: for some reason, the acf is giving what I think are completely
#                                                          bogus plots of auto-correlation.  Looking at lattice plots from
#                                                          earlier in code, doesn't look like that level of autocorrelation
#                                                          could possibly be present.  Maybe something is off with the work
#                                                          we did to employ censtat methods?  Will look into this with non-censtat TKN


M4.final <- lme(data=coc2, Form4a, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
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


#-----------------#
#  Save Formulas  #
#-----------------#

#save important items with phosphorus-specific names
TKN.coc2 <- coc2
TKN.r1X <- r1X
TKN.vf1X <- vf1X
TKN.Form3 <- Form3
TKN.Form4 <- Form4a  #sqrt_traffic + devAge2
TKN.rain <- rain

save(TKN.coc2, TKN.r1X, TKN.vf1X, TKN.Form3, TKN.Form4, TKN.rain, file=here("..", "results", "Frequentist_Total Kjeldahl Nitrogen Models_censtat.RData"))


#-------------------------------#
#  Plot Predictor Coefficients  #
#-------------------------------#

#generate plots that are shown in Rmarkdown script
TKN.null <- gls(data = coc2, result ~ 1, method = "REML") 
TKN.M3 <- lme(data = coc2, Form3, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
TKN.M4 <- lme(data = coc2, Form4, random = r1X, method = "REML", weights = vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

TKN_models <- list(
  Null_Model = TKN.null,
  Categorical_Landuse_Model = TKN.M3,
  Landscape_Predictor_Model = TKN.M4
)

huxtablereg(TKN_models,
            single.row = TRUE, custom.model.names = names(TKN_models)) %>%
  set_bottom_border(1, -1, 0.4) %>%
  set_bold(1, -1, TRUE) 

theme_set(theme_sjplot())
plotreg(TKN_models, custom.title = "Regression Results, Total Kjeldahl Nitrogen", custom.model.names = names(TKN_models))
plot_models(TKN_models, title = "Regression Results, Total Kjeldahl Nitrogen", m.labels = names(TKN_models),legend.title = "Models", show.values = TRUE,show.intercept = TRUE)


# #------------------------------------#
# #  Censored Data Regressions vs LME  #     Is the method of dividing reporting limits by 2 acceptable for TKN?
# #------------------------------------#
# 
# #try a LME model using Regression on Order Statistics (ros) to generate missing values (NADA package)
# tkn.ros <- ros(obs=coc$result, censored=coc$nondetect_flag)  #I think this package fills in for ND values using a distribution 
# tkn.df <- as.data.frame(tkn.ros)  #make a df out of the ROS-generated data
# coc2.sorted <- coc2[order(coc2$result), ]  #sort the coc2 dataframe by result, to match tkn.df
# coc.ros <- cbind(coc2.sorted, tkn.df)
# Form4.ros <- formula(modeled ~ rain + summer + sqrt_CO2_road + devAge2)  #run mle model on the ROS-generated ND values
# M4.ros <- lme(data=coc.ros, Form4.ros, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# 
# #compare ROS model to 0.5*DL model -- both are LME
# summary(M4.ros)$coefficients$fixed
# summary(M4.final)$coefficients$fixed  #censored MLE model (ND values set at 0.5*detection limit)
# #model that utilizes ROS for filling in ND values is close to M4.final (MLE) model using 0.5DL, at least in landscape predictors
# 
# 
# #Try the CenReg function to generate a (linear) model
# coc <- coc %>%
#   mutate(summer = case_when(
#     season=="1" ~ "0",
#     season=="2" ~ "0",
#     season=="3" ~ "1",
#     season=="4" ~ "0")) %>%
#   mutate(summer=as.factor(summer))
# 
# #cenreg is a censored regression function, which accounts for the fact that data are censored; it does NOT use LME model (just linear)
# tkn.reg <- with(coc, cenreg(Cen(oconc, nondetect_flag) ~ daymet_14day_std + summer + sqrt_CO2_road + devAge2, 
#                             dist="lognormal"))
# summary(tkn.reg)$coefficients  ###  NOT A MIXED EFFECTS MODEL!!  ###
# summary(M4.ros)$coefficients$fixed  #mixed effects model using regression on order statistics (ROS)
# summary(M4.final)$coefficients$fixed  #censored MLE model (ND values set at 0.5*detection limit)
# #all of these have similar enough values of landscape predictors and intercept that I think we can safely rely on 
# #   the linear mixed effects model for picking our best fit predictors
# 
# 
# #-----------------------#
# #  Test ROS for Form4a  #     Results are similar to those for Form4b (above)
# #-----------------------#
# 
# #try a LME model using Regression on Order Statistics (ros) to generate missing values (NADA package)
# Form4a.ros <- formula(modeled ~ rain + summer + sqrt_traffic + devAge2)  #run mle model on the ROS-generated ND values
# M4a.ros <- lme(data=coc.ros, Form4a.ros, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# 
# M4a.final <- lme(data=coc2, Form4a, random = r1X, method="REML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
# 
# #compare ROS model to 0.5*DL model -- both are LME
# summary(M4a.ros)$coefficients$fixed
# summary(M4a.final)$coefficients$fixed  #censored MLE model (ND values set at 0.5*detection limit)
# #model that utilizes ROS for filling in ND values is close to M4.final (MLE) model using 0.5DL, at least in landscape predictors
# 
# AIC(M4a.final, M4a.ros, M4.final, M4.ros)  #interesting -- the two ROS models have much lower AIC values than the regular MLE models
# 
# 
# #---------------------------------------------#
# #  Interesting Plotting Package to Play With  #
# #---------------------------------------------#
# 
# install.packages("visreg")
# library(visreg)
# 
# par(mfrow=c(2,2), mar=c(4,4,4,2))
# visreg(Model4, "sqrt_CO2_road", ylab="ln(TKN)")
# visreg(Model4, "devAge2", ylab="ln(TKN)")
# visreg(Model4, "summer", ylab="ln(TKN)")
# visreg(Model4, "rain", ylab="ln(TKN)")
# 
# par(mfrow=c(2,2), mar=c(4,4,4,2))
# visreg(M4.ros, "sqrt_CO2_road", ylab="ln(TKN)", main="ROS model")
# #visreg(tkn.reg, "sqrt_CO2_road", ylab="ln(TKN)", main="Censored Regression model")
# visreg(Model4, "sqrt_CO2_road", ylab="ln(TKN)", main="LME model")
# 
