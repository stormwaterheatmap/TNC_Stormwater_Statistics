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
# Jul 7, 2022 -- sqrt(traffic) used instead of traffic; random effects = agency/location
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
library(stringr)
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
colors_location <- colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)]

if (run_exploratory_code==TRUE) {
  #Cleveland Dotplot - use to detect violations of homogeneity: We are looking for whether the spread of data values
  #   differs between sampling locations (or agencies).  If so, it indicates heterogeneity, and that there may be problems 
  #   with violation ofhomogeneity in a linear regression model applied on these data.  We're also looking for outliers.
  par(mfrow = c(1, 1))
  dotchart(coc$result, groups = coc$loc, pch = 19, col = colors_agency[as.numeric(coc$agency)],
           xlab = "concentration", main = paste("Cleveland Dotplot:", this_param_short))
  
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
  #for copper, AFTER removal of PIE_HDR & PIE_COM, most important predictors are:
  #  intURB, totRES, nodev, sqrt_traffic, pm25_na, CO2_tot, CO2_com, CO2_road, devAge, roof_intURB, roof_intURB_IND
  #  daymet14, 21 or 28day (all -0.2)
}

#-----------------------------------------------------------------------------------------------------#
#  Plot COC vs various predictors: look for potential predictors to model & sources of heterogeneity  #
#-----------------------------------------------------------------------------------------------------#

lp_plots()
pr_plots()
tlp_plots()
mp_plots()

#best predictors for this COC; make sure to only have one version (transformed or not) of each predictor!
best_predictors <- c("intURB", "intURB_IND", "totRES", "grass", "greenery", "impervious", "nodev",
                     "sqrt_traffic", "sqrt_popn", "pm25_na", "sqrt_CO2_tot", "sqrt_CO2_com", "sqrt_CO2_road", 
                     "sqrt_CO2_nonroad", "devAge2", "roof_intURB_IND")

pred_i <- which(predictors %in% best_predictors)
lp_plots(pred_i)

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

  #plot(fitted(myModel), myE, main="residuals", xlab="", ylab="", xaxt="n", yaxt="n", col="gray", pch=16)
  
  
  
  par(mfrow=c(1,1), mar=c(4,4,4,1))
  plot(fitted(M.gls1X), E.1X, main="", xlab="fitted values", xaxt="s", ylab="normalized residuals", col="gray", pch=16)
    

  #Step 3:  Choose a variance structure (if there was heterogeneity)
  #-------  for selecting random structure, use REML to compare - its better at capturing random structure
  #         REML estimates the random effects by considering linear combinations of the data that remove the 
  #         fixed effects. If these fixed effects are changed, the likelihoods of the two models will not be directly comparable.
  #         Compare AIC for the various beyond-optimal models with different variance structures.
  #         Alternately, can use likelihood ratio test (anova) to compare nested models
  
  M.gls2X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|agency))
  M.gls3X <- gls(FormX, data=coc2, method="REML", weights=varIdent(form= ~1|location))

  anova(M.gls1X, M.gls2X, M.gls3X)
  #for copper: AIC is best for location (#3); BIC best for just agency (#2)
  #  try out both in Bayesian model
  
  #residual plots for best fit models -- look for homogeneity of residuals
  E.2X <- resid(object = M.gls2X, type = "normalized")
  boxplots.resids2(M.gls2X, E.2X, "X")

  E.3X <- resid(object = M.gls3X, type = "normalized")
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

#variance function for best fit models so far
vf1X <- varIdent(form= ~1|location)

# Steps 4-6: Find the proper random effects structure; look for temporal and spatial autocorrelation
#-----------  note: look for spatial correlation AFTER setting random effects!

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
  M1.lme2X <- lme(data=coc2, FormX, random = ~1|agency/location, method="REML", weights=varIdent(form= ~1|agency),
                  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  
  anova(M.vf1X, M1.lme1X, M1.lme2X)  #the likelihood ratio test p-value compares the two nested models, and whether they are significantly different
  #For copper, the model with random intercept is significantly better than the model with no random effects
}

#random structure for best fit model so far
r1X <- formula( ~ 1 | agency/location)

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

preds.1 <- best_predictors

if (run_exploratory_code==TRUE) {
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
  gg <- rep("black", length(ee))  ### NOTE: the plot with 1, 2 or 3 predictors doesn't work!  it refers to dd, and ee no longer matches it!
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

  save(my.formulas.m4, my.aics.m4, file=(here("..", "scripts", "formulas_and_aics", "total_copper_m4.RData")))
}

if (run_exploratory_code==FALSE) {
  load(file=here("..", "scripts", "formulas_and_aics", "total_copper_m4.RData"))
}

# #identify which (ranked) formulas have "rain" in them
# aa <- grep("rain", my.formulas.m4)  #ordered list (by AIC value) of all formulae that include rain
# my.aics.m4[aa]
# my.formulas.m4[aa]
# 
# #identify which (ranked) formulas have 2 or less predictors
# bb <- as.character(my.formulas.m4[aa])
# formulas.m4.two <- my.formulas.m4[aa][which(lengths(gregexpr("\\W+", bb))< 4)]
# aics.m4.two <- my.aics.m4[aa][which(lengths(gregexpr("\\W+", bb))< 4)]
# FA.two <- cbind(as.character(formulas.m4.two), round(aics.m4.two, 1))



#--------------------------------------------#
#  Quantiles List for Plotting Interactions  #
#--------------------------------------------#

#list of quantiles for all landscape predictors
qList <- list(
  totRES=quantile(coc2$totRES, probs=c(0,.25,.50,.75,1)),
  sqrt_traffic=quantile(coc2$sqrt_traffic, probs=c(0,.25,.50,.75,1)),
  sqrt_CO2_road=quantile(coc2$sqrt_CO2_road, probs=c(0,.25,.50,.75,1)),
  pm25_na=quantile(coc2$pm25_na, probs=c(0,.25,.50,.75,1)),
  devAge2=quantile(coc2$devAge2, probs=c(0,.25,.50,.75,1)),
  rain=quantile(coc2$rain, probs=c(0,.25,.50,.75,1))
)
  

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

### when sqrt_traffic replaces traffic predictor, Form4b (sqrt_traffic, pm25, devAge2) has an AIC of 787 
#     (compared to 752 for model with sqrt_traffic + totRES). It is the 22nd best model by AIC.  If 
#     pm25:rain interaction is added, it becomes the 12th best model with AIC=781.  Fits to data are good.


#run through the top formulas when only best predictors are on the table;
#  keep only formulas that make sense
myForm <- my.formulas.m4[[6]]
#myForm <- update(Form4b, .~. +rain:pm25_na)

#these lines of code assess fit of this particular model in terms of COC vs. individual predictors, and predictor correlation
myModel <- lme(data=coc2, myForm, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
summary(myModel)
AIC(myModel)
plot.single.preds(myModel)
check.cor(myModel)
boxplots.resids2(myModel, residuals(myModel, type="normalized"), "X")

#formulas that are worth considering (single plots of predictors make sense)
Form4a <- formula(result ~ rain + summer + sqrt_CO2_road + devAge2)  #AIC=699.3; model #2
Form4b <- formula(result ~ rain + summer + sqrt_traffic + devAge2 + pm25_na)  #AIC=700.5; model #4
Form4c <- formula(result ~ rain + summer + sqrt_traffic + devAge2)  #AIC=700.9; model #6

#####  Best fit Model4  ####   NOTE: SELECT ALTERNATE MODEL INSTEAD; SEE NOTE BELOW.
Model4a <- lme(data=coc2, Form4a, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4a <- residuals(object = Model4a, type = "normalized")

#try adding interactions and see if they are significant
M4a.sub1 <- update(Model4a, . ~ . + rain:sqrt_CO2_road) #not significant (p=0.4061)
anova(Model4a, M4.sub1)
M4a.sub2 <- update(Model4a, . ~ . + rain:devAge2)  #not significant (p=0.7118)
anova(Model4a, M4.sub2)
M4a.sub3 <- update(Model4a, . ~ . + sqrt_CO2_road:devAge2)  #not significant (p=0.2972)
anova(Model4a, M4.sub3)


#####  Alternate Model4  ####   SELECT THIS AS BEST MODEL -- COPPER IS MORE CLOSELY RELATED TO TRAFFIC THAN EMISSIONS!
Model4b <- lme(data=coc2, Form4c, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
E4b <- residuals(object = Model4b, type = "normalized")

#try adding interactions and see if they are significant -- NONE significant
M4b.sub1 <- update(Model4b, . ~ . + rain:sqrt_traffic) #not significant (p=0.2367)
anova(Model4b, M4b.sub1)
M4b.sub2 <- update(Model4b, . ~ . + rain:devAge2)  #not significant (p=0.6632)
anova(Model4b, M4b.sub2)
M4b.sub3 <- update(Model4b, . ~ . + sqrt_traffic:devAge2) #not significant (p=0.4392)
anova(Model4b, M4b.sub3)


#----------------------#

#------ Compare models and plot ---------#

#compare all of the models
AIC(Model1, Model3, Model4a, Model4b)
#obtain the delta AIC value
AIC(Model1, Model3, Model4a, Model4b)[, 2] - AIC(Model4a)


#best model #4
Model4 <- Model4b
Form4 <- Form4b

#plot model predictions for each model (above)
if (F) {
  #change from F to T if you want to plot these
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

#check residuals for each model, to see if they tell us anything important
par(mfrow=c(2,2), mar=c(2,4,4,1), oma=c(0,0,0,0))
plot(coc2$location, E1, main="Null Model", ylab="Residuals", xaxt="n", col=colors_location)
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
plot(coc2$location, E3, main="Categorical Landuse Model", ylab="Residuals", xaxt="n", col=colors_location)
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")
# plot(coc2$location, E4a, main="Model 4a", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
# axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
# abline(h=0, col="gray")
plot(coc2$location, E4b, main="Landscape Predictor Model", ylab="Residuals", xaxt="n", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
axis(side=1, at=c(2,3.8,5.2,7,10,13), labels=c("King", "Pie", "POT", "Sea", "Sno", "Tac"))
abline(h=0, col="gray")



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


#-----------------#
#  Save Formulas  #
#-----------------#

#save important items with copper-specific names
Cu.coc2 <- coc2
Cu.r1X <- r1X
Cu.vf1X <- vf1X
Cu.Form3 <- Form3
Cu.Form4 <- Form4
Cu.rain <- rain

save(Cu.coc2, Cu.r1X, Cu.vf1X, Cu.Form3, Cu.Form4, Cu.rain, file=here("..", "results", "Frequentist_Copper Models.RData"))


#-------------------------------#
#  Plot Predictor Coefficients  #
#-------------------------------#

#generate plots that are shown in Rmarkdown script
Cu.null <- gls(data = Cu.coc2, result ~ 1, method = "REML") 
Cu.M3 <- lme(data = Cu.coc2, result ~ landuse + rain + summer, random = Cu.r1X, method = "REML", weights = Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
Cu.M4 <- lme(data = Cu.coc2, Cu.Form4, random = Cu.r1X, method = "REML", weights = Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

Cu_models <- list(
  Null_Model = Cu.null,
  Categorical_Landuse_Model = Cu.M3,
  Landscape_Predictor_Model = Cu.M4)

huxtablereg(Cu_models,
            single.row = TRUE, custom.model.names = names(Cu_models)) %>%
  set_bottom_border(1, -1, 0.4) %>%
  set_bold(1, -1, TRUE) 

#plot parameter estimates for each model
theme_set(theme_sjplot())
plotreg(Cu_models, custom.title = "Regression Results, total Copper", custom.model.names = names(Cu_models))
plot_models(Cu_models,m.labels = names(Cu_models),legend.title = "Models", show.values = TRUE,show.intercept = TRUE)
#plot_models(Cu_models[-c(1)],m.labels = names(Cu_models[-c(1)]),legend.title = "Models", show.values = TRUE,show.intercept = TRUE)


