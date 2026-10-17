# This script calls functions to plot Bayesian model outputs for each COC.  See script "Bayesian
#   Model Plot Functions.R" for functions.

# Author: Eva Dusek Jennings
# Revised: Nov 29, 2022
# Revised: Apr 22, 2025 -- updated Zinc to use greenery model rather than sqrt_traffic + paved
# Revised: Oct 21, 2025 -- updated Zinc to use not_greenBE model
#-------------------------------------------------------------------------------------------


source("scripts/Bayesian Model Plot Functions.R")  #script with plotting functions for Bayesian figures
# source("scripts/Bayesian Credibility Interval Function.R")  #script with credibility interval (& CI plotting) functions

library(cowplot)
library(ggpubr)

theme_set(theme_tidybayes() + panel_border() + grids())

#read in standardized values; this will be used in plots of prediction intervals vs raw data
std.vals <- read.csv(file="processed_data/spatial_predictor_standardization_values.csv")
colnames(std.vals) <- c("predictor", "mean", "sd")

#read in backtransform powers for spatial predictors, used for plotting raw spatial predictor values
#xform.pwr <- read.csv("processed_data/spatial_predictor_backtransform_power.csv", header=TRUE)

# #read (raw) watershed reductions (PSAU version) csv from stormwaterheatmap GitHub page
# wr_psau <- read.csv("https://raw.githubusercontent.com/stormwaterheatmap/TNC_Stormwater_Statistics/main/data/watershed_reductions/psau-predictors.csv")
# wr <- wr_psau[-which(wr_psau$AU_ID==0),] %>%
#   mutate(sqrt_traffic=traffic)   #### get rid of this once we sort out what "traffic" actually means!


#----------------#
#  Total Copper  #
#----------------#

#load data 
load(file="results/Bayesian_Copper.RData")  #Bayesian model
load(file="results/Frequentist_Copper Models.RData")  #frequentist mixed effects model, for coc2 dataframe
#load(file="results/Copper_80_CI_wr_psau.Rdata")  #credibility intervals for all Puget Sound watersheds

###   PLOT 1   ###
# posterior densities and chain traces for important predictors; max of 5 predictors per page
#get_variables(Cu.brm)
plot(Cu.brm, variable=c("b_Intercept", "b_impervious", "b_ln_traffic_k100", "b_rain14:impervious"), bins=50)  #recognizable covariates
plot(Cu.brm, variable=c("b_rain14", "b_summer1", "b_dry"), bins=50)  #recognizable covariates
plot(Cu.brm, variable=c("b_sigma_Intercept", "sd_agency__Intercept", "sd_agency:location__Intercept",
                        "sd_location__sigma_Intercept", "nu"), bins=50)   #covariates dealing with variability

Cu.lim <- c(-1.4,5.1)

###   PLOT 2   ###
# Observed vs Estimated - colored dots (by location) with thin grey lines = 95% interval                 
obs.vs.pred.1(Cu.brm, Cu.coc2, "Cu") + xlim(Cu.lim) + ylim(Cu.lim)

###   PLOT 3   ###
#Observed vs Estimated plots for multiple predictors, with colors showing a selected predictor.  
p1 <- obsPredPlot(Cu.brm, Cu.coc2, "Copper", "impervious", "Reds") + labs(color="imper-\nvious") + xlim(Cu.lim) + ylim(Cu.lim)
p2 <- obsPredPlot(Cu.brm, Cu.coc2, "Copper", "ln_traffic_k100", "Greens") + labs(color="ln_\ntraffic_\nk100") + xlim(Cu.lim) + ylim(Cu.lim)
p3 <- obsPredPlot(Cu.brm, Cu.coc2, "Copper", "rain14", "Blues") + xlim(Cu.lim) + ylim(Cu.lim)
p4 <- obsPredPlot(Cu.brm, Cu.coc2, "Copper", "dry", "Oranges") + xlim(Cu.lim) + ylim(Cu.lim)
#p5 <- obsPredPlot2(Cu.brm, Cu.coc2, "Copper", "summer", paletteCol=c("#fdbe85", "#d94701")) + xlim(Cu.lim) + ylim(Cu.lim)  #oranges for summer palette
p5 <- obsPredPlot2(Cu.brm, Cu.coc2, "Copper", "summer", paletteCol=c("#2d708e", "goldenrod2")) + xlim(Cu.lim) + ylim(Cu.lim)  #oranges for summer palette
# grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
lay <- rbind(c(1,1,2,2,NA,NA),
             c(3,3,4,4,5,5))
grid.arrange(p1, p2, p3, p4, p5, layout_matrix=lay)

###   PLOT 4   ###
#Credible Intervals vs Raw Data plots  
credInt.vs.rawPreds.2preds(Cu.coc2, Cu.brm, "Cu", preds=c("impervious", "ln_traffic_k100"), box.widths=c(0.03, 0.2))

###   PLOT 5   ###
#Global Intercept (+/- 95% CI) with Distributions of Location-Specific Intercepts
plotIntercepts.global.location(Cu.brm, Cu.coc2, xlims=c(0.6, 3.4)) 

# ###   PLOT 6   ###
# #plot transformed & standardized median values (black) with upper/ lower 80% CI's 
# plotBayesianCIs(Cu.CI, "Copper", c("sqrt_traffic", "devAge2"))
# 
# ###   PLOT 7   ###
# #plot raw median values (black) with upper/ lower 80% CI's
# plotBayesianCIs.raw(Cu.CI, "Copper", c("sqrt_traffic", "devAge2"))

#remove large objects from memory
rm(Cu.brm, Cu.CI, Cu.coc2) 


#-------#
#  TSS  #
#-------#

#load data 
load(file="results/Bayesian_TSS.RData")  #Bayesian model
load(file="results/Frequentist_TSS Models.RData")  #frequentist mixed effects model, for coc2 dataframe
#load(file="results/TSS_80_CI_wr_psau.Rdata")  #credibility intervals for all Puget Sound watersheds

###   PLOT 1   ###
# posterior densities and chain traces for important predictors; max of 5 predictors per page
#get_variables(TSS.brm)
plot(TSS.brm, variable=c("b_Intercept", "b_oldDev", "b_sqrt_CO2_road", "b_rain1:oldDev"), bins=50)  #recognizable covariates
plot(TSS.brm, variable=c("b_rain1", "b_rain21", "b_spring1"), bins=50)  #recognizable covariates
plot(TSS.brm, variable=c("b_sigma_Intercept", "sd_agency__Intercept", "sd_agency:location__Intercept",
                        "sd_location__sigma_Intercept", "nu"), bins=50)   #covariates dealing with variability

TSS.lim <- c(5.5, 13.1)

###   PLOT 2   ###
# Observed vs Estimated - colored dots (by location) with thin grey lines = 95% interval                 
obs.vs.pred.1(TSS.brm, TSS.coc2, "TSS") + xlim(TSS.lim) + ylim(TSS.lim)

###   PLOT 3   ###
#Observed vs Estimated plots for multiple predictors, with colors showing a selected predictor.  
p1 <- obsPredPlot(TSS.brm, TSS.coc2, "TSS", "sqrt_CO2_road", "Reds") + labs(color="sqrt_\nCO2_\nroad") + xlim(TSS.lim) + ylim(TSS.lim)
p2 <- obsPredPlot(TSS.brm, TSS.coc2, "TSS", "oldDev", "Greens") + xlim(TSS.lim) + ylim(TSS.lim)
p3 <- obsPredPlot(TSS.brm, TSS.coc2, "TSS", "rain1", "Blues") + xlim(TSS.lim) + ylim(TSS.lim)
p4 <- obsPredPlot(TSS.brm, TSS.coc2, "TSS", "rain21", "Purples") + xlim(TSS.lim) + ylim(TSS.lim)
#p5 <- obsPredPlot2(TSS.brm, TSS.coc2, "TSS", "spring", paletteCol=c("#2d708e", "#c2df23")) + xlim(TSS.lim) + ylim(TSS.lim)  #oranges for summer palette
p5 <- obsPredPlot2(TSS.brm, TSS.coc2, "TSS", "spring", paletteCol=c("#25858e", "#c2df23")) + xlim(TSS.lim) + ylim(TSS.lim)  #oranges for summer palette
#grid.arrange(p1, p2, p3, p4, p5, nrow=2, ncol=3)
lay <- rbind(c(1,1,2,2,NA,NA),
             c(3,3,4,4,5,5))
grid.arrange(p1, p2, p3, p4, p5, layout_matrix=lay)


###   PLOT 4   ###
#Credible Intervals vs Raw Data plots  
credInt.vs.rawPreds.2preds(TSS.coc2, TSS.brm, "TSS", preds=c("sqrt_CO2_road", "oldDev"), box.widths=c(1.4, 0.043))

###   PLOT 5   ###
#Global Intercept (+/- 95% CI) with Distributions of Location-Specific Intercepts
plotIntercepts.global.location(TSS.brm, TSS.coc2, xlims=c(8.8,12))

# ###   PLOT 6   ###
# #plot transformed & standardized median values (black) with upper/ lower 80% CI's 
# plotBayesianCIs(TSS.CI, "TSS", c("sqrt_traffic", "devAge2"))
# 
# ###   PLOT 7   ###
# #plot raw median values (black) with upper/ lower 80% CI's
# plotBayesianCIs.raw(TSS.CI, "TSS", c("sqrt_traffic", "devAge2"))

#remove large objects from memory
rm(TSS.brm, TSS.CI, TSS.coc2) 


#--------------------#
#  Total Phosphorus  #
#--------------------#

#load data 
load(file="results/Bayesian_Phosphorus.RData")  #Bayesian model
load(file="results/Frequentist_Total Phosphorus Models.RData")  #frequentist mixed effects model, for coc2 dataframe
#load(file="../results/Phosphorus_80_CI_wr_psau.Rdata")  #credibility intervals for all Puget Sound watersheds

###   PLOT 1   ###
# posterior densities and chain traces for important predictors; max of 5 predictors per page
#get_variables(P.brm)
plot(P.brm, variable=c("b_Intercept", "b_oldDev", "b_rain14:oldDev"), bins=50)  #recognizable covariates
plot(P.brm, variable=c("b_rain14", "b_rain1", "b_summer1"), bins=50)  #recognizable covariates
plot(P.brm, variable=c("b_sigma_Intercept", "sd_agency__Intercept", "sd_agency:location__Intercept",
                        "sd_location__sigma_Intercept", "nu"), bins=50)   #covariates dealing with variability

P.lim <- c(0.9, 7)

###   PLOT 2   ###
# Observed vs Estimated - colored dots (by location) with thin grey lines = 95% interval                 
obs.vs.pred.1(P.brm, P.coc2, "P") + xlim(P.lim) + ylim(P.lim)

###   PLOT 3   ###
#Observed vs Estimated plots for multiple predictors, with colors showing a selected predictor.  
p1 <- obsPredPlot(P.brm, P.coc2, "P", "oldDev", "Reds") + xlim(P.lim) + ylim(P.lim)
# p.blank <- textGrob("")
p2 <- obsPredPlot(P.brm, P.coc2, "P", "rain1", "Blues") + xlim(P.lim) + ylim(P.lim)
p3 <- obsPredPlot(P.brm, P.coc2, "P", "rain14", "Purples") + xlim(P.lim) + ylim(P.lim)
p4 <- obsPredPlot2(P.brm, P.coc2, "P", "summer", paletteCol=c("#2d708e", "goldenrod2")) + xlim(P.lim) + ylim(P.lim)  #oranges for summer palette
lay <- rbind(c(1,1,2,2),
             c(3,3,4,4))
grid.arrange(p1, p2, p3, p4, layout_matrix=lay)

###   PLOT 4   ###
#Credible Intervals vs Raw Data plots  
credInt.vs.rawPreds.1pred(P.coc2, P.brm, "P", preds="oldDev", box.widths= 0.07)

###   PLOT 5   ###
#Global Intercept (+/- 95% CI) with Distributions of Location-Specific Intercepts
plotIntercepts.global.location(P.brm, P.coc2, xlims=c(2.7,5.9)) 

# ###   PLOT 6   ###
# #plot transformed & standardized median values (black) with upper/ lower 80% CI's 
# plotBayesianCIs(P.CI, "Phosphorus", "sqrt_CO2_road")
# 
# ###   PLOT 7   ###
# #plot raw median values (black) with upper/ lower 80% CI's
# plotBayesianCIs.raw(P.CI, "Phosphorus", "sqrt_CO2_road")

#remove large objects from memory
rm(P.brm, P.CI, P.coc2) 


#--------------#
#  Total Zinc  #     #### NOTE: waiting on update to psau_wr (including sqrt_CO2_tot) so that Zn CI's can be run
#--------------#

#load data 
load(file="results/Bayesian_TotalZinc.RData")  #Bayesian model
load(file="results/Frequentist_Total Zinc Models.RData")  #frequentist mixed effects model, for coc2 dataframe

###   PLOT 1   ###
# posterior densities and chain traces for important predictors; max of 5 predictors per page
#get_variables(totZn.brm)
plot(totZn.brm, variable=c("b_Intercept", "b_sqrt_CO2_transport", "b_oldDev", "b_rain14:oldDev"), bins=50)
plot(totZn.brm, variable=c("b_rain14", "b_summer1"), bins=50)
plot(totZn.brm, variable=c("b_sigma_Intercept", "sd_agency__Intercept", "sd_agency:location__Intercept",
                         "sd_location__sigma_Intercept", "nu"), bins=50)   #covariates dealing with variability

Zn.lim <- c(-1.6, 6.6)

###   PLOT 2   ###
# Observed vs Estimated - colored dots (by location) with thin grey lines = 95% interval                 
obs.vs.pred.1(totZn.brm, totZn.coc2, "Zn") + xlim(Zn.lim) + ylim(Zn.lim)

###   PLOT 3   ###
#Observed vs Predicted plots for multiple predictors, with colors showing a selected predictor.  
p1 <- obsPredPlot(totZn.brm, totZn.coc2, "Zn", "oldDev", "Reds") + xlim(Zn.lim) + ylim(Zn.lim)
p2 <- obsPredPlot(totZn.brm, totZn.coc2, "Zn", "sqrt_CO2_transport", "Greens") + labs(color="sqrt_\nCO2_\ntransport") + xlim(Zn.lim) + ylim(Zn.lim)
p3 <- obsPredPlot(totZn.brm, totZn.coc2, "Zn", "rain14", "Blues") + xlim(Zn.lim) + ylim(Zn.lim)
p4 <- obsPredPlot2(totZn.brm, totZn.coc2, "Zn", "summer", paletteCol=c("#2d708e", "goldenrod2")) + xlim(Zn.lim) + ylim(Zn.lim)  #oranges for summer palette
lay <- rbind(c(1,1,2,2),
             c(3,3,4,4))
grid.arrange(p1, p2, p3, p4, layout_matrix=lay)

###   PLOT 4   ###
#Credible Intervals vs Raw Data plots  
credInt.vs.rawPreds.2preds(totZn.coc2, totZn.brm, "Total Zinc", preds=c("oldDev", "sqrt_CO2_transport"), box.widths=c(.05,1.5))

###   PLOT 5   ###
#Global Intercept (+/- 95% CI) with Distributions of Location-Specific Intercepts
plotIntercepts.global.location(totZn.brm, totZn.coc2, xlims=c(2.3, 5.5)) 

# ###   PLOT 6   ###
# #plot transformed & standardized median values (black) with upper/ lower 80% CI's 
# plotBayesianCIs(totZn.CI, "totZn", c("sqrt_traffic", "devAge2"))
# 

# ###   PLOT 7   ###
# #plot raw median values (black) with upper/ lower 80% CI's
# plotBayesianCIs.raw(totZn.CI, "totZn", c("sqrt_traffic", "devAge2"))

#remove large objects from memory
rm(totZn.brm, #totZn.CI, 
   totZn.coc2)


#---------------------------#
#  Total Kjeldahl Nitrogen  #        #######  Note: the two censored data points (detection limit = 30ppm) seem suspicious, given that
#---------------------------#                       all other samples at those locations have much higher values.  both were collected on 11/29/2010

#best Bayesian model uses censored methods WITHIN Bayesian context, and has Student t-distribution for errors, and
#  agency as a variance covariate

#load data 
load(file="results/Bayesian_TotalKjeldahlNitrogen.RData")  #Bayesian model
load(file="results/Frequentist_Total Kjeldahl Nitrogen Models_censtat.RData")  #frequentist mixed effects model, for coc2 dataframe
#load(file="../results/TotalKjeldahlNitrogen_80_CI_wr_psau.Rdata")  #credibility intervals for all Puget Sound watersheds



# ### this section of code is meant to simulate imputed point values for the censored data (posterior predictions)
# # I'm not sure I want to include the censored data here.  
# this.coc2 <- TKN.coc2
# this.brm <- TKN.brm
# 
# # 1. Filter your dataset to only the censored observations
# censored_data <- this.coc2 %>% 
#   filter(cen1 == "left" | cen1 == 1 | cen1 == TRUE)
# 
# # 2. Predict full posterior paths for these rows (including residual variance)
# simulated_draws <- this.brm %>%
#   predicted_draws(newdata = censored_data) %>%
#   # 3. Truncate the draws at the known detection limit (log of your oconc threshold)
#   filter(.prediction <= log(oconc)) 
# 
# # 4. Summarize the simulated distributions to get the "assigned" point values
# censored_imputed_summary <- simulated_draws %>%
#   group_by(across(everything())) %>% # Keeps your original row identifiers
#   summarize(
#     Imputed_Estimate = median(.prediction),
#     Imputed_Q2.5    = quantile(.prediction, 0.025),
#     Imputed_Q97.5   = quantile(.prediction, 0.975),
#     .groups = "drop"
#   )
# 
# 
# ###



###   PLOT 1   ###
# posterior densities and chain traces for important predictors; max of 5 predictors per page
#get_variables(TKN.brm)
plot(TKN.brm, variable=c("b_Intercept", "b_sqrt_CO2_res", "b_oldDev", "b_rain1:oldDev"), bins=50)  #recognizable covariates
plot(TKN.brm, variable=c("b_rain14", "b_rain1", "b_summer1"), bins=50)  #recognizable covariates
plot(TKN.brm, variable=c("b_sigma_Intercept", "sd_agency__Intercept", "sd_agency:location__Intercept",
                         "sd_agency__sigma_Intercept", "nu"), bins=50)   #covariates dealing with variability

TKN.lim <- c(4.5, 9.2)

###   PLOT 2   ###
# Observed vs Estimated - colored dots (by location) with thin grey lines = 95% interval; censored data shown in triangles
#obs.vs.pred.cen(TKN.brm, TKN.coc2, "TKN")
obs.vs.pred.DL(TKN.brm, TKN.coc2, "TKN") + xlim(TKN.lim) + ylim(TKN.lim)    #this one is probably a better choice, as it doesn't mix frequentist ROS-obtained values with 
                                            #  the Bayesian-generated censored values (which we may not be able to extract?)

###   PLOT 3   ###
#Observed vs Predicted plots for multiple predictors, with colors showing a selected predictor.
#  uncensored points are shown as filled circles, while censored points are shown as filled triangles.
#  censored points are shown at the detection limit
p1 <- obsPredPlot.cen(TKN.brm, TKN.coc2, "TKN", "sqrt_CO2_res", "Reds") + xlim(TKN.lim) + ylim(TKN.lim)+ labs(color="sqrt_\nCO2_\nres") 
p2 <- obsPredPlot.cen(TKN.brm, TKN.coc2, "TKN", "oldDev", "Greens") + xlim(TKN.lim) + ylim(TKN.lim)
p3 <- obsPredPlot.cen(TKN.brm, TKN.coc2, "TKN", "rain1", "Blues") + xlim(TKN.lim) + ylim(TKN.lim)
p4 <- obsPredPlot.cen(TKN.brm, TKN.coc2, "TKN", "rain14", "Purples") + xlim(TKN.lim) + ylim(TKN.lim)
p5 <- obsPredPlot2.cen(TKN.brm, TKN.coc2, "TKN", "summer", paletteCol=c("#2d708e", "goldenrod2")) + xlim(TKN.lim) + ylim(TKN.lim)  #oranges for summer palette
lay <- rbind(c(1,1,2,2,NA,NA),
             c(3,3,4,4,5,5))
grid.arrange(p1, p2, p3, p4, p5, layout_matrix=lay)

###   PLOT 4   ###
#Credible Intervals vs Raw Data plots  
credInt.vs.rawPreds.2preds(TKN.coc2, TKN.brm, "TKN", preds=c("sqrt_CO2_res", "oldDev"), box.widths=c(1, 0.07))

###   PLOT 5   ###
#Global Intercept (+/- 95% CI) with Distributions of Location-Specific Intercepts
plotIntercepts.global.location(TKN.brm, TKN.coc2, xlims=c(5.5,7.5))

# ###   PLOT 6   ###
# #plot transformed & standardized median values (black) with upper/ lower 80% CI's 
# plotBayesianCIs(TKN.CI, "TKN", c("sqrt_traffic", "devAge2"))
# 
# ###   PLOT 7   ###
# #plot raw median values (black) with upper/ lower 80% CI's
# plotBayesianCIs.raw(TKN.CI, "TKN", c("sqrt_traffic", "devAge2"))

#remove large objects from memory
rm(TKN.brm, TKN.CI, TKN.coc2) 







summary(Cu.brm)
summary(TSS.brm)
summary(P.brm)
summary(totZn.brm)
summary(TKN.brm)


