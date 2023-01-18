# This script contains functions for plotting of Bayesian models, which can be applied to each COC.

# Author:  Eva Dusek Jennings
# Date:    July 14, 2022
#-------------------------------------------------------------------------------------------------

library(tidyr)
library(dplyr)
#library(tidyverse)    # ggplot, dplyr, %>%, and friends -- DOESN'T WORK RIGHT NOW...
library(brms)         # Bayesian modeling through Stan
library(tidybayes)    # Manipulate Stan objects in a tidy way
library(broom)        # Convert model objects to data frames
library(forcats)
library(gridExtra)
library(ggplot2)
library(oce)          # approx3d
library(ggdist)


#--------------------------------#
#  Functions for Bayesian Plots  #
#--------------------------------#


#  Observed vs Predicted
#  https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html#spurious-associations

# observed vs predicted with black points and 95% credible interval lines around each point
obs.vs.pred.1 <- function(this.brm, this.coc2, this.chemical) {
  fitted(this.brm) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5
    as_tibble() %>%
    bind_cols(this.coc2) %>%
    ggplot(aes(x = result, y = Estimate)) +        
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
                   size = 1/4, color="gray") +
    geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
                       ymax = Estimate + Est.Error),
                   size = 5/8, color="darkgray") +
    geom_point(size = 1.5, alpha = 3/4) +
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Predicted ln(", this.chemical, ")", sep="") ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}


# observed vs predicted for censored methods models, with black points for uncensored, and pink points for censored.
#   95% credible interval lines shown in grey around each point
obs.vs.pred.cen <- function(this.brm, this.coc2, this.chemical) {
  fitted(this.brm) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5
    as_tibble() %>%
    bind_cols(this.coc2) %>%
    ggplot(aes(x = log(oconc), y = Estimate, group=cen)) +        #plot what came back from the lab; reporting limits for NDs in 2nd color!
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
                   size = 1/4, color="gray") +
    geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
                       ymax = Estimate + Est.Error),
                   size = 5/8, color="darkgray") +
    geom_point(size = 1.5, alpha = 3/4, aes(col=cen)) + 
    scale_color_manual(values=c("black", "deeppink2")) +
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Predicted ln(", this.chemical, ")", sep="") ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}


#function to plot observed vs predicted, with colors showing a selected predictor.  
#   Inputs: brm, coc2, chemical name, predictor name, color brewer palette
obsPredPlot <- function(this.brm, this.coc2, this.chemical, pred, paletteCol) {
  fitted(this.brm) %>%
    as_tibble() %>%
    bind_cols(this.coc2) %>%
    
    ggplot(aes(x = result, y = Estimate, color=get(pred))) +        #here, you can sub out <agency> for other predictors
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_point(size = 1.5, alpha = 3/4) +
    # geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
    #                size = 1/4) +
    # geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
    #                    ymax = Estimate + Est.Error),
    #                size = 1/2) +
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Predicted ln(", this.chemical, ")", sep=""), color=pred) +
    scale_color_distiller(palette = paletteCol, direction=1) +   #colors should be light for low values, dark for high values
    theme_bw() +
    theme(panel.grid = element_blank())
}


#function to plot observed vs predicted, with colors showing a selected discrete predictor.  
#  Inputs: predictor, color brewer palette
obsPredPlot2 <- function(this.brm, this.coc2, this.chemical, pred, paletteCol) {
  fitted(this.brm) %>%
    as_tibble() %>%
    bind_cols(this.coc2) %>%
    
    ggplot(aes(x = result, y = Estimate, color=get(pred))) +        #here, you can sub out <agency> for other predictors
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_point(size = 1.5, alpha = 3/4) +
    # geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
    #                size = 1/4) +
    # geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
    #                    ymax = Estimate + Est.Error),
    #                size = 1/2) +
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Predicted ln(", this.chemical, ")", sep=""), color=pred) +
    scale_color_manual(values=paletteCol) +                     #HEX color values are provided based on # of discrete values for pred 
    theme_bw() +
    theme(panel.grid = element_blank())
}


#function to plot observed vs predicted ("observed" censored points shown as the detection limit), with colors showing a selected predictor.  
#   Inputs: brm, coc2, chemical name, predictor name, color brewer palette
obsPredPlot.cen <- function(this.brm, this.coc2, this.chemical, pred, paletteCol) {

  # this.brm <- TKN.brm
  # this.coc2 <- TKN.coc2
  # this.chemical <- "TKN"
  # pred <- "sqrt_traffic"
  # paletteCol <- "Reds"
  
    fitted(this.brm) %>%
    as_tibble() %>%
    bind_cols(this.coc2) %>%
    
    ggplot(aes(x = log(oconc), y = Estimate, color=get(pred), group=cen)) +        #here, you can sub out <agency> for other predictors
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_point(size = 2, alpha = 3/4, aes(shape=cen)) +
    # geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
    #                size = 1/4) +
    # geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
    #                    ymax = Estimate + Est.Error),
    #                size = 1/2) +
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Predicted ln(", this.chemical, ")", sep=""), color=pred) +
    scale_color_distiller(palette = paletteCol, direction=1) +   #colors should be light for low values, dark for high values
    theme_bw() +
    theme(panel.grid = element_blank())
}


#function to plot observed vs predicted ("observed" censored points shown as the detection limit), with colors showing a selected discrete predictor.  
#  Inputs: predictor, color brewer palette
obsPredPlot2.cen <- function(this.brm, this.coc2, this.chemical, pred, paletteCol) {
  fitted(this.brm) %>%
    as_tibble() %>%
    bind_cols(this.coc2) %>%
    
    ggplot(aes(x = log(oconc), y = Estimate, color=get(pred), group=cen)) +        #here, you can sub out <agency> for other predictors
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_point(size = 2, alpha = 3/4, aes(shape=cen)) +
    # geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
    #                size = 1/4) +
    # geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
    #                    ymax = Estimate + Est.Error),
    #                size = 1/2) +
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Predicted ln(", this.chemical, ")", sep=""), color=pred) +
    scale_color_manual(values=paletteCol) +                     #HEX color values are provided based on # of discrete values for pred 
    theme_bw() +
    theme(panel.grid = element_blank())
}


#function to plot expectation of predicted draws; aka: expected values vs raw predictors, 
#  using the epred_draws function. This one works for TWO landscape predictors
predInt.vs.rawPreds.2preds <- function(this.coc2, this.brm, this.chemical, preds) {

  # ####### FOR TESTING ONLY!!
  # this.coc2 <- Cu.coc2
  # this.brm <- Cu.brm
  # this.chemical <- "Copper"
  # preds <- c("sqrt_traffic", "devAge2")
  # preds.xformPwr <- c(2, 0.5)  #to which power should we transform each predictor, to generate raw values? (0.5, 1, or 2)
  # preds.xlab <- c("Traffic (AADT, raw)", "Development Age (raw)" )  
  # #######  
  
  #throw an error if number of predictors is not 1 or 2
  if(length(preds)!=2) {
    stop("this function only works for 2 landscape predictors")
  }
  
  #epred_draws gives the "expected values of the outcome" - it focuses on the uncertainty in the model parameters, and not
  #   the individual-level residuals.  epred_draws accounts for the uncertainty from two sources:
  #       1.) uncertainty of fixed coefficients
  #       2.) uncertainty of the variance parameters of the groups (eg: sd(Intercept)) (aka: uncertainty in group-level effects);
  #           This does NOT include any elements that deal with heteroskedasticity, which only affect individual-level residuals.
  #  https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/
  #         (see yellow and blue curves ~ 1/4 or 1/3 of the way down, titled: Predicted draws and Expectation of predicted draws)
  
  #take 94 samples total (20% of all data we have) to represent historical (and hence feasible) rain or rain+summer options
  my.sample0 <- this.coc2[which(this.coc2$summer==0)[sample(length(which(this.coc2$summer==0)), 84)], c("rain", "summer")]
  my.sample1 <- this.coc2[which(this.coc2$summer==1)[sample(length(which(this.coc2$summer==1)), 10)], c("rain", "summer")]
  my.sample <- rbind(my.sample0, my.sample1)
  
  #generate simulated data
  newdata0 <- expand_grid(#location=NA, 
    rain=my.sample0$rain, 
    summer=0, 
    pred1=seq(min(this.coc2[[preds[1]]]), max(this.coc2[[preds[1]]]), length.out=8), 
    pred2=seq(min(this.coc2[[preds[2]]]), max(this.coc2[[preds[2]]]), length.out=8) )
  newdata1 <- expand_grid(#location=NA, 
    rain=my.sample1$rain, 
    summer=1, 
    pred1=seq(min(this.coc2[[preds[1]]]), max(this.coc2[[preds[1]]]), length.out=8), 
    pred2=seq(min(this.coc2[[preds[2]]]), max(this.coc2[[preds[2]]]), length.out=8) )
  newdata <- rbind(newdata0, newdata1)

  names(newdata) <- c("rain", "summer", preds)  #this is necessary so that this.brm recognizes the two predictors
    
  #obtain epred_draws
  abc <- epred_draws(this.brm, allow_new_levels=TRUE, newdata=newdata)

  #calculate raw predictor values for this COC, and for the epred draws; these are used in plotting ribbons and points, below
  this.coc2$pred1_raw <- (this.coc2[[preds[1]]] * std.vals$sd[which(std.vals$predictor==preds[1])] + std.vals$mean[which(std.vals$predictor==preds[1])] ) ^ 
    xform.pwr$power[which(xform.pwr$predictor==preds[1])]  
  this.coc2$pred2_raw <- (this.coc2[[preds[2]]] * std.vals$sd[which(std.vals$predictor==preds[2])] + std.vals$mean[which(std.vals$predictor==preds[2])] ) ^ 
    xform.pwr$power[which(xform.pwr$predictor==preds[2])]
  abc$pred1.raw <- (abc[[preds[1]]] * std.vals$sd[which(std.vals$predictor==preds[1])] + std.vals$mean[which(std.vals$predictor==preds[1])] ) ^ 
    xform.pwr$power[which(xform.pwr$predictor==preds[1])]
  abc$pred2.raw <- (abc[[preds[2]]] * std.vals$sd[which(std.vals$predictor==preds[2])] + std.vals$mean[which(std.vals$predictor==preds[2])] ) ^ 
    xform.pwr$power[which(xform.pwr$predictor==preds[2])]

  #plot RAW values for pred1 and pred2; these plots show 50%, 80% and 95% CI's; if median is desired as well, use stat_lineribbon instead of stat_ribbon!
  pred1.raw <- ggplot(abc, aes(x = pred1.raw, y = exp(.epred))) +
    stat_ribbon() +    #NOTE: stat_lineribbon would provide a median line (in black); this is not suitable for our dataset with 2 independent variables
    scale_fill_brewer(palette = "Reds") +
    labs(x = xform.pwr$rawName[which(xform.pwr$predictor==preds[1])], y = paste(this.chemical, "(ppm)"), fill = "Prediction Interval") +
    #  ylim(0, max(abc$.epred)) +
    geom_point(data=this.coc2, mapping=aes(x=pred1_raw, y=exp(result))) +
    theme(legend.position = "bottom")
  
  pred2.raw <- ggplot(abc, aes(x = pred2.raw, y = exp(.epred))) +
    stat_ribbon() +
    scale_fill_brewer(palette = "Purples") +
    labs(x = xform.pwr$rawName[which(xform.pwr$predictor==preds[2])] , y = paste(this.chemical, "(ppm)"), fill = "Prediction Interval") +
    #  ylim(0, max(abc$.epred)) +
    geom_point(data=this.coc2, mapping=aes(x=pred2_raw, y=exp(result))) +
    theme(legend.position = "bottom")
  
  grid.arrange(pred1.raw, pred2.raw, nrow=1, ncol=2)
}


#function to plot expectation of predicted draws; aka: expected values vs raw predictors, 
#  using the epred_draws function. This one works for ONE landscape predictor
predInt.vs.rawPreds.1pred <- function(this.coc2, this.brm, this.chemical, preds) {
  
  # ####### FOR TESTING ONLY!!  Note -- copper won't work here, as it has 2 predictors!
  # this.coc2 <- P.coc2
  # this.brm <- P.brm
  # this.chemical <- "Phosphorus"
  # preds <- "sqrt_CO2_road"
  # #######
  
  #throw an error if number of predictors is not 1 or 2
  if(length(preds)!=1) {
    stop("this function only works for 1 landscape predictor")
  }
  
  #epred_draws gives the "expected values of the outcome" - it focuses on the uncertainty in the model parameters, and not
  #   the individual-level residuals.  epred_draws accounts for the uncertainty from two sources:
  #       1.) uncertainty of fixed coefficients
  #       2.) uncertainty of the variance parameters of the groups (eg: sd(Intercept)) (aka: uncertainty in group-level effects);
  #           This does NOT include any elements that deal with heteroskedasticity, which only affect individual-level residuals.
  #  https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/
  
  #take 94 samples total (20% of all data we have) to represent historical (and hence feasible) rain or rain+summer options
  my.sample0 <- this.coc2[which(this.coc2$summer==0)[sample(length(which(this.coc2$summer==0)), 84)], c("rain", "summer")]
  my.sample1 <- this.coc2[which(this.coc2$summer==1)[sample(length(which(this.coc2$summer==1)), 10)], c("rain", "summer")]
  my.sample <- rbind(my.sample0, my.sample1)
  
  #generate simulated data
  newdata0 <- expand_grid(#location=NA, 
    rain=my.sample0$rain, 
    summer=0, 
    pred1=seq(min(this.coc2[[preds[1]]]), max(this.coc2[[preds[1]]]), length.out=8) )
  newdata1 <- expand_grid(#location=NA, 
    rain=my.sample1$rain, 
    summer=1, 
    pred1=seq(min(this.coc2[[preds[1]]]), max(this.coc2[[preds[1]]]), length.out=8) )
  newdata <- rbind(newdata0, newdata1)
  
  names(newdata) <- c("rain", "summer", preds)  #this is necessary so that this.brm recognizes the two predictors
  
  #obtain epred_draws
  abc <- epred_draws(this.brm, allow_new_levels=TRUE, newdata=newdata)
  
  #calculate raw predictor values for this COC, and for the epred draws; these are used in plotting ribbons and points, below
  this.coc2$pred1_raw <- (this.coc2[[preds[1]]] * std.vals$sd[which(std.vals$predictor==preds[1])] + std.vals$mean[which(std.vals$predictor==preds[1])] ) ^ 
    xform.pwr$power[which(xform.pwr$predictor==preds[1])]
  abc$pred1.raw <- (abc[[preds[1]]] * std.vals$sd[which(std.vals$predictor==preds[1])] + std.vals$mean[which(std.vals$predictor==preds[1])] ) ^ 
    xform.pwr$power[which(xform.pwr$predictor==preds[1])]

  #plot RAW values for pred1 and pred2; these plots show 50%, 80% and 95% CI's; if median is desired as well, use stat_lineribbon instead of stat_ribbon!
  pred1.raw <- ggplot(abc, aes(x = pred1.raw, y = exp(.epred))) +
    stat_ribbon() +    #NOTE: stat_lineribbon would provide a median line (in black); this is not suitable for our dataset with 2 independent variables
    scale_fill_brewer(palette = "Reds") +
    labs(x = xform.pwr$rawName[which(xform.pwr$predictor==preds[1])] , y = paste(this.chemical, "(ppm)"), fill = "Prediction Interval") +
    #  ylim(0, max(abc$.epred)) +
    geom_point(data=this.coc2, mapping=aes(x=pred1_raw, y=exp(result))) +
    theme(legend.position = "bottom")
  
  grid.arrange(pred1.raw, nrow=1, ncol=1)
}


#function to plot location-specific intercepts with densities, plus the population-level intercept and 95% credibility 
# intervals around the population-level intercept
#     Note that the 95% CI does NOT include variability in agencies or in location:agency; it is just the 95% CI
#          for the population-level intercept (ie: 95% of the time, the population-level intercept will fall between 
#          the 95% CI values)
plotIntercepts.global.location <- function(this.brm, xlims=c(0,10)) {
  this.brm %>%
    spread_draws(b_Intercept, r_agency[agency, Intercept], `r_agency:location`[agency_location,Intercept]) %>% 
    mutate(agency_location_mean = b_Intercept + `r_agency` + `r_agency:location`) %>%
    ggplot(aes(y = agency_location, x = agency_location_mean, 
               fill=stat(abs(x - fixef(this.brm)["Intercept", 1]) < (fixef(this.brm)["Intercept", 2]*2) )))  +
    xlim(xlims) +
    stat_halfeye() +
    labs(x="Intercept Value", y="", fill="Population-Level \nIntercept 95% CI") +
    geom_vline(xintercept = fixef(this.brm)["Intercept", 1], linetype = "solid", color="cadet blue4", alpha=0.5, size=2) +
    annotate("text", x=(fixef(this.brm)["Intercept", 1] - 0.1), y=7.2, label="population intercept", color="cadet blue4", angle=90, size=6) +
    geom_vline(xintercept = fixef(this.brm)["Intercept", c(3:4)], linetype = "dashed", color="black", size=1/4) +
    scale_fill_manual(values = c("gray80", "skyblue"))# +
}






