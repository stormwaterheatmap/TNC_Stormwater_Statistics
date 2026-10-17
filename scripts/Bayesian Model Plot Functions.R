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


#---------------------------------------#
#  Set Colors and Shapes for Locations  #
#---------------------------------------#

colors_agency <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
colors_location <- colors_agency[c(1,1,1,2,2,2,3,4,4,4,5,5,5,6,6,6)]

# #custom manual color scale for S8 fitting locations
# colors_location <- c("#E69F00", "#E69F00", "#E69F00",  #king county
#                      "#56B4E9", "#56B4E9", "#56B4E9", #pierce county
#                      "#009E73",  #port of tacoma
#                      "#0072B2", "#0072B2", "#0072B2",   #seattle
#                      "#D55E00", "#D55E00", "#D55E00",   #snohomish county
#                      "#CC79A7", "#CC79A7", "#CC79A7")   #tacoma
# #names(colors_location) <- 

#custom manual symbols for model validation locations
symbols_location <- c(16, 17, 15,  #king co
                      16, 17, 15,          #pierce county
                      16,          #port of tacoma
                      16, 17, 15,  #seattle
                      16, 17, 15,  #snohomish county
                      16, 17, 15)  #tacoma
#names(mySymbols) <- vp.std$Location




#------------------------------#
#  Observed vs Predicted Fxns  #
#------------------------------#


#  Observed vs Predicted
#  https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html#plotting-multivariate-posteriors.
#  go down to 5.1.3.3 POSTERIOR PREDICTION PLOTS

# observed vs predicted with points colored by location and gray 95% credible interval lines around each point
obs.vs.pred.1 <- function(this.brm, this.coc2, this.chemical) {
  fitted(this.brm) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5; note that fitted() ignores residual error
    as_tibble() %>%
    bind_cols(this.coc2) %>%
    ggplot(aes(x = result, y = Estimate, color=location, shape=location)) +       #color points by location (see custom colors below) 
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    # geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
    #                size = 1/4, color="gray") +
    # geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
    #                    ymax = Estimate + Est.Error),
    #                size = 5/8, color="darkgray") +
    geom_point(size = 2.5, alpha = 3/4) +
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Estimated ln(", this.chemical, ")", sep="") ) +
    # scale_color_manual(values=c("#FF9999", "#CC0033", "#993333", "#FF9933", "#FFCC00", "#99FF00", "#66CC33", "#006633", 
    #                                      "#99FFFF", "#0099FF", "#0033FF", "#CC99FF", "#9966FF", "#6600FF")) +  
    scale_shape_manual(name="location", values=symbols_location) +
    scale_color_manual(name="location", values=colors_location) +
    theme_bw() +
    theme(panel.grid = element_blank())
}


#show residuals on a plot, to see which ones are particularly far away from the obs vs pred line
#https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html#Plotting_multivariate_posterior_predictions


# # observed vs predicted for censored methods models, with colors showing location, and triangles showing censored data.
# #   95% credible interval lines shown in grey around each point
# #   "observed" value for censored data points is the Detection Limit! 
# obs.vs.pred.DL <- function(this.brm, this.coc2, this.chemical) {
#   fitted(this.brm) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5; note that fitted() ignores residual error
#     as_tibble() %>%
#     bind_cols(this.coc2) %>%
#     ggplot(aes(x = log(oconc), y = Estimate, color=location, shape=location, group=cen)) +       #color points by location (see custom colors below) 
#     geom_abline(linetype = 2, color = "grey50", size = .5) +
#     geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
#                    size = 1/4, color="gray") +
#     # geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
#     #                    ymax = Estimate + Est.Error),
#     #                size = 5/8, color="darkgray") +
#     geom_point(size = 2.5, alpha = 3/4, aes(shape=cen)) +
#     labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
#          y = paste("Estimated ln(", this.chemical, ")", sep="") ) +
#     # scale_color_manual(values=c("#FF9999", "#CC0033", "#993333", "#FF9933", "#FFCC00", "#99FF00", "#66CC33", "#006633", 
#     #                                      "#99FFFF", "#0099FF", "#0033FF", "#CC99FF", "#9966FF", "#6600FF")) +  
#     scale_shape_manual(name="location", values=symbols_location) +
#     scale_color_manual(name="location", values=colors_location) +
#     scale_linetype_manual(name="Censoring Status", values=c("TRUE"="dashed", "FALSE"="solid")) +
#     theme_bw() +
#     theme(panel.grid = element_blank())
# }



# # observed vs predicted for censored methods models, with colors showing location, and semi-transparent colors showing censored data.
# #   "observed" value for censored data points is the Detection Limit! 
obs.vs.pred.DL <- function(this.brm, this.coc2, this.chemical) {

  # this.brm <- TKN.brm
  # this.coc2 <- TKN.coc2
  # this.chemical <- "TKN"
  

  ### These symbols allow for filled vs unfilled shapes
  #custom manual symbols for model validation locations
  symbols_loc.cen  <- c(21, 22, 24,  #king co
                        21, 22, 24,          #pierce county
                        21,          #port of tacoma
                        21, 22, 24,  #seattle
                        21, 22, 24,  #snohomish county
                        21, 22, 24)  #tacoma
  #names(mySymbols) <- vp.std$Location

  fitted(this.brm) %>% 
    as_tibble() %>% 
    bind_cols(this.coc2) %>% 
    ggplot(aes(x = log(oconc), y = Estimate)) + 
    geom_abline(linetype = 2, color = "grey50", size = .5) + 
    # geom_linerange(aes(ymin = Q2.5, ymax = Q97.5), size = 1/4, color = "gray") + 
    # Use alpha to conditionally change transparency of the fill
    geom_point(aes(shape = location, 
                   color = location, 
                   fill = location, 
                   alpha = cen), 
               size = 2.5) + 
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Estimated ln(", this.chemical, ")", sep="")) + 
    # Standard manual mapping for your custom vectors
    scale_shape_manual(name = "location", values = symbols_loc.cen) + 
    scale_color_manual(name = "location", values = colors_location) + 
    scale_fill_manual(name = "location", values = colors_location) + 
    # Scale alpha to make FALSE completely opaque (1) and TRUE fully hollow (0)
    scale_alpha_manual(name = "Censoring Status", 
                       values = c("FALSE" = 1, "TRUE" = 0.3),
                       labels = c("FALSE" = "Uncensored", "TRUE" = "Censored")) +
    # Guide adjustment prevents the 'location' legend from inheriting the transparency
    # guides(
    #   shape = guide_legend(override.aes = list(alpha = 1)),
    #   color = guide_legend(override.aes = list(alpha = 1)),
    #   fill  = guide_legend(override.aes = list(alpha = 1)),
    #   alpha = guide_legend(override.aes = list(shape = 21, fill = "black", color = "black"))
    # ) +
    theme_bw() + 
    theme(panel.grid = element_blank())
}

# obs.vs.pred.DL <- function(this.brm, this.coc2, this.chemical) {
#   fitted(this.brm) %>% 
#     as_tibble() %>% 
#     bind_cols(this.coc2) %>% 
#     ggplot(aes(x = log(oconc), y = Estimate, color = location, shape = location, group = cen)) + 
#     geom_abline(linetype = 2, color = "grey50", size = .5) + 
#     # Map linetype to cen inside the linerange aesthetics
#     geom_linerange(aes(ymin = Q2.5, ymax = Q97.5, linetype = cen), 
#                    size = 1/4, color = "gray") + 
#     geom_point(size = 2.5, alpha = 3/4) + 
#     labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
#          y = paste("Estimated ln(", this.chemical, ")", sep="")) + 
#     scale_shape_manual(name = "location", values = symbols_location) + 
#     scale_color_manual(name = "location", values = colors_location) + 
#     # Add a scale to name and control your interval line styles for censoring
#     scale_linetype_manual(name = "Censoring Status", values = c(TRUE == "solid", FALSE == "dashed")) +
#     theme_bw() + 
#     theme(panel.grid = element_blank())
# }



# # same plot as above, except for the following:
# #   "observed" value for censored data points is the value that the frequentist ROS method produced 
# #   NOT the one that the Bayesian censored method produced.  :(
# obs.vs.pred.cen <- function(this.brm, this.coc2, this.chemical) {
#   fitted(this.brm) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5; note that fitted() ignores residual error
#     as_tibble() %>%
#     bind_cols(this.coc2) %>%
#     ggplot(aes(x = result, y = Estimate, color=location, group=cen)) +       #color points by location (see custom colors below) 
#     geom_abline(linetype = 2, color = "grey50", size = .5) +
#     geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
#                    size = 1/4, color="gray") +
#     # geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
#     #                    ymax = Estimate + Est.Error),
#     #                size = 5/8, color="darkgray") +
#     geom_point(size = 2.5, alpha = 3/4, aes(shape=cen)) +
#     labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
#          y = paste("Estimated ln(", this.chemical, ")", sep="") ) +
#     scale_color_manual(values=c("#FF9999", "#CC0033", "#993333", "#FF9933", "#FFCC00", "#99FF00", "#66CC33", "#006633", 
#                                          "#99FFFF", "#0099FF", "#0033FF", "#CC99FF", "#9966FF", "#6600FF")) +  
#                                            theme_bw() +
#     theme(panel.grid = element_blank())
# }


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
         y = paste("Estimated ln(", this.chemical, ")", sep=""), color=pred) +
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
         y = paste("Estimated ln(", this.chemical, ")", sep=""), color=pred) +
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
  # pred <- "sqrt_CO2_res"
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
         y = paste("Estimated ln(", this.chemical, ")", sep=""), 
         color=pred,
         shape="censored\ndata") +
    scale_color_distiller(palette = paletteCol, direction=1) +   #colors should be light for low values, dark for high values
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position="right",
      legend.box="vertical"
    ) +
    guides(
      color = guide_colorbar(order=1),
      shape = guide_legend(order=2)
    )
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
         y = paste("Estimated ln(", this.chemical, ")", sep=""), 
         color=pred,
         shape="censored\ndata") +
    scale_color_manual(values=paletteCol) +                     #HEX color values are provided based on # of discrete values for pred 
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position="right",
      legend.box="vertical"
    ) +
    guides(
      color = guide_legend(order=1),
      shape = guide_legend(order=2)
    )
}


#-------------------------------#
#  Credible Interval Plot Fxns  #
#-------------------------------#

#function to plot expectation of predicted draws; aka: expected values vs raw predictors, 
#  using the epred_draws function. This one works for TWO landscape predictors
credInt.vs.rawPreds.2preds <- function(this.coc2, this.brm, this.chemical, preds, box.widths=c(0.9, 0.1)) {

  # ####### FOR TESTING ONLY!!
  # this.coc2 <- TSS.coc2
  # this.brm <- TSS.brm
  # this.chemical <- "TSS"
  # preds <- c("sqrt_CO2_road", "oldDev")
  # box.widths <- c(0.9, 0.1)
  # #######
  
  #throw an error if number of predictors is not 1 or 2
  if(length(preds)!=2) {
    stop("this function only works for 2 landscape predictors")
  }
  
  #names of predictors in this model
  predictors <- variables(this.brm) %>% 
    keep(~ str_detect(.x, "^b_")) %>%     # Keeps only elements starting with "b_"
    str_remove("^b_")                     # Removes the "b_" prefix

  #vector of possible precip-associated predictors -- we'll select from these below  
  p_weather <- c("rain1", "rain14", "rain21", "dry")
  
  #take 94 samples total (10% of all S8 Data we have) to represent feasible rain/ dry/ season combinations. 84 samples from non-summer, 10 from summer; 
  #this matches the summer/ non-summer ratio of 11% summer. If spring is in the model, 75 from non-spring and 19 from spring.
  if("spring1" %in% predictors) {
    my.sample0 <- this.coc2[which(this.coc2$spring==0)[sample(length(which(this.coc2$spring==0)), 75)], c(p_weather[p_weather %in% predictors], "spring")]
    my.sample1 <- this.coc2[which(this.coc2$spring==1)[sample(length(which(this.coc2$spring==1)), 19)], c(p_weather[p_weather %in% predictors], "spring")]
    my.sample <- rbind(my.sample0, my.sample1)

    pred1 <- seq(min(this.coc2[[preds[1]]]), max(this.coc2[[preds[1]]]), length.out=8)
    pred2 <- seq(min(this.coc2[[preds[2]]]), max(this.coc2[[preds[2]]]), length.out=8)
    my.preds <- data.frame(pred1, pred2)
    colnames(my.preds) <- preds
        
    newdata <- crossing(my.preds, my.sample)
    newdata$spring <- as.factor(newdata$spring)
  } else {  #otherwise assume its summer-ish
    my.sample0 <- this.coc2[which(this.coc2$summer==0)[sample(length(which(this.coc2$summer==0)), 84)], c(p_weather[p_weather %in% predictors], "summer")]
    my.sample1 <- this.coc2[which(this.coc2$summer==1)[sample(length(which(this.coc2$summer==1)), 10)], c(p_weather[p_weather %in% predictors], "summer")]
    my.sample <- rbind(my.sample0, my.sample1)
    
    pred1 <- seq(min(this.coc2[[preds[1]]]), max(this.coc2[[preds[1]]]), length.out=8)
    pred2 <- seq(min(this.coc2[[preds[2]]]), max(this.coc2[[preds[2]]]), length.out=8)
    my.preds <- data.frame(pred1, pred2)
    colnames(my.preds) <- preds

    newdata <- crossing(my.preds, my.sample)
    newdata$summer <- as.factor(newdata$summer)
  }
  
  #epred_draws gives the "expected values of the outcome" - it focuses on the uncertainty in the model parameters, and not
  #   the individual-level residuals.  epred_draws accounts for the uncertainty from two sources:
  #       1.) uncertainty of fixed coefficients
  #       2.) uncertainty of the variance parameters of the groups (eg: sd(Intercept)) (aka: uncertainty in group-level effects);
  #           This does NOT include any elements that deal with heteroskedasticity, which only affect individual-level residuals.
  #  https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/
  #         (see yellow and blue curves ~ 1/4 or 1/3 of the way down, titled: Predicted draws and Expectation of predicted draws)
  abc <- epred_draws(this.brm, allow_new_levels=TRUE, newdata=newdata, ndraws=100)

  
  #calculate raw predictor values for this COC, and for the epred draws; these are used in plotting ribbons and points, below
  this.coc2$pred1_raw <- (this.coc2[[preds[1]]] * std.vals$sd[which(std.vals$predictor==preds[1])] + std.vals$mean[which(std.vals$predictor==preds[1])] ) # ^
  this.coc2$pred2_raw <- (this.coc2[[preds[2]]] * std.vals$sd[which(std.vals$predictor==preds[2])] + std.vals$mean[which(std.vals$predictor==preds[2])] ) # ^
  abc$pred1.raw <- (abc[[preds[1]]] * std.vals$sd[which(std.vals$predictor==preds[1])] + std.vals$mean[which(std.vals$predictor==preds[1])] ) # ^
  abc$pred2.raw <- (abc[[preds[2]]] * std.vals$sd[which(std.vals$predictor==preds[2])] + std.vals$mean[which(std.vals$predictor==preds[2])] ) # ^

  #plot RAW values for pred1 and pred2; these plots show 50%, 80% and 95% CI's; if median is desired as well, use stat_lineribbon instead of stat_ribbon!
  pred1.raw <- ggplot(abc, aes(x = pred1.raw, y = exp(.epred))) +
    stat_ribbon() +    #NOTE: stat_lineribbon would provide a median line (in black); this is not suitable for our dataset with 2 independent variables
    scale_fill_brewer(palette = "Reds") +
    # geom_boxplot(data=this.coc2, mapping=aes(x=pred1_raw, y=exp(result), group=location), notch=FALSE, 
    #              width=box.widths[1], position=position_dodge(preserve="single"),      #box width; make boxes uniform in width
    #              outlier.size=0.65, outlier.alpha=0.5,                       #make outlier points smaller and less dominant in the plot
    #              fill="white", alpha=1) +  #add boxplots; width=0.4; TKN width=0.15
    stat_summary(data=this.coc2, mapping=aes(x=pred1_raw, y=exp(result), group=location),
                 fun=median,
                 fun.min=median, fun.max=median,
                 geom="crossbar",
                 width=box.widths[1],
                 color="black",
                 fatten=3 ) + 
    labs(x = preds[1], y = paste(this.chemical, "(ppm)"), fill = "Credible Interval") +
    theme(legend.position = "bottom")
  
  pred2.raw <- ggplot(abc, aes(x = pred2.raw, y = exp(.epred))) +
    stat_ribbon() +
    scale_fill_brewer(palette = "Purples") +
    stat_summary(data=this.coc2, mapping=aes(x=pred2_raw, y=exp(result), group=location),
                 fun=median,
                 fun.min=median, fun.max=median,
                 geom="crossbar",
                 width=box.widths[2],
                 color="black",
                 fatten=3) + 
    # geom_boxplot(data=this.coc2, mapping=aes(x=pred2_raw, y=exp(result), group=location), notch=FALSE, 
    #              width=box.widths[2], position=position_dodge(preserve="single"),      #box width; make boxes uniform in width
    #              outlier.size=0.65, outlier.alpha=0.5,                       #make outlier points smaller and less dominant in the plot
    #              fill="white", alpha=1) +  #add boxplots; width=0.4; TKN width=0.15
    labs(x = preds[2], y = paste(this.chemical, "(ppm)"), fill = "Credible Interval") +
    # geom_point(data=this.coc2, mapping=aes(x=pred2_raw, y=exp(result))) +
    theme(legend.position = "bottom")
  
  grid.arrange(pred1.raw, pred2.raw, nrow=1, ncol=2)

}



#function to plot expectation of predicted draws; aka: expected values vs raw predictors, 
#  using the epred_draws function. This one works for ONE landscape predictor
credInt.vs.rawPreds.1pred <- function(this.coc2, this.brm, this.chemical, preds, box.widths=1) {
  
  # ####### FOR TESTING ONLY!!  Note -- copper won't work here, as it has 2 predictors!
  # this.coc2 <- P.coc2
  # this.brm <- P.brm
  # this.chemical <- "P"
  # preds <- "oldDev"
  # #######

  #throw an error if number of predictors is not 1 or 2
  if(length(preds)!=1) {
    stop("this function only works for 1 landscape predictor")
  }

  #names of predictors in this model
  predictors <- variables(this.brm) %>% 
    keep(~ str_detect(.x, "^b_")) %>%     # Keeps only elements starting with "b_"
    str_remove("^b_")                     # Removes the "b_" prefix
  
  #vector of possible precip-associated predictors -- we'll select from these below  
  p_weather <- c("rain1", "rain14", "rain21", "dry")
  
  #take 94 samples total (10% of all S8 Data we have) to represent feasible rain/ dry/ season combinations. 84 samples from non-summer, 10 from summer; 
  #this matches the summer/ non-summer ratio of 11% summer. If spring is in the model, 75 from non-spring and 19 from spring.
  if("spring1" %in% predictors) {
    my.sample0 <- this.coc2[which(this.coc2$spring==0)[sample(length(which(this.coc2$spring==0)), 75)], c(p_weather[p_weather %in% predictors], "spring")]
    my.sample1 <- this.coc2[which(this.coc2$spring==1)[sample(length(which(this.coc2$spring==1)), 19)], c(p_weather[p_weather %in% predictors], "spring")]
    my.sample <- rbind(my.sample0, my.sample1)
    
    pred1 <- seq(min(this.coc2[[preds[1]]]), max(this.coc2[[preds[1]]]), length.out=8)
    # pred2 <- seq(min(this.coc2[[preds[2]]]), max(this.coc2[[preds[2]]]), length.out=8)
    my.preds <- data.frame(pred1)
    colnames(my.preds) <- preds
    
    newdata <- crossing(my.preds, my.sample)
    newdata$spring <- as.factor(newdata$spring)
  } else {  #otherwise assume its summer-ish
    my.sample0 <- this.coc2[which(this.coc2$summer==0)[sample(length(which(this.coc2$summer==0)), 84)], c(p_weather[p_weather %in% predictors], "summer")]
    my.sample1 <- this.coc2[which(this.coc2$summer==1)[sample(length(which(this.coc2$summer==1)), 10)], c(p_weather[p_weather %in% predictors], "summer")]
    my.sample <- rbind(my.sample0, my.sample1)
    
    pred1 <- seq(min(this.coc2[[preds[1]]]), max(this.coc2[[preds[1]]]), length.out=8)
    # pred2 <- seq(min(this.coc2[[preds[2]]]), max(this.coc2[[preds[2]]]), length.out=8)
    my.preds <- data.frame(pred1)
    colnames(my.preds) <- preds
    
    newdata <- crossing(my.preds, my.sample)
    newdata$summer <- as.factor(newdata$summer)
  }
  #epred_draws gives the "expected values of the outcome" - it focuses on the uncertainty in the model parameters, and not
  #   the individual-level residuals.  epred_draws accounts for the uncertainty from two sources:
  #       1.) uncertainty of fixed coefficients
  #       2.) uncertainty of the variance parameters of the groups (eg: sd(Intercept)) (aka: uncertainty in group-level effects);
  #           This does NOT include any elements that deal with heteroskedasticity, which only affect individual-level residuals.
  #  https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/
  abc <- epred_draws(this.brm, allow_new_levels=TRUE, newdata=newdata)
  
  #calculate raw predictor values for this COC, and for the epred draws; these are used in plotting ribbons and points, below
  this.coc2$pred1_raw <- this.coc2[[preds[1]]] * std.vals$sd[which(std.vals$predictor==preds[1])] + std.vals$mean[which(std.vals$predictor==preds[1])] 
  abc$pred1.raw <- abc[[preds[1]]] * std.vals$sd[which(std.vals$predictor==preds[1])] + std.vals$mean[which(std.vals$predictor==preds[1])]

  pred1.raw <- ggplot(abc, aes(x = pred1.raw, y = exp(.epred))) +
    stat_ribbon() +    #NOTE: stat_lineribbon would provide a median line (in black); this is not suitable for our dataset with 2 independent variables
    scale_fill_brewer(palette = "Reds") +
    # geom_boxplot(data=this.coc2, mapping=aes(x=pred1_raw, y=exp(result), group=location), notch=FALSE, 
    #              width=box.widths[1], position=position_dodge(preserve="single"),      #box width; make boxes uniform in width
    #              outlier.size=0.65, outlier.alpha=0.5,                       #make outlier points smaller and less dominant in the plot
    #              fill="white", alpha=1) +  #add boxplots; width=0.4; TKN width=0.15
    stat_summary(data=this.coc2, mapping=aes(x=pred1_raw, y=exp(result), group=location),
                 fun=median,
                 fun.min=median, fun.max=median,
                 geom="crossbar",
                 width=box.widths,
                 color="black",
                 fatten=3 ) + 
    labs(x = preds[1], y = paste(this.chemical, "(ppm)"), fill = "Credible Interval") +
    theme(legend.position = "bottom")
  
  grid.arrange(pred1.raw, nrow=1, ncol=1)
}




#------------------------------#
#  Intercepts by Location Fxn  #
#------------------------------#

#this version separates out the location-specific variability (loc_draws), agency-specific variability (agency_draws)
#  and then adds them to the population-level intercept to come up with the correct location (nested in agency) intercept
plotIntercepts.global.location <- function(this.brm, this.coc2, xlims=c(0,10)) {
  
  # 1. Clean your original data frame to match how brms formats spaces as periods
  # This acts as our master lookup key
  data_registry <- this.coc2 %>%
    mutate(
      agency_clean = stringr::str_replace_all(agency, " ", "."),
      location_clean = stringr::str_replace_all(location, " ", ".")
    ) %>%
    distinct(agency_clean, location_clean)
  
  # 2. Extract the nested group level effects
  # This extracts columns named 'agency_location' (which brms formatted)
  loc_draws <- this.brm %>%
    spread_draws(`r_agency:location`[agency_location, Intercept])
  
  # 3. Extract the higher level agency effects
  agency_draws <- this.brm %>%
    spread_draws(r_agency[agency, Intercept])
  
  # 4. Combine them cleanly by leveraging our registry lookup
  this.brm %>%
    spread_draws(b_Intercept) %>% 
    left_join(loc_draws, by = ".draw") %>% 
    # Use regular expressions to extract the location code from the index string
    mutate(
      # If brms formatted it as "King.Co_KIC_COM", this separates them cleanly
      location_clean = stringr::str_extract(agency_location, "(KIC|SNO|PIE|SEA|TAC|POT)_[A-Z]+")
    ) %>% 
    # Match the location back to its true agency based on your original data rules
    left_join(data_registry, by = "location_clean") %>% 
    # Now join the agency draws safely by matching the clean agency name and draw id
    left_join(agency_draws, by = c(".draw", "agency_clean" = "agency")) %>% 
    mutate(agency_location_mean = b_Intercept + r_agency + `r_agency:location`) %>%
    
    # --- Your original ggplot engine ---
    ggplot(aes(y = agency_location, x = agency_location_mean, 
               # fill=stat(abs(x - fixef(this.brm)["Intercept", 1]) < (fixef(this.brm)["Intercept", 2]*2) )))  +
               fill = stat(x >= fixef(this.brm)["Intercept", "Q2.5"] & 
                             x <= fixef(this.brm)["Intercept", "Q97.5"]))) +
    xlim(xlims) +
    stat_halfeye() +
    labs(x="Intercept Value", y="", fill="Population-Level \nIntercept 95% CI") +
    geom_vline(xintercept = fixef(this.brm)["Intercept", 1], linetype = "solid", color="cadet blue4", alpha=0.5, size=2) +
    annotate("text", x=(fixef(this.brm)["Intercept", 1] - 0.1), y=7.2, label="population intercept", color="cadet blue4", angle=90, size=6) +
    geom_vline(xintercept = fixef(this.brm)["Intercept", c(3:4)], linetype = "dashed", color="black", size=1/4) +
    scale_fill_manual(values = c("gray80", "skyblue"))
}



### THIS ONE IS WRONG!  It conflates the agency-level predictors with EVERY location level.
#function to plot location-specific intercepts with densities, plus the population-level intercept and 95% credibility 
# intervals around the population-level intercept
#     Note that the 95% CI does NOT include variability in agencies or in location:agency; it is just the 95% CI
#          for the population-level intercept (ie: 95% of the time, the population-level intercept will fall between 
#          the 95% CI values)
# plotIntercepts.global.location <- function(this.brm, xlims=c(0,10)) {
#   this.brm %>%
#     spread_draws(b_Intercept, r_agency[agency, Intercept], `r_agency:location`[agency_location,Intercept]) %>%
#     mutate(agency_location_mean = b_Intercept + `r_agency` + `r_agency:location`) %>%
#     ggplot(aes(y = agency_location, x = agency_location_mean,
#                fill=stat(abs(x - fixef(this.brm)["Intercept", 1]) < (fixef(this.brm)["Intercept", 2]*2) )))  +
#     xlim(xlims) +
#     stat_halfeye() +
#     labs(x="Intercept Value", y="", fill="Population-Level \nIntercept 95% CI") +
#     geom_vline(xintercept = fixef(this.brm)["Intercept", 1], linetype = "solid", color="cadet blue4", alpha=0.5, size=2) +
#     annotate("text", x=(fixef(this.brm)["Intercept", 1] - 0.1), y=7.2, label="population intercept", color="cadet blue4", angle=90, size=6) +
#     geom_vline(xintercept = fixef(this.brm)["Intercept", c(3:4)], linetype = "dashed", color="black", size=1/4) +
#     scale_fill_manual(values = c("gray80", "skyblue"))# +
# }





