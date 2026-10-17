# This script contains functions used in the Model Validation.R script

# Author: Eva Dusek Jennings
# Date: Oct 10, 2026
#---------------------------------------------------------------------

library(brms)
library(purrr)


#--------------------------------------#
#  Data Prep & Manipulation Functions  #
#--------------------------------------#


#provide prediction intervals for a location with particular predictor values.  Prediction intervals predict in what range a future
#   individual observation will fall
interpPI <- function(my.brm, my.coc2, wr_preds) {

  #extract predictor names from the brms model
  predictors <- variables(my.brm) %>%
    keep(~ str_detect(.x, "^b_")) %>%     # Keeps only elements starting with "b_"
    str_remove("^b_")                     # Removes the "b_" prefix

  #vector of possible precip-associated predictors -- we'll select from these below  
  p_weather <- c("rain1", "rain14", "rain21", "dry")
  
  #take 94 samples total (10% of all S8 Data we have) to represent feasible rain/ dry/ season combinations. 84 samples from non-summer, 10 from summer; 
  #this matches the summer/ non-summer ratio of 11% summer. If spring is in the model, 75 from non-spring and 19 from spring.
  if("spring1" %in% predictors) {
    my.sample0 <- my.coc2[which(my.coc2$spring==0)[sample(length(which(my.coc2$spring==0)), 75)], c(p_weather[p_weather %in% predictors], "spring")]
    my.sample1 <- my.coc2[which(my.coc2$spring==1)[sample(length(which(my.coc2$spring==1)), 19)], c(p_weather[p_weather %in% predictors], "spring")]
    my.sample <- rbind(my.sample0, my.sample1)
    
    newdata <- crossing(wr_preds, my.sample)
    newdata$spring <- as.factor(newdata$spring)
  } else {  #otherwise assume its summer-ish
    my.sample0 <- my.coc2[which(my.coc2$summer==0)[sample(length(which(my.coc2$summer==0)), 84)], c(p_weather[p_weather %in% predictors], "summer")]
    my.sample1 <- my.coc2[which(my.coc2$summer==1)[sample(length(which(my.coc2$summer==1)), 10)], c(p_weather[p_weather %in% predictors], "summer")]
    my.sample <- rbind(my.sample0, my.sample1)

    newdata <- crossing(wr_preds, my.sample)
    newdata$summer <- as.factor(newdata$summer)
  }
  
  #prediction interval, which predicts the range in which a future individual observation will fall; these include three levels of 
  #   uncertainty: uncertainty in fixed coefficients, uncertainty in variance parameters for grouping factors, and residual variance
  #   re_formula=NULL indicates that all group-level effects should be considered in the prediction
  #   allow_new_levels=TRUE indicates that new group levels (agency) are fine
  my.site <- posterior_predict(my.brm, newdata = newdata, re_formula = NULL, allow_new_levels=TRUE) %>%
    t() %>%
    as_tibble() %>%
    bind_cols(newdata, .) %>%
    pivot_longer(
      cols=starts_with("V"),
      names_to="pars",
      values_to="est"
    )

  #summarized data for each type of site (specific values of landscape predictors)
  #  this summarizes what might be expected to occur in samples at a site with particular landscape predictor traits 
  #  over a representative set of season (summer/not summer) and rainfall value associated with that season
  target_probs <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
  my.summary_long <- my.site %>%
    group_by(across(all_of(names(wr_preds)))) %>%
    reframe(
      Quantile = target_probs,
      Value = quantile(est, probs = target_probs)
    )
  
  rm(my.site, newdata)
  return(my.summary_long)
}





#----------------------#
#  Plotting Functions  #
#----------------------#

#make a function for scatter plots
scatter_cocs <- function(DF, Title, Subtitle, Caption) {
  p <- ggplot(DF, aes(1, Result)) + geom_jitter() + labs(
    title = Title,
    subtitle = Subtitle,
    caption = Caption,
    x = "Observations"
  )
  p + facet_wrap( ~ Analyte, scales = 'free')+theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(),
                                                    axis.ticks.x=element_blank())
}

# #plot upper and lower limits with data -- all CI's and data are ln-transformed
# plot_PIs <- function(my.PI, my.data) {
#   ggplot() +
#     geom_segment(aes(x=Location, xend=Location, y=lowerCI, yend=upperCI), data=my.PI[[3]], colour=orangePalette[3], linewidth=20) +
#     geom_segment(aes(x=Location, xend=Location, y=lowerCI, yend=upperCI), data=my.PI[[2]], colour=orangePalette[2], linewidth=20) +
#     geom_segment(aes(x=Location, xend=Location, y=lowerCI, yend=upperCI), data=my.PI[[1]], colour=orangePalette[1], linewidth=20) +
#     geom_point(aes(x=Location, y=log(Result)), data=my.data) +  #add ln-transformed data points from this study
#     ggtitle(my.data$Analyte[1]) + 
#     xlab("Location") +
#     ylab(paste("ln-transformed ", my.data$Analyte[1], " (", my.data$Unit[1], ")", sep="")) +
#     theme(plot.title=element_text(face="bold", size=14, hjust=0.5))
# }

# #plot upper and lower limits with data -- all CI's and data are raw values
# plot_PIs_raw <- function(my.PI, my.data) {
# 
#   # my.PI <- comp_Cu_PI %>%
#   #   mutate(location = Location)
#   # my.data <- comp_res_pred[which(comp_res_pred$Analyte=="Total Copper"),]
#   # # my.PI.trim <- my.PI %>%
#   # #   filter(my.PI[[1]]$location %in% unique(my.data$location))
# 
#   ggplot() +
#     geom_segment(aes(x=location, xend=location, y=exp(`0.025`), yend=exp(`0.975`)), data=my.PI, colour=orangePalette[3], linewidth=5) +
#     geom_segment(aes(x=location, xend=location, y=exp(`0.1`), yend=exp(`0.9`)), data=my.PI, colour=orangePalette[2], linewidth=5) +
#     geom_segment(aes(x=location, xend=location, y=exp(`0.25`), yend=exp(`0.75`)), data=my.PI, colour=orangePalette[1], linewidth=5) +
#     geom_point(aes(x=location, y=rawResult), data=my.data) +  #add ln-transformed data points from this study
#     ggtitle(my.data$Analyte[1]) + 
#     xlab("location") +
#     ylab(paste(my.data$Analyte[1], " (", my.data$Unit[1], ")", sep="")) +
#     theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
#           axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
# }

# #plot upper and lower limits with data -- all CI's and data are raw values
# plot_PIs_raw <- function(my.PI, my.data) {
#   
#   # my.PI <- comp_TKN_PI
#   # my.data <- comp_res_pred[which(comp_res_pred$Analyte=="Total Kjeldahl Nitrogen"),]
#   # my.PI.trim <- my.PI %>%
#   #   filter(my.PI[[1]]$Location %in% unique(my.data$Location))
#   # 
#   
#   
#   ggplot() +
#     geom_segment(aes(x=Location, xend=Location, y=exp(lowerCI), yend=exp(upperCI)), data=my.PI[[3]], colour=orangePalette[3], linewidth=5) +
#     geom_segment(aes(x=Location, xend=Location, y=exp(lowerCI), yend=exp(upperCI)), data=my.PI[[2]], colour=orangePalette[2], linewidth=5) +
#     geom_segment(aes(x=Location, xend=Location, y=exp(lowerCI), yend=exp(upperCI)), data=my.PI[[1]], colour=orangePalette[1], linewidth=5) +
#     geom_point(aes(x=Location, y=Result), data=my.data) +  #add ln-transformed data points from this study
#     ggtitle(my.data$Analyte[1]) + 
#     xlab("Location") +
#     ylab(paste(my.data$Analyte[1], " (", my.data$Unit[1], ")", sep="")) +
#     theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
#           axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
# }



#plot upper and lower limits with data -- all CI's and data are raw values
plot_PIs_boxplot_raw <- function(my.PI, my.data, box.width, line.width=7) {

  #generate a wide-format version of the PI dataframe/tibble
  my.PI.wide <- my.PI %>%
    pivot_wider(names_from = Quantile,
                values_from = Value)

  ggplot() +
    # geom_segment(aes(x=location, xend=location, y=exp(`0.025`), yend=exp(`0.975`)), data=my.PI.wide, colour=orangePalette[3], linewidth=line.width) +
    geom_segment(aes(x=location, xend=location, y=exp(`0.1`), yend=exp(`0.9`)), data=my.PI.wide, colour=orangePalette[2], linewidth=line.width) +
    geom_segment(aes(x=location, xend=location, y=exp(`0.25`), yend=exp(`0.75`)), data=my.PI.wide, colour=orangePalette[1], linewidth=line.width) +
    geom_boxplot(aes(x=location, y=rawResult), data=my.data, notch=FALSE, width=box.width, #box.linewidth=1.0,
                 fill="white", alpha=1) +  #add boxplots; width=0.4; TKN width=0.15
    ggtitle(my.data$Analyte[1]) +
    xlab("Location") +
    ylab(paste(my.data$Analyte[1], " (", my.data$Unit[1], ")", sep="")) +
    theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
          axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
}



#  Observed vs Predicted
#  https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html#plotting-multivariate-posteriors.
#  go down to 5.1.3.3 POSTERIOR PREDICTION PLOTS



# observed vs predicted with points colored by location and gray 95% credible interval lines around each point
obs.vs.pred.median.mv <- function(this.brm, this.chemical, mv.dat, myColors, mySymbols) {
  
  # this.brm <- Cu.brm
  # this.chemical <- "Copper"
  # mv.dat <- Cu.mv.dat2
  # myColors <- myColors
  
  fitted(this.brm, newdata=mv.dat, allow_new_levels=TRUE, probs=c(0.25, 0.75)) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5; note that fitted() ignores residual error
  as_tibble() %>%
  bind_cols(mv.dat) %>%
  group_by(location) %>%
  summarise_at(c("Estimate", "result"), median) %>%
  ggplot(aes(x = result, y = Estimate, color=location, shape=location)) +       #color points by location (see custom colors below) 
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(size = 3.5, alpha = 3/4, stroke=3) +
  labs(x = paste("Observed median ln(", this.chemical, ")", sep=""), 
       y = paste("Estimated median ln(", this.chemical, ")", sep="") ) +
  scale_color_manual(name="location", values=myColors) +
  scale_shape_manual(name="location", values=mySymbols) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text=element_text(size=rel(0.8)),
        legend.key.size=unit(0.8, "cm"))
}


# observed vs predicted with points colored by location and gray 95% credible interval lines around each point
obs.vs.pred.median.area.mv <- function(this.brm, this.chemical, mv.dat) {
  
  # this.brm <- totZn.brm
  # this.chemical <- "total Zinc"
  # mv.dat <- Zn.mv.dat2
  # myColors <- myColors
  
  fitted(this.brm, newdata=mv.dat, allow_new_levels=TRUE, probs=c(0.25, 0.75)) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5; note that fitted() ignores residual error
    as_tibble() %>%
    bind_cols(mv.dat) %>%
    group_by(location, sqrt_area_m2) %>%
    summarise_at(c("Estimate", "result"), median) %>%
    ggplot(aes(x = result, y = Estimate, color=log(sqrt_area_m2))) +       #color points by area (log(sqrt(area))) 
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_point(size = 4.5, alpha = 3/4) +
    labs(x = paste("Observed median ln(", this.chemical, ")", sep=""), 
         y = paste("Estimated median ln(", this.chemical, ")", sep="") ) +
    scale_color_distiller(palette = "Purples", direction=-1) +   #colors should be light for low values, dark for high values
    theme_bw() +
    theme(panel.grid = element_blank())
}


# observed vs predicted with points colored by location and gray 95% credible interval lines around each point
obs.vs.pred.mv <- function(this.brm, this.chemical, mv.dat, myColors, mySymbols) {

    fitted(this.brm, newdata=mv.dat, allow_new_levels=TRUE, probs=c(0.25, 0.75)) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5; note that fitted() ignores residual error
    as_tibble() %>%
    bind_cols(mv.dat) %>%
    ggplot(aes(x = result, y = Estimate, color=location, shape=location)) +       #color points by location (see custom colors below) 
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    # geom_linerange(aes(ymin = Q25, ymax = Q75),               #the thin lines are the 95% intervals
    #                size = 1/4, color="gray") +
    # geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
    #                    ymax = Estimate + Est.Error),
    #                size = 5/8, color="darkgray") +
    geom_point(size = 2.5, alpha = 3/4) +
    labs(x = paste("Observed ln(", this.chemical, ")", sep=""), 
         y = paste("Estimated ln(", this.chemical, ")", sep="") ) +
    scale_color_manual(name="location", values=myColors) +
    scale_shape_manual(name="location", values=mySymbols) +
    theme_bw() +
    theme(panel.grid = element_blank())
}


#function to plot observed vs predicted, with colors showing a selected predictor.  
#   Inputs: brm, coc2, chemical name, predictor name, color brewer palette
obsPredPlot.mv <- function(this.brm, mv.dat, this.chemical, pred, paletteCol) {
  
# #testing
# mv.dat <- mv.dat.Cu
# this.brm <- Cu.brm
# this.chemical <-  "Copper"
# pred <- "sqrt_traffic"
# paletteCol <- "Reds"
  
  fitted(this.brm, newdata=mv.dat, allow_new_levels=TRUE, probs=c(0.25, 0.75)) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5; note that fitted() ignores residual error
    as_tibble() %>%
    bind_cols(mv.dat) %>%
  
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
obsPredPlot2.mv <- function(this.brm, mv.dat, this.chemical, pred, paletteCol) {

  fitted(this.brm, newdata=mv.dat, allow_new_levels=TRUE, probs=c(0.25, 0.75)) %>%   #this gives 4 columns: Estimate, Est.Error, Q2.5, and Q97.5; note that fitted() ignores residual error
    as_tibble() %>%
    bind_cols(mv.dat) %>%
    
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


#-----------------------------------------------------#
#  Functions for Single Predictor Plots with S8 Data  #
#-----------------------------------------------------#

#Function to plot result vs one predictor, with location shown in color
Plot.One.Predictor.wAgency.MV <- function(aa, myModel, myCOC, hideLegend=TRUE, xEqn=0.15, yEqn=max(s8$result*1.02, mv$result*1.02) ) {
  coef.int <- myModel$coefficients[[1]]["(Intercept)"]
  coef.slope <- myModel$coefficients[[1]][aa]
  
  #linear equation for the primary predictor, to be printed in the plot
  eqnText <- paste("ln(", myCOC, ") = ", 
                   round(coef.int, 2), " + ",
                   round(coef.slope, 2), "*", aa )
  #plot the primary predictor on x-axis, with secondary predictor shown in color
  ggplot() + 
    geom_point(data=s8, aes(x=s8[,aa], y=result), color="darkgray", size=3) + 
    geom_point(data=mv, aes(x=mv[,aa], y=result, color=Project), size=3) +
    scale_color_manual(values = myColors) +
    geom_abline(intercept=coef.int, slope=coef.slope, col="black", lwd=1.5) +  #, aes(x = aa, y = y.hat), col="darkred", lwd = 1.5)
    ylab(paste("ln", myCOC)) + xlab(aa) +
    geom_text(aes(x=xEqn, y=yEqn, label=eqnText), cex=4.5, color="black")+  #add equation to the top of the plot
    theme_gray() +
    if (hideLegend) {
      theme(legend.position="none")
    }
}    

#function to plot best fit LME model line to a single predictor
plot.single.preds.MV <- function(theModel, myCOC, xEqn=0.15) {
  aa <- summary(theModel)$coefficients$fixed
  aa <- names(aa[names(aa) %in% thePreds])
  if(length(aa)==1) {
    ppp1 <- Plot.One.Predictor.wAgency.MV(aa[1], theModel, myCOC, hideLegend=FALSE, xEqn)
    lay <- matrix(c(1,1,1,1,NA,NA,NA), nrow=1, ncol=7, byrow=TRUE)
    grid.arrange(ppp1, layout_matrix=lay)
  }
  if(length(aa)==2) {
    ppp1 <- Plot.One.Predictor.wAgency.MV(aa[1], theModel, myCOC, hideLegend=TRUE, xEqn)
    ppp2 <- Plot.One.Predictor.wAgency.MV(aa[2], theModel, myCOC, hideLegend=FALSE, xEqn) 
    lay <- matrix(c(1,1,1,2,2,2,2), nrow=1, byrow=TRUE)
    grid.arrange(ppp1, ppp2, layout_matrix=lay)
  }
}
# 
# plot.single.preds.MV(Model, my.COC)
# theModel <- Model
# myCOC <- my.COC
# 
