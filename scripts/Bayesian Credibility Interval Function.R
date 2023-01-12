

# install.packages("tidyverse")
# update.packages("dbplyr")
library(tidyr)
library(dplyr)
#library(tidyverse)    # ggplot, dplyr, %>%, and friends
library(brms)         # Bayesian modeling through Stan
#library(tidybayes)    # Manipulate Stan objects in a tidy way
#library(broom)        # Convert model objects to data frames
#library(forcats)
#library(oce)          # approx3d
#library(pracma)       # interp2


#--------------------------------------------------------------------------------------#
#  Function to Generate Credibility Intervals for each watershed, using interpolation  #
#                        random effect = agency/location                               #
#--------------------------------------------------------------------------------------#

#this function accepts information including a COC's brm model, coc2 dataframe, predictors, and which upper and lower credibility 
#  intervals are desired, then generates a tibble with the lower and upper difference from median values for credibility intervals,
#  along with the associated predictor values.
#  To allow proper sampling of the rain/ summer parameter space, 10% of all data are resampled for rain/ summer, such that 11% are summer
#  rain events (corresponding to the percentage of actual data of summer rain events) and 89% are non-summer rain events.

#INPUTS: LCI = lower credibility interval; defaults to 80% CI's
#        UCI = upper credibility interval; defaults to 80% CI's
#        my.brm = Bayesian model for this COC
#        my.coc2 = coc2 data frame for this COC (saved in the frequentist model output for this COC)
#        wr_pred = df of watershed-reduction landscape predictor values for this COC
interpCI <- function(LCI=0.1, UCI=0.9, my.brm, my.coc2, wr_preds) {
  
  # #testing!
  # LCI <- 0.1
  # UCI <- 0.9
  # my.brm <- Cu.brm
  # my.coc2 <- Cu.coc2
  # wr_preds <- Cu.wr_preds
  # # ######
  

  preds <- colnames(wr_preds)
  
  #take 47 samples total (10% of all data we have) to represent feasible rain or rain/summer options
  #   42 samples from non-summer, 5 from summer; this matches the summer/ non-summer ratio of 11% summer.
  #   NOTE: for the end result, it won't matter whether summer is a predictor or not, as we are drawing the same proportion of 
  #         samples from summertime as are in the data; this should approximately happen during sampling, regardless.
  my.sample0 <- my.coc2[which(my.coc2$summer==0)[sample(length(which(my.coc2$summer==0)), 42)], c("rain", "summer")]
  my.sample1 <- my.coc2[which(my.coc2$summer==1)[sample(length(which(my.coc2$summer==1)), 5)], c("rain", "summer")]
  my.sample <- rbind(my.sample0, my.sample1)
  
  #generate simulated data
  newdata0 <- expand_grid(#location=NA, 
                          rain=my.sample0$rain, 
                          summer=0, 
                          wr_preds)
  newdata1 <- expand_grid(#location=NA, 
                          rain=my.sample1$rain, 
                          summer=1, 
                          wr_preds)
  newdata <- rbind(newdata0, newdata1)

  #conditional mean for an unknown participant, using the population-level intercept plus some variability around where the
  #   location-specific intercept should be.  This is indicated with re_formula=NULL and allow_new_levels=TRUE.  Obtain all estimates  
  #   (not just summarized data) by setting summary=FALSE
  # start_time <- Sys.time()
  my.site <- fitted(my.brm, newdata = newdata, re_formula = NULL, allow_new_levels=TRUE, summary=FALSE) %>%
    t() %>%
    as_tibble() %>%
    bind_cols(newdata, .) %>%
    pivot_longer(
      cols=starts_with("V"),
      names_to="pars",
      values_to="est"
    )
  # end_time <- Sys.time()
  # end_time - start_time
  
  #summarized data for each type of site (specific values of landscape predictors)
  #  this summarizes what might be expected to occur in samples at a site with particular landscape predictor traits 
  #  over a representative set of season (summer/not summer) and rainfall value associated with that season
  my.summary <- my.site %>%
    group_by(across(all_of(names(wr_preds)))) %>%  #group by the location, and include landscape predictors for this COC
    summarise(upperCI = quantile(est, UCI),
              lowerCI = quantile(est, LCI), 
              med=quantile(est, 0.5))

  return(my.summary)
}


#function to plot credibility intervals - transformed & standardized version
#plot median (black) and upper/ lower CI's (1st predictor: oranges; 2nd predictor: purples)
plotBayesianCIs <- function(this.CI, this.chemical, preds) {

  par(mfrow=c(length(preds),1), mar=c(4, 4, 3, 7), xpd=TRUE)
  plot(this.CI$med ~ this.CI[[preds[1]]], type="n", pch=16, col="black", ylim=c(min(this.CI$lowerCI), max(this.CI$upperCI)),
       main="", ylab=paste("ln(", this.chemical, ")", sep=""), xlab=paste("centered", preds[1]) ) 
  points(this.CI$upperCI ~ this.CI[[preds[1]]], pch=16, col="#fd8d3c", cex=0.8)
  points(this.CI$lowerCI ~ this.CI[[preds[1]]], pch=16, col="#d94701", cex=0.8)
  points(this.CI$med ~ this.CI[[preds[1]]], pch=16, col="black", cex=0.8)
  legend("topright", inset=c(-0.18, 0), legend=c("upper CI","median","lower CI"), pch=16, col=c("#fd8d3c", "black", "#d94701"), cex=0.9)

  if (length(preds)==2) {
    plot(this.CI$med ~ this.CI[[preds[2]]], type="n", pch=16, col="black", ylim=c(min(this.CI$lowerCI), max(this.CI$upperCI)),
         main="", ylab=paste("ln(", this.chemical, ")", sep=""), xlab=paste("centered", preds[2]) )
    points(this.CI$upperCI ~ this.CI[[preds[2]]], pch=16, col="#9e9ac8", cex=0.8)
    points(this.CI$lowerCI ~ this.CI[[preds[2]]], pch=16, col="#6a51a3", cex=0.8)
    points(this.CI$med ~ this.CI[[preds[2]]], pch=16, col="black", cex=0.8)
    legend("topright", inset=c(-0.18, 0), legend=c("upper CI","median","lower CI"), pch=16, col=c("#9e9ac8", "black", "#6a51a3"), cex=0.9)
  }
  mtext(paste(this.chemical, "80% Credibility Intervals"),                   # Add main title
        side = 3,
        line = - 2,
        outer = TRUE, cex=1.2, font=2)
  
}


#function to plot credibility intervals - raw version; (median: black; 1st predictor: oranges; 2nd predictor: purples)
plotBayesianCIs.raw <- function(this.CI, this.chemical, preds) {
  
  #calculate raw values for landscape predictors; un-transform CI values and median
  this.CI$pred1_raw <- (this.CI[[preds[1]]] * std.vals$sd[which(std.vals$predictor == preds[1])] +
                          std.vals$mean[which(std.vals$predictor == preds[1])] ) ^ xform.pwr$power[which(xform.pwr$predictor==preds[1])]
  if(length(preds)==2) {
    this.CI$pred2_raw <- (this.CI[[preds[2]]] * std.vals$sd[which(std.vals$predictor == preds[2])] +
                            std.vals$mean[which(std.vals$predictor == preds[2])] ) ^ xform.pwr$power[which(xform.pwr$predictor==preds[2])]
  }
  this.CI$med_raw <- exp(this.CI$med)
  this.CI$upperCI_raw <- exp(this.CI$upperCI)
  this.CI$lowerCI_raw <- exp(this.CI$lowerCI)
  
  #plot median (black) and upper/ lower CI's (red & blue) for untransformed data
  par(mfrow=c(1,length(preds)), mar=c(4,4,3,1))
  plot(this.CI$med_raw ~ this.CI$pred1_raw, type="n", pch=16, col="black", ylim=c(min(this.CI$lowerCI_raw), max(this.CI$upperCI_raw)),
       main="", #paste(this.chemical, "Credibility Intervals"), 
       ylab=paste(this.chemical, "(mg/kg)"), 
       xlab=xform.pwr$rawName[which(xform.pwr$predictor==preds[1])] )
  points(this.CI$upperCI_raw ~ this.CI$pred1_raw, pch=16, col="#fd8d3c", cex=0.8)
  points(this.CI$lowerCI_raw ~ this.CI$pred1_raw, pch=16, col="#d94701", cex=0.8)
  points(this.CI$med_raw ~ this.CI$pred1_raw, pch=16, col="black", cex=0.8)
  legend("topleft", legend=c("upper CI","median","lower CI"), pch=16, col=c("#fd8d3c", "black", "#d94701"), cex=0.9)
  
  
  if(length(preds)==2) {
    plot(this.CI$med_raw ~ this.CI$pred2_raw, pch=16, col="black", ylim=c(min(this.CI$lowerCI_raw), max(this.CI$upperCI_raw)),
       main="", ylab=paste(this.chemical, "(mg/kg)"), 
       xlab=xform.pwr$rawName[which(xform.pwr$predictor==preds[2])] )
    points(this.CI$upperCI_raw ~ this.CI$pred2_raw, pch=16, col="#9e9ac8", cex=0.8)
    points(this.CI$lowerCI_raw ~ this.CI$pred2_raw, pch=16, col="#6a51a3", cex=0.8)
    points(this.CI$med_raw ~ this.CI$pred2_raw, pch=16, col="black", cex=0.8)
    legend("topleft", legend=c("upper CI","median","lower CI"), pch=16, col=c("#9e9ac8", "black", "#6a51a3"), cex=0.9)
  }
  mtext(paste(this.chemical, "80% Credibility Intervals"),                   # Add main title
        side = 3,
        line = - 2,
        outer = TRUE, cex=1.2, font=2)
}



#----------------------------------------------------------------------------------------------#
#  Function to Generate Credibility Intervals, given a min/max and interpolating between them  #  
#                                 random effects = agency                                      #
#----------------------------------------------------------------------------------------------#

#this function accepts information including a COC's brm model, coc2 dataframe, predictors, and which upper and lower credibility 
#  intervals are desired, then generates a tibble with the lower and upper difference from median values for credibility intervals,
#  along with the associated predictor values.
#  To allow proper sampling of the rain/ summer parameter space, 10% of all data are resampled for rain/ summer, such that 11% are summer
#  rain events (corresponding to the percentage of actual data of summer rain events) and 89% are non-summer rain events.

#INPUTS: LCI = lower credibility interval; defaults to 80% CI's
#        UCI = upper credibility interval; defaults to 80% CI's
#        my.brm = Bayesian model for this COC
#        my.coc2 = coc2 data frame for this COC (saved in the frequentist model output for this COC)
#        pred.range = df of min/max landscape predictor values desired for this COC
#        summer = TRUE or FALSE -- was summer used as a predictor for this COC
#        rain = TRUE or FALSE -- was rain used as a predictor for this COC?  Defaults to TRUE
interpCI_old <- function(LCI=0.1, UCI=0.9, my.brm, my.coc2, pred.range, summer, rain=TRUE) {
  
  # #testing!
  # LCI <- 0.1
  # UCI <- 0.9
  # my.brm <- Cu.brm
  # my.coc2 <- Cu.coc2
  # summer <- TRUE
  # rain <- TRUE
  # pred.range <- data.frame(sqrt_traffic=c(-3,5), devAge2=c(-3.19, 3.0), pm25_na=c(-3,3))
  # ######
  
  preds <- colnames(pred.range)
  
  #take 47 samples total (10% of all data we have) to represent feasible rain or rain/summer options
  #   42 samples from non-summer, 5 from summer; this matches the summer/ non-summer ratio of 11% summer.
  #   NOTE: for the end result, it won't matter whether summer is a predictor or not, as we are drawing the same proportion of 
  #         samples from summertime as are in the data; this should approximately happen during sampling, regardless.
  my.sample0 <- my.coc2[which(my.coc2$summer==0)[sample(length(which(my.coc2$summer==0)), 42)], c("rain", "summer")]
  my.sample1 <- my.coc2[which(my.coc2$summer==1)[sample(length(which(my.coc2$summer==1)), 5)], c("rain", "summer")]
  my.sample <- rbind(my.sample0, my.sample1)
  
  #number of samples per landscape predictor for newdata DF
  n.samples <- 12
  p1 <- seq(pred.range[1,1], pred.range[2,1], length.out=n.samples)
  p2 <- seq(pred.range[1,2], pred.range[2,2], length.out=n.samples)
  if (length(preds)==3) {
    p3 <- seq(pred.range[1,3], pred.range[2,3], length.out=n.samples)
  }
  
  #generate simulated data
  if(length(preds)==3) {
    newdata0 <- expand_grid(agency=NA, 
                            rain=my.sample0$rain, 
                            summer=0, 
                            p1, p2, p3)
    newdata1 <- expand_grid(agency=NA, 
                            rain=my.sample1$rain, 
                            summer=1, 
                            p1, p2, p3)
  } else if(length(preds)==2) {
    newdata0 <- expand_grid(agency=NA, 
                            rain=my.sample0$rain, 
                            summer=0, 
                            p1, p2)
    newdata1 <- expand_grid(agency=NA, 
                            rain=my.sample1$rain, 
                            summer=1, 
                            p1, p2)
  }
  newdata <- rbind(newdata0, newdata1)
  names(newdata) <- c("agency", "rain", "summer", preds)  
  
  #conditional mean for an unknown participant, using the population-level intercept plus some variability around where the
  #   agency-specific intercept should be.  This is indicated with re_formula=NULL and allow_new_levels=TRUE.  Obtain all estimates  
  #   (not just summarized data) by setting summary=FALSE
  my.site <- fitted(my.brm, newdata = newdata, re_formula = NULL, allow_new_levels=TRUE, summary=FALSE) %>%
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
  my.summary <- my.site %>%
    group_by(across(all_of( names(my.site)[4:(4+length(preds)-1)]  ))) %>%   #group by the predictors for this COC
    summarise(upperCI = quantile(est, UCI),
              lowerCI = quantile(est, LCI), 
              med=quantile(est, 0.5)) %>%
    mutate(d.upperCI=upperCI-med,
           d.lowerCI=lowerCI-med)
  
  #for COC's with 3 landscape predictors:
  if(length(preds)==3) {
    m <- 200  #final interpolated results will have dimensions [m x m x m] = [xout, yout, zout]
    xout <- rep(seq(pred.range[1,1], pred.range[2,1], length.out=m), each=m^2)
    yout <- rep(rep(seq(pred.range[1,2], pred.range[2,2], length.out=m), each=m), times=m)
    zout <- rep(seq(pred.range[1,3], pred.range[2,3], length.out=m), times=m^2)
    
    x <- p1  #predictor values used for sampling (also in newdata)
    y <- p2
    z <- p3
    
    #three-dimensional interpolation between three predictors
    d.LCI_mat <- array(data = my.summary$d.lowerCI, dim=c(n.samples, n.samples, n.samples))
    xyz_d.LCI <- approx3d(x, y, z, f=d.LCI_mat, xout, yout, zout)
    
    d.UCI_mat <- array(data = my.summary$d.upperCI, dim=c(n.samples, n.samples, n.samples))
    xyz_d.UCI <- approx3d(x, y, z, f=d.UCI_mat, xout, yout, zout)
    
    d.CI_out <- tibble(xout, yout, zout, xyz_d.LCI, xyz_d.UCI)
    names(d.CI_out) <- c(preds, "d.LCI", "d.UCI")
    
  } else if (length(preds)==2) {  #for COC's with 2 landscape predictors
    m <- 200
    xout <- rep(seq(pred.range[1,1], pred.range[2,1], length.out=m), each=m)
    yout <- rep(seq(pred.range[1,2], pred.range[2,2], length.out=m), times=m)
    
    x <- p1
    y <- p2
    
    d.LCI_mat <- array(data = my.summary$d.lowerCI, dim=c(n.samples, n.samples))
    xy_d.LCI <- interp2(x, y, Z=d.LCI_mat, xout, yout, method="linear")
    
    d.UCI_mat <- array(data = my.summary$d.upperCI, dim=c(n.samples, n.samples))
    xy_d.UCI <- interp2(x, y, Z=d.UCI_mat, xout, yout, method="linear")
    
    d.CI_out <- tibble(xout, yout, xy_d.LCI, xy_d.UCI)
    names(d.CI_out) <- c(preds, "d.LCI", "d.UCI")
  }
  
  return(d.CI_out)
}
