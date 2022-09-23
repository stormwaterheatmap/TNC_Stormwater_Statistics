# This script explores possible plots for Bayesian models
# Part 1 contains a function to plot the final plots for presentation in the tech document
# Part 2 has many different exploratory plots

# Author: Eva Dusek Jennings
# Revised: May 3, 2022
#-----------------------------------------------------------------------------------------

library(tidyr)
library(dplyr)
library(tidyverse)    # ggplot, dplyr, %>%, and friends
library(brms)         # Bayesian modeling through Stan
library(tidybayes)    # Manipulate Stan objects in a tidy way
library(broom)        # Convert model objects to data frames
library(forcats)
library(gridExtra)
library(ggplot2)
library(oce)          # approx3d





#-----------------------------------------#
#                                         #
#  PART I -- Function for Bayesian Plots  #
#                                         #
#-----------------------------------------#

# load Bayesian copper model: Cu.brm
load(file="../results/Bayesian_Copper.RData")

# load frequentist mixed effects model
load(file="../results/Copper Models.RData")


# posterior densities and chain traces for important predictors; max of 5 predictors per page
get_variables(Cu.brm)
plot(Cu.brm, variable=c("b_sqrt_traffic", "b_devAge2", "b_pm25_na", "b_rain:pm25_na"))
plot(Cu.brm, variable=c("b_Intercept", "b_rain", "b_summer1"))
#plot(Cu.brm, variable=c("b_sqrt_traffic", "b_devAge2", "b_pm25_na", "b_rain", "b_summer1", "b_rain:pm25_na"))


#parameter values (coefficient plot)
stanplot(Cu.brm)

#how do the predictors relate to each other?  how are their values distributed in relation to each other? 
ggplot(data=Cu.coc2, mapping=aes(x=sqrt_traffic, y=devAge2, color=pm25_na)) +
  geom_point(alpha=0.5, size=5) +
  scale_colour_distiller(palette="RdBu") +
  labs(x="Sqrt(Traffic)", y="(Development Age)^2", color="PM 2.5")


#-------------------------#
#  Observed vs Predicted  #
#-------------------------#

# observed vs predicted with black points and 95% credible interval lines around each point
fitted(Cu.brm) %>%
  as_tibble() %>%
  bind_cols(Cu.coc2) %>%
  
  ggplot(aes(x = result, y = Estimate)) +        #here, you can sub out <agency> for other predictors
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(size = 1.5, alpha = 3/4) +
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
                 size = 1/4) +
  geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
                     ymax = Estimate + Est.Error),
                 size = 1/2) +
  labs(x = "Observed ln(Copper)", 
       y = "Predicted ln(Copper)") +
  theme_bw() +
  theme(panel.grid = element_blank())


#function to plot observed vs predicted, with colors showing a selected predictor.  Inputs: predictor, color brewer palette
obsPredPlot <- function(pred, paletteCol) {
  fitted(Cu.brm) %>%
    as_tibble() %>%
    bind_cols(Cu.coc2) %>%
    
    ggplot(aes(x = result, y = Estimate, color=get(pred))) +        #here, you can sub out <agency> for other predictors
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_point(size = 1.5, alpha = 3/4) +
    geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
                   size = 1/4) +
    geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
                       ymax = Estimate + Est.Error),
                   size = 1/2) +
    labs(x = "Observed ln(Copper)", 
         y = "Predicted ln(Copper)", color=pred) +
    scale_color_distiller(palette = paletteCol, direction=1) +   #colors should be light for low values, dark for high values
    theme_bw() +
    theme(panel.grid = element_blank())
}

p1 <- obsPredPlot("rain", "Blues")
p2 <- obsPredPlot("sqrt_traffic", "Reds")
p3 <- obsPredPlot("pm25_na", "Greens")
p4 <- obsPredPlot("devAge2", "Oranges")
grid.arrange(p1, p2, p4, p3, nrow=2, ncol=2)


#function to plot observed vs predicted, with colors showing a selected predictor.  Inputs: predictor, color brewer palette
obsPredPlot2 <- function(pred, paletteCol) {
  fitted(Cu.brm) %>%
    as_tibble() %>%
    bind_cols(Cu.coc2) %>%
    
    ggplot(aes(x = result, y = Estimate, color=get(pred))) +        #here, you can sub out <agency> for other predictors
    geom_abline(linetype = 2, color = "grey50", size = .5) +
    geom_point(size = 1.5, alpha = 3/4) +
    geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
                   size = 1/4) +
    geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
                       ymax = Estimate + Est.Error),
                   size = 1/2) +
    labs(x = "Observed ln(Copper)", 
         y = "Predicted ln(Copper)", color=pred) +
    scale_color_brewer(palette="Paired") +   #colors should be light for low values, dark for high values
    theme_bw() +
    theme(panel.grid = element_blank())
}


obsPredPlot("dry", "Oranges")
obsPredPlot2("summer", "Oranges")
obsPredPlot2("month", "Oranges")
obsPredPlot2("year", "Oranges")
obsPredPlot("totRES", "Oranges")
obsPredPlot2("season", "Oranges")
obsPredPlot2("landuse", "Oranges")
obsPredPlot2("agency", "Oranges")




#----------------------------------#
#  Plot Results vs Raw Predictors  #
#----------------------------------#

#read in standardized values
stdVals <- read.csv(file="../processed_data/spatial_predictor_standardization_values.csv")
colnames(stdVals) <- c("predictor", "mean", "sd")

Cu.coc2$devAge_raw <- (Cu.coc2$devAge2 * stdVals$sd[which(stdVals$predictor=="devAge2")] + stdVals$mean[which(stdVals$predictor=="devAge2")] ) ^ 0.5
Cu.coc2$pm25_raw <- Cu.coc2$pm25_na * stdVals$sd[which(stdVals$predictor=="pm25_na")] + stdVals$mean[which(stdVals$predictor=="pm25_na")]
Cu.coc2$traffic_raw <- (Cu.coc2$sqrt_traffic * stdVals$sd[which(stdVals$predictor=="sqrt_traffic")] + stdVals$mean[which(stdVals$predictor=="sqrt_traffic")] ) ^ 2  
  
  
  
#RAW pm 2.5 (unstandardized); epred_draws uses the 4000 parameter sets and feeds them into the various permutations of 
#   predictor values
start.time <- Sys.time()
abc <- epred_draws(Cu.brm, newdata=expand_grid(agency=NA, 
                                               rain=seq(min(Cu.coc2$rain), max(Cu.coc2$rain), length.out=4), 
                                               summer=c(0, 1), 
                                               sqrt_traffic=seq(min(Cu.coc2$sqrt_traffic), max(Cu.coc2$sqrt_traffic), length.out=10), 
                                               devAge2=seq(min(Cu.coc2$devAge2), max(Cu.coc2$devAge2), length.out=10),
                                               pm25_na=seq(min(Cu.coc2$pm25_na), max(Cu.coc2$pm25_na), length.out=10)) )
end.time <- Sys.time()
end.time - start.time

#predictors=10 each, rain=3, summer=2:  3.45 sec
#predictors=10 each, rain=5, summer=2:  8.11 sec
#predictors=12 each, rain=4, summer=2:  14.27 sec
#predictors=15 each, rain=4, summer=2:  28.31 sec
#predictors=15 each, rain=5, summer=2:  VECTOR MEMORY EXHAUSTED!
#predictors=20 each, rain=4, summer=2:  VECTOR MEMORY EXHAUSTED!



abc$pm25.raw <- abc$pm25_na * stdVals$sd[which(stdVals$predictor=="pm25_na")] + stdVals$mean[which(stdVals$predictor=="pm25_na")]
abc$traffic.raw <- (abc$sqrt_traffic * stdVals$sd[which(stdVals$predictor=="sqrt_traffic")] + stdVals$mean[which(stdVals$predictor=="sqrt_traffic")] ) ^ 2
abc$devAge.raw <- (abc$devAge2 * stdVals$sd[which(stdVals$predictor=="devAge2")] + stdVals$mean[which(stdVals$predictor=="devAge2")] ) ^ 0.5

pm25.raw.Cu <- ggplot(abc, aes(x = pm25.raw, y = exp(.epred))) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Particulate Matter 2.5 (Raw)", y = "Estimated Copper Concentration", fill = "Credible Interval") +
  ylim(0, y.max) +
#  geom_point(data=Cu.coc2, mapping=aes(x=pm25_raw, y=exp(result))) +
  theme(legend.position = "bottom")
#print(pm25.raw.Cu)

traffic.raw.Cu <- ggplot(abc, aes(x = traffic.raw, y = exp(.epred))) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "raw traffic", y = "Estimated Copper Concentration", fill = "Credible Interval") +
  ylim(0, y.max) +
#  geom_point(data=Cu.coc2, mapping=aes(x=traffic_raw, y=exp(result))) +
  theme(legend.position = "bottom")
#print(traffic.raw.Cu)

devAge.raw.Cu <- ggplot(abc, aes(x = devAge.raw, y = exp(.epred))) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Development Age (Raw)", y = "Estimated Copper Concentration", fill = "Credible Interval") +
  ylim(0, y.max) +
#  geom_point(data=Cu.coc2, mapping=aes(x=devAge_raw, y=exp(result))) +
  theme(legend.position = "bottom")
#print(devAge.raw.Cu)

grid.arrange(traffic.raw.Cu, pm25.raw.Cu, devAge.raw.Cu, nrow=2, ncol=2)



# standardized, transformed plots -- hmmmm...  looks like the higher values have slightly bigger spread for CI's
pm25.Cu <- ggplot(abc, aes(x = pm25_na, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Particulate Matter 2.5 (Std)", y = "Estimated ln(Copper)", fill = "Credible Interval") +
  theme(legend.position = "bottom")
pm25.Cu

traffic.Cu <- ggplot(abc, aes(x = sqrt_traffic, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Sqrt Traffic (Std)", y = "Estimated ln(Copper)", fill = "Credible Interval") +
  theme(legend.position = "bottom")
traffic.Cu

devAge.Cu <- ggplot(abc, aes(x = devAge2, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Dev Age^2 (Std)", y = "Estimated ln(Copper)", fill = "Credible Interval") +
  theme(legend.position = "bottom")
devAge.Cu


#--------------------------------------------#
#  Applying prediction intervals to heatmap  #
#--------------------------------------------#

#plot showing the full range of epred values corresponding to the equally spaced devAge2 values 
ggplot(abc, aes(x = devAge2, y = .epred)) +
  geom_point()

#summarize epred values (as 50, 80 and 95% PIs) for unique combinations of the simulated predictor values
sum_epred <- abc %>%
  group_by(sqrt_traffic, devAge2, pm25_na) %>%
  summarise(med=median(.epred), 
            L95=quantile(.epred, probs=0.025),
            L80=quantile(.epred, probs=0.1),
            L50=quantile(.epred, probs=0.25),
            U50=quantile(.epred, probs=0.75),
            U80=quantile(.epred, probs=0.9),
            U95=quantile(.epred, probs=0.975)) %>%
  mutate(d.L95=L95-med,
         d.L80=L80-med,
         d.L50=L50-med,
         d.U50=U50-med,
         d.U80=U80-med,
         d.U95=U95-med)


#three-dimensional interpolation between three predictors
#approx3d(x, y, z, f, xout, yout, zout)

U50_mat <- array(data = sum_epred$U50, 
                 dim=c(length(unique(sum_epred$sqrt_traffic)), 
                   length(unique(sum_epred$devAge2)), 
                   length(unique(sum_epred$pm25_na))), 
             dimnames=list(unique(sum_epred$sqrt_traffic), unique(sum_epred$devAge2), unique(sum_epred$pm25_na))
)
#U50_mat <- as.matrix(U50_mat)

m <- 500
xout <- rep(seq(-0.75, 1.25, length.out=m), each=m^2)
yout <- rep(rep(seq(-2.0, 1.0, length.out=m), each=m), times=m)
zout <- rep(seq(-1.2, 1.5, length.out=m), times=m^2)

x <- unique(sum_epred$sqrt_traffic)
y <- unique(sum_epred$devAge2)
z <- unique(sum_epred$pm25_na)


start_time <- Sys.time()
xyz_U50 <- approx3d(x, y, z, f=U50_mat, xout, yout, zout)
end_time <- Sys.time()
end_time - start_time

#m=10:  0.001412 sec
#m=100: 0.01526 sec
#m=500: 1.2978 sec

par(mfrow=c(2,2))
plot(xout, xyz_U50, pch=16, col="red")
for (i in 1:length(x)) {
  points(rep(x[i], 100), U50_mat[i, , ])
}

plot(yout, xyz_U50, pch=16, col="blue")
for (i in 1:length(y)) {
  points(rep(y[i], 100), U50_mat[, i, ])
}

plot(zout, xyz_U50, pch=16, col="goldenrod")
for (i in 1:length(y)) {
  points(rep(z[i], 100), U50_mat[, , i])
}




#------------------------------#
#  Agency-Specific Intercepts  #
#------------------------------#

#agency-specific intercepts with densities, plus global intercept and 98% credibility intervals
Cu.brm %>%
  spread_draws(b_Intercept, r_agency[agency,]) %>%
  mutate(agency_mean = b_Intercept + r_agency) %>%
  ggplot(aes(y = agency, x = agency_mean, fill=stat(abs(x-2.21) < 0.87)))  +
  stat_halfeye() +
  labs(x="Agency-Specific Intercepts", y="", fill="95% CI around global intercept") +
  geom_vline(xintercept = c(2.21), linetype = "solid", color="cadet blue4", alpha=0.5, size=2) +
  annotate("text", x=2.18, y=6, label="global intercept", color="cadet blue4", angle=90, size=6) +
  geom_vline(xintercept = c(1.34, 3.08), linetype = "dashed", color="black", size=1/4) +
  scale_fill_manual(values = c("gray80", "skyblue"))# +




#-----------------------------------------#
#                                         #
#  PART II -- Exploratory Bayesian Plots  #
#                                         #
#-----------------------------------------#


#-------------------------------------------------------------------------------------#
# This section from: https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/  #
#-------------------------------------------------------------------------------------#

library(tidyverse)    # ggplot, dplyr, %>%, and friends
library(brms)         # Bayesian modeling through Stan
library(tidybayes)    # Manipulate Stan objects in a tidy way
library(broom)        # Convert model objects to data frames
#install.packages("broom.mixed")
library(broom.mixed)  # Convert brms model objects to data frames
library(emmeans)      # Calculate marginal effects in even fancier ways
# install.packages("vdemdata")
# library(vdemdata)     # Use data from the Varieties of Democracy (V-Dem) project
#install.packages("patchwork")
library(patchwork)    # Combine ggplot objects
#install.packages("ggokabeito")
library(ggokabeito)   # Neat accessible color palette
#install.packages("gghalves")
library(gghalves)     # Special half geoms
#install.packages("ggbeeswarm")
library(ggbeeswarm)   # Special distribution-shaped point jittering

# Custom ggplot theme to make pretty plots
# Get the News Cycle font at https://fonts.google.com/specimen/News+Cycle
theme_clean <- function() {
  theme_minimal(base_family = "News Cycle") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"))
}

# # load Bayesian copper model: Cu.brm
# load(file="../results/Bayesian_Copper.RData")
# plot(Cu.brm)  #parameter estimate plots for all parameters, plus traces of the sampling space for each parameter 

#tidy(Cu.brm)

newdata <- expand_grid(rain=c(0),
                       summer=c(0),
                       sqrt_traffic=seq(-2, 2, by=0.2),
                       devAge2=seq(-2, 2, by=0.2),
                       pm25_na=seq(-2, 2, by=0.2),
                       agency=NA)

# Predicted values (includes all 3 sources of variability)
tidyPred.Cu <- Cu.brm %>%
  predicted_draws(newdata=newdata)

#Expected values (does not include observation-level variance)
tidyEpred.Cu <- Cu.brm %>%
  epred_draws(newdata=newdata)

plot_preds <- bind_rows(
  "Predicted draws" = tidyPred.Cu,
  "Expectation of predicted draws" = rename(tidyEpred.Cu, .prediction=.epred),
  .id="draw_type") %>%
  dplyr::mutate(draw_type = fct_inorder(draw_type))

# ggplot(plot_preds, aes(x = .prediction, fill = summer)) +
#   stat_halfeye() +
#   labs(x = "Predicted media index", y = "Density", fill = "Summer") +
#   facet_wrap(vars(draw_type)) +
#   scale_fill_okabe_ito() +
#   theme_clean() +
#   theme(legend.position = "bottom")


# Posterior predictions across civil_liberties and regions
all_agencies_traffic.Cu <- Cu.brm %>% 
  epred_draws(newdata = expand_grid(rain=c(0),
                                    summer=c(0),
                                    sqrt_traffic=seq(-2, 2, by=0.2),
                                    devAge2=c(0),  #seq(-2, 2, by=0.2),
                                    pm25_na=c(0),  #seq(-2, 2, by=0.2),
                                    agency=levels(Cu.coc2$agency)), 
              re_formula = NULL)

# plot_all_agencies_traffic.Cu <- ggplot(all_agencies_traffic.Cu, 
#                                   aes(x = sqrt_traffic, y = .epred)) +
#   stat_lineribbon() +
#   scale_fill_brewer(palette = "Reds") +
#   labs(x = "Civil liberties index", y = "Predicted media freedom index",
#        fill = "Credible interval") +
#   facet_wrap(vars(agency)) +
#   theme_clean() +
#   theme(legend.position = "bottom")
# 




# Posterior predictions across civil_liberties and regions
all_predictors.Cu <- Cu.brm %>% 
  epred_draws(newdata = expand_grid(rain=c(0),
                                    summer=c(0),
                                    sqrt_traffic=seq(-2, 2, by=0.5),
                                    devAge2=seq(-2, 2, by=0.5),
                                    pm25_na=seq(-2, 2, by=0.5),
                                    agency=NA), 
              re_formula = NULL)

traffic.Cu <- ggplot(all_predictors.Cu, 
                     aes(x = sqrt_traffic, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "standardized traffic", y = "Estimated Copper prediction",
       fill = "Credible interval") +
#  theme_clean() +
  theme(legend.position = "bottom")

devAge.Cu <- ggplot(all_predictors.Cu, 
                     aes(x = devAge2, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "standardized development age", y = "Estimated Copper prediction",
       fill = "Credible interval") +
#  theme_clean() +
  theme(legend.position = "bottom")

pm25.Cu <- ggplot(all_predictors.Cu, 
                  aes(x = pm25_na, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "standardized pm 2.5", y = "Estimated Copper prediction",
       fill = "Credible interval") +
#  theme_clean() +
  theme(legend.position = "bottom")

grid.arrange(traffic.Cu, devAge.Cu, pm25.Cu, nrow=1, ncol=3)
#ggarrange(traffic.Cu, devAge.Cu, pm25.Cu)



#----------------------------------#
#  Read in Standardization Values  #
#----------------------------------#

stdVals <- read.csv(file="../processed_data/spatial_predictor_standardization_values.csv")
colnames(stdVals) <- c("predictor", "mean", "sd")


#vary one predictor at a time (note: this results in some under-estimation of variability!)

#RAW pm 2.5 (unstandardized)
aa <- epred_draws(Cu.brm, newdata=expand_grid(agency=NA, rain=0, summer=0, sqrt_traffic=0, devAge2=0,
                                              pm25_na=seq(min(Cu.coc2$pm25_na), max(Cu.coc2$pm25_na), length.out=1000)) )
aa$pm25.raw <- aa$pm25_na * stdVals$sd[which(stdVals$predictor=="pm25_na")] + stdVals$mean[which(stdVals$predictor=="pm25_na")]

#RAW traffic (unstandardized, back-transformed) 
bb <- epred_draws(Cu.brm, newdata=expand_grid(agency=NA, rain=0, summer=0, pm25_na=0, devAge2=0,
                                              sqrt_traffic=seq(min(Cu.coc2$sqrt_traffic), max(Cu.coc2$sqrt_traffic), length.out=1000)) )
bb$sqrt_traffic.raw <- bb$sqrt_traffic * stdVals$sd[which(stdVals$predictor=="sqrt_traffic")] + stdVals$mean[which(stdVals$predictor=="sqrt_traffic")]
bb$traffic.raw <- bb$sqrt_traffic.raw ^ 2

#RAW devAge (unstandardized, back-transformed)
cc <- epred_draws(Cu.brm, newdata=expand_grid(agency=NA, rain=0, summer=0, pm25_na=0, sqrt_traffic=0,
                                              devAge2=seq(min(Cu.coc2$devAge2), max(Cu.coc2$devAge2), length.out=1000)) )
cc$devAge2.raw <- cc$devAge2 * stdVals$sd[which(stdVals$predictor=="devAge2")] + stdVals$mean[which(stdVals$predictor=="devAge2")]
cc$devAge.raw <- cc$devAge2.raw ^ 0.5

y.epred <- c(exp(aa$.epred), exp(bb$.epred), exp(cc$.epred) )
quantile(y.epred, probs=0.975)
y.max <- max(  quantile(exp(aa$.epred), 0.99),  quantile(exp(bb$.epred), 0.99), quantile(exp(cc$.epred), 0.99)  )

pm25.raw.Cu <- ggplot(aa, aes(x = pm25.raw, y = exp(.epred))) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "raw pm 2.5", y = "Estimated Copper concentration", fill = "") +
  ylim(0, y.max)
#  theme(legend.position = "bottom")
print(pm25.raw.Cu)

traffic.raw.Cu <- ggplot(bb, aes(x = traffic.raw, y = exp(.epred))) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "raw traffic", y = "", fill = "Credible interval") +
  ylim(0, y.max) +
  theme(legend.position = "bottom")
print(traffic.raw.Cu)

devAge.raw.Cu <- ggplot(cc, aes(x = devAge.raw, y = exp(.epred))) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "raw devAge", y = "Estimated raw Copper concentration", fill = "Credible interval") +
  ylim(0, y.max) 
#  theme(legend.position = "bottom") #+
#  geom_point(data=Cu.coc2, mapping=aes(x=devAge_raw, y=exp(result)))
print(devAge.raw.Cu)


grid.arrange(traffic.raw.Cu, pm25.raw.Cu, devAge.raw.Cu, nrow=2, ncol=2)



#--------------------#
#  Some basic plots  #
#--------------------#

plot(Cu.brm)


get_variables(Cu.brm)
plot(Cu.brm, variable=c("b_sqrt_traffic", "b_devAge2", "b_pm25_na"))
plot(Cu.brm, variable=c("b_rain", "b_summer1", "b_rain:pm25_na"))



plot(conditional_effects(Cu.brm, effects = "sqrt_traffic:pm25_na"))


#-------------------------------------------------------------------#
#  From: http://mjskay.github.io/tidybayes/articles/tidy-brms.html  #
#-------------------------------------------------------------------#

library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(posterior)

theme_set(theme_tidybayes() + panel_border())


get_variables(Cu.brm)
# b_intercept is the global mean, and the r_agency[...] variables are offsets from that mean for each agency

#gather variable indices into a separate column in a tidy format
Cu.brm %>%
  spread_draws(r_agency[agency, term]) %>%
  head(10)

Cu.brm %>%
  spread_draws(b_Intercept, b_sigma_Intercept) %>%
  head(10)

Cu.brm %>%
  spread_draws(b_Intercept, b_sigma_Intercept) %>%
  median_qi(b_Intercept, b_sigma_Intercept)

Cu.brm %>%
  spread_draws(b_Intercept, b_sigma_Intercept) %>%
  median_qi()

Cu.brm %>%
  gather_draws(b_Intercept, b_sigma_Intercept) %>%
  median_qi()

#agency-specific summary information (intercept is offset [difference] from mean intercept)
Cu.brm %>%
  spread_draws(r_agency[agency,]) %>%
  median_qi()

#intercepts for each agency (added to the overall mean intercept)
Cu.brm %>%
  spread_draws(b_Intercept, r_agency[agency,]) %>%
  median_qi(agency_mean = b_Intercept + r_agency)

#plot agency-specific intercept values (mean, 66th percentile, 95th percentile)
Cu.brm %>%
  spread_draws(b_Intercept, r_agency[agency,]) %>%
  median_qi(agency_mean = b_Intercept + r_agency, .width = c(.95, .66)) %>%  #width sets the percentiles around the mean
  ggplot(aes(y = agency, x = agency_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() 

#intervals with densities
Cu.brm %>%
  spread_draws(b_Intercept, r_agency[agency,]) %>%
  mutate(agency_mean = b_Intercept + r_agency) %>%
  ggplot(aes(y = agency, x = agency_mean)) +
  stat_halfeye() +
  labs(x="Agency-Specific Intercepts", y="") 



#-----------------------------------------------------------------------------------------------------------------------#
#  Statistical Rethinking - translated to brms by A. J. Kurz, from: 
#     https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/multivariate-linear-models.html#spurious-associations
#-----------------------------------------------------------------------------------------------------------------------#

# Observed vs Predicted
fitted(Cu.brm) %>%
  as_tibble() %>%
  bind_cols(Cu.coc2) %>%
  
  ggplot(aes(x = result, y = Estimate)) +
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(size = 1.5, color = "firebrick4", alpha = 3/4) +
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),               #the thin lines are the 95% intervals
                 size = 1/4, color = "firebrick4") +
  geom_linerange(aes(ymin = Estimate - Est.Error,              #the thicker lines are +/- the posterior SD
                     ymax = Estimate + Est.Error),
                 size = 1/2, color = "firebrick4") +
  # Note our use of the dot placeholder, here: https://magrittr.tidyverse.org/reference/pipe.html
  # geom_text(data = . %>% filter(Loc %in% c("ID", "UT")),   #this can be used to label certain points, if desired...
  #           aes(label = Loc), 
  #           hjust = 0, nudge_x = - 0.65) +
  labs(x = "Observed ln(Copper)", 
       y = "Predicted ln(Copper)") +
  theme_bw() +
  theme(panel.grid = element_blank())




# ggdist::stat_halfeye(Cu.brm)
# 
# conditional_effects(Cu.brm, surface=TRUE)