# This script reads in saved Bayesian models for various COCs, and allows credibility 
#  intervals to be generated for a given dataframe of landscape predictor values, using
#  the function interpCI from Bayesian Credibility Interval Function.R

# NOTE: See https://janhove.github.io/visualise_uncertainty/#462_brm() 
#       for definitions of various credibility intervals.
#       For this study, we are using the "Conditional Mean for an Unknown Participant".
#       This predicts the mean COC concentration for a new, unknown location, assuming that
#       the new location is sampled many times.  Since location is nested in agency, we need to 
#       include the random effects, acknowledging that we don't know the intercept.  

# Author: Eva Dusek Jennings
# Revised: Jan 18, 2023
#------------------------------------------------------------------------------------

rm(list=ls(all=T))

#read watershed reductions
wr_psau <- read.csv("../data/watershed_reductions/psau-predictors.csv")
wr_huc12 <- read.csv("../data/watershed_reductions/huc12-predictors.csv")
wr_nhdplus <- read.csv("../data/watershed_reductions/nhdplus-predictors.csv")
# #read (raw) watershed reductions (PSAU version) csv from stormwaterheatmap GitHub page
# wr_psau <- read.csv("https://raw.githubusercontent.com/stormwaterheatmap/TNC_Stormwater_Statistics/main/data/watershed_reductions/psau-predictors.csv")

#add a column to each watershed reduction, that indicates watershed ID; each watershed reduction has a diff't column name for ID
wr_psau$watershed_ID <- wr_psau$AU_ID
wr_huc12$watershed_ID <- wr_huc12$system.index
wr_nhdplus$watershed_ID <- wr_nhdplus$Watershed.Name

# psau - are there any rows of redundant watersheds?  If so, remove them!
n_occur <- data.frame(table(wr_psau$watershed_ID))    
n_occur[n_occur$Freq > 1,]
wr_psau[wr_psau$watershed_ID %in% n_occur$Var1[n_occur$Freq > 1],]  #no repeat occurrences in wr_psau
nrow(wr_psau)  #check if the number of rows = number of unique rows
length(unique(wr_psau$watershed_ID))  

# huc12 - are there any rows of redundant watersheds?  If so, remove them!
n_occur <- data.frame(table(wr_huc12$watershed_ID))   
n_occur[n_occur$Freq > 1,]
wr_huc12[wr_huc12$watershed_ID %in% n_occur$Var1[n_occur$Freq > 1],]  #there are 10 entries with "0.00E+00" in system.index; each has diff't predictor values
wr_huc12 <- wr_huc12[-which(wr_huc12$watershed_ID=="0.00E+00"),]  #remove rows that have repeat watershed_ID values
nrow(wr_huc12)  #check if the number of rows = number of unique rows
length(unique(wr_huc12$watershed_ID))  

# nhdplus - are there any rows of redundant watersheds?  If so, remove them!
n_occur <- data.frame(table(wr_nhdplus$watershed_ID))  
n_occur[n_occur$Freq > 1,]
wr_nhdplus[wr_nhdplus$watershed_ID %in% n_occur$Var1[n_occur$Freq > 1],]  #there are 10 entries with "0.00E+00" in system.index; each has diff't predictor values
nrow(wr_nhdplus)  #check if the number of rows = number of unique rows
length(unique(wr_nhdplus$watershed_ID))  

#combine watershed reductions into a list; 1=psau, 2=huc12, 3=nhdplus
wr.list <- list("psau"=wr_psau, 
                "huc12"=wr_huc12,
                "nhdplus"=wr_nhdplus)
wr.names <- names(wr.list)

rm(wr_psau, wr_huc12, wr_nhdplus)

#read in standardization values for spatial predictors, for unstandardizing and untransforming purposes
std.vals <- read.csv("../processed_data/spatial_predictor_standardization_values.csv", header=TRUE)
names(std.vals)[1] <- "predictor"

#read in backtransform powers for spatial predictors, used for plotting raw spatial predictor values
xform.pwr <- read.csv("../processed_data/spatial_predictor_backtransform_power.csv", header=TRUE)

#read in credibility interval function
source("Bayesian Credibility Interval Function.R")


#----------#
#  Copper  #
#----------#

# load Bayesian copper model: Cu.brm
load(file="../results/Bayesian_Copper.RData")
load(file="../results/Frequentist_Copper Models.RData")

#summary(Cu.brm)$fixed  #provides predictors for this COC

#select which watershed to use;  1=psau; 2=huc12; 3=nhdplus
my.wr <- 2
wr <- wr.list[[my.wr]]
this.wr <- wr.names[my.wr]  #which watershed reduction are we using?

#predictor values for this COC, by watershed
Cu.wr_preds <- wr[, c("watershed_ID", "sqrt_traffic", "devAge2")]
names(Cu.wr_preds) <- c("location", "sqrt_traffic", "devAge2")

#run the 80% credibility interval function for all watersheds
Cu.CI <- interpCI(LCI=0.1, UCI=0.9, Cu.brm, Cu.coc2, Cu.wr_preds)

# if(F) {
#   load(file="../results/CI_80_huc12_Copper.Rdata")   #make sure to adjust my.wr above to match!
# }

#plot transformed & standardized median values (black) with upper/ lower 80% CI's
plotBayesianCIs(Cu.CI, "Copper", c("sqrt_traffic", "devAge2"))

#plot raw median values (black) with upper/ lower 80% CI's
plotBayesianCIs.raw(Cu.CI, "Copper", c("sqrt_traffic", "devAge2"))

#save the credibility intervals for this particular watershed reduction, using the watershed reduction acronym
save(Cu.CI, file= paste("../results/CI_80_", this.wr, "_Copper.Rdata", sep="") )

rm(Cu.brm, Cu.CI)  #remove large items associated with previous calculations


#--------------#
#  Phosphorus  #
#--------------#

# load Bayesian phosphorus model: P.brm and phosphorus coc2 data frame
load(file="../results/Bayesian_Phosphorus.RData")
load(file="../results/Frequentist_Total Phosphorus Models.RData")

summary(P.brm)$fixed

#select which watershed to use;  1=psau; 2=huc12; 3=nhdplus
my.wr <- 3
wr <- wr.list[[my.wr]]
this.wr <- wr.names[my.wr]  #which watershed reduction are we using?

#predictor values for this COC, by watershed
P.wr_preds <- wr[, c("watershed_ID", "sqrt_CO2_road")]
names(P.wr_preds) <- c("location", "sqrt_CO2_road")

#run the 80% credibility interval function for all watersheds
P.CI <- interpCI(LCI=0.1, UCI=0.9, P.brm, P.coc2, P.wr_preds)

#plot transformed & standardized median values (black) with upper/ lower 80% CI's
plotBayesianCIs(P.CI, "Phosphorus", c("sqrt_CO2_road"))

#plot raw median values (black) with upper/ lower 80% CI's
plotBayesianCIs.raw(P.CI, "Phosphorus", c("sqrt_CO2_road"))

#save the credibility intervals for this particular watershed reduction, using the watershed reduction acronym
save(P.CI, file= paste("../results/CI_80_", this.wr, "_Phosphorus.Rdata", sep="") )

rm(P.brm, P.CI)   #remove large items associated with previous calculations


#--------------#
#  Total Zinc  #
#--------------#

# load Bayesian total Zinc model: Zn.brm and phosphorus coc2 data frame
load(file="../results/Bayesian_TotalZinc.RData")
load(file="../results/Frequentist_Total Zinc Models.RData")

summary(totZn.brm)$fixed

#select which watershed to use;  1=psau; 2=huc12; 3=nhdplus
my.wr <- 3
wr <- wr.list[[my.wr]]
this.wr <- wr.names[my.wr]  #which watershed reduction are we using?

#predictor values for this COC, by watershed
totZn.wr_preds <- wr[, c("watershed_ID", "sqrt_traffic", "paved")]
names(totZn.wr_preds) <- c("location", "sqrt_traffic", "paved")

#run the 80% credibility interval function for all watersheds
totZn.CI <- interpCI(LCI=0.1, UCI=0.9, totZn.brm, totZn.coc2, totZn.wr_preds)

#plot transformed & standardized median values (black) with upper/ lower 80% CI's
plotBayesianCIs(totZn.CI, "Total Zinc", c("sqrt_traffic", "paved"))

#plot raw median values (black) with upper/ lower 80% CI's
plotBayesianCIs.raw(totZn.CI, "Total Zinc", c("sqrt_traffic", "paved"))

#save the credibility intervals for this particular watershed reduction, using the watershed reduction acronym
save(totZn.CI, file= paste("../results/CI_80_", this.wr, "_TotalZinc.Rdata", sep="") )

rm(totZn.brm, totZn.CI)   #remove large items associated with previous calculations


#---------------------------#
#  Total Kjeldahl Nitrogen  #
#---------------------------#

# load Bayesian TKN model: TKN.brm and TKN coc2 data frame
load(file="../results/Bayesian_TotalKjeldahlNitrogen.RData")
load(file="../results/Frequentist_Total Kjeldahl Nitrogen Models_censtat.RData")

summary(TKN.brm)$fixed

#select which watershed to use;  1=psau; 2=huc12; 3=nhdplus
my.wr <- 3
wr <- wr.list[[my.wr]]
this.wr <- wr.names[my.wr]  #which watershed reduction are we using?

#predictor values for this COC, by watershed
TKN.wr_preds <- wr[, c("watershed_ID", "sqrt_traffic", "devAge2")]
names(TKN.wr_preds) <- c("location", "sqrt_traffic", "devAge2")

#run the 80% credibility interval function for all watersheds
TKN.CI <- interpCI(LCI=0.1, UCI=0.9, TKN.brm, TKN.coc2, TKN.wr_preds)

#plot transformed & standardized median values (black) with upper/ lower 80% CI's
plotBayesianCIs(TKN.CI, "Total Kjeldahl Nitrogen", c("sqrt_traffic", "devAge2"))

#plot raw median values (black) with upper/ lower 80% CI's
plotBayesianCIs.raw(TKN.CI, "Total Kjeldahl Nitrogen", c("sqrt_traffic", "devAge2"))

#save the credibility intervals for this particular watershed reduction, using the watershed reduction acronym
save(TKN.CI, file= paste("../results/CI_80_", this.wr, "_TotalKjeldahlNitrogen.Rdata", sep="") )

rm(TKN.brm, TKN.CI)   #remove large items associated with previous calculations


#--------------------------#
#  Total Suspended Solids  #
#--------------------------#

# load Bayesian copper model: P.brm and phosphorus coc2 data frame
load(file="../results/Bayesian_TSS.RData")
load(file="../results/Frequentist_TSS Models.RData")

summary(TSS.brm)$fixed

#select which watershed to use;  1=psau; 2=huc12; 3=nhdplus
my.wr <- 1
wr <- wr.list[[my.wr]]
this.wr <- wr.names[my.wr]  #which watershed reduction are we using?

#predictor values for this COC, by watershed
TSS.wr_preds <- wr[, c("watershed_ID", "sqrt_traffic", "devAge2")]
names(TSS.wr_preds) <- c("location", "sqrt_traffic", "devAge2")

#run the 80% credibility interval function for all watersheds
TSS.CI <- interpCI(LCI=0.1, UCI=0.9, TSS.brm, TSS.coc2, TSS.wr_preds)

#plot transformed & standardized median values (black) with upper/ lower 80% CI's
plotBayesianCIs(TSS.CI, "Total Suspended Solids", c("sqrt_traffic", "devAge2"))

#plot raw median values (black) with upper/ lower 80% CI's
plotBayesianCIs.raw(TSS.CI, "Total Suspended Solids", c("sqrt_traffic", "devAge2"))

#save the credibility intervals for this particular watershed reduction, using the watershed reduction acronym
save(TSS.CI, file= paste("../results/CI_80_", this.wr, "_TSS.Rdata", sep="") )

rm(TSS.brm, TSS.CI)   #remove large items associated with previous calculations








#--------------------------  May Not Need To Run Anything Beyond This Point  ------------------------------#


#---------------------------------------------#
#  Truncate CIs Function -- Doesn't Work Yet  #
#---------------------------------------------#

# truncate_CI <- function(CIs, trunc.df, m=200) {
# 
#   #testing
#   CIs <- Cu.CI
#   trunc.df <- data.frame(sqrt_traffic=200, pm25_na=50, devAge2=13)
#   m <- 200
#   ####
# 
#   preds <- names(trunc.df)
# 
#   #weed through to eliminate some predictor values
#   T1a <- unique(CIs[, preds[1]])
#   T1b <- T1a[seq(1, nrow(T1a), (m/trunc.df[,1]) )]
# 
#   T2a <- unique(CIs[, preds[2]])
#   T2b <- T2a[seq(1, nrow(T2a), (m/trunc.df[,2]) )]
# 
#   T3a <- unique(CIs[, preds[3]])
#   T3b <- T3a[seq(1, nrow(T3a), (m/trunc.df[,3]) )]
# 
#   Cu.CI_trunc <- Cu.CI %>%
#     dplyr::filter(devAge2 %in% d) %>%
#     dplyr::filter(pm25_na %in% b)
# }






#--------------------------------#
#  Obtain Credibility Intervals  #
#--------------------------------#

#look at a histogram of rain and summer, to see how they look/ how to divide up our simulated data
hist(Cu.coc2$rain)
hist(as.numeric(Cu.coc2$summer))
#  It doesn't look like there is a simple way to parse this in the newdata frame below.  Particularly considering how
#  much data we will be generating if we do summary=FALSE...

#generate simulated data.
#  NOTE: here, we will want to sample a number of times from the distribution of "rain" values in our data, and have a 
#        more representative perspective on summer vs not summer.  Maybe just sample 20+ times from the data for both rain 
#        AND summer, to capture the possibilities of what can happen?  Find out which number of times to sample by testing
#        several options and seeing how they pan out on a scatterplot and how representative they look of the full data
newdata <- expand_grid(agency=NA, 
                       rain=seq(min(Cu.coc2$rain), max(Cu.coc2$rain), length.out=3), 
                       summer=c(0, 1), 
                       sqrt_traffic=seq(min(Cu.coc2$sqrt_traffic), max(Cu.coc2$sqrt_traffic), length.out=7), 
                       devAge2=seq(min(Cu.coc2$devAge2), max(Cu.coc2$devAge2), length.out=7),
                       pm25_na=seq(min(Cu.coc2$pm25_na), max(Cu.coc2$pm25_na), length.out=7))

# # conditional mean for the population, using the population-level intercept -- this is indicated with re_formula=NA
# #   NOTE: black dot is the mean; if median is desired (and median absolute deviation instead of standard deviation), 
# #         use robust=TRUE
# avg.Cu.pop <- fitted(Cu.brm, newdata = newdata, re_formula = NA, robust=TRUE, probs=c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)) %>%
#   as_tibble() %>%
#   bind_cols(newdata, .) %>%
#   mutate(label=paste0("tr=", round(sqrt_traffic, 2), ", dA=", round(devAge2, 2), ", pm=", round(pm25_na, 2)))
# avg.Cu.pop
# 
# avg.Cu.site <- fitted(Cu.brm, newdata = newdata, re_formula = NULL, robust=TRUE, probs=c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)) %>%
#   as_tibble() %>%
#   bind_cols(newdata, .) %>%
#   mutate(label=paste0("tr=", round(sqrt_traffic, 2), ", dA=", round(devAge2, 2), ", pm=", round(pm25_na, 2)))
# avg.Cu.site

#obtain all estimates (not just summaries) by setting summary=FALSE
Cu.site <- fitted(Cu.brm, newdata = newdata, re_formula = NULL, allow_new_levels=TRUE, summary=FALSE) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(newdata, .) %>%
  pivot_longer(
    cols=starts_with("V"),
    names_to="pars",
    values_to="est"
  )
Cu.site

#summarized data for each type of site (specific values of landscape predictors)
#  this summarizes what might be expected to occur in samples at a site with particular landscape predictor traits 
#  over a representative set of season (summer/not summer) and rainfall value associated with that season
Cu.summ <- Cu.site %>%
  group_by(sqrt_traffic, devAge2, pm25_na) %>%
  summarise(U95=quantile(est, 0.975),
            U80=quantile(est, 0.9),
            U50=quantile(est, 0.75),
            med=quantile(est, 0.5),
            L50=quantile(est, 0.25),
            L80=quantile(est, 0.1),
            L95=quantile(est, 0.025)) %>%
  mutate(d.U95=U95-med,
         d.U80=U80-med,
         d.U50=U50-med,
         d.L50=L50-med,
         d.L80=L80-med,
         d.L95=L95-med)



#---------------------------------------------#
#  Applying credibility intervals to heatmap  #
#---------------------------------------------#

# How do we want to deal with summer and rain?  Summer could easily be drawn from, to represent percent of the year  in summer
# vs not summer.  Alternately, we could assume that the proportion of sapmpling events in summer vs not summer is representative
# of the schedule others might employ for summer sampling vs non-summer sampling.
# Rain is a little more complicated, whereby we would need to know the distribution in amount of rainfall over
# a 3-week period (for copper), then sample from that distribution for the data acquisition.  In fitted.brmsfit, summary can be
# set to FALSE to obtain all rows of data (similar to epred).  These can then be averaged over all parameter estimates/ rainfall/
# summer values for a given set of landscape predictors, to generate the 50%, 80% and 95% credibility intervals.

# June 3, 2022:  We decided to draw from sampled data for rainfall values and for summer/ not summer values, so that we have
#                a representative set of rainfall values to capture the range and frequency of certain rain events in our
#                generated predictions that are then summarized, above.

#three-dimensional interpolation between three predictors
d.U80_mat <- array(data = Cu.summ$d.U80,
                 dim=c(length(unique(Cu.summ$sqrt_traffic)),
                       length(unique(Cu.summ$devAge2)),
                       length(unique(Cu.summ$pm25_na))),
                 dimnames=list(unique(Cu.summ$sqrt_traffic), unique(Cu.summ$devAge2), unique(Cu.summ$pm25_na))
)

m <- 200
xout <- rep(seq(-1.2, 1.3, length.out=m), each=m^2)
yout <- rep(rep(seq(-2.0, 1.0, length.out=m), each=m), times=m)
zout <- rep(seq(-1.2, 1.5, length.out=m), times=m^2)

x <- unique(newdata$sqrt_traffic)
y <- unique(newdata$devAge2)
z <- unique(newdata$pm25_na)


start_time <- Sys.time()
xyz_d.U80 <- approx3d(x, y, z, f=d.U80_mat, xout, yout, zout)
end_time <- Sys.time()
end_time - start_time

Cu.d.U80 <- tibble(sqrt_traffic=xout, devAge2=yout, pm25_na=zout, d.U80=xyz_d.U80)
write.csv(Cu.d.U80, "../results/Copper_diff_U80.csv")
save(Cu.d.U80, file="../results/Copper_diff_U80.Rdata")


#m=10:  0.001412 sec
#m=100: 0.01526 sec
#m=500: 1.2978 sec

par(mfrow=c(2,2))
plot(xout, xyz_U80, pch=16, col="red")
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


#-------------------------------------------------------------------------#
#  How many samples of rain & summer are needed for good representation?  #
#-------------------------------------------------------------------------#

my.sample <- Cu.coc2[sample(nrow(Cu.coc2), 40),]

ggplot(data=Cu.coc2, mapping=aes(x=rain, y=summer)) +
  geom_point(color="black") +
  geom_point(data=my.sample, color="red")

#summer samples make up 11% of total samples
length(which(Cu.coc2$summer==1))/length(which(Cu.coc2$summer==0))

#try taking 30 samples total; 27 from non-summer, 3 from summer
my.sample0 <- Cu.coc2[which(Cu.coc2$summer==0)[sample(length(which(Cu.coc2$summer==0)), 27)], c("rain", "summer")]
my.sample1 <- Cu.coc2[which(Cu.coc2$summer==1)[sample(length(which(Cu.coc2$summer==1)), 3)], c("rain", "summer")]
my.sample <- rbind(my.sample0, my.sample1)

ggplot(data=Cu.coc2, mapping=aes(x=rain, y=summer)) +
  geom_point(color="black") +
  geom_point(data=my.sample0, color="red") +
  geom_point(data=my.sample1, color="orange")
#this looks pretty good as far as diversity of situations goes.  Could also see how various draws affect the outcome, or jsut go with
#  what we have.  Can also get 10% of all samples - that would be 47 total, with 5 from summer and 42 from rainy season

##--------------------------------------------#
##  Applying PREDICTION INTERVALS to heatmap  #
##--------------------------------------------#
# 
# #plot showing the full range of epred values corresponding to the equally spaced devAge2 values 
# ggplot(abc, aes(x = devAge2, y = .epred)) +
#   geom_point()
# 
# #summarize epred values (as 50, 80 and 95% PIs) for unique combinations of the simulated predictor values
# sum_epred <- abc %>%
#   group_by(sqrt_traffic, devAge2, pm25_na) %>%
#   summarise(med=median(.epred), 
#             L95=quantile(.epred, probs=0.025),
#             L80=quantile(.epred, probs=0.1),
#             L50=quantile(.epred, probs=0.25),
#             U50=quantile(.epred, probs=0.75),
#             U80=quantile(.epred, probs=0.9),
#             U95=quantile(.epred, probs=0.975)) %>%
#   mutate(d.L95=L95-med,
#          d.L80=L80-med,
#          d.L50=L50-med,
#          d.U50=U50-med,
#          d.U80=U80-med,
#          d.U95=U95-med)
# 
# 
# #three-dimensional interpolation between three predictors
# #approx3d(x, y, z, f, xout, yout, zout)
# 
# U50_mat <- array(data = sum_epred$U50, 
#                  dim=c(length(unique(sum_epred$sqrt_traffic)), 
#                    length(unique(sum_epred$devAge2)), 
#                    length(unique(sum_epred$pm25_na))), 
#              dimnames=list(unique(sum_epred$sqrt_traffic), unique(sum_epred$devAge2), unique(sum_epred$pm25_na))
# )
# #U50_mat <- as.matrix(U50_mat)
# 
# m <- 500
# xout <- rep(seq(-0.75, 1.25, length.out=m), each=m^2)
# yout <- rep(rep(seq(-2.0, 1.0, length.out=m), each=m), times=m)
# zout <- rep(seq(-1.2, 1.5, length.out=m), times=m^2)
# 
# x <- unique(sum_epred$sqrt_traffic)
# y <- unique(sum_epred$devAge2)
# z <- unique(sum_epred$pm25_na)
# 
# 
# start_time <- Sys.time()
# xyz_U50 <- approx3d(x, y, z, f=U50_mat, xout, yout, zout)
# end_time <- Sys.time()
# end_time - start_time
# 
# #m=10:  0.001412 sec
# #m=100: 0.01526 sec
# #m=500: 1.2978 sec
# 
# par(mfrow=c(2,2))
# plot(xout, xyz_U50, pch=16, col="red")
# for (i in 1:length(x)) {
#   points(rep(x[i], 100), U50_mat[i, , ])
# }
# 
# plot(yout, xyz_U50, pch=16, col="blue")
# for (i in 1:length(y)) {
#   points(rep(y[i], 100), U50_mat[, i, ])
# }
# 
# plot(zout, xyz_U50, pch=16, col="goldenrod")
# for (i in 1:length(y)) {
#   points(rep(z[i], 100), U50_mat[, , i])
# }
# 



#------------------------------------------------------------#
#  Components of Function to Generate Credibility Intervals  #
#------------------------------------------------------------#


#take 47 samples total (10% of all data we have); 42 from non-summer, 5 from summer; this matches the summer/ non-summer ratio of 11% summer
my.sample0 <- Cu.coc2[which(Cu.coc2$summer==0)[sample(length(which(Cu.coc2$summer==0)), 42)], c("rain", "summer")]
my.sample1 <- Cu.coc2[which(Cu.coc2$summer==1)[sample(length(which(Cu.coc2$summer==1)), 5)], c("rain", "summer")]
my.sample <- rbind(my.sample0, my.sample1)


#generate simulated data.
#  NOTE: here, we will want to sample a number of times from the distribution of "rain" values in our data, and have a 
#        more representative perspective on summer vs not summer.  Maybe just sample 20+ times from the data for both rain 
#        AND summer, to capture the possibilities of what can happen?  Find out which number of times to sample by testing
#        several options and seeing how they pan out on a scatterplot and how representative they look of the full data
newdata0 <- expand_grid(agency=NA, 
                       rain=my.sample0$rain, 
                       summer=0, 
                       sqrt_traffic=seq(min(Cu.coc2$sqrt_traffic), max(Cu.coc2$sqrt_traffic), length.out=12), 
                       devAge2=seq(min(Cu.coc2$devAge2), max(Cu.coc2$devAge2), length.out=12),
                       pm25_na=seq(min(Cu.coc2$pm25_na), max(Cu.coc2$pm25_na), length.out=12))

newdata1 <- expand_grid(agency=NA, 
                        rain=my.sample1$rain, 
                        summer=1, 
                        sqrt_traffic=seq(min(Cu.coc2$sqrt_traffic), max(Cu.coc2$sqrt_traffic), length.out=12), 
                        devAge2=seq(min(Cu.coc2$devAge2), max(Cu.coc2$devAge2), length.out=12),
                        pm25_na=seq(min(Cu.coc2$pm25_na), max(Cu.coc2$pm25_na), length.out=12))

newdata <- rbind(newdata0, newdata1)

#conditional mean for an unknown participant, using the population-level intercept plus some variability around where the
#   agency-specific intercept should be.  This is indicated with re_formula=NULL and allow_new_levels=TRUE.  Obtain all estimates  
#   (not just summarized data) by setting summary=FALSE

start_time <- Sys.time()
Cu.site <- fitted(Cu.brm, newdata = newdata, re_formula = NULL, allow_new_levels=TRUE, summary=FALSE) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(newdata, .) %>%
  pivot_longer(
    cols=starts_with("V"),
    names_to="pars",
    values_to="est"
  )
end_time <- Sys.time()
end_time - start_time

#summarized data for each type of site (specific values of landscape predictors)
#  this summarizes what might be expected to occur in samples at a site with particular landscape predictor traits 
#  over a representative set of season (summer/not summer) and rainfall value associated with that season
Cu.summ <- Cu.site %>%
  group_by(sqrt_traffic, devAge2, pm25_na) %>%
  summarise(U95=quantile(est, 0.975),
            U80=quantile(est, 0.9),
            U50=quantile(est, 0.75),
            med=quantile(est, 0.5),
            L50=quantile(est, 0.25),
            L80=quantile(est, 0.1),
            L95=quantile(est, 0.025)) %>%
  mutate(d.U95=U95-med,
         d.U80=U80-med,
         d.U50=U50-med,
         d.L50=L50-med,
         d.L80=L80-med,
         d.L95=L95-med)


#three-dimensional interpolation between three predictors
d.U80_mat <- array(data = Cu.summ$d.U80,
                   dim=c(length(unique(Cu.summ$sqrt_traffic)),
                         length(unique(Cu.summ$devAge2)),
                         length(unique(Cu.summ$pm25_na))),
                   dimnames=list(unique(Cu.summ$sqrt_traffic), unique(Cu.summ$devAge2), unique(Cu.summ$pm25_na))
)

m <- 200
xout <- rep(seq(-1.2, 1.3, length.out=m), each=m^2)
yout <- rep(rep(seq(-2.0, 1.0, length.out=m), each=m), times=m)
zout <- rep(seq(-1.2, 1.5, length.out=m), times=m^2)

x <- unique(newdata$sqrt_traffic)
y <- unique(newdata$devAge2)
z <- unique(newdata$pm25_na)


start_time <- Sys.time()
xyz_d.U80 <- approx3d(x, y, z, f=d.U80_mat, xout, yout, zout)
end_time <- Sys.time()
end_time - start_time

Cu.d.U80 <- tibble(sqrt_traffic=xout, devAge2=yout, pm25_na=zout, d.U80=xyz_d.U80)
#write.csv(Cu.d.U80, "../results/Copper_diff_U80.csv")
#save(Cu.d.U80, file="../results/Copper_diff_U80.Rdata")






# #---------------------------------------#
# #  Generalized Additive Model for CI's  #
# #---------------------------------------#
# 
# # -- NOTE:  using GAM, we don't get an equation -- we get a plot and can make predictions.  Probably not what we want??  Christian
# #           is after an equation, I think...
# 
# install.packages("gam")
# library(gam)
# 
# gam(d.upperCI ~ s(sqrt_traffic, ))
# 
# 
# 
# fit1 <- lm(my.summary$d.upperCI ~ poly(cbind(my.summary$sqrt_traffic, my.summary$devAge2, my.summary$pm25_na), degree=3) )
# 
# #create a scatterplot of x vs. y
# plot(my.summary$sqrt_traffic, my.summary$d.upperCI, pch=19, xlab='sqrt_traffic', ylab='d.upperCI')
# 
# #define x-axis values
# x_axis <- seq(1, 15, length=15)
# 
# #add curve of each model to plot
# lines(x_axis, predict(fit1, data.frame(x=x_axis)), col='green')
# 
# 
