# This script takes model validation test data in the form of dataframes of composite 
#   samples and of grab samples collected on multiple occasions over time, and compares 
#   the test data to model predictions fitted to the S8 data.

# In line 196 (building mv.full_res_pred) we remove the data points where 1-day precip was 0 for the day

# Author: Eva Dusek Jennings
# Date: Oct 10, 2026
#------------------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(tidyr)
library(brms)
library(dplyr)
library(cowplot)
library(stringr)
library(gridExtra)


#load prepared model validation data; run script with functions required for this project
load("data/model_validation_prepped_data.RData")
source("scripts/Model Validation Functions.R")  #run script with functions that are used in this script

#source("Model Validation_EIM Data.R")  #run script that preps EIM data
#source("Model Validation_Other Data.R")  #run script that preps other data (including MLK street sweeping data)
#source("Model Validation_Data Prep.R")  #run script that generates comp_res and grab_res dataframes from composite data and EIM grab data

#load RData files with info from Bayesian analyses (including coc2 objects), using the brms models with Location Names from Model Validation
load(file="results/Best Fit Bayesian Models Using Model Validation Location Names.RData")

#read in model validation predictor data, and scaling/centering values. validation predictor data are raw, and require transformation and scaling/centering
vp <- read.csv(file="data/validation_watersheds_mean_predictors_all_raw_2026-09-29-0649.csv")  #in the csv, make sure columns are numbers, with no commas to denote thousands!
stdvals <- read.csv(file="processed_data/spatial_predictor_standardization_values.csv")

#
s8p <- read.csv(file="data/s8_mean_predictors_all_raw_2026-09-29-0649.csv")  #in the csv, make sure columns are numbers, with no commas to denote thousands!


#----------------------------------#
#  Prep Validation Predictor Data  #
#----------------------------------#

head(vp)
head(stdvals)

#these situations don't logically make sense.  Fix them!
vp$dev_1975_1990[which(vp$Name=="WSDOT_I5_MP210.85")] <- 1   #this square of highway is fully paved, but image file indicates that dev_1975_1990=0.36 and all other dev's=0.  Fix it!
vp$dev_pre_1975[which(vp$Name=="WSDOT_I5_MP197.35")] <- 1   #this square of highway is fully paved, but image file indicates that dev_pre_1975=0.89 and all other dev's=0.  Fix it!
vp$impervious[which(vp$Name=="MLK Street Sweep SS5")] <- mean(vp$impervious[which(vp$Name %in% c("MLK Street Sweep SS2", "MLK Street Sweep SS3", "MLK Street Sweep SS4"))])  
   #set SS5 impervious to mean of SS2, SS3, SS4; all should be 100%, but are not...  Maybe tree cover is obscuring it in the satellite photo?

vp <- vp %>%
  mutate(Project = case_when(str_starts(Name, "MLK") ~ "MLK Street Sweep",
                             str_starts(Name, "Coupeville") ~ "Coupeville",
                             str_starts(Name, "GardenEd") ~ "GardenEd",
                             str_starts(Name, "I5") ~ "I5 Bridge",
                             str_starts(Name, "Lahti") ~ "LahtiDr",
                             str_starts(Name, "TacS8C") ~ "Tacoma S8C",
                             str_starts(Name, "WSDOT") ~ "WSDOT"))

#prep the validation predictors by calculating devAge and transforming where applicable.  Also, remove locations where predictors make no sense
vp.mod <- vp %>%
  mutate(devAge = 4*dev_pre_1975 + 3*dev_1975_1990 + 2*dev_1990_2000 + 1*dev_2000_2014,
         oldDev = 2*dev_pre_1975 + 1*dev_1975_1990,
         totRES = urbRES + ruRES,
         roof_totRES = roof_urbRES + roof_ruRES) %>%
  rename(paved = imperv_ground,
         sqrt_popn = sqrt_pop_per_ha,
         sqrt_CO2_res = sqrt_CO2_residential,
         trees = tnc_trees) %>%
  #remove predictors that seem highly out of whack with the actual situation
  # filter(! Location=="520-W") %>%  #520-W should be fully paved and was built in 1963. The image file indicates paved=0 and devAge2=0 (unpaved, undeveloped)
  # filter(! Name== "I5-Bridge_v2") %>%   #remove the wrong version of the I5 Bridge polygon (the right one is a sub-section of the v2 one, per e-mail from Dylan Ahearn Nov 15, 2024)
  # filter(! Project=="Hwy DB") %>%    #remove the Highway DB sites.  I wasn't able to find actual locations for these, so guessed...
  mutate(Location = sub("^[^_]*_", "", Name) ) %>%  #only remove the project name (before the first "_" in the Name column))
  mutate(Location = case_when(Project=="Tacoma S8C" ~ paste("OF", Location, sep=""),  #add "OF" (outfall) in front of Tac S8C location names
                              Project=="MLK Street Sweep" ~ str_remove(Location, "MLK Street Sweep "),  #MLK samples should just be SS2, SS3, etc (NOT MLK Street Sweep SS2)
                              TRUE ~ Location)) %>%
  select(Project, Location, #area_m2, totRES, ln_traffic_k50, sqrt_popn, devAge, paved, trees, roof_totRES,
         ln_traffic_k100, oldDev, sqrt_CO2_res, sqrt_CO2_road, sqrt_CO2_transport, impervious)


#set order of levels of Location for the validation predictors
Loc_levels <- vp.mod$Location[c(5:8, 9:10, 12:19, 1:4, 20:24, 11)]  #coupeville, GardenEd, Lahti, Tacoma, MLK, WSDOT, I5 Bridge
#Loc_levels <- vp.mod$Location[c(5:8, 17:18, 9:10, 12:16, 19:25, 1:4, 26:30, 11)]
vp.mod$Location <- factor(vp.mod$Location, levels=Loc_levels)  #make Location a factor with a set order of levels

#comp_res projects: "MLK Street Sweep", "Tacoma S8C", "Lake WA PFAS", "Hwy DB", "WSDOT", "I5 Bridge"
#grab_res projects: "Coupeville", "LahtiDr", "GardenEd", "MasonCo"

vp.std <- vp.mod %>%
  mutate(
    # totRES=(totRES - stdvals$mean[which(stdvals$X %in% "totRES")])/stdvals$sd[which(stdvals$X %in% "totRES")], 
    # ln_traffic_k50=(ln_traffic_k50 - stdvals$mean[which(stdvals$X %in% "ln_traffic_k50")])/stdvals$sd[which(stdvals$X %in% "ln_traffic_k50")], 
    ln_traffic_k100=(ln_traffic_k100 - stdvals$mean[which(stdvals$X %in% "ln_traffic_k100")])/stdvals$sd[which(stdvals$X %in% "ln_traffic_k100")], 
    # sqrt_popn=(sqrt_popn - stdvals$mean[which(stdvals$X %in% "sqrt_popn")])/stdvals$sd[which(stdvals$X %in% "sqrt_popn")], 
    oldDev=(oldDev - stdvals$mean[which(stdvals$X %in% "oldDev")])/stdvals$sd[which(stdvals$X %in% "oldDev")], 
    # devAge=(devAge - stdvals$mean[which(stdvals$X %in% "devAge")])/stdvals$sd[which(stdvals$X %in% "devAge")],
    sqrt_CO2_res=(sqrt_CO2_res - stdvals$mean[which(stdvals$X %in% "sqrt_CO2_res")])/stdvals$sd[which(stdvals$X %in% "sqrt_CO2_res")], 
    sqrt_CO2_road=(sqrt_CO2_road - stdvals$mean[which(stdvals$X %in% "sqrt_CO2_road")])/stdvals$sd[which(stdvals$X %in% "sqrt_CO2_road")],
    sqrt_CO2_transport=(sqrt_CO2_transport - stdvals$mean[which(stdvals$X %in% "sqrt_CO2_transport")])/stdvals$sd[which(stdvals$X %in% "sqrt_CO2_transport")],
    # paved=(paved - stdvals$mean[which(stdvals$X %in% "paved")])/stdvals$sd[which(stdvals$X %in% "paved")],
    # trees=(trees - stdvals$mean[which(stdvals$X %in% "trees")])/stdvals$sd[which(stdvals$X %in% "trees")], 
    # roof_totRES=(roof_totRES - stdvals$mean[which(stdvals$X %in% "roof_totRES")])/stdvals$sd[which(stdvals$X %in% "roof_totRES")],
    impervious=(impervious - stdvals$mean[which(stdvals$X %in% "impervious")])/stdvals$sd[which(stdvals$X %in% "impervious")]
  )

#summary table of min/ max for each predictor
# vp.std %>%
#   summarize(across(-c(Project, Location), range)) %>%
#   mutate(value = c("min", "max"), .before = 1)

vp.std %>%
  summarize(across(-c(Project, Location), min))

vp.std %>%
  summarize(across(-c(Project, Location), max))




#s8 predictor min/max

#prep the validation predictors by calculating devAge and transforming where applicable.  Also, remove locations where predictors make no sense
s8p.mod <- s8p %>%
  mutate(oldDev = 2*dev_pre_1975 + 1*dev_1975_1990) %>%
  rename(sqrt_CO2_res = sqrt_CO2_residential) %>%
  select(Name, area_m2, ln_traffic_k100, oldDev, sqrt_CO2_res, sqrt_CO2_road, sqrt_CO2_transport, impervious)

s8p.std <- s8p.mod %>%
  mutate(
    ln_traffic_k100=(ln_traffic_k100 - stdvals$mean[which(stdvals$X %in% "ln_traffic_k100")])/stdvals$sd[which(stdvals$X %in% "ln_traffic_k100")], 
    oldDev=(oldDev - stdvals$mean[which(stdvals$X %in% "oldDev")])/stdvals$sd[which(stdvals$X %in% "oldDev")], 
    sqrt_CO2_res=(sqrt_CO2_res - stdvals$mean[which(stdvals$X %in% "sqrt_CO2_res")])/stdvals$sd[which(stdvals$X %in% "sqrt_CO2_res")], 
    sqrt_CO2_road=(sqrt_CO2_road - stdvals$mean[which(stdvals$X %in% "sqrt_CO2_road")])/stdvals$sd[which(stdvals$X %in% "sqrt_CO2_road")],
    sqrt_CO2_transport=(sqrt_CO2_transport - stdvals$mean[which(stdvals$X %in% "sqrt_CO2_transport")])/stdvals$sd[which(stdvals$X %in% "sqrt_CO2_transport")],
    impervious=(impervious - stdvals$mean[which(stdvals$X %in% "impervious")])/stdvals$sd[which(stdvals$X %in% "impervious")] )
    
s8p.std %>%
  summarize(across(-c(Name, area_m2), min))

s8p.std %>%
  summarize(across(-c(Name, area_m2), max))




#make a new object with capped traffic and paved values.  Traffic max is 14.75, paved max is 3.078


# for now, ignore cap level...  use it when we got to the next stage AFTER we've selected best fit frequentist models!
# vp.cap <- vp.std


cap_level <- 3
cap_level_traffic <- 3
#cap_level_paved <- 1
vp.cap <- vp.std %>%
  mutate(ln_traffic_k100 = case_when(ln_traffic_k100 > cap_level_traffic ~ cap_level_traffic,
                                  TRUE ~ ln_traffic_k100),
         impervious = case_when(impervious > cap_level ~ cap_level,
                           TRUE ~ impervious),
         sqrt_CO2_road = case_when(sqrt_CO2_road > cap_level ~ cap_level,
                                   TRUE ~ sqrt_CO2_road),
         sqrt_CO2_res = case_when(sqrt_CO2_res > cap_level ~ cap_level,
                           TRUE ~ sqrt_CO2_res),
         sqrt_CO2_transport = case_when(sqrt_CO2_transport > cap_level ~ cap_level,
                                       TRUE ~ sqrt_CO2_transport),
         oldDev = case_when(oldDev > cap_level ~ cap_level,
                           TRUE ~ oldDev))


#-----------------------------------#
#  Landscape Predictor Value Plots  #
#-----------------------------------#

# plot(vp.std$sqrt_traffic, ylim=c(-5, 15), ylab="standardized sqrt_traffic", xlab="", main="sqrt_traffic values at model validation locations", xaxt="n")
# abline(h=c(-3, 3), col="red")
# abline(h=c(-2, 2), col="orange")
# 
# plot(vp.std$paved, ylim=c(-4, 4), ylab="standardized paved", xlab="", main="paved values at model validation locations", xaxt="n")
# abline(h=c(-3, 3), col="red")
# abline(h=c(-2, 2), col="orange")
# 
# plot(vp.std$sqrt_CO2_road, ylim=c(-4, 4), ylab="standardized sqrt_CO2_road", xlab="", main="sqrt_CO2_road values at model validation locations", xaxt="n")
# abline(h=c(-3, 3), col="red")
# abline(h=c(-2, 2), col="orange")

colors_preds <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
vp.ylim <- c(-3.6, 3.6)
vp.yint <- c(-3,3)

#plot validation predictor values for each standardized predictor (sqrt_traffic, devAge2, paved, sqrt_CO2_road)
ggplot() + 
  geom_point(aes(x=Location, y=sqrt_CO2_res, colour="sqrt_CO2_res"), data=vp.std) +  #add ln-transformed data points from this study
  geom_point(aes(x=Location, y=sqrt_CO2_road, colour="sqrt_CO2_road"), data=vp.std) +  #add ln-transformed data points from this study
  geom_point(aes(x=Location, y=sqrt_CO2_transport, colour="sqrt_CO2_transport"), data=vp.std) +  #add ln-transformed data points from this study
  geom_point(aes(x=Location, y=ln_traffic_k100, colour="ln_traffic_k100"), data=vp.std) +  #add ln-transformed data points from this study
  geom_point(aes(x=Location, y=impervious, colour="impervious"), data=vp.std) +  #add ln-transformed data points from this study
  geom_point(aes(x=Location, y=oldDev, colour="oldDev"), data=vp.std) +  #add ln-transformed data points from this study
  scale_color_manual(values = c(sqrt_CO2_res="#E69F00", sqrt_CO2_road="#56B4E9", sqrt_CO2_transport="#009E73", ln_traffic_k100="#0072B2", impervious="#D55E00", oldDev="#CC79A7")) +
  ylab("standardized predictor values") +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

vp1 <- ggplot() + 
  geom_point(aes(x=Location, y=sqrt_CO2_res), data=vp.std, colour=colors_preds[1], size=3) +  #add ln-transformed data points from this study
  geom_hline(yintercept=vp.yint, color=colors_preds[1], linetype="dashed", size=0.8) +
  xlab("") +   #("Location") +
  ylab("std sqrt_CO2_res") + ylim(vp.ylim) +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        panel.background = element_rect(fill = "gray95", color = "black") )

vp2 <- ggplot() + 
  geom_point(aes(x=Location, y=sqrt_CO2_road), data=vp.std, colour=colors_preds[2], size=3) +  #add ln-transformed data points from this study
  geom_hline(yintercept=vp.yint, color=colors_preds[2], linetype="dashed", size=0.8) +
  xlab("") +   #("Location") +
  ylab("std sqrt_CO2_road") + ylim(vp.ylim) +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        panel.background = element_rect(fill = "gray95", color = "black") )

vp3 <- ggplot() + 
  geom_point(aes(x=Location, y=sqrt_CO2_transport), data=vp.std, colour=colors_preds[3], size=3) +  #add ln-transformed data points from this study
  geom_hline(yintercept=vp.yint, color=colors_preds[3], linetype="dashed", size=0.8) +
  xlab("") +   #("Location") +
  ylab("std sqrt_CO2_transport") + ylim(vp.ylim) +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        panel.background = element_rect(fill = "gray95", color = "black") )

vp4 <- ggplot() + 
  geom_point(aes(x=Location, y=ln_traffic_k100), data=vp.std, colour=colors_preds[4], size=3) +  #add ln-transformed data points from this study
  geom_hline(yintercept=vp.yint, colour=colors_preds[4], linetype="dashed", size=0.8) +
  xlab("") +   #("Location") +
  ylab("std ln_traffic_k100") +  ylim(vp.ylim) +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        panel.background = element_rect(fill = "gray95", color = "black") )

vp5 <- ggplot() + 
  geom_point(aes(x=Location, y=impervious), data=vp.std, colour=colors_preds[5], size=3) +  #add ln-transformed data points from this study
  geom_hline(yintercept=vp.yint, colour=colors_preds[5], linetype="dashed", size=0.8) +
  xlab("") +   #("Location") +
  ylab("std impervious") + ylim(vp.ylim) +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        panel.background = element_rect(fill = "gray95", color = "black") )

vp6 <- ggplot() + 
  geom_point(aes(x=Location, y=oldDev), data=vp.std, colour=colors_preds[6], size=3) +  #add ln-transformed data points from this study
  geom_hline(yintercept=vp.yint, colour=colors_preds[6], linetype="dashed", size=0.8) +
  xlab("") +   #("Location") +
  ylab("std oldDev") + ylim(vp.ylim) +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        panel.background = element_rect(fill = "gray95", color = "black") )

ggarrange(vp1, vp2, vp3, vp4, vp5, vp6, nrow=3, ncol=2)


#-------------------------------------------------#
#  Add Predictor Values to Model Validation Data  #
#-------------------------------------------------#

comp_res.mod <- comp_res %>%
  mutate(Location = case_when(Project=="I5 Bridge" ~ "I5 Bridge",
                              TRUE ~ Location)) %>%
  filter(! is.na(Result),  #remove all results that are "NA"
         ! Project=="Hwy DB")


#add predictor values to the comp_res and grab_res dataframes
comp_res_pred <- left_join(x=comp_res.mod, y=vp.cap, join_by(Project, Location) )
grab_res_pred <- left_join(x=grab_res, y=vp.cap, join_by(Project, Location) )
grab_res_4plus_pred <- left_join(x=grab_res_4plus, y=vp.cap, join_by(Project, Location) )
grab_res_3_pred <- left_join(x=grab_res_3, y=vp.cap, join_by(Project, Location) )

#which locations are missing some predictor data? (Usually b/c the polygon was wrong, or s'thing wasn't working. Paved is arbitrary here)
unique(comp_res_pred$Location[which(is.na(comp_res_pred$impervious))] )
unique(grab_res_pred$Location[which(is.na(grab_res_pred$impervious))] )

#get rid of any rows that don't have predictor data
comp_res_pred <- comp_res_pred[!is.na(comp_res_pred$impervious),]
grab_res_pred <- grab_res_pred[!is.na(grab_res_pred$impervious),]
grab_res_4plus_pred <- grab_res_4plus_pred[!is.na(grab_res_4plus_pred$impervious),]
comp_res_3_pred <- grab_res_3_pred[!is.na(grab_res_3_pred$impervious),]


# #all model validation data, including both composite samples AND discrete samples
# mv.dat.full <- rbind(comp_res_pred, grab_res_pred) %>%
#   mutate(agency=NA,
#          result=log(Result)) %>%
#   rename(location=Location)

#dataframe of the predictors for each location; does not include stormwater outfall data (but, unlike vp.cap, DOES include analyte!)
mv.full_res_pred <- rbind(comp_res_pred, grab_res_pred) %>%
  mutate(Project = as.factor(Project),
         Location = factor(Location, levels=Loc_levels),  #set order of levels of Location for the validation predictors
         location = Location,
         result = log(Result)) %>%
  mutate(agency = case_when(Project=="Tacoma S8C" ~ "Tacoma",
                            TRUE ~ Project),
         agency = as.factor(agency)) %>%
  mutate(spring = case_when(month %in% 4:6 ~ 1,
                            TRUE ~ 0)) %>%
  filter(!daymet_precip_std < -1.48) %>% #For model validation, remove the (Tacoma) samples that were collected when standardized 1-day precip < -1.48 (0mm rainfall)
  filter(daymet_precip > 0) %>%
  rename(rain1 = daymet_precip_std,
         rain14 = daymet_14day_std,
         rain21 = daymet_21day_std,
         dry = antecedant_dry_days_std)
  

# #what does the s8 fitting data look like in terms of daily rainfall (based on Daymet estimates)
# s8 <- read.csv(file="../processed_data/s8data_with_spatial_predictors.csv")
# s8 %>% 
#   filter(daymet_precip < 1) %>% 
#   group_by(loc) %>% 
#   count()  #how many <1mm precip samples at each location?
# 
# # how does that compare to model validation data?
# mv.full_res_pred %>%
#   filter(daymet_precip < 1) %>%
#   group_by(Location) %>%
#   count()  #how many <1mm precip samples at each location?
# ## In comparing the S8 fitting data to the model validation data, there are similar percentages of samples at rainfall levels 
# #    between 1mm and 6mm in both data sets.  Conclude that we should just get rid of the data where daymet_precip <= 0.

mv.full_preds <- mv.full_res_pred %>%
  group_by(Project, agency, Location, location, Analyte) %>%
  select(Project, agency, Location, location, Analyte, sqrt_CO2_res, sqrt_CO2_road, sqrt_CO2_transport, oldDev, ln_traffic_k100, impervious) %>%
  ungroup() %>%
  distinct()

# #dataframe of the predictors for each location; does not include stormwater outfall data (but, unlike vp.cap, DOES include analyte!)
# mv.full_preds <- mv.full_res_pred %>%
#   group_by(Project, Location, Analyte) %>%
#   select(Project, Location, Analyte, sqrt_CO2_road, devAge2, sqrt_traffic, paved) %>%
#   ungroup() %>%
#   distinct()

# #delete later
# Zn.mv.full <- mv.full_res_pred %>%
#   filter(Analyte=="Total Zinc") %>%
#   filter(!Result > 600)
# write.csv(Zn.mv.full, "Zinc model validation data for Teddy.csv")

#---------------------------------------#
#  Model Validation Data Frames by COC  #
#---------------------------------------#

Cu.mv.dat <- mv.full_res_pred %>% 
  filter(Analyte=="Total Copper") %>% 
  mutate(summer = as.factor(summer),
         rawResult = Result) %>%
  select(Project, agency, Location, location, Analyte, Unit, rawResult, result, summer, rain14, dry, ln_traffic_k100, impervious)

Zn.mv.dat <- mv.full_res_pred %>% 
  filter(Analyte=="Total Zinc") %>% 
  mutate(summer = as.factor(summer),
         rawResult = Result) %>%
  select(Project, agency, Location, location, Analyte, Unit, rawResult, result, summer, rain14, sqrt_CO2_transport, oldDev)

TSS.mv.dat <- mv.full_res_pred %>% 
  filter(Analyte=="Total Suspended Solids") %>% 
  mutate(spring = as.factor(spring),
         rawResult = Result) %>%
  select(Project, agency, Location, location, Analyte, Unit, rawResult, result, rain1, rain21, spring, sqrt_CO2_road, oldDev)
  
P.mv.dat <- mv.full_res_pred %>% 
  filter(Analyte=="Total Phosphorus") %>% 
  mutate(summer = as.factor(summer),
         rawResult = Result) %>%
  select(Project, agency, Location, location, Analyte, Unit, rawResult, result, summer, rain14, rain1, oldDev)

TKN.mv.dat <- mv.full_res_pred %>% 
  filter(Analyte=="Total Kjeldahl Nitrogen") %>% 
  mutate(summer = as.factor(summer),
         rawResult = Result) %>%
  select(Project, agency, Location, location, Analyte, Unit, rawResult, result, summer, rain14, rain1, sqrt_CO2_res, oldDev)



#-------------------------------------------#
#  Fitting vs. Model Validation Data Plots  #
#-------------------------------------------#

#zinc - compare range for fitting data and for model validation data
ggplot() + 
  geom_point(aes(x=Location, y=rawResult, colour="ModelValidationData"), data=Zn.mv.dat) +  #add raw model validation data 
  geom_point(aes(x=location, y=exp(result), colour="FittingData"), data=totZn.coc2) +  #add raw fitting data
  scale_color_manual(values = c(ModelValidationData = "purple", FittingData = "orange")) +
  ylab("Total Zinc (ug/kg)") +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))






#----------------------------#
#  Find and Remove Outliers  #
#----------------------------#

colors <- c("red", "orange", "yellow", "green", "darkgreen", "blue", "purple", "magenta", "brown")

#Cleveland Dot Plots - looking for outliers.
dotchart(Cu.mv.dat$result, groups = Cu.mv.dat$Project, 
         pch = 19, color = colors[as.numeric(Cu.mv.dat$Project)],
         xlab = "concentration", main = "Cleveland Dotplot: Copper", ann=TRUE)
dotchart(Zn.mv.dat$result, groups = Zn.mv.dat$Project, pch = 19, col = colors[as.numeric(Zn.mv.dat$Project)],
         xlab = "concentration", main = "Cleveland Dotplot: Zinc")
dotchart(TSS.mv.dat$result, groups = TSS.mv.dat$Project, pch = 19, col = colors[as.numeric(TSS.mv.dat$Project)],
         xlab = "concentration", main = "Cleveland Dotplot: TSS")
dotchart(TKN.mv.dat$result, groups = TKN.mv.dat$Project, pch = 19, col = colors[as.numeric(TKN.mv.dat$Project)],
         xlab = "concentration", main = "Cleveland Dotplot: TKN")
dotchart(P.mv.dat$result, groups = P.mv.dat$Project, pch = 19, col = colors[as.numeric(P.mv.dat$Project)],
         xlab = "concentration", main = "Cleveland Dotplot: Phosphorus")

#Remove outliers from each COC
Cu.mv.dat2 <- Cu.mv.dat %>% filter(!result < -2)  #one super low outlier at Tac OF237B -- ALREADY REMOVED ABOVE w/NO RAIN SAMPLES
Zn.mv.dat2 <- Zn.mv.dat %>% filter(!result > 8)  #two super high outliers at Tac OF235
TSS.mv.dat2 <- TSS.mv.dat %>% filter(!result > 13.5) %>% filter(!result < 7)  #one high outliers at MLK SS4; two really low outliers where daymet_precip <=1mm
TKN.mv.dat2 <- TKN.mv.dat %>% filter(!result > 10)  #one super high outlier at Hwy SR9_MP17.92
P.mv.dat2 <- P.mv.dat %>% filter(!result > 9)  #one super high outlier at Hwy SR9_MP17.92



#---------------------------------------#
#  Color Palettes / Plotting Specifics  #
#---------------------------------------#

orangePalette=c(
# "#d94701", 
  "#fd8d3c", "#fdbe85", "#ffdab9")
#orangePalette=c("#d94701", "#fd8d3c", "#fdbe85", "#feedde")


myColors <- c("#E69F00", "#E69F00", "#E69F00", "#E69F00",  #MLK sites
              "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9",  #Coupeville sites
              "#009E73", "#009E73",  #Garden Ed sites
              "#D55E00",  #I-5 bridge
              "#009E73",  #Lahti Drive
              "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2",  #Tacoma S8C
              "#D55E00", "#D55E00", "#D55E00", "#D55E00", "#D55E00")  #Highway sites

names(myColors) <- vp.std$Location

# [1] SS2             SS3             SS4             SS5             CoupevilleWharf FrontAndMain    BoatLaunch      PerkinsBluff   
# [9] W1              W2              I5 Bridge       LahtiDrive      OF230           OF235           OF237A          OF237B         
# [17] OF243           OF245           OF254           I5_MP197.27     I5_MP197.35     I5_MP210.71     I5_MP210.85     SR9_MP17.92    


#custom manual symbols for model validation locations
mySymbols <- c(15, 16, 17, 18,  #MLK way
               16, 15, 17, 18,  #coupeville
               15, 16,   #garden ed (1-2), I5-bridge, lahti swale
               15, #I5-bridge, 
               17, #lahti swale
               15, 16, 17, 18, 3, 8, 6,  #tacoma
               16, 17, 18, 8, 3)  #WSDOT
names(mySymbols) <- vp.std$Location


#----------------------------------------------------------------#
#  Bar charts with Prediction Intervals & Data Plotted Atop PIs  #
#----------------------------------------------------------------#

theme_set(theme_bw())
#theme_set(theme_gray())  #undo above


#the my.coc2 is only to obtain summer and rainfall information (nothing else, such as location, is used)
mv.full_Cu_PI <- interpPI(my.brm=Cu.brm, my.coc2=Cu.coc2, wr_preds=mv.full_preds[which(mv.full_preds$Analyte=="Total Copper"), c("agency", "location", "ln_traffic_k100", "impervious")])
mv.full_Zn_PI <- interpPI(my.brm=totZn.brm, my.coc2=totZn.coc2, wr_preds=mv.full_preds[which(mv.full_preds$Analyte=="Total Zinc"), c("agency", "location", "oldDev", "sqrt_CO2_transport")])
mv.full_TSS_PI <- interpPI(my.brm=TSS.brm, my.coc2=TSS.coc2, wr_preds=mv.full_preds[which(mv.full_preds$Analyte=="Total Suspended Solids"), c("agency", "location", "oldDev", "sqrt_CO2_road")])
mv.full_TKN_PI <- interpPI(my.brm=TKN.brm, my.coc2=TKN.coc2, wr_preds=mv.full_preds[which(mv.full_preds$Analyte=="Total Kjeldahl Nitrogen"), c("agency", "location", "oldDev", "sqrt_CO2_res")])
mv.full_P_PI <- interpPI(my.brm=P.brm, my.coc2=P.coc2, wr_preds=mv.full_preds[which(mv.full_preds$Analyte=="Total Phosphorus"), c("agency", "location", "oldDev")])


mv.box1 <- plot_PIs_boxplot_raw(mv.full_Cu_PI, Cu.mv.dat2, box.width=0.3, line.width=12) + theme(plot.title=element_blank())
mv.box2 <- plot_PIs_boxplot_raw(mv.full_Zn_PI, Zn.mv.dat2, box.width=0.3, line.width=12) + theme(plot.title=element_blank())
mv.box3 <- plot_PIs_boxplot_raw(mv.full_TSS_PI, TSS.mv.dat2, box.width=0.3, line.width=10) + theme(plot.title=element_blank())
mv.box4 <- plot_PIs_boxplot_raw(mv.full_TKN_PI, TKN.mv.dat2, box.width=0.15, line.width=18) + ylim(-50, 10000) + theme(plot.title=element_blank())
mv.box5 <- plot_PIs_boxplot_raw(mv.full_P_PI, P.mv.dat2, box.width=0.2, line.width=13.5) + theme(plot.title=element_blank()) + ylim(0, 1520)
ggarrange(mv.box1, mv.box2, mv.box3, mv.box4, mv.box5)

mv.box1
mv.box2 
mv.box3 + ylim(0, 350000)
mv.box4
mv.box5


#-----------------------------------#
#  Median per Location Obs vs Pred  #
#-----------------------------------#

# obs.vs.pred.median.mv(Cu.brm, "Copper", Cu.mv.dat2, myColors)
p1.Cu <- obs.vs.pred.median.mv(Cu.brm, "Cu", Cu.mv.dat2, myColors, mySymbols) + xlim(0.7, 4.7) + ylim(0.7, 4.7)  #remove Perkins Bluff
p2.Zn <- obs.vs.pred.median.mv(totZn.brm, "Zn", Zn.mv.dat2, myColors, mySymbols) + xlim(2, 6) + ylim(2, 6) #remove Perkins Bluff
p3.TSS <- obs.vs.pred.median.mv(TSS.brm, "TSS", TSS.mv.dat2, myColors, mySymbols) + xlim(8.7, 11.7) + ylim(8.7, 11.7)
p4.TKN <- obs.vs.pred.median.mv(TKN.brm, "TKN", TKN.mv.dat2, myColors, mySymbols) + xlim(6, 8) + ylim(6, 8)
p5.P <- obs.vs.pred.median.mv(P.brm, "P", P.mv.dat2, myColors, mySymbols) + xlim(3.5, 5.7) + ylim(3.5, 5.7)

p1.Cu
p2.Zn
p3.TSS
p4.TKN
p5.P

# ggarrange(p1.Cu, p2.Zn, p3.TSS, p4.TKN, p5.P)


#--------------------------------#
#  Median per Area, Obs vs Pred  #
#--------------------------------#

# # obs.vs.pred.median.mv(Cu.brm, "Copper", Cu.mv.dat2, myColors)
# p1.Cu.area <- obs.vs.pred.median.area.mv(Cu.brm, "Copper", Cu.mv.dat2[-which(Cu.mv.dat2$location=="PerkinsBluff"),]) + xlim(0,4) + ylim(0,4)  #remove Perkins Bluff
# p2.Zn.area <- obs.vs.pred.median.area.mv(totZn.brm, "Zinc", Zn.mv.dat2[-which(Zn.mv.dat2$location=="PerkinsBluff"),]) + xlim(2.5, 6.2) + ylim(2.5, 6.2) #remove Perkins Bluff
# p3.TSS.area <- obs.vs.pred.median.area.mv(TSS.brm, "TSS", TSS.mv.dat2) + xlim(8, 12) + ylim(8, 12)
# p4.TKN.area <- obs.vs.pred.median.area.mv(TKN.brm, "TKN", TKN.mv.dat2) + xlim(6, 8) + ylim(6, 8)
# p5.P.area <- obs.vs.pred.median.area.mv(P.brm, "Phosphorus", P.mv.dat2) + xlim(3.5, 6) + ylim(3.5, 6)
# 
# ggarrange(p1.Cu.area, p2.Zn.area, p3.TSS.area, p4.TKN.area, p5.P.area)


#---------------#
#  Obs vs Pred  #
#---------------#

# obs.vs.pred.mv(Cu.brm, "Copper", Cu.mv.dat2, myColors)
obs.vs.pred.mv(Cu.brm, "Cu", Cu.mv.dat2, myColors, mySymbols) + xlim(-1.5,5.1) + ylim(-1.5,5.1)
obs.vs.pred.mv(totZn.brm, "Zn", Zn.mv.dat2, myColors, mySymbols) + xlim(0.6, 6.5) + ylim(0.6, 6.5) 
obs.vs.pred.mv(TSS.brm, "TSS", TSS.mv.dat2, myColors, mySymbols) + xlim(7.2, 13.5) + ylim(7.2, 13.5)
obs.vs.pred.mv(TKN.brm, "TKN", TKN.mv.dat2, myColors, mySymbols) + xlim(5.8, 9.2) + ylim(5.8, 9.2)
obs.vs.pred.mv(P.brm, "P", P.mv.dat2, myColors, mySymbols) + xlim(2.8, 8.3) + ylim(2.8, 8.3)


#------------------------------------------#
#  Obs vs Pred showing a select Predictor  #
#------------------------------------------#

#Observed vs Estimated plots for multiple predictors, with colors showing a selected predictor.  
Cu.lim <- c(-1.5, 5.1)
p1 <- obsPredPlot.mv(Cu.brm, Cu.mv.dat2, "Copper", "impervious", "Reds") + labs(color="imper-\nvious") + xlim(Cu.lim) + ylim(Cu.lim)
p2 <- obsPredPlot.mv(Cu.brm, Cu.mv.dat2, "Copper", "ln_traffic_k100", "Greens") + labs(color="ln_\ntraffic_\nk100") + xlim(Cu.lim) + ylim(Cu.lim)
p3 <- obsPredPlot.mv(Cu.brm, Cu.mv.dat2, "Copper", "rain14", "Blues") + xlim(Cu.lim) + ylim(Cu.lim)
p4 <- obsPredPlot.mv(Cu.brm, Cu.mv.dat2, "Copper", "dry", "Oranges") + xlim(Cu.lim) + ylim(Cu.lim)
p5 <- obsPredPlot2.mv(Cu.brm, Cu.mv.dat2, "Copper", "summer", paletteCol=c("#2d708e", "goldenrod2")) + xlim(Cu.lim) + ylim(Cu.lim)  #oranges for summer palette
# grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
lay <- rbind(c(1,1,2,2,NA,NA),
             c(3,3,4,4,5,5))
grid.arrange(p1, p2, p3, p4, p5, layout_matrix=lay)


#Observed vs Predicted plots for multiple predictors, with colors showing a selected predictor.  
Zn.lim <- c(0.6, 6.5)
p1 <- obsPredPlot.mv(totZn.brm, Zn.mv.dat2, "Zn", "oldDev", "Reds") + xlim(Zn.lim) + ylim(Zn.lim)
p2 <- obsPredPlot.mv(totZn.brm, Zn.mv.dat2, "Zn", "sqrt_CO2_transport", "Greens") + labs(color="sqrt_\nCO2_\ntransport") + xlim(Zn.lim) + ylim(Zn.lim)
p3 <- obsPredPlot.mv(totZn.brm, Zn.mv.dat2, "Zn", "rain14", "Blues") + xlim(Zn.lim) + ylim(Zn.lim)
p4 <- obsPredPlot2.mv(totZn.brm, Zn.mv.dat2, "Zn", "summer", paletteCol=c("#2d708e", "goldenrod2")) + xlim(Zn.lim) + ylim(Zn.lim)  #oranges for summer palette
lay <- rbind(c(1,1,2,2),
             c(3,3,4,4))
grid.arrange(p1, p2, p3, p4, layout_matrix=lay)


#Observed vs Predicted plots for multiple predictors, with colors showing a selected predictor.  
TSS.lim <- c(7.2, 13.5)
p1 <- obsPredPlot.mv(TSS.brm, TSS.mv.dat2, "TSS", "sqrt_CO2_road", "Reds") + labs(color="sqrt_\nCO2_\nroad") + xlim(TSS.lim) + ylim(TSS.lim)
p2 <- obsPredPlot.mv(TSS.brm, TSS.mv.dat2, "TSS", "oldDev", "Greens") + xlim(TSS.lim) + ylim(TSS.lim)
p3 <- obsPredPlot.mv(TSS.brm, TSS.mv.dat2, "TSS", "rain1", "Blues") + xlim(TSS.lim) + ylim(TSS.lim)
p4 <- obsPredPlot.mv(TSS.brm, TSS.mv.dat2, "TSS", "rain21", "Purples") + xlim(TSS.lim) + ylim(TSS.lim)
p5 <- obsPredPlot2.mv(TSS.brm, TSS.mv.dat2, "TSS", "spring", paletteCol=c("#25858e", "#c2df23")) + xlim(TSS.lim) + ylim(TSS.lim)  #oranges for summer palette
lay <- rbind(c(1,1,2,2,NA,NA),
             c(3,3,4,4,5,5))
grid.arrange(p1, p2, p3, p4, p5, layout_matrix=lay)


#Observed vs Predicted plots for multiple predictors, with colors showing a selected predictor.  
TKN.lim <- c(5.8, 9.2)
p1 <- obsPredPlot.mv(TKN.brm, TKN.mv.dat2, "TKN", "oldDev", "Reds") + xlim(TKN.lim) + ylim(TKN.lim)
p2 <- obsPredPlot.mv(TKN.brm, TKN.mv.dat2, "TKN", "sqrt_CO2_res", "Greens") + labs(color="sqrt_\nCO2_\nres") + xlim(TKN.lim) + ylim(TKN.lim)
p3 <- obsPredPlot.mv(TKN.brm, TKN.mv.dat2, "TKN", "rain1", "Blues") + xlim(TKN.lim) + ylim(TKN.lim)
p4 <- obsPredPlot.mv(TKN.brm, TKN.mv.dat2, "TKN", "rain14", "Purples") + xlim(TKN.lim) + ylim(TKN.lim)
p5 <- obsPredPlot2.mv(TKN.brm, TKN.mv.dat2, "TKN", "summer", paletteCol=c("#2d708e", "goldenrod2")) + xlim(TKN.lim) + ylim(TKN.lim)  #oranges for summer palette
lay <- rbind(c(1,1,2,2,NA,NA),
             c(3,3,4,4,5,5))
grid.arrange(p1, p2, p3, p4, p5, layout_matrix=lay)


#Observed vs Predicted plots for multiple predictors, with colors showing a selected predictor.  
P.lim <- c(2.8, 8.3)
p1 <- obsPredPlot.mv(P.brm, P.mv.dat2, "P", "oldDev", "Reds") + xlim(P.lim) + ylim(P.lim)
p2 <- obsPredPlot.mv(P.brm, P.mv.dat2, "P", "rain1", "Blues") + xlim(P.lim) + ylim(P.lim)
p3 <- obsPredPlot.mv(P.brm, P.mv.dat2, "P", "rain14", "Purples") + xlim(P.lim) + ylim(P.lim)
p4 <- obsPredPlot2.mv(P.brm, P.mv.dat2, "P", "summer", paletteCol=c("#2d708e", "goldenrod2")) + xlim(P.lim) + ylim(P.lim)  #oranges for summer palette
lay <- rbind(c(1,1,2,2),
             c(3,3,4,4))
grid.arrange(p1, p2, p3, p4, layout_matrix=lay)



#------------------------------------#
#  Single Predictor Plots w/S8 Data  #
#------------------------------------#

myColors <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

thePreds <- c("totRES", "roof_totRES", "sqrt_popn", "ln_traffic_k50", "ln_traffic_k100", "oldDev", "devAge", "trees", "impervious", "paved",
              "sqrt_CO2_road", "sqrt_CO2_res", "sqrt_CO2_transport")

load(file="results/Bayesian_Copper.RData")  #Bayesian model
load(file="results/Frequentist_Copper Models.RData")  #frequentist mixed effects model, for coc2 dataframe

s8 <- Cu.coc2
mv <- Cu.mv.dat2
my.COC <- "Cu"
Form <- Cu.Form4  
Model <- lme(data=s8, Form, method="REML", random=Cu.r1X, weights=Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC, xEqn=0.4)


load(file="results/Bayesian_TSS.RData")  #Bayesian model
load(file="results/Frequentist_TSS Models.RData")  #frequentist mixed effects model, for coc2 dataframe

s8 <- TSS.coc2
mv <- TSS.mv.dat2
my.COC <- "TSS"
Form <- TSS.Form4  
Model <- lme(data=s8, Form, method="REML", random=TSS.r1X, weights=TSS.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC, xEqn=0)


load(file="results/Bayesian_Phosphorus.RData")  #Bayesian model
load(file="results/Frequentist_Total Phosphorus Models.RData")  #frequentist mixed effects model, for coc2 dataframe

s8 <- P.coc2
mv <- P.mv.dat2
my.COC <- "P"
Form <- P.Form4  
Model <- lme(data=s8, Form, method="REML", random=P.r1X, weights=P.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC, xEqn=-0.4)
summary(P.brm)


load(file="results/Bayesian_TotalZinc.RData")  #Bayesian model
load(file="results/Frequentist_Total Zinc Models.RData")  #frequentist mixed effects model, for coc2 dataframe

s8 <- totZn.coc2
mv <- Zn.mv.dat2
my.COC <- "Zn"
Form <- totZn.Form4  
Model <- lme(data=s8, Form, method="REML", random=totZn.r1X, weights=totZn.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC, xEqn=0)
summary(totZn.brm)


load(file="results/Bayesian_TotalKjeldahlNitrogen.RData")  #Bayesian model
load(file="results/Frequentist_Total Kjeldahl Nitrogen Models_censtat.RData")  #frequentist mixed effects model, for coc2 dataframe

s8 <- TKN.coc2
mv <- TKN.mv.dat2
my.COC <- "TKN"
Form <- TKN.Form4 
# TKN.vf1X.final <- varIdent(form= ~1|agency)  #this seems to yield big differences btwn Bayesian values & lme.  Stick with ~1|location here
Model <- lme(data=s8, Form, method="REML", random=TKN.r1X, weights=TKN.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC, xEqn=0.05)
summary(TKN.brm)




#----------------------------   OLD CODE BELOW, JUST FOR REFERENCE  -----------------------------------


# 
# #-------------------------#
# #  Composite Sample Data  #
# #-------------------------#
# 
# #dataframe of the predictors for each location; does not include stormwater outfall data (but, unlike vp.cap, DOES include analyte!)
# comp_preds <- comp_res_pred %>%
#   group_by(Project, Location, Analyte) %>%
#   select(Project, Location, Analyte, sqrt_CO2_road, devAge2, sqrt_traffic, paved) %>%
#   ungroup() %>%
#   distinct()
# 
# comp_Cu_PI <- make_PIs(my.brm=Cu.brm, my.coc2=Cu.coc2, wr.preds=comp_preds[which(comp_preds$Analyte=="Total Copper"), c("Location", "sqrt_traffic", "devAge2")])
# comp_Zn_PI <- make_PIs(my.brm=totZn.brm, my.coc2=totZn.coc2, wr.preds=comp_preds[which(comp_preds$Analyte=="Total Zinc"), c("Location", "sqrt_traffic", "paved")])
# comp_TSS_PI <- make_PIs(my.brm=TSS.brm, my.coc2=TSS.coc2, wr.preds=comp_preds[which(comp_preds$Analyte=="Total Suspended Solids"), c("Location", "sqrt_traffic", "devAge2")])
# comp_TKN_PI <- make_PIs(my.brm=TKN.brm, my.coc2=TKN.coc2, wr.preds=comp_preds[which(comp_preds$Analyte=="Total Kjeldahl Nitrogen"), c("Location", "sqrt_traffic", "devAge2")])
# comp_P_PI <- make_PIs(my.brm=P.brm, my.coc2=P.coc2, wr.preds=comp_preds[which(comp_preds$Analyte=="Total Phosphorus"), c("Location", "sqrt_CO2_road")])
# 
# comp1 <- plot_PIs_raw(comp_Cu_PI, comp_res_pred[which(comp_res_pred$Analyte=="Total Copper"),])
# comp2 <- plot_PIs_raw(comp_Zn_PI, comp_res_pred[which(comp_res_pred$Analyte=="Total Zinc"),]) + ylim(-50, 1500)
# comp3 <- plot_PIs_raw(comp_TSS_PI, comp_res_pred[which(comp_res_pred$Analyte=="Total Suspended Solids"),])
# comp4 <- plot_PIs_raw(comp_TKN_PI, comp_res_pred[which(comp_res_pred$Analyte=="Total Kjeldahl Nitrogen"),]) + ylim(-50, 10000)
# comp5 <- plot_PIs_raw(comp_P_PI, comp_res_pred[which(comp_res_pred$Analyte=="Total Phosphorus"),]) + ylim(-50, 2500)
# ggarrange(comp1, comp2, comp3, comp4, comp5)
# 
# comp1
# comp2
# comp3
# comp4
# comp5
# 
# comp1 + ylim(-10, 150)   #cut out high copper PIs & data
# comp2 + ylim(-50, 3200)  #cut out high zinc PIs & data
# comp3 + ylim(-50, 1e+06)
# comp4 + ylim(-50, 10000)
# comp5 + ylim(-50, 2500)  #cut out the super high phosphorus data points
# 
# 
# 
# #######  ALL ITEMS BELOW STILL INCLUDE THE OUTLIERS!  #########
# 
# #--------------------#
# #  Grab Sample Data  #   ##### NEED TO THINK ABOUT -- Since these are grab samples, which metric should we use?  How should we approach the fact that
# #--------------------#       # we would expect to see much higher variability here than with composite samples?
# 
# #dataframe of the predictors for each location; does not include stormwater outfall data (but, unlike vp.cap, DOES include analyte!)
# grab_preds <- grab_res_pred %>%
#   group_by(Project, Location, Analyte) %>%
#   select(Project, Location, Analyte, sqrt_CO2_road, devAge2, sqrt_traffic, paved) %>%
#   ungroup() %>%
#   distinct()
# 
# grab_preds[which(grab_preds$Analyte=="Total Kjeldahl Nitrogen"),]
# 
# grab_Cu_PI <- make_PIs(my.brm=Cu.brm, my.coc2=Cu.coc2, wr.preds=grab_preds[which(grab_preds$Analyte=="Total Copper"), c("Location", "sqrt_traffic", "devAge2")])
# grab_Zn_PI <- make_PIs(my.brm=totZn.brm, my.coc2=totZn.coc2, wr.preds=grab_preds[which(grab_preds$Analyte=="Total Zinc"), c("Location", "sqrt_traffic", "paved")])
# grab_TSS_PI <- make_PIs(my.brm=TSS.brm, my.coc2=TSS.coc2, wr.preds=grab_preds[which(grab_preds$Analyte=="Total Suspended Solids"), c("Location", "sqrt_traffic", "devAge2")])
# #grab_TKN_PI <- make_PIs(my.brm=TKN.brm, my.coc2=TKN.coc2, wr.preds=grab_preds[which(grab_preds$Analyte=="Total Kjeldahl Nitrogen"), c("Location", "sqrt_traffic", "devAge2")])
# grab_P_PI <- make_PIs(my.brm=P.brm, my.coc2=P.coc2, wr.preds=grab_preds[which(grab_preds$Analyte=="Total Phosphorus"), c("Location", "sqrt_CO2_road")])
# 
# grab1 <- plot_PIs_raw(grab_Cu_PI, grab_res_pred[which(grab_res_pred$Analyte=="Total Copper"),])
# grab2 <- plot_PIs_raw(grab_Zn_PI, grab_res_pred[which(grab_res_pred$Analyte=="Total Zinc"),])
# grab3 <- plot_PIs_raw(grab_TSS_PI, grab_res_pred[which(grab_res_pred$Analyte=="Total Suspended Solids"),])
# #grab4 <- plot_PIs_raw(grab_TKN_PI, grab_res_pred[which(grab_res_pred$Analyte=="Total Kjeldahl Nitrogen"),])
# grab5 <- plot_PIs_raw(grab_P_PI, grab_res_pred[which(grab_res_pred$Analyte=="Total Phosphorus"),])
# ggarrange(grab1, grab2, grab3, grab5)
# 
# grab1
# grab2
# grab3
# #grab4
# grab5
# 
# 
# #----------------------------#
# #  S8 Data Used for Fitting  #
# #----------------------------#
# 
# #dataframe of the predictors for each location; does not include stormwater outfall data (but, unlike vp.cap, DOES include analyte!)
# s8fit_preds_A <- Cu.coc2 %>%  #any COC can be used -- its just for predictors (not for analyte data)
#   group_by(location, agency) %>%
#   select(location, agency, sqrt_CO2_road, devAge2, sqrt_traffic) %>%
#   ungroup() %>%
#   distinct()
# 
# s8fit_preds_B <- totZn.coc2 %>%
#   group_by(location, agency) %>%
#   select(location, agency, paved) %>%
#   ungroup() %>%
#   distinct()
# 
# s8fit_preds <- full_join(s8fit_preds_A, s8fit_preds_B, by=c("location", "agency")) %>%
#   rename(Location = location,
#          Agency = agency)
# 
# 
# s8fit_Cu_PI <- make_PIs(my.brm=Cu.brm, my.coc2=Cu.coc2, wr.preds=s8fit_preds[, c("Location", "sqrt_traffic", "devAge2")])
# s8fit_Zn_PI <- make_PIs(my.brm=totZn.brm, my.coc2=totZn.coc2, wr.preds=s8fit_preds[, c("Location", "sqrt_traffic", "paved")])
# s8fit_TSS_PI <- make_PIs(my.brm=TSS.brm, my.coc2=TSS.coc2, wr.preds=s8fit_preds[, c("Location", "sqrt_traffic", "devAge2")])
# s8fit_TKN_PI <- make_PIs(my.brm=TKN.brm, my.coc2=TKN.coc2, wr.preds=s8fit_preds[, c("Location", "sqrt_traffic", "devAge2")])
# s8fit_P_PI <- make_PIs(my.brm=P.brm, my.coc2=P.coc2, wr.preds=s8fit_preds[, c("Location", "sqrt_CO2_road")])
# 
# s8fit1 <- plot_PIs_raw(s8fit_Cu_PI, Cu.coc2 %>% mutate(Result=exp(result)) %>% rename(Location=location)) + ggtitle("Total Copper")
# s8fit2 <- plot_PIs_raw(s8fit_Zn_PI, totZn.coc2 %>% mutate(Result=exp(result)) %>% rename(Location=location)) + ggtitle("Total Zinc")
# s8fit3 <- plot_PIs_raw(s8fit_TSS_PI, TSS.coc2 %>% mutate(Result=exp(result)) %>% rename(Location=location)) + ggtitle("Total Suspended Solids")
# s8fit4 <- plot_PIs_raw(s8fit_TKN_PI, TKN.coc2 %>% mutate(Result=exp(result)) %>% rename(Location=location)) + ggtitle("Total Kjeldahl Nitrogen")
# s8fit5 <- plot_PIs_raw(s8fit_P_PI, P.coc2 %>% mutate(Result=exp(result)) %>% rename(Location=location)) + ggtitle("Total Phosphorus")
# ggarrange(s8fit1, s8fit2, s8fit3, s8fit4, s8fit5)
# 
# s8fit1
# s8fit2
# s8fit3
# s8fit4
# s8fit5
# 
# #summary table of min/ max for each predictor
# s8fit_preds %>%
#   summarize(across(-c(Agency, Location), range)) %>%
#   mutate(value = c("min", "max"), .before = 1)
# 
# 
# 
# 
# #-------------------------------------------------------------#
# #  Compare Tacoma Rainfall: S8 Data vs Model Validation Data  #
# #-------------------------------------------------------------#
# 
# rain.tac.s8 <- TSS.coc2 %>%
#   filter(agency=="Tacoma") %>%
#   select(location, rain, result)
# 
# rain.tac.mv <- TSS.mv.dat2 %>%
#   filter(Project=="Tacoma S8C") %>%
#   select(location, result, rain)
#   
# par(mfrow=c(1,1))
# plot(rain.tac.s8$rain ~ rain.tac.s8$result, col="purple", pch=16, xlab="ln(TSS)", ylab="std 1-day rainfall", main="Tacoma Sites", xlim=c(6,13), ylim=c(-2, 4))
# points(rain.tac.mv$rain ~ rain.tac.mv$result, col="orange", pch=16)
# legend("topleft", c("s8 fitting data", "model validation data"), col=c("purple", "orange"), pch=16)
# 
# rain.tac.mv[which(rain.tac.mv$rain < -1.3),] #these are the samples with precip=0mm; its about 43 samples. 
# #standardization values for 1-day precip are: mean=10.9374 sd=7.375
# -1.483009*7.37515+10.9374  #rain = -1.483009 is the equivalent of precip = 0mm.  Maybe remove all samples that are < -1.48 for 1-day rain?
# # For model validation, try removing the Tacoma data that was collected when 1-day precip < -1.48
# 
# 
# 



