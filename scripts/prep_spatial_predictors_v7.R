# Get raw and transformed spatial predictors and prepare them by standardizing.
# Note that this version reads in a csv file of spatial predictors from Austin Orr (Geosyntec).

# Landuse percentage predictors are based on WA commerce landuses

# avg_AADT = average annual daily traffic

# Author: Eva Dusek Jennings
# Date: Sep 29, 2026
#---------------------------------------------

library(tidyverse)
library("lattice")

sp <- read.csv(file="data/s8_mean_predictors_all_raw_2026-09-29-0649.csv") %>%            ### Make sure its the right version here!
  dplyr::rename(location=Name) %>%   #rename column "Location_N" to "location"
  dplyr::mutate(loc = case_when(  #add shorthand location names
    location == "KIC_LDR_new" ~ "KIC_LDR",
    location == "SNO_COM_new" ~ "SNO_COM",
    location == "SNO_LDR_new" ~ "SNO_LDR",
    location == "SEA_COM_new" ~ "SEA_COM",
    location == "SEA_HDR_new" ~ "SEA_HDR",
    location == "SEA_IND_new" ~ "SEA_IND",
    location == "PIE_COM_new" ~ "PIE_COM",
    location == "PIE_HDR_new" ~ "PIE_HDR",
    location == "PIE_LDR_new" ~ "PIE_LDR",
    location == "POT_COM_new" ~ "POT_COM",
    TRUE ~ location
    )) %>%
  filter(! (loc %in% c("SEA_HDR_old", "SEA_COM_old")))
 

#create a new column for proportion of development by age as a continuous variable, where the larger
#  the number, the larger the proportion of old development.  Numbers close to 0 either have very new
#  development, or very little development.
sp$devAge <- 4*sp$dev_pre_1975 + 3*sp$dev_1975_1990 + 2*sp$dev_1990_2000 + 1*sp$dev_2000_2014
sp$devAge2 <- 16*sp$dev_pre_1975 + 9*sp$dev_1975_1990 + 4*sp$dev_1990_2000 + 1*sp$dev_2000_2014
sp$oldDev <- 2*sp$dev_pre_1975 + 1*sp$dev_1975_1990
sp$oldDev2 <- 4*sp$dev_pre_1975 + 1*sp$dev_1975_1990

#combine some similar categories
sp$totRES <- sp$ruRES + sp$urbRES  #urban + rural residential
sp$intURB_IND <- sp$intURB + sp$IND  #intensive urban + industrial
sp$roof_intURB_IND <- sp$roof_intURB + sp$roof_IND  #roofs: intensive urban + industrial
sp$roof_totRES <- sp$roof_ruRES + sp$roof_urbRES  #roofs: intensive urban + industrial
sp$greenery <- sp$grass_low_veg + sp$shrub_med_veg + sp$tnc_trees  #low, med and high veg
sp$greenery_bareEarth <- sp$greenery + sp$bare_earth
sp$not_greenBE <- 1-sp$greenery_bareEarth
sp$not_green <- 1-sp$greenery


#look at histograms of the data; determine which columns require log-transformation;
#  look for explanatory variables that have extreme observations; these will likely need
#  a transformation applied to them
ggplot(gather(subset(sp, select=-c(system.index, Description, location, loc, .geo))), aes(value)) +
  geom_histogram(bins=10, fill = "cadetBlue") +
  facet_wrap(~key, scales = "free")

sp$loc <- as.factor(sp$loc)

#use cleveland dotcharts to determine which predictors have extreme values and thus require xforms
op <- par(mfrow=c(3,3), mar=c(3,3,3,1))
dotchart(sp$intURB, main="intURB", group=sp$loc)  #intensive urban
dotchart(sp$IND, main="IND", group=sp$loc)  #industrial
dotchart(sp$urbRES, main="urbRES", group=sp$loc)  #urban residential
dotchart(sp$ruRES, main="ruRES", group=sp$loc)  #rural residential
dotchart(sp$totRES, main="totRES", group=sp$loc)  #total residential (urban + rural res)
dotchart(sp$intURB_IND, main="intURB + IND", group=sp$loc)
#intURB, IND, urbRES and ruRES have too many locations with 0's.  get rid of these as individual predictors and use the aggregated ones instead (totRES, intURB_IND)

op <- par(mfrow=c(3,4), mar=c(3,3,3,1))
dotchart(sp$grass_low_veg, main="grass", group=sp$loc)  
dotchart(sp$shrub_med_veg, main="shrub", group=sp$loc) 
dotchart(sp$greenery, main="greenery", group=sp$loc)
dotchart(sp$imperv_ground, main="paved", group=sp$loc)  
dotchart(sp$imperv_roofs, main="roofs", group=sp$loc)  
dotchart(sp$impervious, main="impervious", group=sp$loc)
dotchart(sp$no_dev, main="nodev", group=sp$loc)  
dotchart(sp$tnc_trees, main="trees", group=sp$loc)  
dotchart(sp$pop_per_ha, main="popn", group=sp$loc) 
dotchart(sp$sqrt_pop_per_ha, main="sqrt popn", group=sp$loc) 
#pop'n should be sqrt_transformed; get rid of shrub - too few shrubs ;)


op <- par(mfrow=c(3,3), mar=c(3,3,3,1))
dotchart(sp$traffic, main="traffic", group=sp$loc)  
dotchart(sp$traffic_k50, main="traffic_k50", group=sp$loc)  
dotchart(sp$traffic_k100, main="traffic_k100", group=sp$loc)  
dotchart(sp$sqrt_traffic, main="sqrt_traffic", group=sp$loc)  
dotchart(sp$sqrt_traffic_k50, main="sqrt_traffic_k50", group=sp$loc)  
dotchart(sp$sqrt_traffic_k100, main="sqrt_traffic_k100", group=sp$loc)  
dotchart(sp$ln_traffic, main="ln_traffic", group=sp$loc)  
dotchart(sp$ln_traffic_k50, main="ln_traffic_k50", group=sp$loc)  
dotchart(sp$ln_traffic_k100, main="ln_traffic_k100", group=sp$loc)  
#try both ln- and sqrt- transformed traffic.  I suspect the ln-transformed will work best.  probably use the k_50 or k_100.


op <- par(mfrow=c(2,2), mar=c(3,3,3,1))
dotchart(sp$greenery, main="greenery", group=sp$loc)
dotchart(sp$greenBE, main="greenery/ bare earth", group=sp$loc)  
dotchart(sp$not_greenBE, main="not green/ bare", group=sp$loc)  
dotchart(sp$impervious, main="impervious", group=sp$loc)
#not_greenBE is identical to impervious.  Remove not_greenBE as well as greenBE (maybe also greenery?)

op <- par(mfrow=c(2,3), mar=c(3,3,3,1))
dotchart(sp$NO_2, main="NO_2", group=sp$loc)
dotchart(sp$pm25, main="pm25", group=sp$loc)  
dotchart(sp$PM25_NA, main="PM25_NA", group=sp$loc)  
dotchart(sp$particulate_surface_area, main="partSA", group=sp$loc)
dotchart(sp$slope, main="slope", group=sp$loc)
#Christian found PM25_NA and preferred it over pm25.  They are somewhat similar; remove pm25, as I don't know where it is from

par(mfrow=c(3,4))
dotchart(sp$CO2_airport, main="CO2_airport", group=sp$loc)
# dotchart(sp$CO2_almostTotal, main="CO2_almostTotal", group=sp$loc)
dotchart(sp$CO2_cmv, main="CO2_cmv", group=sp$loc)
dotchart(sp$CO2_commercial, main="CO2_com", group=sp$loc)
dotchart(sp$CO2_elec_prod, main="CO2_elec", group=sp$loc)
dotchart(sp$CO2_industrial, main="CO2_ind", group=sp$loc)
dotchart(sp$CO2_nonroad, main="CO2_nonroad", group=sp$loc)
dotchart(sp$CO2_rail, main="CO2_rail", group=sp$loc)
dotchart(sp$CO2_residential, main="CO2_res", group=sp$loc)
dotchart(sp$CO2_road, main="CO2_road", group=sp$loc)
dotchart(sp$CO2_total, main="CO2_total", group=sp$loc)
dotchart(sp$CO2_transport, main="CO2_transport", group=sp$loc)
#CO2 airport, cmv, electrical, industrial, and rail should not be individual predictors -- too few sites have them

par(mfrow=c(3,4))
dotchart(sp$sqrt_CO2_airport, main="sqrt_CO2_airport", group=sp$loc)
# dotchart(sp$sqrt_CO2_almostTotal, main="sqrt_CO2_almostTotal", group=sp$loc)
dotchart(sp$sqrt_CO2_cmv, main="sqrt_CO2_cmv", group=sp$loc)
dotchart(sp$sqrt_CO2_commercial, main="sqrt_CO2_com", group=sp$loc)
dotchart(sp$sqrt_CO2_elec_prod, main="sqrt_CO2_elec", group=sp$loc)
dotchart(sp$sqrt_CO2_industrial, main="sqrt_CO2_ind", group=sp$loc)
dotchart(sp$sqrt_CO2_nonroad, main="sqrt_CO2_nonroad", group=sp$loc)
dotchart(sp$sqrt_CO2_rail, main="sqrt_CO2_rail", group=sp$loc)
dotchart(sp$sqrt_CO2_residential, main="sqrt_CO2_res", group=sp$loc)
dotchart(sp$sqrt_CO2_road, main="sqrt_CO2_road", group=sp$loc)
dotchart(sp$sqrt_CO2_total, main="sqrt_CO2_total", group=sp$loc)
dotchart(sp$sqrt_CO2_transport, main="sqrt_CO2_transport", group=sp$loc)
#sqrt_CO2 airport, cmv, electrical, industrial, and rail should not be individual predictors -- too few sites have them
#also, CO2_almostTotal is a weird one, now that additional sources have been added.  Remove it, it doesn't make sense.


par(mfrow=c(3,4))
dotchart(sp$CO2_commercial, main="CO2_com", group=sp$loc)
dotchart(sp$sqrt_CO2_commercial, main="sqrt_CO2_com", group=sp$loc)
dotchart(sp$CO2_nonroad, main="CO2_nonroad", group=sp$loc)
dotchart(sp$sqrt_CO2_nonroad, main="sqrt_CO2_nonroad", group=sp$loc)
dotchart(sp$CO2_residential, main="CO2_res", group=sp$loc)
dotchart(sp$sqrt_CO2_residential, main="sqrt_CO2_res", group=sp$loc)
dotchart(sp$CO2_road, main="CO2_road", group=sp$loc)
dotchart(sp$sqrt_CO2_road, main="sqrt_CO2_road", group=sp$loc)
dotchart(sp$CO2_total, main="CO2_total", group=sp$loc)
dotchart(sp$sqrt_CO2_total, main="sqrt_CO2_total", group=sp$loc)
dotchart(sp$CO2_transport, main="CO2_transport", group=sp$loc)
dotchart(sp$sqrt_CO2_transport, main="sqrt_CO2_transport", group=sp$loc)
#you could really go either way here, wrt transforming or not the CO2 predictors.  Try both sets in the models to see what fits well, then decide?
# this is especially true of CO2_total and CO2_transport.  The others benefit from sqrt_transform

par(mfrow=c(3,3))
dotchart(sp$no_dev, main="nodev", group=sp$loc)  
dotchart(sp$dev_pre_1975, main="dev pre 1975", group=sp$loc)  
dotchart(sp$dev_1975_1990, main="dev 1975-1990", group=sp$loc)  
dotchart(sp$dev_1990_2000, main="dev 1990-2000", group=sp$loc)  
dotchart(sp$dev_2000_2014, main="dev 2000-2014", group=sp$loc)  
dotchart(sp$devAge, main="dev age", group=sp$loc)  
dotchart(sp$oldDev, main="old dev", group=sp$loc)  
dotchart(sp$devAge2, main="dev age2", group=sp$loc)  
dotchart(sp$oldDev2, main="old dev2", group=sp$loc)  
#oldDev looks very similar to devAge2.  dev1990_2000 and dev2000_2014 do not need to be predictors on their own.  Remove devAge2.  keep devAge and oldDev
# try oldDev2 just for curiosity's sake, but probably remove it, as oldDev is just as spread out and no need to multiply oldest development by 4.

par(mfrow=c(2,3))
dotchart(sp$roof_intURB, main="roof_intURB", group=sp$loc)  
dotchart(sp$roof_IND, main="roof_IND", group=sp$loc)  
dotchart(sp$roof_urbRES, main="roof_urbRES", group=sp$loc) 
dotchart(sp$roof_ruRES, main="roof_ruRES", group=sp$loc) 
dotchart(sp$roof_intURB_IND, main="roof_intURB + IND", group=sp$loc)
dotchart(sp$roof_totRES, main="roof_totRES", group=sp$loc)
#get rid of roof_urbRES, roof_ruRes, roof_IND and roof_intURB. Just stick with the summed roofs (roof_totRES and roof_intURB_IND).  Individual types have too many 0's

#create new columns for transformed data
sp_t <- sp %>%
  dplyr::relocate(location, loc) %>%  #relocate "location" column to first column
  dplyr::rename(roofs=imperv_roofs, 
                paved=imperv_ground,
                trees=tnc_trees,
                nodev=no_dev,
                grass=grass_low_veg,
                shrub=shrub_med_veg,
                popn=pop_per_ha,
                sqrt_popn = sqrt_pop_per_ha,
                no2=NO_2,
                CO2_res=CO2_residential,
                CO2_com=CO2_commercial,
                pm25_na=PM25_NA,
                partSA=particulate_surface_area,
                sqrt_CO2_res=sqrt_CO2_residential,
                sqrt_CO2_com=sqrt_CO2_commercial
                ) %>%
  dplyr::select(location,
                loc,
                totRES, 
                intURB_IND,
                grass,
                trees,
                paved,
                roofs,
                impervious,
                nodev, 
                sqrt_traffic_k100,
                ln_traffic_k50,
                ln_traffic_k100,
                sqrt_popn,
                no2,
                pm25_na,
                partSA,
                # CO2_res,
                # CO2_com,
                # CO2_road,
                # CO2_nonroad,
                # CO2_transport,
                # CO2_total,
                sqrt_CO2_res,
                sqrt_CO2_com,
                sqrt_CO2_road,
                sqrt_CO2_nonroad,
                sqrt_CO2_transport,
                sqrt_CO2_total,
                dev_pre_1975,
                dev_1975_1990,
                devAge,
                oldDev,
                # oldDev2,
                roof_intURB_IND, #sqrt_roof_intURB_IND,
                roof_totRES #sqrt_roof_totRES,
  )
                
              

sp_means <- sp_sds <- rep(NA, ncol(sp_t))
sp_means[3:ncol(sp_t)] <- colMeans(sp_t[, 3:ncol(sp_t)])
sp_sds[3:ncol(sp_t)] <- apply(sp_t[,3:ncol(sp_t)], 2, sd)

sp_std <- t((t(sp_t[, 3:ncol(sp_t)])-sp_means[3:ncol(sp_t)])/sp_sds[3:ncol(sp_t)])
sp_std <- data.frame(location=as.factor(sp_t[,c("location")]), loc=as.factor(sp_t[,c("loc")]), sp_std)
sp_final <- sp_std

write.csv(sp_final, file="processed_data/spatial_predictors_standardized.csv", row.names=FALSE)


#save mean and sd for standardization as a csv file
sp_standardization_values <- data.frame(mean=sp_means[3:length(sp_means)],
                                        sd=sp_sds[3:length(sp_sds)])
rownames(sp_standardization_values) <- colnames(sp_t[3:ncol(sp_t)])

write.csv(sp_standardization_values, file="processed_data/spatial_predictor_standardization_values.csv", row.names=TRUE)










#---------#
#  plots  #
#---------#


sp_test <- sp_final %>%
  select(!c("location", "sqrt_traffic_k100", "ln_traffic_k50", "dev_1975_1990", "dev_pre_1975") )

#sp_test <- sp_final[-1]
sp_2 <- pivot_longer(sp_test,
             cols=c(2:ncol(sp_test)))
dotplot(loc~value | name,  data=sp_2)

#dotchart(as.matrix(sp_test[, c(2:ncol(sp_test))]), groups=sp_test$loc)


#create new columns for transformed data
sp_1 <- sp %>%
  dplyr::relocate(loc) %>%  #relocate "loc" column to first column
  dplyr::rename(roofs=imperv_roofs, 
                paved=imperv_ground,
                trees=percent_tree_cover,
                nodev=no_dev,
                grass=grass_low_veg,
                popn=pop_per_ha,
                no2=NO_2,
                CO2_res=CO_emissions_residential,
                CO2_com=CO_emissions_commercial,
                CO2_road=CO_emissions_onroad,
                CO2_nonroad=CO_emissions_nonroad,
                CO2_almostTotal=CO_almostTotal,
                CO2_transport=CO_transport,
                pm25_na=PM25_NA,
                partSA=particulate_surface_area
  ) %>%
  dplyr::select(loc,
                ruRES,
                urbRES,
                totRES, #sqrt_totRES,
                intURB, #sqrt_intURB,
                intURB_IND,
                grass,
                greenery,
                trees,
                greenery_bareEarth,
                not_greenBE,
                paved,
                roofs,
                impervious,
                nodev, #sqrt_nodev,
                traffic, #traffic,
                popn,
                no2,
                #pm25,
                pm25_na,
                partSA,
                slope, #slope,
                CO2_res,
                CO2_com,
                CO2_road,
                CO2_nonroad,
                CO2_almostTotal,
                CO2_transport,
                #dev_pre_1975,
                #dev_1975_1990,
                #dev_1990_2000,
                #dev_2000_2014,
                #devAge, 
                devAge,
                roof_intURB, #sqrt_roof_intURB,
                roof_urbRES, #sqrt_roof_urbRES,
                roof_intURB_IND, #sqrt_roof_intURB_IND,
                roof_totRES #sqrt_roof_totRES,
                #AG
  )

#standardize the sp_1 (untransformed) predictors
sp_1_means <- sp_1_sds <- rep(NA, ncol(sp_1))
sp_1_means[2:ncol(sp_1)] <- colMeans(sp_1[, 2:ncol(sp_1)])
sp_1_sds[2:ncol(sp_1)] <- apply(sp_1[,2:ncol(sp_1)], 2, sd)
sp_1_std <- t((t(sp_1[, 2:ncol(sp_1)])-sp_1_means[2:ncol(sp_1)])/sp_1_sds[2:ncol(sp_1)])
sp_1_std <- data.frame(loc=as.factor(sp_1[,c("loc")]), sp_1_std)

#make the dataframe longer -- combine columns by predictor name
sp_1_std.a <- pivot_longer(sp_1_std,
                           cols=c(2:ncol(sp_1_std)))
names(sp_1_std.a) <- c("loc", "predictor", "standardized_value")




#create new columns for transformed data
sp_2 <- sp_1 %>%
  dplyr::mutate(sqrt_popn=sqrt(popn),
                sqrt_traffic=sqrt(traffic),
                sqrt_slope=sqrt(slope),
                sqrt_CO2_res=sqrt(CO2_res),
                sqrt_CO2_com=sqrt(CO2_com),
                sqrt_CO2_road=sqrt(CO2_road),
                sqrt_CO2_nonroad=sqrt(CO2_nonroad),
                sqrt_CO2_almostTotal=sqrt(CO2_almostTotal),
                sqrt_CO2_transport=sqrt(CO2_transport),
                devAge2=devAge^2,
  ) %>%
  dplyr::select(loc,
                sqrt_traffic, #traffic,
                sqrt_popn,
                sqrt_slope, #slope,
                sqrt_CO2_res,
                sqrt_CO2_com,
                sqrt_CO2_road,
                sqrt_CO2_nonroad,
                sqrt_CO2_almostTotal,
                sqrt_CO2_transport,
                devAge2
  )


#standardize the transformed sp_2 predictors
sp_2_means <- sp_2_sds <- rep(NA, ncol(sp_2))
sp_2_means[2:ncol(sp_2)] <- colMeans(sp_2[, 2:ncol(sp_2)])
sp_2_sds[2:ncol(sp_2)] <- apply(sp_2[,2:ncol(sp_2)], 2, sd)
sp_2_std <- t((t(sp_2[, 2:ncol(sp_2)])-sp_2_means[2:ncol(sp_2)])/sp_2_sds[2:ncol(sp_2)])
sp_2_std <- data.frame(loc=as.factor(sp_2[,c("loc")]), sp_2_std)

#make the dataframe longer -- combine columns by predictor name
sp_2_std.a <- pivot_longer(sp_2_std,
                           cols=c(2:ncol(sp_2)))
names(sp_2_std.a) <- c("loc", "predictor", "standardized_value")



#----------------------------#
#     PLOTS FOR TECH DOC     #
#----------------------------#

#plot of untransformed predictors
dotplot(loc~standardized_value | predictor,  data=sp_1_std.a, col="blue", layout=c(6,6))

#plot of untransformed predictors
dotplot(loc~standardized_value | predictor,  data=sp_2_std.a, col="red", layout=c(6,2))  #layout=(column, rows)





 
