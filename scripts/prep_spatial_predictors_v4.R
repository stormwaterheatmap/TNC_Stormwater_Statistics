# Get raw spatial predictors and prepare them by standardizing as needed.
# Note that this version obtains the raw spatial predictors extracted directly from Earth Engine,
# rather than the spatial predictors in the python notebook, that were copied and pasted into a csv file.

# Additionally, v4 also processes landuse percentage predictors, with landuse based on WA commerce landuses

# avg_AADT = average annual daily traffic

# Author: Eva Dusek Jennings
# Date: July 29, 2021
#---------------------------------------------

library(tidyverse)
library(here)
#library(magrittr)

sp <- read.csv(file=here("data", "spatial_predictors_raw.csv")) %>%
  dplyr::rename(location=Location_N) %>%   #rename column "Location_N" to "location"
  dplyr::filter(! (location %in% c("POSOUTFALL_60",  #remove the following sites: Port of Seattle outfall (unrepresentative of other watersheds)
                            "PIEHIRES_OUT",   #                            Pierce County High Residential (data were collected in the middle of a stream)
                            "PIELORES_OUT"))) %>%  #                            Pierce County Low Residential (data collected in the middle of a stream)  
                    #NOTE: Centering and scaling happen in this script, so the mean & sd for standardization will NOT include these sites.
  dplyr::mutate(loc = case_when(  #add shorthand location names
    location == "KICLDRS8D_OUT" ~ "KIC_LDR",
    location == "KICHDRS8D_OUT" ~ "KIC_HDR",
    location == "KICCOMS8D_OUT" ~ "KIC_COM",
    location == "SEAI1S8D_OUT" ~ "SEA_IND",
    location == "SEAC1S8D_OUT" ~ "SEA_COM",
    location == "SEAR1S8D_OUT" ~ "SEA_HDR",
    location == "SNO_COM" ~ "SNO_COM",
    location == "SNO_LDR" ~ "SNO_LDR",
    location == "SNO_HDR" ~ "SNO_HDR",
    location == "TAC001S8D_OF235" ~ "TAC_COM",
    location == "TAC003S8D_OF245" ~ "TAC_IND",
    location == "TFWFD1" ~ "TAC_HDR",
    location == "PIECOMM_OUT" ~ "PIE_COM",
    location == "POT564S8D_OUT" ~ "POT_COM"
    ))
 

#create a new column for proportion of development by age as a continuous variable, where the larger
#  the number, the larger the proportion of old development.  Numbers close to 0 either have very new
#  development, or very little development.
sp$devAge <- 4*sp$dev_pre_1975 + 3*sp$dev_1975_1990 + 2*sp$dev_1990_2000 + 1*sp$dev_2000_2014

#combine some similar categories
sp$totRES <- sp$ruRES + sp$urbRES  #urban + rural residential
sp$intURB_IND <- sp$intURB + sp$IND  #intensive urban + industrial
sp$roof_intURB_IND <- sp$roof_intURB + sp$roof_IND  #roofs: intensive urban + industrial
sp$roof_totRES <- sp$roof_ruRES + sp$roof_urbRES  #roofs: intensive urban + industrial
sp$greenery <- sp$grass_low_veg + sp$shrub_med_veg + sp$percent_tree_cover/100  #low, med and high veg

#look at histograms of the data; determine which columns require log-transformation;
#  look for explanatory variables that have extreme observations; these will likely need
#  a transformation applied to them
ggplot(gather(subset(sp, select=-c(location, loc))), aes(value)) +
  geom_histogram(bins=10, fill = "cadetBlue") +
  facet_wrap(~key, scales = "free")

sp$loc <- as.factor(sp$loc)

#use cleveland dotcharts to determine which predictors have extreme values and thus require xforms
op <- par(mfrow=c(3,3), mar=c(3,3,3,1))
dotchart(sp$AG, main="AG", group=sp$loc)  #agriculture and silviculture
#dotchart(sqrt(sp$AG), main="sqrtAG", group=sp$loc)  #agriculture
dotchart(sp$intURB, main="intURB", group=sp$loc)  #intensive urban
dotchart(sqrt(sp$intURB), main="sqrt intURB", group=sp$loc)  #intensive urban
#dotchart(log(sp$intURB), main="logIntURB", group=sp$loc)  #intensive urban
dotchart(sp$IND, main="IND", group=sp$loc)  #industrial
dotchart(sp$urbRES, main="urbRES", group=sp$loc)  #urban residential
dotchart(sp$ruRES, main="ruRES", group=sp$loc)  #rural residential
dotchart(sp$totRES, main="totRES", group=sp$loc)  #total residential (urban + rural res)
dotchart(sqrt(sp$totRES), main="sqrt_totRES", group=sp$loc)  #total residential (urban + rural res)
dotchart(sp$intURB_IND, main="intURB + IND", group=sp$loc)
#looks like the following require strong transformations: 
#the following may require weak transformations: intURB
#the following are uninformative b/c of low data quantities: AG, ruRES, IND
#transformations applied: sqrt(intURB), sqrt(totRES); keep intURB and totRES as well
#1. predictors to keep: intURB, sqrt(intURB), urbRES, totRES, sqrt(totRES), intURB+IND

op <- par(mfrow=c(3,4), mar=c(3,3,3,1))
dotchart(sp$grass_low_veg, main="grass", group=sp$loc)  
dotchart(sp$shrub_med_veg, main="shrub", group=sp$loc) 
dotchart(sp$greenery, main="greenery", group=sp$loc)
dotchart(sp$imperv_ground, main="paved", group=sp$loc)  
dotchart(sp$imperv_roofs, main="roofs", group=sp$loc)  
dotchart(sp$impervious, main="impervious", group=sp$loc)
#dotchart(log(sp$no_dev*100+1), main="log_nodev", group=sp$loc)  
dotchart(sp$no_dev, main="nodev", group=sp$loc)  
dotchart(sqrt(sp$no_dev), main="sqrt nodev", group=sp$loc)  
dotchart(sp$percent_tree_cover, main="trees", group=sp$loc)  
dotchart(sp$traffic, main="traffic", group=sp$loc)  
dotchart(sp$pop_per_ha, main="popn", group=sp$loc) 
dotchart(sqrt(sp$pop_per_ha), main="sqrt popn", group=sp$loc) 
#dotchart(log(sp$pop_per_ha), main="log popn", group=sp$loc)  #too much
#looks like the following require strong transformations:
#the following may require weak transformations: nodev, popn
#the following are uninformative b/c of low data quantities: shrub
#transformations: sqrt(popn), sqrt(nodev); keep nodev as well...
#2. predictors to keep: grass, greenery, paved, roofs, impervious, sqrt(nodev), nodev,
#   trees, traffic, sqrt(popn)

op <- par(mfrow=c(2,3), mar=c(3,3,3,1))
dotchart(sp$NO_2, main="NO_2", group=sp$loc)
dotchart(sp$pm25, main="pm25", group=sp$loc)  
dotchart(sp$PM25_NA, main="PM25_NA", group=sp$loc)  
dotchart(sp$particulate_surface_area, main="partSA", group=sp$loc)
#dotchart(sqrt(sp$particulate_surface_area), main="sqrt partSA", group=sp$loc)  #no improvement on partSA
dotchart(sp$slope, main="slope", group=sp$loc)
dotchart(sqrt(sp$slope), main="sqrt slope", group=sp$loc)
#looks like the following require strong transformations:
#the following may require weak transformations: nodev, slope, partSA
#the following are uninformative b/c of low data quantities: shrub
#transformations: sqrt(slope); keep slope as well...
#3. predictors to keep: NO_2, pm25, PM25_NA, partSA, slope, sqrt(slope)

par(mfrow=c(3,4))
dotchart(sp$CO_emissions_residential, main="CO2_res", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_residential), main="sqrt CO2_res", group=sp$loc)
dotchart(sp$CO_emissions_total, main="CO2_tot", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_total), main="sqrt CO2_tot", group=sp$loc)
dotchart(sp$CO_emissions_commercial, main="CO2_com", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_commercial), main="sqrt CO2_com", group=sp$loc)
dotchart(sp$CO_emissions_nonroad, main="CO2_nonroad", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_nonroad), main="sqrt CO2_nonroad", group=sp$loc)
dotchart(sp$CO_emissions_onroad, main="CO2_road", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_onroad), main="sqrt CO2_road", group=sp$loc)
#the following are uninformative b/c of low data quantities: CO_emissions_industrial
#transformations: sqrt(CO2_res), sqrt(CO2_tot); sqrt(CO2_road), sqrt(CO2_nonroad), sqrt(CO2_com); 
#     keep CO2_tot as well...
#4. predictors to keep: sqrt(CO2_res), CO2_tot, sqrt(CO2_tot), sqrt(CO2_road), sqrt(CO2_nonroad), 
#     sqrt(CO2_com)

par(mfrow=c(3,2))
#dotchart(sp$no_dev, main="nodev", group=sp$loc)  
dotchart(sp$dev_pre_1975, main="dev pre 1975", group=sp$loc)  
dotchart(sp$dev_1975_1990, main="dev 1975-1990", group=sp$loc)  
dotchart(sp$dev_1990_2000, main="dev 1990-2000", group=sp$loc)  
dotchart(sp$dev_2000_2014, main="dev 2000-2014", group=sp$loc)  
dotchart(sp$devAge, main="dev age", group=sp$loc)  
# dotchart(sqrt(sp$devAge), main="sqrt dev age", group=sp$loc)
# dotchart(log(sp$devAge), main="ln dev age", group=sp$loc)
# dotchart(exp(sp$devAge), main="exp dev age", group=sp$loc)
dotchart((sp$devAge)^2, main="dev age squared", group=sp$loc)
#transformations: devAge^2; keep devAge as well...
#5. predictors to keep: dev_pre_1975, dev_1975_1990, dev_1990_2000, dev_2000_2014,
#   devAge, devAge2

par(mfrow=c(3,3))
#dotchart(sp$roof_AG, main="roof_AG", group=sp$loc)  #only one site, with very low proportion; don't use!
#dotchart(sp$roof_IND, main="roof_IND", group=sp$loc)  
dotchart(sp$roof_intURB, main="roof_intURB", group=sp$loc)  
dotchart(sqrt(sp$roof_intURB), main="sqrt roof_intURB", group=sp$loc)  
#dotchart(sp$roof_ruRES, main="roof_ruRES", group=sp$loc)  
dotchart(sp$roof_urbRES, main="roof_urbRES", group=sp$loc) 
dotchart(sqrt(sp$roof_urbRES), main="sqrt roof_urbRES", group=sp$loc) 
dotchart(sp$roof_intURB_IND, main="roof_intURB + IND", group=sp$loc)
dotchart(sqrt(sp$roof_intURB_IND), main="sqrt roof_intURB + IND", group=sp$loc)
dotchart(sp$roof_totRES, main="roof_totRES", group=sp$loc)
dotchart(sqrt(sp$roof_totRES), main="sqrt roof_totRES", group=sp$loc)
#looks like the following require strong transformations:
#the following may require weak transformations: roof_intURB, roof_urbRES, roof_intURB+IND, roof_totRES
#the following are uninformative b/c of low data quantities: roof_AG, roof_IND, roof_ruRES
#transformations: sqrt(roof_intURB), sqrt(roof_urbRES), sqrt(roof_intURB+IND); sqrt(roof_totRES);
#  keep untransformed versions as well, to see if untransformed data fit well for various COCs
#6. predictors to keep: roof_intURB, sqrt(roof_intURB), roof_urbRES, sqrt(roof_urbRES),
#   roof_intURB+IND, sqrt(roof_intURB+IND), roof_totRES, sqrt(roof_totRES)




#create new columns for transformed data
sp_t <- sp %>%
  dplyr::relocate(location, loc) %>%  #relocate "location" column to first column
  # dplyr::select(-c(loc)) %>%  #remove temporary column "loc"
  dplyr::rename(roofs=imperv_roofs, 
                paved=imperv_ground,
                trees=percent_tree_cover,
                nodev=no_dev,
                grass=grass_low_veg,
                popn=pop_per_ha,
                no2=NO_2,
                CO2_res=CO_emissions_residential,
                CO2_tot=CO_emissions_total,
                CO2_com=CO_emissions_commercial,
                CO2_road=CO_emissions_onroad,
                CO2_nonroad=CO_emissions_nonroad,
                pm25_na=PM25_NA,
                partSA=particulate_surface_area
                ) %>%
  dplyr::mutate(sqrt_intURB=sqrt(intURB),
                sqrt_totRES=sqrt(totRES),
                sqrt_nodev=sqrt(nodev),
                sqrt_popn=sqrt(popn),
                sqrt_slope=sqrt(slope),
                sqrt_CO2_res=sqrt(CO2_res),
                sqrt_CO2_tot=sqrt(CO2_tot),
                sqrt_CO2_com=sqrt(CO2_com),
                sqrt_CO2_road=sqrt(CO2_road),
                sqrt_CO2_nonroad=sqrt(CO2_nonroad),
                devAge2=devAge^2,
                sqrt_roof_intURB=sqrt(roof_intURB),
                sqrt_roof_urbRES=sqrt(roof_urbRES),
                sqrt_roof_intURB_IND=sqrt(roof_intURB_IND),
                sqrt_roof_totRES=sqrt(roof_totRES)
  ) %>%
  dplyr::select(location,
                loc,
                ruRES,
                urbRES,
                totRES, #sqrt_totRES,
                intURB, #sqrt_intURB,
                intURB_IND,
                grass,
                greenery,
                paved,
                roofs,
                impervious,
                nodev, #sqrt_nodev,
                trees,
                traffic,
                sqrt_popn,
                no2,
                #pm25,
                pm25_na,
                partSA,
                sqrt_slope, #slope,
                sqrt_CO2_res,
                #CO2_tot, 
                sqrt_CO2_tot,
                sqrt_CO2_com,
                sqrt_CO2_road,
                sqrt_CO2_nonroad,
                #dev_pre_1975,
                #dev_1975_1990,
                #dev_1990_2000,
                #dev_2000_2014,
                #devAge, 
                devAge2,
                roof_intURB, #sqrt_roof_intURB,
                roof_urbRES, #sqrt_roof_urbRES,
                roof_intURB_IND, #sqrt_roof_intURB_IND,
                roof_totRES #sqrt_roof_totRES,
                #AG
  )
                
              

  #1. predictors to keep: intURB, sqrt(intURB), urbRES, totRES, sqrt(totRES), intURB+IND
  #2. predictors to keep: grass, greenery, paved, roofs, impervious, sqrt(nodev), nodev,
  #   trees, traffic, sqrt(popn)
  #3. predictors to keep: NO_2, pm25, PM25_NA, partSA, slope, sqrt(slope)
  #4. predictors to keep: sqrt(CO2_res), CO2_tot, sqrt(CO2_tot), sqrt(CO2_road), sqrt(CO2_nonroad), 
  #     sqrt(CO2_com)
  #5. predictors to keep: dev_pre_1975, dev_1975_1990, dev_1990_2000, dev_2000_2014,
  #   devAge, devAge2
  #6. predictors to keep: roof_intURB, sqrt(roof_intURB), roof_urbRES, sqrt(roof_urbRES),
  #   roof_intURB_IND, sqrt(roof_intURB_IND), roof_totRES, sqrt(roof_totRES)
  
  
  
sp_means <- sp_sds <- rep(NA, ncol(sp_t))
sp_means[3:ncol(sp_t)] <- colMeans(sp_t[, 3:ncol(sp_t)])
sp_sds[3:ncol(sp_t)] <- apply(sp_t[,3:ncol(sp_t)], 2, sd)

sp_std <- t((t(sp_t[, 3:ncol(sp_t)])-sp_means[3:ncol(sp_t)])/sp_sds[3:ncol(sp_t)])
# sp_std <- data.frame(location=as.factor(sp_t[,c("location")]), sp_std) %>%
#   dplyr::select(c(location,
# #                  loc,
#                   sqrtAG_std=sqrtAG, 
#                   logCOM_std=logCOM, 
#                   IND_std=IND, 
#                   OPEN_std=OPEN, 
#                   sqrtRES_std=sqrtRES, 
#                   TRANS_std=TRANS, 
#                   logCOM_IND_std=logCOM_IND,
#                   dev_pre_1975_std=dev_pre_1975, 
#                   dev_1975_1990_std=dev_1975_1990, 
#                   dev_1990_2000_std=dev_1990_2000, 
#                   dev_2000_2014_std=dev_2000_2014, 
#                   grass_std=grass, 
#                   paved_std=paved, 
#                   roofs_std=roofs, 
#                   impervious_std=impervious,
#                   no2_std=no2,
#                   sqrtNodev_std=sqrtNodev, 
#                   logTrees_std=logTrees,
#                   pm25_std=pm25, 
#                   sqrtSlope_std=sqrtSlope, 
#                   sqrtPopn_std=sqrtPopn, 
#                   logTraffic_std=logTraffic,
#                   roof_AG_std=roof_AG,
#                   roof_COM_std=roof_COM,
#                   roof_IND_std=roof_IND,
#                   roof_RES_std=roof_RES,
#                   roof_TRANS_std=roof_TRANS,
#                   roof_nonRES_std=roof_nonRES,
#                   roof_COM_IND_std=roof_COM_IND,
#                   devAge2_std=devAge2
#                   ))#%>%
#   # dplyr::relocate(dev_pre_1975_std, .before=dev_1975_1990_std) %>%
#   # dplyr::relocate(logCOM_std, .after=AG_std)

sp_std <- data.frame(location=as.factor(sp_t[,c("location")]), loc=as.factor(sp_t[,c("loc")]), sp_std)
  
# sp_final <- cbind(sp_std, select(sp_t, -c("location"))) %>%
#   dplyr::relocate(loc, location)

sp_final <- sp_std

write.csv(sp_final, here("processed_data", "spatial_predictors_standardized.csv"), row.names=FALSE)


#---------#
#  plots  #
#---------#

#untransformed predictors
op <- par(mfrow=c(6,5), mar=c(1,1,1,1), oma=c(0,0,0,0), xaxt="n", cex.main=1.5)
dotchart(sp$urbRES, main="urbRES", group=sp$loc)  #urban residential
dotchart(sp$ruRES, main="ruRES", group=sp$loc)  #rural residential
dotchart(sp$totRES, main="totRES", group=sp$loc)  #total residential (urban + rural res)
dotchart(sp$intURB, main="intURB", group=sp$loc, yaxt="n")  #intensive urban
dotchart(sp$IND, main="IND", group=sp$loc)  #industrial
dotchart(sp$intURB_IND, main="intURB_IND", group=sp$loc)
dotchart(sp$grass_low_veg, main="grass", group=sp$loc)  
dotchart(sp$greenery, main="greenery", group=sp$loc)
dotchart(sp$percent_tree_cover, main="trees", group=sp$loc)  
dotchart(sp$no_dev, main="nodev", group=sp$loc)  
dotchart(sp$impervious, main="impervious", group=sp$loc)
dotchart(sp$imperv_ground, main="paved", group=sp$loc)  
dotchart(sp$imperv_roofs, main="roofs", group=sp$loc)  
dotchart(sp$traffic, main="traffic", group=sp$loc)  
dotchart(sp$pop_per_ha, main="popn", group=sp$loc) 
dotchart(sp$NO_2, main="NO_2", group=sp$loc)
dotchart(sp$PM25_NA, main="PM25_NA", group=sp$loc)  
dotchart(sp$particulate_surface_area, main="partSA", group=sp$loc)
dotchart(sp$slope, main="slope", group=sp$loc)
dotchart(sp$CO_emissions_residential, main="CO2_res", group=sp$loc)
dotchart(sp$CO_emissions_total, main="CO2_tot", group=sp$loc)
dotchart(sp$CO_emissions_commercial, main="CO2_com", group=sp$loc)
dotchart(sp$CO_emissions_nonroad, main="CO2_nonroad", group=sp$loc)
dotchart(sp$CO_emissions_onroad, main="CO2_road", group=sp$loc)
dotchart(sp$devAge, main="dev age", group=sp$loc)  
dotchart(sp$roof_intURB, main="roof_intURB", group=sp$loc)  
dotchart(sp$roof_ruRES, main="roof_ruRES", group=sp$loc)  
dotchart(sp$roof_urbRES, main="roof_urbRES", group=sp$loc) 
dotchart(sp$roof_intURB_IND, main="roof_intURB_IND", group=sp$loc)
dotchart(sp$roof_totRES, main="roof_totRES", group=sp$loc)


#transformed predictors
op <- par(mfrow=c(6,5), mar=c(1,1,1,1), oma=c(0,0,0,0), xaxt="n", cex.main=1.5)
#dotchart(sqrt(sp$totRES), main="sqrt_totRES", group=sp$loc)  #total residential (urban + rural res)
#dotchart(sqrt(sp$intURB), main="sqrt_intURB", group=sp$loc, yaxt="n")  #intensive urban
#dotchart(sqrt(sp$no_dev), main="sqrt_nodev", group=sp$loc)  
dotchart(sqrt(sp$pop_per_ha), main="sqrt_popn", group=sp$loc) 
dotchart(sqrt(sp$slope), main="sqrt_slope", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_residential), main="sqrt_CO2_res", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_total), main="sqrt_CO2_tot", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_commercial), main="sqrt_CO2_com", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_nonroad), main="sqrt_CO2_nonroad", group=sp$loc)
dotchart(sqrt(sp$CO_emissions_onroad), main="sqrt_CO2_road", group=sp$loc)
dotchart((sp$devAge)^2, main="devAge^2", group=sp$loc)  
#dotchart(sqrt(sp$roof_intURB), main="sqrt_roof_intURB", group=sp$loc)  
#dotchart(sqrt(sp$roof_urbRES), main="sqrt_roof_urbRES", group=sp$loc) 
#dotchart(sqrt(sp$roof_intURB_IND), main="sqrt_roof_intURB_IND", group=sp$loc)
#dotchart(sqrt(sp$roof_totRES), main="sqrt_roof_totRES", group=sp$loc)
              
