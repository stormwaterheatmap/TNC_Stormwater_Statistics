# Get raw spatial predictors and prepare them by standardizing as needed.
# Note that this version obtains the raw spatial predictors extracted directly from Earth Engine,
# rather than the spatial predictors in the python notebook, that were copied and pasted into a csv file.

# Additionally, v3 also processes landuse percentage predictors

# avg_AADT = average annual daily traffic

# Author: Eva Dusek Jennings
# Date: May 27, 2021
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
 

# #fix location names for TAC001S8D and TAC003S8D so that they match with the COC data   ###THIS WAS DONE ALREADY IN NEW WATERSHED EE DATA!
# sp$location <- as.character(sp$location)
# sp$location[which(sp$location=="TAC001S8D_OF2")] <- "TAC001S8D_OF235"
# sp$location[which(sp$location=="TAC003S8D_OF2")] <- "TAC003S8D_OF245"

#create a new column for proportion of development by age as a continuous variable, where the larger
#  the number, the larger the proportion of old development.  Numbers close to 0 either have very new
#  development, or very little development.
sp$devAge <- 4*sp$dev_pre_1975 + 3*sp$dev_1975_1990 + 2*sp$dev_1990_2000 + 1*sp$dev_2000_2014

#create a new column for landuse as a continuous variable, where the larger the number,
#  the more commercially/ industrially developed the area is.
sp$LU <- 1*sp$AG + 2*sp$RES + 3*sp$COM + 4*sp$TRANS + 5*sp$IND

#look at histograms of the data; determine which columns require log-transformation;
#  look for explanatory variables that have extreme observations; these will likely need
#  a transformation applied to them
ggplot(gather(subset(sp, select=-c(location, loc))), aes(value)) +
  geom_histogram(bins=10, fill = "cadetBlue") +
  facet_wrap(~key, scales = "free")

#use cleveland dotcharts to determine which predictors have extreme values and thus require xforms
sp$loc <- as.factor(sp$loc)
op <- par(mfrow=c(4,4), mar=c(3,3,3,1))
dotchart(sp$AG, main="AG", group=sp$loc)  #agriculture and silviculture
#dotchart(sqrt(sp$AG), main="sqrtAG", group=sp$loc)  #agriculture and silviculture
dotchart(sp$COM, main="COM", group=sp$loc)  #commercial
#dotchart(sqrt(sp$COM), main="sqrtCOM", group=sp$loc)  #commercial
#dotchart(log(sp$COM), main="logCOM", group=sp$loc)  #commercial
dotchart(sp$IND, main="IND", group=sp$loc)  #industrial
dotchart(sp$RES, main="RES", group=sp$loc)  #residential
dotchart(sp$OPEN, main="OPEN", group=sp$loc)  #open space
#dotchart(sqrt(sp$RES), main="sqrtRES", group=sp$loc)  #residential
dotchart(sp$TRANS, main="TRANS", group=sp$loc)  #transportation (proportion of landuse)
dotchart(sp$grass_low_veg, main="grass", group=sp$loc)  
dotchart(sp$imperv_ground, main="paved", group=sp$loc)  
dotchart(sp$imperv_roofs, main="roofs", group=sp$loc)  
dotchart(sp$impervious, main="impervious", group=sp$loc)
dotchart(sp$NO_2, main="NO_2", group=sp$loc)
#dotchart(sqrt(sp$no_dev), main="sqrt_nodev", group=sp$loc)  
#dotchart(log(sp$no_dev*100+1), main="log_nodev", group=sp$loc)  
dotchart(sp$no_dev, main="nodev", group=sp$loc)  
dotchart(sp$percent_tree_cover, main="trees", group=sp$loc)  
dotchart(sp$pm25, main="pm25", group=sp$loc)  
dotchart(sp$slope, main="slope", group=sp$loc)  
dotchart(sp$traffic, main="traffic", group=sp$loc)  
dotchart(sp$pop_per_ha, main="popn", group=sp$loc)  
#looks like the following require strong transformations: COM, trees, traffic, popn
#  these might need weak transformations: AG, IND, RES, nodev, slope

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


par(mfrow=c(3,3))
dotchart(sp$roof_AG, main="roof_AG", group=sp$loc)  #AG and TIMBER combined
dotchart(sp$roof_IND, main="roof_IND", group=sp$loc)  
dotchart(sp$roof_COM, main="roof_COM", group=sp$loc)  
dotchart(sqrt(sp$roof_COM), main="sqrt_roof_COM", group=sp$loc)  
dotchart(log(sp$roof_COM), main="log_roof_COM", group=sp$loc)  
dotchart(sp$roof_RES, main="roof_RES", group=sp$loc)  
dotchart(sp$roof_TRANS, main="roof_TRANS", group=sp$loc) 
dotchart(sp$roof_COM + sp$roof_IND, main="roof_COM_IND", group=sp$loc)  
dotchart(sqrt(sp$roof_COM + sp$roof_IND), main="sqrt_roof_COM_IND", group=sp$loc)  
dotchart(sp$roof_AG + sp$roof_COM + sp$roof_IND + sp$roof_TRANS, main="roof_nonRES", group=sp$loc)
dotchart(sqrt(sp$roof_AG + sp$roof_COM + sp$roof_IND + sp$roof_TRANS), main="sqrt_roof_nonRES", group=sp$loc)
#The various roof types should either ALL be transformed, or none should be transformed,
#  as they might be combined together in various ways.  Choose to NOT transform any of them.
#  Initial thought was to not transform roofs.  However, for Zinc, the relationship really should be
#  with non-RES roofs, and the nonRES roofs are too clustered on left side of the plot.  Change plan.
#  Transform roof_nonRES with sqrt xform


#create a new column for landuse as a continuous variable, where the larger the number,
#  the more commercially/ industrially developed the area is.
#sp$LU <- 1*sp$AG + 2*sp$RES + 3*sp$TRANS + 4*sp$COM + 5*sp$IND
#### This gets into all sorts of machinations for how to apply multipliers, and I can't seem to come
#       up with something that really makes sense.  Best to skip this, as it will be hard to explain
#       exactly how a certain value was chosen.

#look at various combinations of land use types (COM + IND)
par(mfrow=c(3,2))
dotchart(sp$AG, main="AG", group=sp$loc)
dotchart(sp$RES, main="RES", group=sp$loc)
dotchart(sp$COM, main="COM", group=sp$loc)  
dotchart(sp$TRANS, main="TRANS", group=sp$loc)  
dotchart(sp$IND, main="IND", group=sp$loc)  
# dotchart((sp$COM + sp$IND), main="COM+IND", group=sp$loc)  
# dotchart(log((sp$COM + sp$IND)*100+1), main="logCOM+IND", group=sp$loc)
#dotchart(sp$LU, main="land use", group=sp$loc)
#log+1 transform the combination of COM+IND


#try various transformations for COM, trees, traffic, population
#  note: log+1 transformation is chosen when some values in the dataset are 0
par(mfrow=c(4,3), mar=c(4,4,4,1))
dotchart(sp$COM, main="COM", group=sp$loc)  #commercial
dotchart(sqrt(sp$COM), main="sqrtCOM", group=sp$loc)  
dotchart(log(sp$COM*100+1), main="logCOM", group=sp$loc)  #log+1 transform is still not quite right, but slightly better than sqrt
dotchart(sp$percent_tree_cover, main="trees", group=sp$loc)  
dotchart(sqrt(sp$percent_tree_cover), main="sqrt_trees", group=sp$loc)  
dotchart(log(sp$percent_tree_cover+1), main="log_trees", group=sp$loc)  #definitely log+1-transform tree cover!
dotchart(sp$traffic, main="traffic", group=sp$loc)  
dotchart(sqrt(sp$traffic), main="sqrt_traffic", group=sp$loc)  
dotchart(log(sp$traffic), main="log_traffic", group=sp$loc)  #definitely log-transform traffic
dotchart(sp$pop_per_ha, main="popn", group=sp$loc)  
dotchart(sqrt(sp$pop_per_ha), main="sqrt_popn", group=sp$loc)  #sqrt transform is better than log for pop'n
dotchart(log(sp$pop_per_ha), main="log_popn", group=sp$loc) 
#conclusion: log-transform traffic; log+1 transform COM and trees; sqrt transform pop'n

#try sqrt transformations for AG, IND, RES, nodev and slope
par(mfrow=c(3,4), mar=c(4,4,4,1))
dotchart(sp$AG, main="AG", group=sp$loc)  
dotchart(sqrt(sp$AG), main="sqrtAG", group=sp$loc)  #choose sqrt(AG)
dotchart(sp$IND, main="IND", group=sp$loc)  
dotchart(sqrt(sp$IND), main="sqrtIND", group=sp$loc)  #this xform doesn't do much; stick with regular IND
dotchart(sp$RES, main="RES", group=sp$loc)  
dotchart(sqrt(sp$RES), main="sqrtRES", group=sp$loc)  #choose sqrt(RES)
dotchart(sp$no_dev, main="nodev", group=sp$loc)  
dotchart(sqrt(sp$no_dev), main="sqrt_nodev", group=sp$loc)  #choose sqrt(nodev)
dotchart(sp$slope, main="slope", group=sp$loc)  
dotchart(sqrt(sp$slope), main="sqrt_slope", group=sp$loc)  #choose sqrt(slope)
#conclusion: sqrt transform AG, RES, nodev, slope; no transform on IND

#create new columns for transformed data
sp_t <- sp %>%
  dplyr::relocate(location, loc) %>%  #relocate "location" column to first column
  # dplyr::select(-c(loc)) %>%  #remove temporary column "loc"
  dplyr::mutate(logCOM=log(COM*100+1),  #COM: log+1 xform
                logCOM_IND=log((COM + IND)*100+1),  #COM + IND: log+1 xform
                logTrees=log(percent_tree_cover+1),  #tree cover: log+1 xform
                logTraffic=log(traffic),  #traffic: log xform
                sqrtPopn=sqrt(pop_per_ha),  #population: sqrt xform
                sqrtAG=sqrt(AG),  #AG (agriculture/tree farming): sqrt xform
                sqrtRES=sqrt(RES),  #RES (residential): sqrt xform
                sqrtNodev=sqrt(no_dev),  #nodev: sqrt xform
                sqrtSlope=sqrt(slope),  #slope: sqrt xform
                roof_COM_IND_orig=(roof_COM + roof_IND),  #commercial + industrial landuse combined
                roof_nonRES_orig = (roof_COM + roof_IND + roof_AG + roof_TRANS + roof_OPEN),
                roof_nonRES=sqrt(roof_nonRES_orig),  #non-residential roofs
                roof_COM_IND=sqrt(roof_COM_IND_orig),
                devAge2=(devAge)^2) %>%   #development age, as a calculated continuous variable
  dplyr::rename(roofs=imperv_roofs, 
                paved=imperv_ground,
                trees=percent_tree_cover,
                nodev=no_dev,
                grass=grass_low_veg,
                popn=pop_per_ha,
                no2=NO_2) 

sp_means <- sp_sds <- rep(NA, ncol(sp_t))
sp_means[3:ncol(sp_t)] <- colMeans(sp_t[, 3:ncol(sp_t)])
sp_sds[3:ncol(sp_t)] <- apply(sp_t[,3:ncol(sp_t)], 2, sd)

sp_std <- t((t(sp_t[, 3:ncol(sp_t)])-sp_means[3:ncol(sp_t)])/sp_sds[3:ncol(sp_t)])
sp_std <- data.frame(location=as.factor(sp_t[,c("location")]), sp_std) %>%
  dplyr::select(c(location,
#                  loc,
                  sqrtAG_std=sqrtAG, 
                  logCOM_std=logCOM, 
                  IND_std=IND, 
                  OPEN_std=OPEN, 
                  sqrtRES_std=sqrtRES, 
                  TRANS_std=TRANS, 
                  logCOM_IND_std=logCOM_IND,
                  dev_pre_1975_std=dev_pre_1975, 
                  dev_1975_1990_std=dev_1975_1990, 
                  dev_1990_2000_std=dev_1990_2000, 
                  dev_2000_2014_std=dev_2000_2014, 
                  grass_std=grass, 
                  paved_std=paved, 
                  roofs_std=roofs, 
                  impervious_std=impervious,
                  no2_std=no2,
                  sqrtNodev_std=sqrtNodev, 
                  logTrees_std=logTrees,
                  pm25_std=pm25, 
                  sqrtSlope_std=sqrtSlope, 
                  sqrtPopn_std=sqrtPopn, 
                  logTraffic_std=logTraffic,
                  roof_AG_std=roof_AG,
                  roof_COM_std=roof_COM,
                  roof_IND_std=roof_IND,
                  roof_RES_std=roof_RES,
                  roof_TRANS_std=roof_TRANS,
                  roof_nonRES_std=roof_nonRES,
                  roof_COM_IND_std=roof_COM_IND,
                  devAge2_std=devAge2
                  ))#%>%
  # dplyr::relocate(dev_pre_1975_std, .before=dev_1975_1990_std) %>%
  # dplyr::relocate(logCOM_std, .after=AG_std)

sp_final <- cbind(sp_std, select(sp_t, -c("location"))) %>%
  dplyr::relocate(loc, location)

write.csv(sp_final, here("processed_data", "spatial_predictors_standardized.csv"), row.names=FALSE)


