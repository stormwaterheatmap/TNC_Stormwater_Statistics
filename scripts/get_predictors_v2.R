## Header ---------------------------
##
## Script name:get_predictors.R
##
## Abstract: This script is used to get landscape predictors from Google Earth Engine
## for use in stormwaterheatmap regressions.
##
## Author: Christian Nilsen, Geosyntec Consultants
## Email: cnilsen@geosyntec.com
##
## Date Created: 2021-01-24, Christian Nilsen
## Date Modified: 2021-01-31, Eva Dusek Jennings
##
## Copyright (c) Geosyntec Consultants, 2021
##
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, version 3.0
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## For a copy of the GNU General Public License
##  see <https://www.gnu.org/licenses/>.
##


# Installation ------------------------------------------------------------

# # Uncomment to install the first time (only need to do this once)
# #
# # Install the rgee package from GitHub
# install.packages("devtools")
# devtools::install_github("r-spatial/rgee")
# 
# #
# # rgee depends on reticulate because it has some Python dependencies (i.e. numpy and ee), run as follows to install them:
# rgee::ee_install()
# # If you are a Windows user reticulate requires miniconda/anaconda.
# # The use of rgee::ee_install() is not mandatory, you can count on with your own custom installation.
# #
# # After install rgee, you might use the function below for checking the status of rgee.
ee_check() # Check non-R dependencies

rgee::ee_install_upgrade()



# Libraries ---------------------------------------------------------------
library(rgee)
library(mapview)
library(purrr)
library(tidyverse)
library(data.table)


ee$Authenticate()  #the dollar sign on this function tells Python to do this.  ee is a wrapper around python commands

# initialize earth engine api:
ee$Initialize()



# Data  -------------------------------------------------------------------
watersheds <- ee$FeatureCollection("users/cnilsen/s8_watersheds")$select(c("Location_N"), c("Location")) # $filter(ee$Filter$neq('Location_N', 'POSOUTFALL_60')

## Get Predictor images from earth engine
## trees
tree_cover <- ee$Image("USGS/NLCD/NLCD2016")$select("percent_tree_cover")

## traffic
traffic <- ee$Image(0)$blend(ee$Image("users/cnilsen/traffic_raw"))$rename("traffic")

## population density
population <- ee$Image("users/stormwaterheatmap/population_per_ha")

## pm 2.5
pm25 <- (ee$Image("users/cnilsen/pm25clipped")$rename("pm25"))

## imperviousness:
tnc_landcover <- ee$Image("users/jrobertson2000/psLandCover_1m_finPS_roofs")#get land cover:
impervious <- tnc_landcover$eq(6)$Or(tnc_landcover$eq(7))$rename("impervious")# impervious (excl roofs) is coded as 6; roofs are coded as 7; this imperviousness is both combined
imp_ground <- tnc_landcover$eq(6)$rename("imperv_ground")  #roads & parking areas & sidewalks; Christian wouldn't do roads separately-- too correlated with impervious & traffic
imp_roofs <- tnc_landcover$eq(7)$rename("imperv_roofs")
grass_low_veg <- tnc_landcover$eq(1)$rename("grass_low_veg")
#The tnc_landcover image can also be used for the following:
# values: [0, 1, 2, 3, 4,5, 6, 7]
# labels: [ 'No data',
#           'Grass/Low Vegetation',
#           'Shrub/Medium Vegetation',
#           'Trees/Forest',
#           'Bare soil',
#           'Water',
#           'Impervious – Except Roofs',
#           'Impervious - Roofs'
# ]

# age of development
#age_of_development <- ee$Image("JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1")$select("built")$rename("mean_dev_age") #mean value; this one is confusing b/c value of 2=no development
no_dev <- ee$Image("JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1")$select("built")$eq(2)$rename("no_dev") #proportion that has no development
age_2000_2014 <- ee$Image("JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1")$select("built")$eq(3)$rename("dev_2000_2014") #proportion that was built in 2000-2014
age_1990_2000 <- ee$Image("JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1")$select("built")$eq(4)$rename("dev_1990_2000") #proportion that was built in 1990-2000
age_1975_1990 <- ee$Image("JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1")$select("built")$eq(5)$rename("dev_1975_1990") #proportion that was built in 1975-1990
age_pre_1975 <- ee$Image("JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1")$select("built")$eq(6)$rename("dev_pre_1975") #proportion that was built pre-1975

# slope
elevation = ee$Image("USGS/NED")
slope = ee$Terrain$slope(elevation)


## Other available data:
# landuse <- ee$Image("users/stormwaterheatmap/gen_landuse")

# Make Map of One Landscape Predictor ------------------------------------------
map_image <- pm25  #can change this to other variables to see other things (this predictor is used in map)
map_viz <- list(min = 4.5, max = 6.8, palette = list("blue", "green", "yellow", "orange", "red"), opacity = 0.5) #0 to 100 for trees; 0 to 1 for proportions; 4.5 to 6.8 for pm25
Map$centerObject(eeObject = watersheds, zoom = 7)
Map$addLayer(map_image, visParams = map_viz) +
  Map$addLayer(watersheds)

# Reduce Predictors -------------------------------------------------------
## combine predictors in one dataset (one band each) - this is an image
predictors <- ee$Image(0)$blend(
  ee$Image$cat(impervious, imp_ground, imp_roofs, grass_low_veg, tree_cover, traffic, population, pm25, slope,
               no_dev, age_2000_2014, age_1990_2000, age_1975_1990, age_pre_1975)
)

## calculate mean stats from earth engine
ee_stats <- predictors$reduceRegions(
  collection = watersheds,  #this is what we gave it for 16 watersheds; Christian has one for Puget Lowlands; will try in Sandbox to see if it works
  scale = 30,  #30m per pixel
  reducer = ee$Reducer$mean()
)

# evaluate ee object - pulls data from server to client
ee_df <- ee_stats$getInfo()

# wrangle the data
df_predictors <- ee_df$features %>%  #band names diff't than image names; constant is traffic
  map("properties") %>%
  rbindlist()
#built refers to buildup index.  Those are categorical data - maybe shouldn't be averaged?

# Results -----------------------------------------------------------------
View(df_predictors)

# pivot for charting
df_long <- df_predictors %>%
  pivot_longer(cols = -c(Location_N))

# plot value of each predictor for the 16 watersheds
ggplot(df_long) +
  geom_col(aes(x = Location_N, y = value), fill = "cadetBlue", position = "dodge") +
  facet_wrap(~name, scales = "free")


# save csv file of predictors
write.csv(df_predictors, "../data/spatial_predictors_raw.csv", row.names=FALSE)


#---------------------------------------------
# To extract what was in other watersheds?

#whole model domain: we want to predict what is goign on in urbanized areas.  Constrain just to puget lowlands?
#probably 2000 watersheds.  EE might crash?  Or do at a larger scale?  



