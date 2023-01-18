# Get and prepare Stormwater Permit Outfall data from WA Department of Ecology

# Note: these data are available from two locations in Dept of Ecology:
# 1.) https://data.wa.gov/Natural-Resources-Environment/Municipal-Stormwater-Permit-Outfall-Data/d958-q2ci.
# 2.) These can be obtained directly using the fxn read.socrata("https://data.wa.gov/resource/rc6b-fvgb.json")

# The first location is from EIM - Ecology's information management.  The second is the raw data.
# Only the EIM data have a working timestamp associated with field sampling, which is vital for linking
# sample data to precipitation data from BigQuery.  

# A .csv file is saved in ```../data/S8_data_cleaned.csv```, which includes precip & antecedant dry days

# also note: this version differs from v2 in that the daymet data is added to the joined data rather than the eim data,
#   to ensure that all raw data have a precip value associated with them; it differs from v3 in that antecedant dry day
#   data are calculated from daymet data, to generate complete ADD data.
# 2nd note: daymet data were updated on Dec 15, 2020 to V4 (see daymet website for documentation)

# Author: Eva Dusek Jennings  (Based on scripts from Christian Nilsen)
# Date: Feb 25, 2021
#---------------------------------------------

library(tidyverse)
library(magrittr)
library(lubridate)
library(RSocrata)
library(readr)
library(fuzzyjoin)  #not sure I ended up using this one...
library(daymetr) #for DayMet data


#------------------------------------------------------------------#
#  Download EIM Data; Use for Sampling Start/End Times, Lat, Long  #
#------------------------------------------------------------------#

# load s8 data downloaded from eim.  This source includes sampling start and end times
eim_data <- readRDS("../data/S8_data_eim_raw.rds")
colnames(eim_data) %<>% tolower  # all lowercase col names
# format dates
eim_data$field_collection_start_date <- as.Date.character(eim_data$field_collection_start_date,format="%m/%d/%Y")
eim_data$field_collection_end_date <- as.Date(eim_data$field_collection_end_date,format="%m/%d/%Y")

# eim_data <- eim_data %>%
#   dplyr::mutate(location_date=paste(location_id, "_", field_collection_start_date, "_", field_collection_end_date, sep=""))

#locations of outfalls of interest (note: TFWFD6 is not included here b/c it is a sediment trap only - no stormwater!)
locs <- c("KICCOMS8D_OUT", "KICHDRS8D_OUT", "KICLDRS8D_OUT",
          "PIECOMM_OUT", "PIEHIRES_OUT", "PIELORES_OUT",
          "SEAC1S8D_OUT", "SEAI1S8D_OUT", "SEAR1S8D_OUT",
          "SNO_COM", "SNO_HDR", "SNO_LDR",
          "TAC001S8D_OF235", "TAC003S8D_OF245", "TFWFD1",
          "POT564S8D_OUT")

eim_params <- c("Copper", "Lead", "Total Suspended Solids", "Zinc", "Total Phosphorus", "Nitrite-Nitrate", "Total Kjeldahl Nitrogen as N")

# filter eim_data to obtain start and end time of sampling, along with latitude & longitude
eim_time_location <- eim_data %>%
  dplyr::mutate(location_date=paste(location_id, "_", field_collection_start_date, "_", field_collection_end_date, sep="")) %>%
  dplyr::filter(result_parameter_name!="Storm Event Flow Volume") %>%
  dplyr::filter(result_parameter_name!="Sample Event Flow Volume") %>%
  dplyr::filter(result_parameter_name!="Precipitation") %>%
  #  dplyr::filter(result_parameter_name %in% eim_params) %>%
  dplyr::mutate(start_date_time=lubridate::mdy_hm(field_collection_start_date_time)) %>%
  dplyr::mutate(end_date_time=lubridate::mdy_hm(field_collection_end_date_time)) %>%
  dplyr::select(location_date, location_id,
                field_collection_start_date, field_collection_start_time, start_date_time,
                field_collection_end_date, field_collection_end_time, end_date_time,
                #                result_parameter_name,
                #                result_value,
                calculated_latitude_decimal_degrees_nad83harn,
                calculated_longitude_decimal_degrees_nad83harn) %>%
  dplyr::filter(location_id %in% locs) %>%
  dplyr::filter(field_collection_end_time - field_collection_start_time != 0) %>%
  dplyr::mutate(time_elapsed=end_date_time - start_date_time) %>%
  dplyr::mutate(rain_gage=case_when(
    location_id %in% c("PIECOMM_OUT", "PIEHIRES_OUT", "PIELORES_OUT",
                       "SEAC1S8D_OUT", "SEAI1S8D_OUT", "SEAR1S8D_OUT",
                       "SNO_COM", "SNO_HDR", "SNO_LDR", "POT564S8D_OUT") ~ location_id,
    location_id %in% c("TAC001S8D_OF235", "TAC003S8D_OF245", "TFWFD1") ~ as.factor("TAC004S8D_RAIN"),
    location_id %in% c("KICCOMS8D_OUT", "KICHDRS8D_OUT") ~ as.factor("KIC19US8D_RAIN"),
    location_id %in% "KICLDRS8D_OUT" ~ as.factor("KICRRRS8D_RAIN"))) %>%
  dplyr::mutate(rain_gage_start_date=paste(rain_gage, "_", field_collection_start_date, sep="")) %>%
  dplyr::mutate(rain_gage_end_date=paste(rain_gage, "_", field_collection_end_date, sep="")) %>%
  unique()


#which(eim_time_location$field_collection_end_time - eim_time_location$field_collection_start_time == 0)
#abc <- eim_data[which(eim_data$location_id=="TAC003S8D_OF245" & eim_data$field_collection_start_date=="2011-05-02"),]
# TAC003S8D_OF245 has two sampling times: 8:30 and 10:00.  The 8:30 corresponds to when most samples were collected, and should be kept.
# remove the row with the 10:00 sample
abc <- which(eim_time_location$location_date=="TAC003S8D_OF245_2011-05-02_2011-05-02" & as.character(eim_time_location$field_collection_start_time)=="10:00:00")
eim_time_location <- eim_time_location[-abc,]

eim_dups <- which(duplicated(eim_time_location$location_date))
eim_time_location[c(eim_dups-1, eim_dups),]


#delete the first 6 TFWFD1 samples
#the sampling date 2012-01-09 has two start times: 9:45 and 20:30.  The COC results in raw_data matches those associated with the 20:30 sampling time.
#the other TFWFD1 samples from 2011 were collected 1-3 minutes apart from their duplicates; keep the second set.  
eim_time_location <- eim_time_location[-(eim_dups-1),]  #Note: since we're only using this for start/end time of precip, doesn't matter which we use

# #based on the analysis below, the TFWFD1_2012-01-09_2012-01-10 samples collected starting at 20:30 match those in raw_data
# aa <- s8data %>%
#   dplyr::filter(Location=="TFWFD1") %>%
#   dplyr::filter(field_collection_start_date=="1/9/2012") %>%
# bb <- eim_time_location %>%
#   dplyr::filter(location_date=="TFWFD1_2012-01-09_2012-01-10") %>%
#   dplyr::select(field_collection_start_time, result_parameter_name, result_value)
# cc <- eim_data %>%
#   dplyr::filter(location_date=="TFWFD1_2012-01-09_2012-01-10")
# dd <- all.S8.data %>%
#   dplyr::mutate(location_date=paste(location_id, "_", field_collection_start_date, "_", field_collection_end_date, sep="")) %>%
#   dplyr::filter(location_date=="TFWFD1_1/9/2012_1/10/2012")


## ALSO NOTE: I suspect that all CPAH, Total Phthalate, Total PAH and HPAH samples have no field collection end date!
#length(which(is.na(s8data$field_collection_end_date)))

#how many sampling dates for each location?
n_loc <- rep(0, length(locs))
for(i in 1:length(locs)) {
  n_loc[i] <- length(which(eim_time_location$location_id==locs[i]))
}
names(n_loc) <- locs

# eim_time_location[which(eim_time_location$location_id=="TFWFD6"),]  #sediment trap only!
# 

# save information about flow and preciption
eim_storm_event_flow_vol <- eim_data %>%
  dplyr::filter(result_parameter_name=="Storm Event Flow Volume") %>%
  dplyr::filter(location_id %in% locs)
eim_sample_event_flow_vol <- eim_data %>%
  dplyr::filter(result_parameter_name=="Sample Event Flow Volume") %>%
  dplyr::filter(location_id %in% locs)
eim_precipitation <- eim_data %>%
  dplyr::filter(result_parameter_name=="Precipitation") %>%
  dplyr::filter(location_id %in% c(locs, "KIC19US8D_RAIN", "KICRRRS8D_RAIN", "TAC004S8D_RAIN")) %>%
  dplyr::mutate(start_date_time=lubridate::mdy_hm(field_collection_start_date_time)) %>%
  dplyr::mutate(end_date_time=lubridate::mdy_hm(field_collection_end_date_time)) %>%
  dplyr::select(location_id, field_collection_start_date, field_collection_end_date,
                field_collection_start_time, field_collection_end_time, start_date_time, end_date_time,
                result_value, result_value_units) %>%
  dplyr::mutate(location_date=paste(location_id, "_", field_collection_start_date, "_", field_collection_end_date, sep="")) %>%
  dplyr::mutate(time_elapsed=end_date_time - start_date_time) %>%
  dplyr::mutate(inches_rain_per_hour=result_value/as.numeric(time_elapsed)) %>%
  dplyr::mutate(rain_gage_start_date=paste(location_id, "_", field_collection_start_date, sep="")) %>%
  dplyr::mutate(rain_gage_end_date=paste(location_id, "_", field_collection_end_date, sep=""))


###generate start_date_time and end_date_time.  Subtract end from start, then divide total rain by hour.  Apply to each site this way!
### for sites where the rainfall collection times match the sampling times, it'll line up perfectly.
### for sites where it doesn't match exactly, it won't line up perfectly, but will provide good data.
### alternately, could get rain off of daily precip gages, then match up amnt of rain per hour with length of sampling time
### either option is defensible.

# ggplot(eim_precipitation, aes(field_collection_start_date, result_value)) +
#   geom_point(aes(color=location_id))
# 

#------------------------------------------------------------------#
#  Join Precip Per Hour Over Sampling Period w/ eim_time_location  #
#------------------------------------------------------------------#

eim_time_loc_precip <- 
  left_join(eim_time_location, 
            eim_precipitation %>% 
              dplyr::select(-c(location_date, location_id, field_collection_start_date, field_collection_end_date,
                        field_collection_start_time, field_collection_end_time)),
            by="rain_gage_start_date")                    

#### Try to figure out why there's duplicate observations for some locations/times
length(which(duplicated(eim_time_loc_precip$location_date)))  #16 dups here
precip_dups <- which(duplicated(eim_time_loc_precip$location_date))
aa <- eim_time_loc_precip[c(precip_dups-1, precip_dups),] %>%  #data frame of the duplicate samples
  dplyr::select(location_date, start_date_time.x, end_date_time.x, time_elapsed.x, start_date_time.y, end_date_time.y,
         result_value, time_elapsed.y, inches_rain_per_hour) %>%
  dplyr::mutate(index=c(precip_dups-1, precip_dups)) %>%  #index keeps track of the row's index within eim_time_loc_precip
  dplyr::mutate(start_dif = as.numeric(start_date_time.x - start_date_time.y, units="hours")) %>%  #COC sample start time - precip start time
  dplyr::mutate(end_dif = as.numeric(end_date_time.x - end_date_time.y, units="hours")) %>%  #COC sample end time - precip end time
  dplyr::mutate(sum_difs = abs(start_dif) + abs(end_dif))  #sum of the differences betweens COC and precip start and end times

discard_i <- rep(NA, length(precip_dups))
average_i <- rep(NA, 2*length(precip_dups))
aa <- aa[order(aa$location_date, aa$sum_difs), ]  #rearrange aa so that duplicates are next to each other, and smaller sum_difs is first
for (i in seq(1, nrow(aa), by=2)) { #go through dataframe aa, by every other row (since duplicate rows are next to each other)
  #compare sum_difs for i and i-1
  if (aa$sum_difs[i]==aa$sum_difs[i+1] & aa$inches_rain_per_hour[i]==aa$inches_rain_per_hour[i+1]) {  #if the two rows are identical, discard second occurrence
    discard_i[(i+1)/2] <- aa$index[i+1]
  } else if (aa$sum_difs[i]/aa$sum_difs[i+1] < 0.75) {  #if there is a big difference between the sum_difs, discard the row with the larger sum_dif
    discard_i[(i+1)/2] <- aa$index[i+1]    
  } else {  #if the sum_difs between duplicates are close to each other, rainfall data needs to be averaged; save index in average_i
    average_i[i] <- aa$index[i]
    average_i[i+1] <- aa$index[i+1]
  }
}

#remove NA values from index vectors average_i and discard_i
average_i <- average_i[-which(is.na(average_i))]
discard_i <- discard_i[-which(is.na(discard_i))]


#average the precip for location_date duplicates that have similar sum_difs
for (i in seq(1, length(average_i), by=2)) {
  eim_time_loc_precip$start_date_time.y[average_i[i]] <- NA
  eim_time_loc_precip$end_date_time.y[average_i[i]] <- NA
  eim_time_loc_precip$result_value[average_i[i]] <- mean(eim_time_loc_precip$result_value[average_i[i]:average_i[i+1]])
  eim_time_loc_precip$time_elapsed.y[average_i[i]] <- NA
  eim_time_loc_precip$inches_rain_per_hour[average_i[i]] <- mean(eim_time_loc_precip$inches_rain_per_hour[average_i[i]:average_i[i+1]])
  eim_time_loc_precip$rain_gage_end_date.y[average_i[i]] <- NA
  discard_i <- c(discard_i, average_i[i+1])
}

#discard entries in eim_time_loc_precip that are duplicates identified for discarding
eim_time_loc_precip <- eim_time_loc_precip[-discard_i, ]
which(duplicated(eim_time_loc_precip)==TRUE)  #verify that all duplicates have been removed

#add year and Julian day
eim_time_loc_precip <- eim_time_loc_precip %>%
  add_column(
    yday = yday(.$field_collection_start_date),
    year = year(.$field_collection_start_date)
  )

# Clean up the data - keep only columns we need
eim_time_loc_precip <- eim_time_loc_precip %>%
  dplyr::select(location_id, eim_col_start_date=field_collection_start_date, field_collection_start_time, field_collection_end_date, field_collection_end_time,
         latitude_decDeg=calculated_latitude_decimal_degrees_nad83harn, longitude_decDeg=calculated_longitude_decimal_degrees_nad83harn,
         start_date_time=start_date_time.x, end_date_time=end_date_time.x,
         precip=result_value, precip_units=result_value_units, inches_rain_per_hour,
         yday, year) %>%
  add_column(eim_data=TRUE)


#-----------------------------------------------------#
#  Download Raw Data; Use for COC Values & Access_ID  #
#-----------------------------------------------------#

# use the socrata REST api to download raw data from Ecology (not EIM).  This source includes access_id
raw_data <- read.socrata("https://data.wa.gov/resource/rc6b-fvgb.json") 

params <- c(  # list of parameters in raw_data
  "Zinc - Water - Total",
  "Copper - Water - Total",
  "Nitrite-Nitrate - Water - Dissolved",
  "Lead - Water - Total",
  "Total Phosphorus - Water - Total",
  "Total Suspended Solids - Water - Total",
  "Total Phthalate - Water - Total",
  "Total PAH - Water - Total",
  "CPAH - Water - Total",
  "HPAH - Water - Total",
  "Total Kjeldahl Nitrogen - Water - Total",
  "Zinc - Water - Dissolved"
)

raw_precipitation <- raw_data %>%  ##has the same number of rows as eim_precipitation
  dplyr::filter(parameter=="Precipitation - Water - Total") %>%
  dplyr::filter(location_id %in% locs)

raw_storm_event_flow <- raw_data %>%  ###has the same number of row as eim_storm_event_flow_vol
  dplyr::filter(parameter=="Storm Event Flow Volume - Water - Total") %>%
  dplyr::filter(location_id %in% locs)

# clean up raw_data to remove extraneous locations and parameters; add year and Julian day columns
raw_data <- raw_data %>%
  dplyr::filter(location_id %in% locs) %>%  #keep only locations we want
  dplyr::filter(parameter %in% params) %>%  #keep only parameters we want
  dplyr::select(-c(field_collection_start_time, field_collection_end_time, field_collection_end_date)) %>%  #remove uninformative columns for start/end time and end date
  #  dplyr::mutate(location_date=paste(location_id, "_", field_collection_start_date, "_", field_collection_end_date, sep="")) #new column for joining with eim data
  add_column(
    yday = yday(.$field_collection_start_date),
    year = year(.$field_collection_start_date),
    raw_data=TRUE
  )


#-----------------------------------------------------------#
#  Obtain DayMet Data for Daily Precip & Antecedant Precip  #
#-----------------------------------------------------------#

# save eim data's site list and xy data to file
siteList <- unique(
  tibble(
    site = eim_time_loc_precip$location_id,
    lat = eim_time_loc_precip$latitude_decDeg,
    lon = eim_time_loc_precip$longitude_decDeg
  )
)
write.csv(siteList, "../data/siteList.csv", row.names = FALSE)

# Download the nearest rainfall gage from daymet
daymet_p <- download_daymet_batch(file_location = "../data/siteList.csv", start = 2009,
                                  end = 2013, simplify = T)  #note that first collection day is 2009-02-16, and last day is 2013-04-18

# rename columns
colnames(daymet_p)[colnames(daymet_p) == "site"] <- "location_id"
colnames(daymet_p)[colnames(daymet_p) == "value"] <- "daymet_precip"
colnames(daymet_p)[colnames(daymet_p) == "measurement"] <- "daymet_units"

# get just the precip data
p <- daymet_p %>%
  filter(daymet_units == "prcp..mm.day.") %>%
  select(-c(tile,altitude))

#calculate cummulative two day, three day, and seven day antecedant precip  
p$daymet_2day <- zoo::rollsum(p$daymet_precip, 2, fill=0, align="right")  #these are antecedant 2, 3, 5 and 7 day precip (mm)
p$daymet_3day <- zoo::rollsum(p$daymet_precip, 3, fill=0, align="right")
p$daymet_5day <- zoo::rollsum(p$daymet_precip, 5, fill=0, align="right")
p$daymet_7day <- zoo::rollsum(p$daymet_precip, 7, fill=0, align="right")
p$daymet_10day <- zoo::rollsum(p$daymet_precip, 10, fill=0, align="right")
p$daymet_14day <- zoo::rollsum(p$daymet_precip, 14, fill=0, align="right")
p$daymet_21day <- zoo::rollsum(p$daymet_precip, 21, fill=0, align="right")
p$daymet_28day <- zoo::rollsum(p$daymet_precip, 28, fill=0, align="right")

#calculate Antecedant Dry Days
p$ADD <- 0
#Work backward through all rows of p. Since first date was 2/25, no need to separate out by location, as there was definitely rain between Jan 1, 2009 and Feb 25, 2009
for (i in nrow(p):2) {
  j <- 1
  while(p$daymet_precip[i-j]==0) {
    j <- j+1
  }
  p$ADD[i] <- j-1
}

#obtain month for each year/ yday combination
p$origin <- as.Date(paste0(p$year, "-01-01"), tz="UTC") - days(1)
p$date <- as.Date(p$yday, origin=p$origin, tz="UTC")
p$month <- month(p$date)

#create tibble with monthly cumulative precip for each location & year
p2 <- p %>%
  group_by(location_id, year, month) %>%
  summarize(mPrecip=sum(daymet_precip))


# #-----------------------#
# #  Get Antecedant Days  #  #this dataset is missing data for all PIE and KIC sites
# #-----------------------#
# 
# Antecedant_days <- read_csv("../data/Antecedant_days.csv") %>% 
#   dplyr::filter(location_id %in% locs) %>%  #keep only locations we're interested in
#   dplyr::select(c(access_id, antecedant_dry_days=ADD)) %>%
#   dplyr::distinct()
# 
# add.antecedant.days <- function(df) {
#   df_add <- dplyr::left_join(df, Antecedant_days, by="access_id")
#   df_add$antecedant_dry_days <- as.numeric(df_add$antecedant_dry_days)
#   return(df_add)
# }


#----------------------#
#  Cleaning Functions  #
#----------------------#

clean_data <- function(df) {
  filtered.df <- (filter(df, !result_data_qualifier %in% "REJ")) %>%  #filter out rejected data
    filter(!sample_replicate_flag %in% "Y") %>%    #filter out replicates
    dplyr::mutate(month = month(start_date)) %>%   #add month column
    dplyr::mutate(agency=case_when(  #add column for "agency", with short names for the long study_name
      study_name %in% "City of Seattle Phase I Municipal Stormwater Permit" ~ "Seattle",
      study_name %in% "City of Tacoma Phase I Municipal Stormwater Permit" ~ "Tacoma",
      study_name %in% "King County Phase I Municipal Stormwater Permit" ~ "King Co",
      study_name %in% "Snohomish County Phase I Municipal Stormwater Permit" ~ "Snohomish Co",
      study_name %in% "Pierce County Phase I Municipal Stormwater Permit" ~ "Pierce Co",
      study_name %in% "Port of Tacoma Phase I Municipal Stormwater Permit" ~ "Port of Tacoma"
    )) %>%
    dplyr::select(-c(study_name))  #remove cumbersome study_name column
    
  # change nondetect warnings to detects
  warnings <- filtered.df$nondetect_flag == "WARNING"
  filtered.df$nondetect_flag[warnings] <- FALSE
  
  # Change NA to detect
  filtered.df$nondetect_flag[is.na(filtered.df$nondetect_flag)] <- FALSE
  
  # Change season and month to factors
  filtered.df$season <- as.factor(filtered.df$season)
  filtered.df$month <- as.factor(filtered.df$month)
  
  # Change access id to numeric
  filtered.df$access_id <- as.numeric(filtered.df$access_id)
  return(filtered.df)
}

select_columns <- function(df) {
  # select columns used in analysis
  df <- df %>%
    dplyr::select(
      location_id,
      parameter, result, units, nondetect_flag,
      start_date_time, start_date, start_time,
      end_date_time, end_date, end_time,
      precip, precip_units, inches_rain_per_hour,
      yday, month, year, season,
      latitude, longitude, land_use=type,
      daymet_units, daymet_precip, daymet_2day, daymet_3day, daymet_5day, daymet_7day,
      daymet_10day, daymet_14day, daymet_21day, daymet_28day, mPrecip,
      antecedant_dry_days=ADD,
      access_id, agency
    )

  df$nondetect_flag <- as.logical(df$nondetect_flag)
  df$result <- as.numeric(df$result)
  return(df)
}

#-------------------------------------#
#  Join All Data Together & Clean Up  # 
#-------------------------------------#

# join raw data (for all COCs) with eim data (no COC's, only time, lat/long, collected precip data)
joined_data <- plyr::join_all(list(raw_data, eim_time_loc_precip), by = c("location_id", "year", "yday"),
                              type = "left", match="all") %>%
  dplyr::select(location_id, parameter, result=new_result_value, units=new_result_units, start_date_time, end_date_time, 
                start_date=field_collection_start_date, start_time=field_collection_start_time, 
                end_date=field_collection_end_date, end_time=field_collection_end_time,
                precip, precip_units, inches_rain_per_hour, eim_data, raw_data,
                yday, year, latitude_decDeg, longitude_decDeg, 
                type, nondetect_flag, access_id, sample_replicate_flag, result_data_qualifier, study_name, season
  )

# Join raw/eim data with daymet precip (including cummulative precip over 2,3,5 and 7 antecedant days, AND antecedant dry days (ADD))
joined_data <- plyr::join_all(list(joined_data, p), by = c("location_id", "year", "yday"),
                              type = "left",match="all")

#add monthly precip to the joined data
joined_data <- plyr::join_all(list(joined_data, p2), by=c("location_id", "year", "month"))

#run final cleaning functions, add antecedant dry days
s8data <- joined_data %>%
  clean_data() %>%     #4848 observations to here, 35 variables
#  add.antecedant.days() %>%   #the csv version is NOT complete!  use calculated version from daymet data
  select_columns()  #4848 observations to here, 30 variables 

write.csv(s8data, "../data/s8data.csv", row.names = FALSE)

#NOTE: there are many dates for which sampled precip data are not available.  I wanted to match up daymet precip for the
#  sampling duration, but sample start time/ end time and end date are not available for these instances.
