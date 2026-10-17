# Convert .kml and .kmz files to .geojson files


library(sf)

# Set your working directory where the files are located
setwd("data/model validation data/model validation catchment polygons")

# 1. Find all input files in the working directory
kmz_files <- list.files(pattern = "\\.kmz$", ignore.case = TRUE)
kml_files <- list.files(pattern = "\\.kml$", ignore.case = TRUE)

# 2. Extract each KMZ into its own unique subfolder to prevent overwriting doc.kml
base_temp_dir <- file.path(tempdir(), "extracted_kmz")
all_kml_paths <- kml_files

for (kmz in kmz_files) {
  kmz_base_name <- tools::file_path_sans_ext(kmz)
  unique_subfolder <- file.path(base_temp_dir, kmz_base_name)
  if (!dir.exists(unique_subfolder)) dir.create(unique_subfolder, recursive = TRUE)
  
  unzip(kmz, exdir = unique_subfolder)
  
  extracted_kmls <- list.files(unique_subfolder, pattern = "\\.kml$", full.names = TRUE, ignore.case = TRUE)
  all_kml_paths <- c(all_kml_paths, extracted_kmls)
}

# 3. Initialize an empty list to store spatial data frames
sf_list <- list()

# 4. Loop through files and layers to collect all polygon data
for (kml_path in all_kml_paths) {
  
  # Identify the clean original filename (.kmz name or original .kml name)
  if (grepl("extracted_kmz", kml_path)) {
    base_name <- basename(dirname(kml_path)) 
  } else {
    base_name <- tools::file_path_sans_ext(basename(kml_path)) 
  }
  
  # Check for available layers
  layer_info <- st_layers(kml_path)
  layer_names <- layer_info$name
  total_layers <- length(layer_names)
  
  for (layer in layer_names) {
    sf_data <- st_read(kml_path, layer = layer, quiet = TRUE)
    if (nrow(sf_data) == 0) next
    
    # Track down the existing 'Name' attribute inside the spatial data frame
    col_names <- colnames(sf_data)
    name_col <- col_names[tolower(col_names) == "name"]
    
    # Initialize a clean, standardized data frame containing only Geometry
    sf_clean <- sf_data[, "geometry", drop = FALSE]
    
    # Assign the Polygon Name based on your exact requirements
    if (length(name_col) > 0 && !all(is.na(sf_data[[name_col[1]]]))) {
      # Use the internal Name attribute if it exists
      sf_clean$polygon_name <- as.character(sf_data[[name_col[1]]])
    } else if (total_layers > 1) {
      # Use the layer name if it's a multi-polygon KML missing internal names
      sf_clean$polygon_name <- paste0(base_name, "_", layer)
    } else {
      # Fall back to the original file name
      sf_clean$polygon_name <- base_name
    }
    
    # Append to our merge list
    sf_list[[length(sf_list) + 1]] <- sf_clean
  }
}

# 5. Merge all layers together, drop Z-dimension, and export to one final GeoJSON file
if (length(sf_list) > 0) {
  combined_sf <- do.call(rbind, sf_list)
  
  # Drop the Z (and M if present) dimensions to force strict XY coordinates
  combined_sf <- st_zm(combined_sf, drop = TRUE, what = "ZM")
  
  # Rename the standardized column to exactly "Name" for your GeoJSON
  names(combined_sf)[names(combined_sf) == "polygon_name"] <- "Name"
  
  # Define your single output filename
  output_geojson <- "model_validation_polygons_2026_09_28.geojson"
  
  # Write the unified dataset
  st_write(combined_sf, output_geojson, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  
  cat("\nSuccess! Stripped Z-dimensions and combined all polygons into a single file:", output_geojson, "\n")
} else {
  cat("\nNo valid spatial data found to combine.\n")
}

setwd("../../..") 

# # 5. Merge all layers together and export to one final GeoJSON file
# if (length(sf_list) > 0) {
#   combined_sf <- do.call(rbind, sf_list)
#   
#   # Rename the standardized column to exactly "Name" for your GeoJSON
#   names(combined_sf)[names(combined_sf) == "polygon_name"] <- "Name"
#   
#   # Define your single output filename
#   output_geojson <- "model_validation_polygons_2026_09_25.geojson"
# 
#   # Write the unified dataset
#   st_write(combined_sf, output_geojson, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
#   
#   cat("\nSuccess! Combined all polygons into a single file with names intact:", output_geojson, "\n")
# } else {
#   cat("\nNo valid spatial data found to combine.\n")
# }



