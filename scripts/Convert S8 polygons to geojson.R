# Install and load the required package
if (!requireNamespace("sf", quietly = TRUE)) {
  install.packages("sf")
}
library(sf)

# Define file paths (Update these with your actual paths)
kmz_file <- "data/S8 polygons/S8 polygons 2026_09_16.kmz"
geojson_file <- "data/S8 polygons/S8 polygons 2026_09_16.geojson"

# Create a temporary directory to extract the KMZ contents
temp_dir <- tempdir()
unzip(kmz_file, exdir = temp_dir)

# Find the unzipped KML file(s)
kml_files <- list.files(temp_dir, pattern = "\\.kml$", full.names = TRUE)

if (length(kml_files) == 0) {
  stop("No KML file found inside the KMZ archive.")
}

# Assume we are working with the first KML file found
target_kml <- kml_files[1]

# List layers to ensure we capture all polygon names/attributes properly
layers <- st_layers(target_kml)

# Read the KML layers into an sf object
# (KML attributes like "Name" are automatically preserved by sf)
spatial_data <- st_read(target_kml, layer = layers$name[1])

# Save the spatial data as a GeoJSON file
# This preserves the feature names and attributes in the GeoJSON properties
st_write(spatial_data, geojson_file, driver = "GeoJSON", delete_dsn = TRUE)

cat("Success! Converted KMZ to GeoJSON saved at:", geojson_file, "\n")
