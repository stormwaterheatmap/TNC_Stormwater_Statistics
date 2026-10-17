library(sf)
library(ggplot2)
library(leaflet)
library(viridis)

# 1. Read the GeoJSON file
#geojson_path <- "data/s8_watersheds_v5.geojson"
#geojson_path <- "data/S8 polygons/S8 polygons 2026_09_16.geojson"
geojson_path <- "data/model validation data/model validation catchment polygons/combined_polygons.geojson"

geo_data <- st_read(geojson_path)



# 
# # 2. Define your target column name
# # Replace "location_attribute" with the actual column name in your GeoJSON
# location_col <- "Name"
# 
# # 3. Create the static map
# ggplot(data = geo_data) +
#   # Use aes_string to dynamically select the column by its text name
#   geom_sf(aes_string(fill = location_col), color = "white", size = 0.2) +
#   # Apply a high-quality color palette (discrete scales for categories)
#   scale_fill_viridis_d(option = "viridis", name = "Location") +
#   theme_minimal() +
#   labs(
#     title = "GeoJSON Polygons Color-Coded by Location",
#     caption = "Generated via sf and ggplot2"
#   ) +
#   theme(legend.position = "right")
# 







# Ensure the data matches leaflet's expected map projection (WGS 84)
geo_data <- st_transform(geo_data, crs = 4326)


geo_size <- data.frame(locName=geo_data$Location_N, m2=geo_data$area)
geo_size$acres <- geo_size$m2 * 0.000247105
# geo_size$loc <- c("KIC_LDR", "KIC_HDR", "SEA_IND", "SNO_COM", "SNO_LDR", "SEA_COM", "SEA_HDR", "POT_COM", 
#              "PIE_COM", "SNO_HDR", "KIC_COM", "TAC_COM", "TAC_IND", "TAC_HDR")
# geo_size <- geo_size %>%
#   arrange(loc)


# 2. Extract your location factor variable
# Replace 'Location_N' with your exact column name if different
#location_values <- as.factor(geo_data$Location_N)
location_values <- as.factor(geo_data$Name)
num_unique_colors <- length(levels(location_values))

# 3. Create a highly distinct qualitative color palette
# Using topo.colors or rainbow provides sharp contrast instead of soft gradations
distinct_colors <- topo.colors(num_unique_colors, alpha = 1)

pal <- colorFactor(
  palette = distinct_colors, 
  domain = location_values
)

# 4. Build the map centered on Puget Sound
leaflet(data = geo_data) %>%
  # Use a neutral, light background map to make the distinct colors pop
  addProviderTiles(providers$CartoDB.Positron) %>% 
  
  # Set the focus area to the Puget Sound region
  setView(lng = -122.45, lat = 47.60, zoom = 9) %>% 
  
  # Overlay your GeoJSON polygons
  addPolygons(
#    fillColor = ~pal(Location_N),
    fillColor = ~pal(Name),
    
    weight = 1.5,
    opacity = 1,
    color = "white",       # Clear borders separating the colored shapes
    fillOpacity = 0.75,    # High opacity to emphasize distinct colors
    highlightOptions = highlightOptions(
      weight = 3,
      color = "#333",
      fillOpacity = 0.9,
      bringToFront = TRUE
    ),
#    label = ~as.character(Location_N)
    label = ~as.character(Name)
  ) %>%
  
  # Add a legend map key
  addLegend(
    pal = pal, 
    values = location_values, 
    opacity = 0.8, 
    title = "Locations",
    position = "bottomright"
  )