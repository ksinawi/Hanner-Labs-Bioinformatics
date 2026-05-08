# Created By: Faisal Isam
# Last Updated: May 4, 2026

# Check if the 'plotly' package is installed; if not, install it
if(!require(plotly)) install.packages("plotly")

# Load libraries used for data handling, mapping, and visualization
library(tidyverse)   # For data wrangling and manipulation
library(sf)          # For handling spatial data (simple features)
library(plotly)      # For creating interactive maps
library(readr)       # For reading CSV files efficiently
library(htmlwidgets) # For saving interactive visualizations as HTML
library(htmltools)   # For working with HTML content
library(base64enc)   # For working with images


# Directories containing taxa (biological) data and mine data
taxa_dir <- "~/Desktop/Bioinformatics/taxa_dir"
mine_dir <- "~/Desktop/Bioinformatics/mine_dir"

# GeoJSON file containing Canada’s provincial boundaries
geojson_file <- "~/Desktop/Bioinformatics/scripts/python_scripts/ca.json"

logo_png <- "~/Desktop/Bioinformatics/assets/Acknowledgments_Team_Funders.jpg"


# Get a list of all CSV files in the taxa directory
taxa_csv_files <- list.files(path = taxa_dir, pattern = ".csv", full.names = TRUE)
if (length(taxa_csv_files) == 0) stop("Error: Directory <taxa_dir> does not contain any files")

# Get a list of all CSV files in the mine directory
mine_csv_files <- list.files(path = mine_dir, pattern = ".csv", full.names = TRUE)
if (length(mine_csv_files) == 0) stop("Error: Directory <mine_dir> does not contain any files")

# Read and combine all taxa CSV files into one dataframe
# Convert lat/lon to numeric and remove rows with missing coordinates
taxa_data <- map_df(taxa_csv_files, ~ read_csv(.x, col_types = cols(.default = col_character()), show_col_types = FALSE)) %>%
  mutate(lat = as.numeric(lat),
         lon = as.numeric(lon)) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  distinct(family_name, lat, lon, .keep_all = TRUE)  # Remove duplicate locations per family

# Assign each taxa family a unique colour
taxa_families <- unique(taxa_data$family_name)
taxa_colours <- setNames(rainbow(length(taxa_families)), taxa_families)
taxa_data <- taxa_data %>% mutate(colour = taxa_colours[family_name])


# Read and combine all mine CSV files
mine_data <- map_df(mine_csv_files, ~ read_csv(.x, col_types = cols(.default = col_character()), show_col_types = FALSE))

# Create empty data frames to store different mine datasets
cs <- tibble()
mine <- tibble()
pper <- tibble()

# Extract EcoReg DNA Research Sites (if columns exist)
if (all(c("Station_latitude", "Station_longitude") %in% colnames(mine_data))) {
  cs <- mine_data %>%
    select(Station_latitude, Station_longitude, Metal, Mine, everything()) %>%
    rename(lat = Station_latitude, lon = Station_longitude, metal = Metal, site_name = Mine) %>%
    mutate(
      # Clean coordinates to ensure numeric format
      lat = as.numeric(str_replace_all(lat, "[^0-9.-]", "")),
      lon = as.numeric(str_replace_all(lon, "[^0-9.-]", ""))
    ) %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    filter(!str_detect(tolower(metal), "n/a"))%>%
    distinct(site_name, .keep_all = TRUE)
}

# Extract MDMER Mine data (alternative column format)
if (all(c("Latitude", "Longitude") %in% colnames(mine_data))) {
  mine <- mine_data %>%
    select(Latitude, Longitude, Metal, `Organization name/Nom de l'organisation`, everything()) %>%
    rename(lat = Latitude, lon = Longitude, metal = Metal, site_name = `Organization name/Nom de l'organisation`) %>%
    mutate(
      lat = as.numeric(str_replace_all(lat, "[^0-9.-]", "")),
      lon = as.numeric(str_replace_all(lon, "[^0-9.-]", ""))
    ) %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    filter(!str_detect(tolower(metal), "n/a"))
}

# Extract PPER Mill data (coordinates stored as a single text field)
if ("Coordinates" %in% colnames(mine_data)) {
  pper <- mine_data %>%
    select(Coordinates, Product, Company, everything()) %>%
    filter(!str_detect(tolower(Product), "n/a")) %>%
    separate(Coordinates, into = c("lat", "lon"), sep = ",\\s*", convert = TRUE) %>%
    filter(!is.na(lat), !is.na(lon))
}


# Read the Canada GeoJSON file
canada <- st_read(geojson_file, quiet = TRUE)

# Merge duplicate province polygons into single shapes
canada <- canada %>%
  group_by(name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# Project Canada shapefile into Lambert Conformal Conic projection (EPSG:3347)
canada_proj <- st_transform(canada, 3347)

# Generate province label coordinates
province_labels <- suppressWarnings(
  canada_proj %>%
    st_point_on_surface() %>%
    st_transform(4326) %>%
    st_coordinates() %>%
    as_tibble() %>%
    bind_cols(canada %>% st_drop_geometry()) %>%
    rename(lon = X, lat = Y)
)


taxa_sf <- st_as_sf(taxa_data, coords = c("lon", "lat"), crs = 4326)
taxa_sf <- taxa_sf[st_within(taxa_sf, st_union(canada), sparse = FALSE), ]

# Prepare taxa points for plotting
taxa_plot <- taxa_sf %>%
  st_coordinates() %>%
  as_tibble() %>%
  bind_cols(taxa_sf %>% st_drop_geometry())


# Initialize Plotly map
fig <- plot_ly(height = 900, width = 1795) %>% 
  layout(
    title = "EEM Benthic Research Sampling",
    mapbox = list(
      style = "open-street-map",
      center = list(lat = 56, lon = -106),
      zoom = 3
    ),
    showlegend = TRUE,
    title = "EEM Benthic Research Sampling"
  )

fig <- fig %>% layout(
  images = list(
    list(
      source = base64enc::dataURI(file = logo_png, mime = "image/png"),
      xref = "paper",  
      yref = "paper",
      x = 0.00,         
      y = 1.00,         
      sizex = 0.3,      
      sizey = 0.3,      
      xanchor = "left",
      yanchor = "top",
      opacity = 1,
      layer = "above"   
    )
  )
)

# Add provincial boundary lines
for (i in seq_len(nrow(canada))) {
  lines_i <- st_cast(canada[i, ], "MULTILINESTRING")
  coords <- st_coordinates(lines_i) %>% as_tibble()
  
  for (g in unique(coords$L1)) {
    seg <- coords %>% filter(L1 == g)
    fig <- fig %>% add_trace(
      type = "scattermapbox",
      lon = seg$X,
      lat = seg$Y,
      mode = "lines",
      line = list(color = "black", width = 1),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
}

# Add EcoReg DNA Research Sites to map
fig <- fig %>% add_trace(
  type = "scattermapbox", 
  lat = cs$lat,
  lon = cs$lon,
  mode = "markers",
  marker = list(size = 12, color = "#000000", symbol = "circle"),
  text = paste("<br>Site: ", cs$site_name, "<br>Metal:", cs$metal),
  hoverinfo = "text",
  name = "EcoReg DNA Research Sites"
)

# Add MDMER Mines
fig <- fig %>% add_trace(
  type = 'scattermapbox', 
  lat = mine$lat,
  lon = mine$lon,
  mode = "markers",
  marker = list(size = 12, color = "#012345", symbol = "circle"),
  text = paste("<br>Site: ", mine$site_name, "<br>Metal:", mine$metal),
  hoverinfo = "text",
  name = "MDMER Mines"
)

# Add PPER Mills
fig <- fig %>% add_trace(
  type = "scattermapbox", 
  lat = pper$lat,
  lon = pper$lon,
  mode = "markers",
  marker = list(size = 12, color = "#808080", symbol = "circle"),
  text = paste("<br>Site: ", pper$Company, "<br>Product:", pper$Product),
  hoverinfo = "text",
  name = "PPER Mills"
)

# Add taxa (biological sampling) points
fig <- fig %>% add_trace(
  type = "scattermapbox", 
  lat = taxa_plot$Y,
  lon = taxa_plot$X,
  mode = "markers",
  marker = list(size = 8, color = taxa_plot$colour, symbol = "circle"),
  text = taxa_plot$family_name,
  hoverinfo = "text",
  name = taxa_plot$family_name
)

# Display the interactive map
fig

# Save HTML with embedded image + plot
htmlwidgets::saveWidget(
  widget = fig,
  file = "dynamic_map_may_4_2026.html",
  selfcontained = TRUE
)

browseURL("dynamic_map_may_4_2026.html")
