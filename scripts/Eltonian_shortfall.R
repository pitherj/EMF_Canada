# Author: Jason Pither
# June 1, 2025
# This script pulls code from various places (primarily latest lab manuscript
# from MSc graduate student Hannah Pilat)
# it focuses on downloading widely-used abiotic geospatial data from public
# repos, and processing those to use with North American extent (north of Mexico)

# increase file download timeout options for downloading large files
options(timeout = max(6000, getOption("timeout")))

# assumed directory structure:

# EMF_Canada/
# ├── data/
# │   ├── raw
# │       ├── EMF_sequence_data_GlobalFungi_genbank_Canada_2025_06_01.csv
# ├── scripts
# ├── figures
# └── _README.md

# increase file download timeout options for downloading large files
options(timeout = max(6000, getOption("timeout")))

# Get session info

sessionInfo()

library(here)
library(tidyverse)
library(geodata) 
library(terra)
library(sf)
library(BIEN)
library(purrr)
library(mapview)

#####
## Admin boundaries
#####

# Download and load GADM admin1 boundaries for each country, SpatVector format
canada_bound <- geodata::gadm(country = "CAN", level = 1, path = here::here("data", "raw"))
usa_bound    <- geodata::gadm(country = "USA", level = 1, path = here::here("data", "raw"))
mexico_bound <- geodata::gadm(country = "MEX", level = 1, path = here::here("data", "raw"))

# Filter USA to include:
# - All states except Hawaii, Puerto Rico, and minor territories
# - Retain Alaska and the Lower 48
usa_filtered <- usa_bound[!usa_bound$NAME_1 %in% c(
  "Hawaii", "Puerto Rico", "Guam", "American Samoa",
  "Northern Mariana Islands", "United States Virgin Islands"
), ]

# crop out anything south of 24 degrees North in the USA,
# and Aleutian islands
usa_cont_extent <- terra::ext(-170, -65, 24, 75)
usa_contiguous <- terra::crop(usa_filtered, usa_cont_extent)

canada_unioned <- sf::st_union(sf::st_as_sf(canada_bound))
usa_unioned <- sf::st_union(sf::st_as_sf(usa_contiguous))
mexico_unioned <- sf::st_union(sf::st_as_sf(mexico_bound))

# merge individual country polygons:
north_america <- c(canada_unioned, usa_unioned, mexico_unioned)

# reproject to North America Albers equal-area conic
# https://spatialreference.org/ref/esri/102008/
# define CRS
#new_crs <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs"
#north_america <- terra::project(north_america, new_crs)

## convert to sf object and simplify vectors
north_america_simplified <- sf::st_simplify(north_america, dTolerance = 1000)

#######
## Canada tree species range maps
#######


