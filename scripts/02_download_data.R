# =============================================================================
# EMF Canada Analysis - Data Download Script
# =============================================================================
# Purpose: Download and prepare all external datasets
# Author: Jason Pither, Clara Qin, Isaac Eckert
# =============================================================================

# Required packages
required_packages <- c("here", "geodata", "terra", "rnaturalearth", "sf", "tidyverse")

# Load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
}

# Load coordinate systems
if (!exists("crdref")) {
  source(here::here("scripts", "01_setup.R"))
}

# Load utility functions for temp directory management
source(here::here("scripts", "utils.R"))

#' Download administrative boundaries for North America
#' 
#' @param rawdata_dir Character. Directory to save data
#' @param outdata_dir Character. Directory for output data
#' @return List of file paths to saved boundary files
download_admin_boundaries <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Checking administrative boundaries...")
  
  # Create temporary directory for intermediate files
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Define expected output files in rawdata admin_boundaries subdirectory
  admin_dir <- file.path(rawdata_dir, "admin_boundaries")
  canada_path <- file.path(admin_dir, "canada_bound_raw.gpkg")
  usa_path <- file.path(admin_dir, "usa_bound_raw.gpkg")
  mexico_path <- file.path(admin_dir, "mexico_bound_raw.gpkg")
  
  # Check if files already exist
  if (file.exists(canada_path) && file.exists(usa_path) && file.exists(mexico_path)) {
    message("✓ Administrative boundaries already downloaded")
    return(list(
      canada = canada_path,
      usa = usa_path,
      mexico = mexico_path
    ))
  }
  
  message("Downloading administrative boundaries...")
  
  tryCatch({
    # Create admin directory
    if (!dir.exists(admin_dir)) {
      dir.create(admin_dir, recursive = TRUE)
    }
    
    # Download GADM boundaries (geodata will create temp files in temp_dir)
    canada_bound <- geodata::gadm(country = "CAN", level = 1, path = temp_dir)
    usa_bound <- geodata::gadm(country = "USA", level = 1, path = temp_dir)
    mexico_bound <- geodata::gadm(country = "MEX", level = 1, path = temp_dir)
    
    # Save immediately to final locations to avoid pointer issues
    terra::writeVector(canada_bound, canada_path, overwrite = TRUE)
    terra::writeVector(usa_bound, usa_path, overwrite = TRUE)
    terra::writeVector(mexico_bound, mexico_path, overwrite = TRUE)
    
    message("✓ Administrative boundaries downloaded and saved")
    
    return(list(
      canada = canada_path,
      usa = usa_path,
      mexico = mexico_path
    ))
    
  }, error = function(e) {
    stop(paste("Error downloading administrative boundaries:", e$message))
  })
}

#' Download and extract ecoregions data
#' 
#' @param rawdata_dir Character. Directory to save data
#' @param outdata_dir Character. Directory for output data
#' @return Character. Path to ecoregions shapefile
download_ecoregions <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Checking ecoregions data...")
  
  # Create temporary directory for downloads
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Check if ecoregions already exist
  ecoregions_dir <- file.path(rawdata_dir, "ecoregions")
  existing_shp <- list.files(ecoregions_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
  
  if (length(existing_shp) > 0 && file.exists(existing_shp[1])) {
    message("✓ Ecoregions already downloaded")
    return(existing_shp[1])
  }
  
  message("Downloading ecoregions data...")
  
  tryCatch({
    # Create ecoregions directory
    if (!dir.exists(ecoregions_dir)) {
      dir.create(ecoregions_dir, recursive = TRUE)
    }
    
    # Download ecoregions to temporary directory first
    ecoregions_url <- 'https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip'
    ecoregions_zip <- file.path(temp_dir, "ecoregions.zip")
    
    download.file(ecoregions_url, destfile = ecoregions_zip)
    unzip(ecoregions_zip, exdir = ecoregions_dir, overwrite = TRUE)
    
    # Discover the actual structure after unzipping
    # Look for the shapefile in subdirectories
    shp_files <- list.files(ecoregions_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
    
    if (length(shp_files) == 0) {
      stop("No shapefile found after extracting ecoregions")
    }
    
    # Use the first (and likely only) shapefile found
    ecoregions_path <- shp_files[1]
    message(paste("Found ecoregions shapefile at:", ecoregions_path))
    
    # Download ecozone names to the same directory as the shapefile
    ecozone_names_url <- 'https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/zn_names.dbf'
    ecozone_names_path <- file.path(dirname(ecoregions_path), "ecozone_names.dbf")
    
    if (!file.exists(ecozone_names_path)) {
      download.file(ecozone_names_url, destfile = ecozone_names_path)
    }
    
    message("✓ Ecoregions data downloaded")
    
    return(ecoregions_path)
    
  }, error = function(e) {
    stop(paste("Error downloading ecoregions:", e$message))
  })
}

#' Download USTreeAtlas from GitHub
#' 
#' @param rawdata_dir Character. Directory to save data
#' @param outdata_dir Character. Directory for output data
#' @return Character. Path to USTreeAtlas shapefiles directory
download_ustreeatlas <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Checking USTreeAtlas...")
  
  # Create temporary directory for downloads
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Check if USTreeAtlas already exists
  shp_dir <- file.path(rawdata_dir, "USTreeAtlas", "USTreeAtlas-main", "shp")
  
  if (dir.exists(shp_dir) && length(list.dirs(shp_dir, recursive = FALSE)) > 0) {
    message("✓ USTreeAtlas already downloaded")
    return(shp_dir)
  }
  
  message("Downloading USTreeAtlas from GitHub...")
  
  tryCatch({
    # Create USTreeAtlas directory
    ustreeatlas_dir <- file.path(rawdata_dir, "USTreeAtlas")
    if (!dir.exists(ustreeatlas_dir)) {
      dir.create(ustreeatlas_dir, recursive = TRUE)
    }
    
    # Download USTreeAtlas zip to temporary directory first
    ustreeatlas_url <- "https://github.com/wpetry/USTreeAtlas/archive/refs/heads/main.zip"
    ustreeatlas_zip <- file.path(temp_dir, "USTreeAtlas-main.zip")
    
    download.file(ustreeatlas_url, destfile = ustreeatlas_zip, mode = "wb")
    unzip(ustreeatlas_zip, exdir = ustreeatlas_dir, overwrite = TRUE)
    
    # Check that the shapefiles directory exists
    if (!dir.exists(shp_dir)) {
      stop("USTreeAtlas shapefiles directory not found after download")
    }
    
    message("✓ USTreeAtlas downloaded and extracted")
    
    return(shp_dir)
    
  }, error = function(e) {
    stop(paste("Error downloading USTreeAtlas:", e$message))
  })
}

#' Download Natural Earth data for mapping
#' 
#' @param rawdata_dir Character. Directory to save data
#' @param outdata_dir Character. Directory for output data
#' @return List of file paths to Natural Earth data
download_natural_earth <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Checking Natural Earth data...")
  
  # Create temporary directory
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Define expected output files
  ne_dir <- file.path(rawdata_dir, "natural_earth")
  canada_path <- file.path(ne_dir, "canada_ne.gpkg")
  lakes_path <- file.path(ne_dir, "lakes_ne.gpkg")
  
  # Check if files already exist
  if (file.exists(canada_path) && file.exists(lakes_path)) {
    message("✓ Natural Earth data already downloaded")
    return(list(
      canada = canada_path,
      lakes = lakes_path
    ))
  }
  
  message("Downloading Natural Earth data...")
  
  tryCatch({
    # Create natural_earth directory
    if (!dir.exists(ne_dir)) {
      dir.create(ne_dir, recursive = TRUE)
    }
    
    # Download Natural Earth data - using only medium scale to avoid rnaturalearthhires dependency
    canada_ne <- rnaturalearth::ne_countries(scale = "medium", country = "Canada", returnclass = "sf") 
    canada_ne <-  terra::vect(canada_ne)
    
    # For lakes, use medium scale and avoid any high-res dependencies
    lakes_ne <- rnaturalearth::ne_download(scale = "medium", type = "lakes", 
                                           category = "physical", returnclass = "sf") 
    lakes_ne <-  terra::vect(lakes_ne)
    
    # Save immediately
    terra::writeVector(canada_ne, canada_path, overwrite = TRUE)
    terra::writeVector(lakes_ne, lakes_path, overwrite = TRUE)
    
    message("✓ Natural Earth data downloaded and saved")
    
    return(list(
      canada = canada_path,
      lakes = lakes_path
    ))
    
  }, error = function(e) {
    stop(paste("Error downloading Natural Earth data:", e$message))
  })
}

#' Check and validate required input files
#' 
#' @param rawdata_dir Character. Directory containing data files
#' @return Logical. TRUE if all required files exist
validate_input_files <- function(rawdata_dir = here::here("rawdata")) {
  
  required_files <- c(
    file.path(rawdata_dir, "species_ranges.csv"),
    file.path(rawdata_dir, "fungalroot", "clean_fungalroot_genus.csv"),
    file.path(rawdata_dir, "fungalroot", "clean_fungalroot_species.csv"),
    file.path(rawdata_dir, "Stephanie", "EMFungiCanada_2025_06_03.csv")
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  
  if (length(missing_files) > 0) {
    warning("Missing required files:")
    for (file in missing_files) {
      warning(paste(" -", file))
    }
    return(FALSE)
  }
  
  message("✓ All required input files found")
  return(TRUE)
}

#' Main function to download all required data
#' 
#' @param rawdata_dir Character. Directory to save data
#' @param outdata_dir Character. Directory for output data
#' @return List of file paths to downloaded data
download_all_data <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Starting data download process...")
  
  # Validate input files first
  if (!validate_input_files(rawdata_dir)) {
    stop("Required input files are missing. Please check rawdata directory.")
  }
  
  # Download external data
  admin_paths <- download_admin_boundaries(rawdata_dir, outdata_dir)
  ecoregions_path <- download_ecoregions(rawdata_dir, outdata_dir)
  ustreeatlas_path <- download_ustreeatlas(rawdata_dir, outdata_dir)
  ne_paths <- download_natural_earth(rawdata_dir, outdata_dir)
  
  # Create summary of downloaded files
  download_summary <- list(
    admin_boundaries = admin_paths,
    ecoregions = ecoregions_path,
    ustreeatlas = ustreeatlas_path,
    natural_earth = ne_paths,
    timestamp = Sys.time()
  )
  
  # Save download summary
  saveRDS(download_summary, file.path(rawdata_dir, "download_summary.rds"))
  
  message("✓ All data download completed successfully")
  
  return(download_summary)
}