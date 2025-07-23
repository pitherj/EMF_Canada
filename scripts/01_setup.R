# =============================================================================
# EMF Canada Analysis - Setup Script
# =============================================================================
# Purpose: Load libraries, set options, create directories
# Author: Jason Pither, Clara Qin, Isaac Eckert
# =============================================================================

# Required packages for setup
required_packages <- c("here")

# Load packages with error handling
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
}

#' Setup the EMF Canada analysis environment
#' 
#' This function sets global options and creates necessary directories.
#' 
#' @param project_root Character. Path to the project root directory
#' @return NULL (invisible)
setup_emf_analysis <- function(project_root = here::here()) {
  
  message("Setting up EMF Canada analysis environment...")
  
  # Set global options
  options(
    timeout = max(6000, getOption("timeout")),  # For large file downloads
    stringsAsFactors = FALSE,
    scipen = 999
  )
  
  # Set encoding
  Sys.setlocale("LC_CTYPE", "en_CA.UTF-8")
  
  # Create directory structure
  dirs_to_create <- c(
    file.path(project_root, "rawdata", "vanGalen_Data"),
    file.path(project_root, "rawdata", "fungalroot"),
    file.path(project_root, "rawdata", "Stephanie"),
    file.path(project_root, "rawdata", "ecoregions"),
    file.path(project_root, "rawdata", "USTreeAtlas"),
    file.path(project_root, "outdata"),
    file.path(project_root, "figures"),
    file.path(project_root, "scripts")
  )
  
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      message(paste("✓ Created directory:", dir))
    }
  }
  
  # Set coordinate reference system (make them available globally)
  crdref <<- "+proj=longlat +datum=WGS84"
  canada_albers <<- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  
  message("✓ Coordinate systems defined globally")
  message("✓ EMF analysis environment setup complete")
  
  return(invisible(NULL))
}

#' Get session information and save to file
#' 
#' @param output_dir Character. Directory to save session info
#' @return Data frame with session information
save_session_info <- function(output_dir = here::here("outdata")) {
  
  session_info <- sessionInfo()
  
  # Create a formatted session info table
  session_df <- data.frame(
    Component = c("R Version", "Platform", "Running under", "Locale"),
    Details = c(
      paste(session_info$R.version$major, session_info$R.version$minor, sep = "."),
      session_info$platform,
      session_info$running,
      session_info$locale
    ),
    stringsAsFactors = FALSE
  )
  
  # Save as CSV
  write.csv(session_df, file.path(output_dir, "session_info.csv"), row.names = FALSE)
  
  # Save full session info as RDS
  saveRDS(session_info, file.path(output_dir, "session_info.rds"))
  
  message("✓ Session information saved")
  
  return(session_df)
}

# Run setup when sourced
setup_emf_analysis()