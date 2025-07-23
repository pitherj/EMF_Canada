# =============================================================================
# EMF Canada Analysis - Spatial Data Processing Script
# =============================================================================
# Purpose: Process spatial data, create boundaries, handle projections
# Author: Jason Pither, Clara Qin, Isaac Eckert
# =============================================================================

# Required packages
required_packages <- c("here", "terra", "sf", "readr", "tidyterra", "dplyr")

# Load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
}

# Load coordinate systems
if (!exists("canada_albers")) {
  source(here::here("scripts", "01_setup.R"))
}

# Load utility functions for temp directory management
source(here::here("scripts", "utils.R"))

#' Safely load terra objects with error handling
#' 
#' @param filepath Character. Path to the spatial file
#' @return SpatVector object
safe_load_terra <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  
  tryCatch({
    terra::vect(filepath)
  }, error = function(e) {
    stop(paste("Error loading spatial file", filepath, ":", e$message))
  })
}

#' Process administrative boundaries
#' 
#' @param rawdata_dir Character. Directory containing boundary files
#' @param outdata_dir Character. Directory to save processed files
#' @return Character. Path to simplified Canada boundary file
process_admin_boundaries <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Processing administrative boundaries...")
  
  # Create temporary directory for intermediate processing
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Load raw boundaries
  canada_bound <- safe_load_terra(file.path(rawdata_dir, "admin_boundaries", "canada_bound_raw.gpkg"))
  usa_bound <- safe_load_terra(file.path(rawdata_dir, "admin_boundaries", "usa_bound_raw.gpkg"))
  mexico_bound <- safe_load_terra(file.path(rawdata_dir, "admin_boundaries", "mexico_bound_raw.gpkg"))
  
  # Filter USA to exclude Hawaii and territories
  usa_filtered <- usa_bound[!usa_bound$NAME_1 %in% c(
    "Hawaii", "Puerto Rico", "Guam", "American Samoa",
    "Northern Mariana Islands", "United States Virgin Islands"
  ), ]
  
  # Crop USA to continental extent
  usa_cont_extent <- terra::ext(-170, -65, 24, 75)
  usa_contiguous <- terra::crop(usa_filtered, usa_cont_extent)
  
  # Union individual country polygons (save intermediate files to temp)
  canada_unioned <- terra::aggregate(canada_bound, dissolve = TRUE)
  usa_unioned <- terra::aggregate(usa_contiguous, dissolve = TRUE)
  mexico_unioned <- terra::aggregate(mexico_bound, dissolve = TRUE)
  
  # Save intermediate files to temp directory
  terra::writeVector(canada_unioned, file.path(temp_dir, "canada_unioned_temp.gpkg"), overwrite = TRUE)
  terra::writeVector(usa_unioned, file.path(temp_dir, "usa_unioned_temp.gpkg"), overwrite = TRUE)
  terra::writeVector(mexico_unioned, file.path(temp_dir, "mexico_unioned_temp.gpkg"), overwrite = TRUE)
  
  # Merge countries
  north_america <- terra::vect(c(canada_unioned, usa_unioned, mexico_unioned))
  
  # Simplify geometries
  north_america_simplified <- terra::simplifyGeom(north_america, tolerance = 0.1)
  canada_simple <- terra::simplifyGeom(canada_unioned, tolerance = 0.1)
  
  # Save processed boundaries to final locations
  terra::writeVector(canada_simple, file.path(outdata_dir, "canada_simple.gpkg"), overwrite = TRUE)
  terra::writeVector(north_america_simplified, file.path(outdata_dir, "north_america_simple.gpkg"), overwrite = TRUE)
  
  message("✓ Administrative boundaries processed")
  
  return(file.path(outdata_dir, "canada_simple.gpkg"))
}

#' Process tree species range maps from USTreeAtlas
#' 
#' @param rawdata_dir Character. Directory containing data
#' @param outdata_dir Character. Directory to save processed files
#' @return Character. Path to processed tree ranges file
process_tree_ranges <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Processing tree species range maps from USTreeAtlas...")
  
  # Create temporary directory for intermediate processing
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Load range metadata
  range_metadata <- readr::read_csv(file.path(rawdata_dir, "species_ranges.csv"), show_col_types = FALSE)
  
  # Path to USTreeAtlas shapefiles
  ustreeatlas_shp_dir <- file.path(rawdata_dir, "USTreeAtlas", "USTreeAtlas-main", "shp")
  
  if (!dir.exists(ustreeatlas_shp_dir)) {
    stop("USTreeAtlas shapefiles directory not found. Run download script first.")
  }
  
  # Create vector of filenames using USTreeAtlas structure
  range_filenames <- file.path(ustreeatlas_shp_dir, range_metadata$code, range_metadata$shapefile)
  
  # Check if range files exist
  missing_ranges <- range_filenames[!file.exists(range_filenames)]
  if (length(missing_ranges) > 0) {
    warning("Missing range map files:")
    for (file in missing_ranges) {
      warning(paste(" -", file))
    }
  }
  
  existing_ranges <- range_filenames[file.exists(range_filenames)]
  existing_metadata <- range_metadata[file.exists(range_filenames), ]
  
  if (length(existing_ranges) == 0) {
    stop("No range map files found in USTreeAtlas")
  }
  
  # Process existing range files
  message(paste("Processing", length(existing_ranges), "range map files..."))
  
  all_ranges_list <- list()
  
  for (i in seq_along(existing_ranges)) {
    tryCatch({
      # Read the shapefile
      species_shp <- terra::vect(existing_ranges[i])
      
      # Union all polygons for this species
      species_union <- terra::aggregate(species_shp, dissolve = TRUE)
      
      # Add species identifier
      species_union$species_name <- existing_metadata$code[i]
      
      all_ranges_list[[i]] <- species_union
      
      # Save intermediate result every 50 species to temp directory
      if (i %% 50 == 0) {
        temp_file <- file.path(temp_dir, paste0("ranges_batch_", i, ".gpkg"))
        if (length(all_ranges_list) > 0) {
          valid_ranges <- all_ranges_list[!sapply(all_ranges_list, is.null)]
          if (length(valid_ranges) > 0) {
            batch_ranges <- do.call(rbind, valid_ranges)
            terra::writeVector(batch_ranges, temp_file, overwrite = TRUE)
          }
        }
        message(paste("Processed", i, "of", length(existing_ranges), "species"))
      }
      
    }, error = function(e) {
      warning(paste("Error processing", existing_ranges[i], ":", e$message))
    })
  }
  
  # Combine all valid ranges
  valid_ranges <- all_ranges_list[!sapply(all_ranges_list, is.null)]
  
  if (length(valid_ranges) == 0) {
    stop("No valid range maps could be processed")
  }
  
  all_tree_species_ranges <- do.call(rbind, valid_ranges)
  
  # Load Canada boundary and filter ranges
  canada_simple <- safe_load_terra(file.path(outdata_dir, "canada_simple.gpkg"))
  
  # Find ranges that intersect with Canada
  canada_tree_species_ranges <- all_tree_species_ranges[canada_simple, ]
  
  # Add taxonomic info
  canada_tree_species_ranges <- tidyterra::left_join(
    canada_tree_species_ranges,
    existing_metadata, 
    by = c("species_name" = "code")
  )
  
  # Save processed ranges
  output_path <- file.path(outdata_dir, "tree_ranges_canada.gpkg")
  terra::writeVector(canada_tree_species_ranges, output_path, overwrite = TRUE)
  
  message(paste("✓ Tree ranges processed:", nrow(canada_tree_species_ranges), "species in Canada"))
  
  return(output_path)
}

#' Process ecoregions data
#' 
#' @param rawdata_dir Character. Directory containing data
#' @param outdata_dir Character. Directory to save processed files
#' @param ecoregions_path Character. Path to ecoregions shapefile (from download)
#' @return Character. Path to processed ecoregions file
process_ecoregions <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata"), ecoregions_path = NULL) {
  
  message("Processing ecoregions...")
  
  # Create temporary directory for intermediate processing
  temp_dir <- create_temp_dir(outdata_dir)
  
  # If path not provided, try to find it
  if (is.null(ecoregions_path)) {
    ecoregions_dir <- file.path(rawdata_dir, "ecoregions")
    shp_files <- list.files(ecoregions_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
    
    if (length(shp_files) == 0) {
      stop("No ecoregions shapefile found. Run download script first.")
    }
    
    ecoregions_path <- shp_files[1]
  }
  
  # Load ecoregions
  ecoregions <- safe_load_terra(ecoregions_path)
  
  # Make valid geometries and save intermediate result
  ecoregions <- terra::makeValid(ecoregions)
  terra::writeVector(ecoregions, file.path(temp_dir, "ecoregions_valid_temp.gpkg"), overwrite = TRUE)
  
  # Project to WGS84
  terra::crs(ecoregions) <- "EPSG:4326"
  
  # Load ecozone names from the same directory as the shapefile
  ecozone_names_path <- file.path(dirname(ecoregions_path), "ecozone_names.dbf")
  if (file.exists(ecozone_names_path)) {
    ecozone_names <- sf::st_read(ecozone_names_path, quiet = TRUE)
    ecozone_names <- as.data.frame(ecozone_names)  # Convert to data frame
    
    # Save ecozone names
    write.csv(ecozone_names, file.path(outdata_dir, "ecozone_names.csv"), row.names = FALSE)
  }
  
  # Save processed ecoregions
  output_path <- file.path(outdata_dir, "ecoregions_processed.gpkg")
  terra::writeVector(ecoregions, output_path, overwrite = TRUE)
  
  message("✓ Ecoregions processed")
  
  return(output_path)
}

#' Process Natural Earth data and project to Canada Albers
#' 
#' @param rawdata_dir Character. Directory containing data
#' @param outdata_dir Character. Directory to save processed files
#' @return List of paths to projected Natural Earth data
process_natural_earth <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Processing Natural Earth data...")
  
  # Create temporary directory for intermediate processing
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Load Natural Earth data
  ne_dir <- file.path(rawdata_dir, "natural_earth")
  canada_ne <- safe_load_terra(file.path(ne_dir, "canada_ne.gpkg"))
  lakes_ne <- safe_load_terra(file.path(ne_dir, "lakes_ne.gpkg"))
  
  # Project to Canada Albers using the CRS string
  canada_ne_albers <- terra::project(canada_ne, canada_albers)
  lakes_ne_albers <- terra::project(lakes_ne, canada_albers)
  
  # Save intermediate projected files to temp
  terra::writeVector(canada_ne_albers, file.path(temp_dir, "canada_ne_albers_temp.gpkg"), overwrite = TRUE)
  terra::writeVector(lakes_ne_albers, file.path(temp_dir, "lakes_ne_albers_temp.gpkg"), overwrite = TRUE)
  
  # Clip lakes to Canada boundaries
  lakes_canada <- terra::intersect(lakes_ne_albers, canada_ne_albers)
  
  # Save projected data to final locations
  paths <- list(
    canada_albers = file.path(outdata_dir, "canada_ne_albers.gpkg"),
    lakes_canada = file.path(outdata_dir, "lakes_canada_albers.gpkg")
  )
  
  terra::writeVector(canada_ne_albers, paths$canada_albers, overwrite = TRUE)
  terra::writeVector(lakes_canada, paths$lakes_canada, overwrite = TRUE)
  
  message("✓ Natural Earth data processed and projected")
  
  return(paths)
}

#' Main function to process all spatial data
#' 
#' @param rawdata_dir Character. Directory containing data
#' @param outdata_dir Character. Directory to save processed files
#' @param download_summary List. Summary from download operations (optional)
#' @return List of paths to processed spatial files
process_all_spatial <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata"), download_summary = NULL) {
  
  message("Starting spatial data processing...")
  
  # Create temporary directory
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Check for existing checkpoint in temp directory
  checkpoint_path <- file.path(temp_dir, "spatial_processing_checkpoint.rds")
  if (file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    message("Found existing spatial processing checkpoint...")
  } else {
    checkpoint <- list(
      canada_boundary = FALSE,
      tree_ranges = FALSE,
      ecoregions = FALSE,
      natural_earth = FALSE
    )
  }
  
  # Process admin boundaries
  canada_boundary_path <- file.path(outdata_dir, "canada_simple.gpkg")
  if (file.exists(canada_boundary_path)) {
    message("✓ Administrative boundaries already processed")
    canada_boundary <- canada_boundary_path
    checkpoint$canada_boundary <- TRUE
  } else {
    message("Processing administrative boundaries...")
    canada_boundary <- process_admin_boundaries(rawdata_dir, outdata_dir)
    checkpoint$canada_boundary <- TRUE
  }
  saveRDS(checkpoint, checkpoint_path)
  
  # Process tree ranges
  tree_ranges_path <- file.path(outdata_dir, "tree_ranges_canada.gpkg")
  if (file.exists(tree_ranges_path)) {
    message("✓ Tree ranges already processed")
    tree_ranges <- tree_ranges_path
    checkpoint$tree_ranges <- TRUE
  } else {
    message("Processing tree species range maps...")
    tree_ranges <- process_tree_ranges(rawdata_dir, outdata_dir)
    checkpoint$tree_ranges <- TRUE
  }
  saveRDS(checkpoint, checkpoint_path)
  
  # Process ecoregions
  ecoregions_path <- file.path(outdata_dir, "ecoregions_processed.gpkg")
  if (file.exists(ecoregions_path)) {
    message("✓ Ecoregions already processed")
    ecoregions <- ecoregions_path
    checkpoint$ecoregions <- TRUE
  } else {
    message("Processing ecoregions...")
    # Get ecoregions path from download summary if available, otherwise discover it
    ecoregions_source_path <- NULL
    if (!is.null(download_summary) && "ecoregions" %in% names(download_summary)) {
      ecoregions_source_path <- download_summary$ecoregions
    }
    ecoregions <- process_ecoregions(rawdata_dir, outdata_dir, ecoregions_source_path)
    checkpoint$ecoregions <- TRUE
  }
  saveRDS(checkpoint, checkpoint_path)
  
  # Process Natural Earth
  natural_earth_path <- file.path(outdata_dir, "canada_ne_albers.gpkg")
  if (file.exists(natural_earth_path)) {
    message("✓ Natural Earth data already processed")
    natural_earth <- list(
      canada_albers = file.path(outdata_dir, "canada_ne_albers.gpkg"),
      lakes_canada = file.path(outdata_dir, "lakes_canada_albers.gpkg")
    )
    checkpoint$natural_earth <- TRUE
  } else {
    message("Processing Natural Earth data...")
    natural_earth <- process_natural_earth(rawdata_dir, outdata_dir)
    checkpoint$natural_earth <- TRUE
  }
  saveRDS(checkpoint, checkpoint_path)
  
  # Create summary
  spatial_summary <- list(
    canada_boundary = canada_boundary,
    tree_ranges = tree_ranges,
    ecoregions = ecoregions,
    natural_earth = natural_earth,
    timestamp = Sys.time(),
    checkpoint_completed = TRUE
  )
  
  # Save summary to main outdata directory
  saveRDS(spatial_summary, file.path(outdata_dir, "spatial_processing_summary.rds"))
  
  message("✓ All spatial data processing completed")
  
  return(spatial_summary)
}