# =============================================================================
# EMF Canada Analysis - Metrics Calculation Script
# =============================================================================
# Purpose: Calculate richness, coverage, and spatial metrics
# Author: Jason Pither, Clara Qin, Isaac Eckert
# =============================================================================

# Required packages
required_packages <- c("here", "terra")

# Load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
}

# Load coordinate systems and spatial functions
if (!exists("canada_albers")) {
  source(here::here("scripts", "01_setup.R"))
}
if (!exists("safe_load_terra")) {
  source(here::here("scripts", "03_process_spatial.R"))
}

#' Calculate ecoregion overlap with EMF habitat
#' 
#' @param outdata_dir Character. Directory containing data files
#' @return List containing overlap statistics and file paths
calculate_ecoregion_overlap <- function(outdata_dir = here::here("outdata")) {
  
  message("Calculating ecoregion overlap with EMF habitat...")
  
  # Load required data
  emf_hosts <- safe_load_terra(file.path(outdata_dir, "emf_host_tree_ranges.gpkg"))
  ecoregions <- safe_load_terra(file.path(outdata_dir, "ecoregions_processed.gpkg"))
  
  # Ensure consistent CRS
  terra::crs(emf_hosts) <- "EPSG:4326"
  terra::crs(ecoregions) <- "EPSG:4326"
  
  # Union all EMF host ranges
  emf_habitat_union <- terra::aggregate(emf_hosts, dissolve = TRUE)
  
  # Make geometries valid
  ecoregions <- terra::makeValid(ecoregions)
  emf_habitat_union <- terra::makeValid(emf_habitat_union)
  
  # Find ecoregions that intersect with EMF habitat
  overlapping_indices <- terra::relate(ecoregions, emf_habitat_union, "intersects")
  ecoregions_with_habitat <- ecoregions[overlapping_indices, ]
  
  # Calculate statistics
  overlap_stats <- list(
    total_ecoregions = nrow(ecoregions),
    ecoregions_with_habitat = nrow(ecoregions_with_habitat),
    percent_with_habitat = round(100 * nrow(ecoregions_with_habitat) / nrow(ecoregions), 2)
  )
  
  # Save ecoregions with habitat
  habitat_ecoregions_path <- file.path(outdata_dir, "ecoregions_with_emf_habitat.gpkg")
  terra::writeVector(ecoregions_with_habitat, habitat_ecoregions_path, overwrite = TRUE)
  
  message(paste("✓ Ecoregion overlap calculated:", overlap_stats$ecoregions_with_habitat, 
                "of", overlap_stats$total_ecoregions, "ecoregions have EMF habitat"))
  
  return(list(
    stats = overlap_stats,
    habitat_ecoregions_path = habitat_ecoregions_path
  ))
}

#' Calculate sampling coverage by ecoregion
#' 
#' @param outdata_dir Character. Directory containing data files
#' @return List containing sampling coverage statistics
calculate_sampling_coverage <- function(outdata_dir = here::here("outdata")) {
  
  message("Calculating sampling coverage by ecoregion...")
  
  # Load required data
  aggregated_samples <- safe_load_terra(file.path(outdata_dir, "emf_samples_aggregated_1000m.gpkg"))
  habitat_ecoregions <- safe_load_terra(file.path(outdata_dir, "ecoregions_with_emf_habitat.gpkg"))
  
  # Ensure consistent CRS
  terra::crs(aggregated_samples) <- "EPSG:4326"
  terra::crs(habitat_ecoregions) <- "EPSG:4326"
  
  # Intersect samples with habitat ecoregions
  sampled_ecoregions <- terra::intersect(aggregated_samples, habitat_ecoregions)
  
  # Calculate coverage statistics
  coverage_stats <- list(
    habitat_ecoregions_total = length(unique(habitat_ecoregions$ECOREGION)),
    habitat_ecoregions_sampled = length(unique(sampled_ecoregions$ECOREGION)),
    sampling_locations_total = nrow(aggregated_samples),
    sampling_locations_in_habitat = nrow(sampled_ecoregions)
  )
  
  coverage_stats$percent_habitat_ecoregions_sampled <- round(
    100 * coverage_stats$habitat_ecoregions_sampled / coverage_stats$habitat_ecoregions_total, 2
  )
  
  coverage_stats$percent_habitat_ecoregions_unsampled <- round(
    100 * (coverage_stats$habitat_ecoregions_total - coverage_stats$habitat_ecoregions_sampled) / 
    coverage_stats$habitat_ecoregions_total, 2
  )
  
  # Save sampled ecoregions
  sampled_ecoregions_path <- file.path(outdata_dir, "sampled_habitat_ecoregions.gpkg")
  terra::writeVector(sampled_ecoregions, sampled_ecoregions_path, overwrite = TRUE)
  
  message(paste("✓ Sampling coverage calculated:", coverage_stats$habitat_ecoregions_sampled,
                "of", coverage_stats$habitat_ecoregions_total, "habitat ecoregions sampled"))
  
  return(list(
    stats = coverage_stats,
    sampled_ecoregions_path = sampled_ecoregions_path
  ))
}

#' Create raster grid for richness analysis
#' 
#' @param outdata_dir Character. Directory containing data files
#' @param resolution_deg Numeric. Resolution in degrees (default 1)
#' @return Character. Path to raster template file
create_raster_grid <- function(outdata_dir = here::here("outdata"), resolution_deg = 1) {
  
  message(paste("Creating", resolution_deg, "degree raster grid..."))
  
  # Load Canada boundary
  canada_boundary <- safe_load_terra(file.path(outdata_dir, "canada_simple.gpkg"))
  
  # Get extent and round to nearest degree
  canada_extent <- terra::ext(canada_boundary)
  
  xmin_grid <- floor(canada_extent[1] / resolution_deg) * resolution_deg
  xmax_grid <- ceiling(canada_extent[2] / resolution_deg) * resolution_deg
  ymin_grid <- floor(canada_extent[3] / resolution_deg) * resolution_deg
  ymax_grid <- ceiling(canada_extent[4] / resolution_deg) * resolution_deg
  
  # Create raster template
  raster_template <- terra::rast(
    xmin = xmin_grid, xmax = xmax_grid,
    ymin = ymin_grid, ymax = ymax_grid,
    resolution = resolution_deg, crs = "EPSG:4326"
  )
  
  # Initialize with values first
  terra::values(raster_template) <- 1
  
  # Then mask to Canada
  raster_template <- terra::mask(raster_template, canada_boundary)
  
  # Save raster template
  template_path <- file.path(outdata_dir, paste0("raster_template_", resolution_deg, "deg.tif"))
  terra::writeRaster(raster_template, template_path, overwrite = TRUE)
  
  message(paste("✓ Raster grid created:", dim(raster_template)[1], "rows x", 
                dim(raster_template)[2], "columns"))
  
  return(template_path)
}

#' Calculate species richness per raster cell
#' 
#' @param outdata_dir Character. Directory containing data files
#' @param resolution_deg Numeric. Resolution in degrees
#' @return Character. Path to richness raster file
calculate_species_richness <- function(outdata_dir = here::here("outdata"), resolution_deg = 1) {
  
  message("Calculating species richness per raster cell...")
  
  # Load data
  emf_hosts <- safe_load_terra(file.path(outdata_dir, "emf_host_tree_ranges.gpkg"))
  template_path <- file.path(outdata_dir, paste0("raster_template_", resolution_deg, "deg.tif"))
  
  if (!file.exists(template_path)) {
    template_path <- create_raster_grid(outdata_dir, resolution_deg)
  }
  
  raster_template <- terra::rast(template_path)
  
  # Initialize richness raster
  richness_raster <- raster_template
  terra::values(richness_raster) <- 0
  
  # Calculate richness by rasterizing each species and summing
  n_species <- nrow(emf_hosts)
  message(paste("Processing", n_species, "species..."))
  
  for (i in 1:n_species) {
    species_polygon <- emf_hosts[i, ]
    
    # Rasterize species presence
    species_raster <- terra::rasterize(species_polygon, raster_template, field = 1)
    
    # Add to richness (convert NA to 0)
    species_binary <- terra::ifel(is.na(species_raster), 0, 1)
    richness_raster <- richness_raster + species_binary
    
    if (i %% 20 == 0) {
      message(paste("Processed", i, "of", n_species, "species"))
    }
  }
  
  # Save richness raster
  richness_path <- file.path(outdata_dir, paste0("species_richness_", resolution_deg, "deg.tif"))
  terra::writeRaster(richness_raster, richness_path, overwrite = TRUE)
  
  max_richness <- max(terra::values(richness_raster), na.rm = TRUE)
  cells_with_species <- sum(terra::values(richness_raster) > 0, na.rm = TRUE)
  
  message(paste("✓ Species richness calculated - Max richness:", max_richness,
                "Cells with species:", cells_with_species))
  
  return(richness_path)
}

#' Calculate EMF sequence data richness per raster cell
#' 
#' @param outdata_dir Character. Directory containing data files
#' @param resolution_deg Numeric. Resolution in degrees
#' @return Character. Path to EMF richness raster file
calculate_emf_richness <- function(outdata_dir = here::here("outdata"), resolution_deg = 1) {
  
  message("Calculating EMF sequence data richness per raster cell...")
  
  # Load data
  emf_hosts <- safe_load_terra(file.path(outdata_dir, "emf_host_tree_ranges.gpkg"))
  basic_stats <- readRDS(file.path(outdata_dir, "basic_emf_stats.rds"))
  
  template_path <- file.path(outdata_dir, paste0("raster_template_", resolution_deg, "deg.tif"))
  raster_template <- terra::rast(template_path)
  
  # Initialize EMF richness raster
  emf_richness_raster <- raster_template
  terra::values(emf_richness_raster) <- 0
  
  # Get species with EMF sequence data
  if (!"native_hosts_with_sequences" %in% names(basic_stats)) {
    warning("No native hosts with sequences found in basic stats")
    return(NULL)
  }
  
  # Load sample data to get species with sequences
  sample_points <- safe_load_terra(file.path(outdata_dir, "emf_samples_canada.gpkg"))
  stephanie_samples <- sample_points[sample_points$dataset == "Stephanie", ]
  
  if (nrow(stephanie_samples) == 0) {
    warning("No Stephanie samples found")
    return(NULL)
  }
  
  stephanie_df <- as.data.frame(stephanie_samples)
  
  if (!"Plant_Species" %in% names(stephanie_df)) {
    warning("Plant_Species column not found in sample data")
    return(NULL)
  }
  
  # Get unique host species with sequence data
  unique_hosts <- unique(stephanie_df$Plant_Species)
  unique_hosts <- unique_hosts[!is.na(unique_hosts)]
  hosts_formatted <- gsub(" ", "_", unique_hosts)
  
  # Find matching species in EMF host ranges
  species_with_data <- emf_hosts$species %in% gsub("_", " ", hosts_formatted)
  emf_hosts_with_data <- emf_hosts[species_with_data, ]
  
  if (nrow(emf_hosts_with_data) == 0) {
    warning("No matching EMF host species found")
    return(NULL)
  }
  
  # Calculate EMF richness
  n_emf_species <- nrow(emf_hosts_with_data)
  message(paste("Processing", n_emf_species, "EMF species with sequence data..."))
  
  for (i in 1:n_emf_species) {
    species_polygon <- emf_hosts_with_data[i, ]
    
    # Rasterize species presence
    species_raster <- terra::rasterize(species_polygon, raster_template, field = 1)
    
    # Add to EMF richness
    species_binary <- terra::ifel(is.na(species_raster), 0, 1)
    emf_richness_raster <- emf_richness_raster + species_binary
  }
  
  # Save EMF richness raster
  emf_richness_path <- file.path(outdata_dir, paste0("emf_richness_", resolution_deg, "deg.tif"))
  terra::writeRaster(emf_richness_raster, emf_richness_path, overwrite = TRUE)
  
  max_emf_richness <- max(terra::values(emf_richness_raster), na.rm = TRUE)
  cells_with_emf <- sum(terra::values(emf_richness_raster) > 0, na.rm = TRUE)
  
  message(paste("✓ EMF richness calculated - Max EMF richness:", max_emf_richness,
                "Cells with EMF data:", cells_with_emf))
  
  return(emf_richness_path)
}

#' Calculate proportion of species with EMF data per cell
#' 
#' @param outdata_dir Character. Directory containing data files
#' @param resolution_deg Numeric. Resolution in degrees
#' @return List containing proportion raster path and summary statistics
calculate_emf_proportion <- function(outdata_dir = here::here("outdata"), resolution_deg = 1) {
  
  message("Calculating proportion of species with EMF data per cell...")
  
  # Load richness rasters
  richness_path <- file.path(outdata_dir, paste0("species_richness_", resolution_deg, "deg.tif"))
  emf_richness_path <- file.path(outdata_dir, paste0("emf_richness_", resolution_deg, "deg.tif"))
  
  if (!file.exists(richness_path) || !file.exists(emf_richness_path)) {
    stop("Richness rasters not found. Calculate richness first.")
  }
  
  richness_raster <- terra::rast(richness_path)
  emf_richness_raster <- terra::rast(emf_richness_path)
  
  # Calculate proportion (avoid division by zero)
  proportion_raster <- terra::ifel(
    richness_raster > 0,
    emf_richness_raster / richness_raster,
    NA
  )
  
  # Save proportion raster
  proportion_path <- file.path(outdata_dir, paste0("emf_proportion_", resolution_deg, "deg.tif"))
  terra::writeRaster(proportion_raster, proportion_path, overwrite = TRUE)
  
  # Calculate summary statistics
  prop_values <- terra::values(proportion_raster)
  prop_values_clean <- prop_values[!is.na(prop_values)]
  
  proportion_stats <- list(
    cells_with_data = length(prop_values_clean),
    cells_zero_proportion = sum(prop_values_clean == 0, na.rm = TRUE),
    cells_nonzero_proportion = sum(prop_values_clean > 0, na.rm = TRUE),
    mean_proportion = round(mean(prop_values_clean, na.rm = TRUE), 3),
    max_proportion = round(max(prop_values_clean, na.rm = TRUE), 3),
    min_proportion = round(min(prop_values_clean, na.rm = TRUE), 3)
  )
  
  # Save statistics
  saveRDS(proportion_stats, file.path(outdata_dir, "proportion_stats.rds"))
  
  message(paste("✓ EMF proportion calculated - Mean proportion:", proportion_stats$mean_proportion,
                "Max proportion:", proportion_stats$max_proportion))
  
  return(list(
    proportion_path = proportion_path,
    stats = proportion_stats
  ))
}

#' Main function to calculate all metrics
#' 
#' @param outdata_dir Character. Directory containing data
#' @param resolution_deg Numeric. Resolution in degrees for raster analysis
#' @return List containing all calculated metrics and file paths
calculate_all_metrics <- function(outdata_dir = here::here("outdata"), resolution_deg = 1) {
  
  message("Starting comprehensive metrics calculation...")
  
  # Calculate ecoregion overlap
  overlap_results <- calculate_ecoregion_overlap(outdata_dir)
  
  # Calculate sampling coverage
  coverage_results <- calculate_sampling_coverage(outdata_dir)
  
  # Create raster grid
  template_path <- create_raster_grid(outdata_dir, resolution_deg)
  
  # Calculate species richness
  richness_path <- calculate_species_richness(outdata_dir, resolution_deg)
  
  # Calculate EMF richness
  emf_richness_path <- calculate_emf_richness(outdata_dir, resolution_deg)
  
  # Calculate EMF proportion
  proportion_results <- calculate_emf_proportion(outdata_dir, resolution_deg)
  
  # Compile all results
  metrics_summary <- list(
    ecoregion_overlap = overlap_results,
    sampling_coverage = coverage_results,
    raster_template = template_path,
    species_richness = richness_path,
    emf_richness = emf_richness_path,
    emf_proportion = proportion_results,
    resolution_deg = resolution_deg,
    timestamp = Sys.time()
  )
  
  # Save comprehensive summary
  saveRDS(metrics_summary, file.path(outdata_dir, "metrics_summary.rds"))
  
  message("✓ All metrics calculation completed")
  
  return(metrics_summary)
}