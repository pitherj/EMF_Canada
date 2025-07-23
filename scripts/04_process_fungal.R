# =============================================================================
# EMF Canada Analysis - Fungal Data Processing Script
# =============================================================================
# Purpose: Process FungalRoot database and EMF sample data
# Author: Jason Pither, Clara Qin, Isaac Eckert
# =============================================================================

# Required packages
required_packages <- c("here", "terra", "readr", "tidyterra", "dplyr", "magrittr")

# Load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
}

# Load coordinate systems and spatial functions
if (!exists("crdref")) {
  source(here::here("scripts", "01_setup.R"))
}
if (!exists("safe_load_terra")) {
  source(here::here("scripts", "03_process_spatial.R"))
}

# Load utility functions for temp directory management
source(here::here("scripts", "utils.R"))

#' Load and validate FungalRoot database files
#' 
#' @param rawdata_dir Character. Directory containing data files
#' @return List containing species and genus level mycorrhizal data
load_fungalroot_data <- function(rawdata_dir = here::here("rawdata")) {
  
  message("Loading FungalRoot database...")
  
  # File paths
  species_path <- file.path(rawdata_dir, "fungalroot", "clean_fungalroot_species.csv")
  genus_path <- file.path(rawdata_dir, "fungalroot", "clean_fungalroot_genus.csv")
  
  # Check file existence
  if (!file.exists(species_path)) {
    stop(paste("Species-level FungalRoot file not found:", species_path))
  }
  if (!file.exists(genus_path)) {
    stop(paste("Genus-level FungalRoot file not found:", genus_path))
  }
  
  # Load FungalRoot database
  fungalroot_species <- readr::read_csv(species_path, show_col_types = FALSE)
  fungalroot_genus <- readr::read_csv(genus_path, show_col_types = FALSE)
  
  # Validate required columns
  required_species_cols <- c("UpdatedPlantBinomial", "MycorrhizalType")
  required_genus_cols <- c("UpdatedGenus", "MycorrhizalType")
  
  if (!all(required_species_cols %in% names(fungalroot_species))) {
    stop("Missing required columns in species FungalRoot data")
  }
  if (!all(required_genus_cols %in% names(fungalroot_genus))) {
    stop("Missing required columns in genus FungalRoot data")
  }
  
  message(paste("✓ FungalRoot data loaded:", 
                nrow(fungalroot_species), "species records,",
                nrow(fungalroot_genus), "genus records"))
  
  return(list(
    species = fungalroot_species,
    genus = fungalroot_genus
  ))
}

#' Assign mycorrhizal types to tree species
#' 
#' @param tree_ranges_path Character. Path to tree ranges spatial file
#' @param fungalroot_data List. FungalRoot database from load_fungalroot_data()
#' @param outdata_dir Character. Directory to save output
#' @return Character. Path to EMF host tree species file
assign_mycorrhizal_types <- function(tree_ranges_path, fungalroot_data, outdata_dir = here::here("outdata")) {
  
  message("Assigning mycorrhizal types to tree species...")
  
  # Create temporary directory for intermediate processing
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Load tree ranges
  canada_tree_ranges <- safe_load_terra(tree_ranges_path)
  
  # Join with species-level mycorrhizal info
  canada_myco_host_range_files <- tidyterra::left_join(canada_tree_ranges, 
                                                       fungalroot_data$species, by = c("species" = "UpdatedPlantBinomial"))
  
  # Rename variable to denote species-level myco designation
  canada_myco_host_range_files <- canada_myco_host_range_files %>%
    rename(MycorrhizalType_SP = MycorrhizalType)
  
  # Save intermediate result to temp
  terra::writeVector(canada_myco_host_range_files, 
                     file.path(temp_dir, "myco_host_species_temp.gpkg"), 
                     overwrite = TRUE)
  
  # Get genus-level myco assignment
  canada_myco_host_range_files <- tidyterra::left_join(canada_myco_host_range_files, 
                                                       fungalroot_data$genus, by = c("genus" = "UpdatedGenus"))
  
  # Rename variable to denote genus-level myco designation
  canada_myco_host_range_files <- canada_myco_host_range_files %>%
    rename(MycorrhizalType_GN = MycorrhizalType)
  
  # Filter to EM hosts
  canada_emf_host_range_files <- canada_myco_host_range_files %>%
    filter(grepl("EM", MycorrhizalType_SP) | grepl("EM", MycorrhizalType_GN))
  
  # Clean up column names to avoid duplicates and issues
  # Get column names and remove duplicates
  col_names <- names(canada_emf_host_range_files)
  
  # Remove any duplicate or problematic columns
  # Keep only the columns we need
  keep_cols <- c("species_name", "species", "genus", "common", "sp", "shapefile", 
                 "MycorrhizalType_SP", "MycorrhizalType_GN")
  
  # Find which columns actually exist
  existing_cols <- intersect(keep_cols, col_names)
  
  # Select only the columns we want to keep
  canada_emf_host_range_files <- canada_emf_host_range_files[, existing_cols]
  
  # Save EMF host ranges
  output_path <- file.path(outdata_dir, "emf_host_tree_ranges.gpkg")
  
  # Remove file if it exists to avoid conflicts
  if (file.exists(output_path)) {
    file.remove(output_path)
  }
  
  terra::writeVector(canada_emf_host_range_files, output_path, overwrite = TRUE)
  
  message(paste("✓ Mycorrhizal types assigned:", 
                nrow(canada_emf_host_range_files), "EMF host species identified"))
  
  return(output_path)
}

#' Load and process EMF sample data
#' 
#' @param rawdata_dir Character. Directory containing data files
#' @param outdata_dir Character. Directory to save output
#' @return Character. Path to processed sample points file
process_emf_samples <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Processing EMF sample data...")
  
  # Create temporary directory for intermediate processing
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Load Stephanie's data
  stephanie_path <- file.path(rawdata_dir, "Stephanie", "EMFungiCanada_2025_06_03.csv")
  if (!file.exists(stephanie_path)) {
    stop(paste("Stephanie's EMF data not found:", stephanie_path))
  }
  
  stephanie_data <- readr::read_csv(stephanie_path, show_col_types = FALSE)
  
  # Create spatial points for Stephanie's data
  stephanie_pts <- terra::vect(stephanie_data, 
                               geom = c("Longitude", "Latitude"), 
                               crs = crdref)
  stephanie_pts$dataset <- "Stephanie"
  
  # Save intermediate result to temp
  terra::writeVector(stephanie_pts, 
                     file.path(temp_dir, "stephanie_pts_temp.gpkg"), 
                     overwrite = TRUE)
  
  # Check for van Galen data
  vangalen_path <- file.path(rawdata_dir, "vanGalen_Data", "GFv5_Ecm_unassigned_per_sample.csv")
  
  sample_points_list <- list(stephanie_pts)
  
  if (file.exists(vangalen_path)) {
    message("Loading van Galen data...")
    
    vangalen_data <- readr::read_csv(vangalen_path, show_col_types = FALSE)
    vangalen_canada <- vangalen_data %>% dplyr::filter(Country == "Canada")
    
    if (nrow(vangalen_canada) > 0) {
      vangalen_pts <- terra::vect(vangalen_canada,
                                  geom = c("Longitude", "Latitude"),
                                  crs = crdref)
      vangalen_pts$dataset <- "vanGalen"
      
      # Save intermediate result to temp
      terra::writeVector(vangalen_pts, 
                         file.path(temp_dir, "vangalen_pts_temp.gpkg"), 
                         overwrite = TRUE)
      
      sample_points_list <- list(stephanie_pts, vangalen_pts)
    }
  }
  
  # Combine all sample points
  all_sample_points <- do.call(rbind, sample_points_list)
  
  # Save combined points to temp
  terra::writeVector(all_sample_points, 
                     file.path(temp_dir, "all_sample_points_temp.gpkg"), 
                     overwrite = TRUE)
  
  # Intersect with Canada boundary to ensure points are within Canada
  canada_boundary <- safe_load_terra(file.path(outdata_dir, "canada_simple.gpkg"))
  terra::crs(canada_boundary) <- "EPSG:4326"
  terra::crs(all_sample_points) <- "EPSG:4326"
  
  canada_samples <- terra::intersect(all_sample_points, canada_boundary)
  
  # Save processed samples
  output_path <- file.path(outdata_dir, "emf_samples_canada.gpkg")
  terra::writeVector(canada_samples, output_path, overwrite = TRUE)
  
  message(paste("✓ EMF samples processed:", nrow(canada_samples), "sample points"))
  
  return(output_path)
}

#' Aggregate sample points by proximity
#' 
#' @param sample_points_path Character. Path to sample points file
#' @param distance_threshold Numeric. Distance threshold in meters for aggregation
#' @param outdata_dir Character. Directory to save output
#' @return Character. Path to aggregated sample points file
aggregate_sample_points <- function(sample_points_path, distance_threshold = 1000, outdata_dir = here::here("outdata")) {
  
  message(paste("Aggregating sample points within", distance_threshold, "meters..."))
  
  # Create temporary directory for intermediate processing
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Load sample points
  sample_points <- safe_load_terra(sample_points_path)
  
  # Function to aggregate by proximity
  aggregate_proximity <- function(spatvect, dist = distance_threshold) {
    # Calculate distance matrix
    dist_matrix <- as.matrix(terra::distance(spatvect))
    
    # Create clusters based on distance threshold
    clusters <- rep(NA, nrow(spatvect))
    cluster_id <- 1
    
    for (i in 1:nrow(spatvect)) {
      if (is.na(clusters[i])) {
        # Find all points within distance threshold
        nearby_points <- which(dist_matrix[i, ] <= dist)
        clusters[nearby_points] <- cluster_id
        cluster_id <- cluster_id + 1
      }
    }
    
    # Create result dataframe
    result_data <- data.frame(
      cluster = clusters,
      dataset = spatvect$dataset,
      stringsAsFactors = FALSE
    )
    
    # Determine dataset value for each cluster
    cluster_datasets <- aggregate(dataset ~ cluster, result_data, 
                                  function(x) {
                                    unique_datasets <- unique(x)
                                    if (length(unique_datasets) == 1) {
                                      return(unique_datasets[1])
                                    } else {
                                      return("Both")
                                    }
                                  })
    
    # Get representative point for each cluster (centroid)
    coords <- terra::crds(spatvect)
    cluster_coords <- aggregate(coords, by = list(cluster = clusters), mean)
    
    # Create new SpatVector with aggregated data
    new_coords <- as.matrix(cluster_coords[, c("x", "y")])
    result_spatvect <- terra::vect(new_coords, crs = terra::crs(spatvect))
    
    # Add dataset information
    result_spatvect$dataset <- cluster_datasets$dataset[match(cluster_coords$cluster, cluster_datasets$cluster)]
    
    return(result_spatvect)
  }
  
  # Aggregate points
  aggregated_points <- aggregate_proximity(sample_points, distance_threshold)
  
  # Save intermediate result to temp
  terra::writeVector(aggregated_points, 
                     file.path(temp_dir, paste0("aggregated_points_", distance_threshold, "m_temp.gpkg")), 
                     overwrite = TRUE)
  
  # Save aggregated points
  output_path <- file.path(outdata_dir, paste0("emf_samples_aggregated_", distance_threshold, "m.gpkg"))
  terra::writeVector(aggregated_points, output_path, overwrite = TRUE)
  
  message(paste("✓ Sample points aggregated:", 
                nrow(sample_points), "→", nrow(aggregated_points), "points"))
  
  return(output_path)
}

#' Calculate basic EMF data statistics
#' 
#' @param sample_points_path Character. Path to sample points file
#' @param emf_hosts_path Character. Path to EMF host ranges file
#' @param outdata_dir Character. Directory to save results
#' @return List containing basic statistics
calculate_basic_stats <- function(sample_points_path, emf_hosts_path, outdata_dir = here::here("outdata")) {
  
  message("Calculating basic EMF statistics...")
  
  # Load data
  sample_points <- safe_load_terra(sample_points_path)
  emf_hosts <- safe_load_terra(emf_hosts_path)
  
  # Calculate unique sampling locations for ALL data (Stephanie + van Galen)
  all_samples_df <- as.data.frame(sample_points, geom = "xy")
  all_unique_locations <- all_samples_df %>% 
    select(x, y) %>%
    distinct(x, y)
  total_unique_locations_all_data <- nrow(all_unique_locations)
  
  # Filter to Stephanie's data for detailed analysis
  stephanie_samples <- sample_points[sample_points$dataset == "Stephanie", ]
  
  # Get unique sampling locations for Stephanie's data (with EMF information)
  if (nrow(stephanie_samples) > 0) {
    stephanie_samples_unique_samples <- as.data.frame(stephanie_samples, geom = "xy") %>% 
      select(dataset, x, y) %>%
      distinct(x, y)
    unique_locations_stephanie <- nrow(stephanie_samples_unique_samples)
  } else {
    unique_locations_stephanie <- 0
  }
  
  # Basic counts
  stats <- list(
    total_sequence_records = nrow(stephanie_samples),
    total_sampling_locations_all_data = total_unique_locations_all_data,
    unique_locations_with_emf_data = unique_locations_stephanie,
    emf_host_species_count = nrow(emf_hosts)
  )
  
  # Extract unique values from Stephanie's data if available
  if (nrow(stephanie_samples) > 0) {
    stephanie_df <- as.data.frame(stephanie_samples)
    
    if ("EM_Species" %in% names(stephanie_df)) {
      stats$unique_emf_taxa <- length(unique(stephanie_df$EM_Species))
    }
    
    if ("Plant_Species" %in% names(stephanie_df)) {
      unique_hosts <- unique(stephanie_df$Plant_Species)
      unique_hosts <- unique_hosts[!is.na(unique_hosts)]
      stats$unique_host_species <- length(unique_hosts)
      
      # Compare with Canadian tree species
      host_species_formatted <- gsub(" ", "_", unique_hosts)
      canadian_tree_species <- gsub(" ", "_", emf_hosts$species)
      
      native_hosts_with_sequences <- intersect(host_species_formatted, canadian_tree_species)
      stats$native_hosts_with_sequences <- length(native_hosts_with_sequences)
      stats$percent_hosts_with_data <- round(100 * length(native_hosts_with_sequences) / nrow(emf_hosts), 2)
      
      # Genus-level analysis
      unique_genera <- unique(sapply(strsplit(native_hosts_with_sequences, "_"), function(x) x[1]))
      canadian_genera <- unique(emf_hosts$genus)
      stats$unique_host_genera <- length(unique_genera)
      stats$percent_genera_with_data <- round(100 * length(unique_genera) / length(canadian_genera), 2)
    }
  }
  
  # Save statistics
  saveRDS(stats, file.path(outdata_dir, "basic_emf_stats.rds"))
  write.csv(data.frame(Metric = names(stats), Value = unlist(stats)), 
            file.path(outdata_dir, "basic_emf_stats.csv"), row.names = FALSE)
  
  message("✓ Basic statistics calculated and saved")
  
  return(stats)
}

#' Main function to process all fungal data
#' 
#' @param rawdata_dir Character. Directory containing data
#' @param outdata_dir Character. Directory to save processed files
#' @return List of paths to processed fungal data files
process_all_fungal <- function(rawdata_dir = here::here("rawdata"), outdata_dir = here::here("outdata")) {
  
  message("Starting fungal data processing...")
  
  # Create temporary directory
  temp_dir <- create_temp_dir(outdata_dir)
  
  # Load FungalRoot database
  fungalroot_data <- load_fungalroot_data(rawdata_dir)
  
  # Get tree ranges path (should exist from spatial processing)
  tree_ranges_path <- file.path(outdata_dir, "tree_ranges_canada.gpkg")
  if (!file.exists(tree_ranges_path)) {
    stop("Tree ranges file not found. Run spatial processing first.")
  }
  
  # Assign mycorrhizal types
  emf_hosts_path <- assign_mycorrhizal_types(tree_ranges_path, fungalroot_data, outdata_dir)
  
  # Process sample data
  sample_points_path <- process_emf_samples(rawdata_dir, outdata_dir)
  
  # Aggregate sample points
  aggregated_points_path <- aggregate_sample_points(sample_points_path, 1000, outdata_dir)
  
  # Calculate basic statistics
  basic_stats <- calculate_basic_stats(sample_points_path, emf_hosts_path, outdata_dir)
  
  # Calculate comprehensive EMF taxon analysis (includes legacy EMF taxon sampling)
  comprehensive_results <- calculate_comprehensive_emf_analysis(sample_points_path, outdata_dir)
  
  # Create summary
  fungal_summary <- list(
    emf_hosts = emf_hosts_path,
    sample_points = sample_points_path,
    aggregated_points = aggregated_points_path,
    basic_stats = basic_stats,
    comprehensive_analysis = comprehensive_results,
    timestamp = Sys.time()
  )
  
  # Save summary
  saveRDS(fungal_summary, file.path(outdata_dir, "fungal_processing_summary.rds"))
  
  message("✓ All fungal data processing completed")
  
  return(fungal_summary)
}