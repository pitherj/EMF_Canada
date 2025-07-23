# =============================================================================
# EMF Canada Analysis - Utility Functions
# =============================================================================
# Purpose: Helper functions for data processing and table formatting
# Author: Jason Pither, Clara Qin, Isaac Eckert
# =============================================================================

# Required packages
required_packages <- c("here", "dplyr", "knitr", "terra", "tidyr")

# Load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    # kableExtra and DT are optional - load if available
    if (pkg %in% c("kableExtra", "DT")) {
      next
    }
    stop(paste("Package", pkg, "is required but not installed."))
  }
}

# Try to load optional packages
tryCatch({
  library(kableExtra)
}, error = function(e) {
  message("kableExtra not available - using basic kable")
})

tryCatch({
  library(DT)
}, error = function(e) {
  message("DT not available - interactive tables disabled")
})

#' Create and ensure temporary directory exists
#' 
#' @param outdata_dir Character. Base output directory
#' @return Character. Path to temporary directory
create_temp_dir <- function(outdata_dir = here::here("outdata")) {
  temp_dir <- file.path(outdata_dir, "temp")
  
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
    message(paste("✓ Created temporary directory:", temp_dir))
  }
  
  return(temp_dir)
}

#' Extract genus from species name
#' 
#' @param species_names Character vector of species names
#' @return Character vector of genus names
extract_genus <- function(species_names) {
  # Handle different separators and extract first part
  genus_names <- sapply(species_names, function(x) {
    if (is.na(x) || x == "") return(NA)
    # Split by underscore or space and take first part
    parts <- strsplit(as.character(x), "[_ ]")[[1]]
    return(parts[1])
  })
  return(unname(genus_names))
}

#' Comprehensive EMF taxon analysis
#' 
#' @param sample_points_path Character. Path to sample points file
#' @param outdata_dir Character. Directory to save output
#' @return List containing comprehensive taxon statistics
calculate_comprehensive_emf_analysis <- function(sample_points_path, outdata_dir = here::here("outdata")) {
  
  message("Calculating comprehensive EMF taxon analysis...")
  
  # Load sample points (should be in main outdata directory)
  if (!file.exists(sample_points_path)) {
    stop(paste("Sample points file not found:", sample_points_path))
  }
  
  # Load coordinate systems if not available
  if (!exists("safe_load_terra")) {
    source(here::here("scripts", "03_process_spatial.R"))
  }
  
  sample_points <- safe_load_terra(sample_points_path)
  
  # Filter to Stephanie's data for detailed EMF analysis
  stephanie_samples <- sample_points[sample_points$dataset == "Stephanie", ]
  
  if (nrow(stephanie_samples) == 0) {
    warning("No Stephanie samples found for comprehensive EMF analysis")
    return(NULL)
  }
  
  # Convert to dataframe for analysis
  stephanie_df <- as.data.frame(stephanie_samples, geom = "xy")
  
  # Check for required columns
  required_cols <- c("Other_ID", "UNITE_ID", "EM_Species", "Plant_Species", "x", "y")
  missing_cols <- setdiff(required_cols, names(stephanie_df))
  if (length(missing_cols) > 0) {
    warning(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # Clean the data - remove rows where all taxon identifiers are NA
  stephanie_clean <- stephanie_df %>%
    filter(!is.na(Other_ID) | !is.na(UNITE_ID) | !is.na(EM_Species)) %>%
    mutate(
      # Extract genus from EM_Species
      EM_Genus = extract_genus(EM_Species),
      # Create location identifier
      location_id = paste(round(x, 6), round(y, 6), sep = "_")
    )
  
  if (nrow(stephanie_clean) == 0) {
    warning("No valid EMF data found after cleaning")
    return(NULL)
  }
  
  message(paste("Processing", nrow(stephanie_clean), "EMF records..."))
  
  # =============================================================================
  # 1. BASIC TAXON COUNTS AND DIVERSITY
  # =============================================================================
  
  # Count unique values for each taxonomic level
  unique_counts <- list(
    unique_other_id = length(unique(stephanie_clean$Other_ID[!is.na(stephanie_clean$Other_ID)])),
    unique_unite_id = length(unique(stephanie_clean$UNITE_ID[!is.na(stephanie_clean$UNITE_ID)])),
    unique_em_genera = length(unique(stephanie_clean$EM_Genus[!is.na(stephanie_clean$EM_Genus)])),
    unique_em_species = length(unique(stephanie_clean$EM_Species[!is.na(stephanie_clean$EM_Species)])),
    total_records = nrow(stephanie_clean),
    unique_locations = length(unique(stephanie_clean$location_id))
  )
  
  # =============================================================================
  # 2. SAMPLING LOCATIONS PER TAXON (WALLACEAN SHORTFALL ANALYSIS)
  # =============================================================================
  
  message("Calculating sampling locations per taxon...")
  
  # Function to calculate locations per taxon for a given taxonomic level
  calc_locations_per_taxon <- function(data, taxon_col) {
    data %>%
      filter(!is.na(.data[[taxon_col]])) %>%
      group_by(.data[[taxon_col]]) %>%
      summarise(
        unique_locations = n_distinct(location_id),
        total_records = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(unique_locations)) %>%
      rename(taxon = 1)
  }
  
  # Calculate for each taxonomic level
  locations_per_other_id <- calc_locations_per_taxon(stephanie_clean, "Other_ID")
  locations_per_unite_id <- calc_locations_per_taxon(stephanie_clean, "UNITE_ID")
  locations_per_genus <- calc_locations_per_taxon(stephanie_clean, "EM_Genus")
  locations_per_species <- calc_locations_per_taxon(stephanie_clean, "EM_Species")
  
  # Save detailed location statistics
  write.csv(locations_per_other_id, file.path(outdata_dir, "locations_per_other_id.csv"), row.names = FALSE)
  write.csv(locations_per_unite_id, file.path(outdata_dir, "locations_per_unite_id.csv"), row.names = FALSE)
  write.csv(locations_per_genus, file.path(outdata_dir, "locations_per_genus.csv"), row.names = FALSE)
  write.csv(locations_per_species, file.path(outdata_dir, "locations_per_species.csv"), row.names = FALSE)
  
  # =============================================================================
  # 3. HOST-TAXON ASSOCIATION MATRICES (ELTONIAN SHORTFALL ANALYSIS)
  # =============================================================================
  
  message("Creating host-taxon association matrices...")
  
  # Filter data to only include records with host information
  host_data <- stephanie_clean %>%
    filter(!is.na(Plant_Species) & Plant_Species != "")
  
  if (nrow(host_data) == 0) {
    warning("No host association data found")
    host_associations <- NULL
  } else {
    
    # Function to create host-taxon matrix
    create_host_taxon_matrix <- function(data, taxon_col, matrix_name) {
      # Create presence/absence matrix
      matrix_data <- data %>%
        filter(!is.na(.data[[taxon_col]])) %>%
        select(Plant_Species, !!sym(taxon_col)) %>%
        distinct() %>%
        mutate(presence = 1) %>%
        pivot_wider(
          names_from = !!sym(taxon_col),
          values_from = presence,
          values_fill = 0
        )
      
      # Save the matrix
      write.csv(matrix_data, file.path(outdata_dir, paste0("host_taxon_matrix_", matrix_name, ".csv")), row.names = FALSE)
      
      return(matrix_data)
    }
    
    # Create matrices for each taxonomic level
    host_other_id_matrix <- create_host_taxon_matrix(host_data, "Other_ID", "other_id")
    host_unite_id_matrix <- create_host_taxon_matrix(host_data, "UNITE_ID", "unite_id")
    host_genus_matrix <- create_host_taxon_matrix(host_data, "EM_Genus", "genus")
    host_species_matrix <- create_host_taxon_matrix(host_data, "EM_Species", "species")
    
    # =============================================================================
    # 4. HOST-TAXON ASSOCIATION STATISTICS
    # =============================================================================
    
    message("Calculating host-taxon association statistics...")
    
    # Function to calculate association statistics
    calc_association_stats <- function(data, taxon_col, level_name) {
      clean_data <- data %>%
        filter(!is.na(.data[[taxon_col]]) & !is.na(Plant_Species) & Plant_Species != "")
      
      # EMF taxa per host
      taxa_per_host <- clean_data %>%
        group_by(Plant_Species) %>%
        summarise(
          unique_emf_taxa = n_distinct(.data[[taxon_col]]),
          total_records = n(),
          .groups = "drop"
        ) %>%
        arrange(desc(unique_emf_taxa))
      
      # Hosts per EMF taxon
      hosts_per_taxon <- clean_data %>%
        group_by(.data[[taxon_col]]) %>%
        summarise(
          unique_hosts = n_distinct(Plant_Species),
          total_records = n(),
          .groups = "drop"
        ) %>%
        arrange(desc(unique_hosts)) %>%
        rename(emf_taxon = 1)
      
      # Save results
      write.csv(taxa_per_host, file.path(outdata_dir, paste0("emf_taxa_per_host_", level_name, ".csv")), row.names = FALSE)
      write.csv(hosts_per_taxon, file.path(outdata_dir, paste0("hosts_per_emf_taxon_", level_name, ".csv")), row.names = FALSE)
      
      # Summary statistics
      summary_stats <- list(
        total_host_species = n_distinct(clean_data$Plant_Species),
        total_emf_taxa = n_distinct(clean_data[[taxon_col]]),
        mean_taxa_per_host = round(mean(taxa_per_host$unique_emf_taxa), 2),
        max_taxa_per_host = max(taxa_per_host$unique_emf_taxa),
        mean_hosts_per_taxon = round(mean(hosts_per_taxon$unique_hosts), 2),
        max_hosts_per_taxon = max(hosts_per_taxon$unique_hosts)
      )
      
      return(list(
        taxa_per_host = taxa_per_host,
        hosts_per_taxon = hosts_per_taxon,
        summary = summary_stats
      ))
    }
    
    # Calculate for each taxonomic level
    other_id_associations <- calc_association_stats(host_data, "Other_ID", "other_id")
    unite_id_associations <- calc_association_stats(host_data, "UNITE_ID", "unite_id")
    genus_associations <- calc_association_stats(host_data, "EM_Genus", "genus")
    species_associations <- calc_association_stats(host_data, "EM_Species", "species")
    
    host_associations <- list(
      other_id = other_id_associations,
      unite_id = unite_id_associations,
      genus = genus_associations,
      species = species_associations
    )
  }
  
  # =============================================================================
  # 5. COMPREHENSIVE SUMMARY STATISTICS
  # =============================================================================
  
  # Calculate summary statistics for sampling intensity
  sampling_intensity <- list(
    other_id = if(nrow(locations_per_other_id) > 0) {
      list(
        mean_locations = round(mean(locations_per_other_id$unique_locations), 2),
        median_locations = median(locations_per_other_id$unique_locations),
        max_locations = max(locations_per_other_id$unique_locations),
        min_locations = min(locations_per_other_id$unique_locations)
      )
    } else NULL,
    
    unite_id = if(nrow(locations_per_unite_id) > 0) {
      list(
        mean_locations = round(mean(locations_per_unite_id$unique_locations), 2),
        median_locations = median(locations_per_unite_id$unique_locations),
        max_locations = max(locations_per_unite_id$unique_locations),
        min_locations = min(locations_per_unite_id$unique_locations)
      )
    } else NULL,
    
    genus = if(nrow(locations_per_genus) > 0) {
      list(
        mean_locations = round(mean(locations_per_genus$unique_locations), 2),
        median_locations = median(locations_per_genus$unique_locations),
        max_locations = max(locations_per_genus$unique_locations),
        min_locations = min(locations_per_genus$unique_locations)
      )
    } else NULL,
    
    species = if(nrow(locations_per_species) > 0) {
      list(
        mean_locations = round(mean(locations_per_species$unique_locations), 2),
        median_locations = median(locations_per_species$unique_locations),
        max_locations = max(locations_per_species$unique_locations),
        min_locations = min(locations_per_species$unique_locations)
      )
    } else NULL
  )
  
  # Combine all results
  comprehensive_results <- list(
    unique_counts = unique_counts,
    sampling_intensity = sampling_intensity,
    host_associations = host_associations,
    location_data = list(
      other_id = locations_per_other_id,
      unite_id = locations_per_unite_id,
      genus = locations_per_genus,
      species = locations_per_species
    )
  )
  
  # Save comprehensive summary
  saveRDS(comprehensive_results, file.path(outdata_dir, "comprehensive_emf_analysis.rds"))
  
  # Create summary table
  summary_table <- data.frame(
    Taxonomic_Level = c("Other_ID", "UNITE_ID", "Genus", "Species"),
    Unique_Taxa = c(
      unique_counts$unique_other_id,
      unique_counts$unique_unite_id, 
      unique_counts$unique_em_genera,
      unique_counts$unique_em_species
    ),
    Mean_Locations_Per_Taxon = c(
      ifelse(is.null(sampling_intensity$other_id), NA, sampling_intensity$other_id$mean_locations),
      ifelse(is.null(sampling_intensity$unite_id), NA, sampling_intensity$unite_id$mean_locations),
      ifelse(is.null(sampling_intensity$genus), NA, sampling_intensity$genus$mean_locations),
      ifelse(is.null(sampling_intensity$species), NA, sampling_intensity$species$mean_locations)
    ),
    Max_Locations_Per_Taxon = c(
      ifelse(is.null(sampling_intensity$other_id), NA, sampling_intensity$other_id$max_locations),
      ifelse(is.null(sampling_intensity$unite_id), NA, sampling_intensity$unite_id$max_locations),
      ifelse(is.null(sampling_intensity$genus), NA, sampling_intensity$genus$max_locations),
      ifelse(is.null(sampling_intensity$species), NA, sampling_intensity$species$max_locations)
    ),
    Host_Species_Count = if(!is.null(host_associations)) {
      c(
        host_associations$other_id$summary$total_host_species,
        host_associations$unite_id$summary$total_host_species,
        host_associations$genus$summary$total_host_species,
        host_associations$species$summary$total_host_species
      )
    } else rep(NA, 4),
    Mean_Taxa_Per_Host = if(!is.null(host_associations)) {
      c(
        host_associations$other_id$summary$mean_taxa_per_host,
        host_associations$unite_id$summary$mean_taxa_per_host,
        host_associations$genus$summary$mean_taxa_per_host,
        host_associations$species$summary$mean_taxa_per_host
      )
    } else rep(NA, 4),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_table, file.path(outdata_dir, "comprehensive_emf_summary.csv"), row.names = FALSE)
  
  message("✓ Comprehensive EMF analysis completed")
  message(paste("✓ Unique Other_ID values:", unique_counts$unique_other_id))
  message(paste("✓ Unique UNITE_ID values:", unique_counts$unique_unite_id))
  message(paste("✓ Unique EMF genera:", unique_counts$unique_em_genera))
  message(paste("✓ Unique EMF species:", unique_counts$unique_em_species))
  if(!is.null(host_associations)) {
    message(paste("✓ Host-taxon associations analyzed for", 
                  host_associations$species$summary$total_host_species, "host species"))
  }
  
  return(comprehensive_results)
}

#' Legacy function - now calls comprehensive analysis
#' 
#' @param sample_points_path Character. Path to sample points file
#' @param outdata_dir Character. Directory to save output
#' @return Data frame with EMF taxon sampling statistics
calculate_emf_taxon_sampling <- function(sample_points_path, outdata_dir = here::here("outdata")) {
  
  message("Running comprehensive EMF analysis (legacy function called)...")
  
  # Run comprehensive analysis
  comprehensive_results <- calculate_comprehensive_emf_analysis(sample_points_path, outdata_dir)
  
  if (is.null(comprehensive_results)) {
    return(NULL)
  }
  
  # Extract legacy format for backward compatibility
  emf_taxon_stats <- comprehensive_results$location_data$species
  
  # Create legacy summary stats
  summary_stats <- list(
    total_unique_emf_taxa = comprehensive_results$unique_counts$unique_em_species,
    total_unique_locations_all_taxa = comprehensive_results$unique_counts$unique_locations,
    mean_locations_per_taxon = comprehensive_results$sampling_intensity$species$mean_locations,
    median_locations_per_taxon = comprehensive_results$sampling_intensity$species$median_locations,
    max_locations_per_taxon = comprehensive_results$sampling_intensity$species$max_locations,
    min_locations_per_taxon = comprehensive_results$sampling_intensity$species$min_locations
  )
  
  # Save legacy format files for backward compatibility
  write.csv(emf_taxon_stats, file.path(outdata_dir, "emf_taxon_sampling_stats.csv"), row.names = FALSE)
  saveRDS(summary_stats, file.path(outdata_dir, "emf_taxon_summary_stats.rds"))
  write.csv(t(summary_stats), file.path(outdata_dir, "emf_taxon_summary_stats.csv"), 
            col.names = c("Statistic", "Value"))
  
  return(list(
    taxon_stats = emf_taxon_stats,
    summary_stats = summary_stats
  ))
}

#' Create formatted summary tables
#' 
#' @param outdata_dir Character. Directory containing data files
#' @param output_dir Character. Directory to save tables
#' @return List of formatted tables
create_summary_tables <- function(outdata_dir = here::here("outdata"), output_dir = here::here("outdata")) {
  
  message("Creating summary tables...")
  
  # Load basic statistics
  basic_stats_path <- file.path(outdata_dir, "basic_emf_stats.rds")
  if (!file.exists(basic_stats_path)) {
    stop("Basic statistics not found. Run fungal data processing first.")
  }
  
  basic_stats <- readRDS(basic_stats_path)
  
  # Load comprehensive EMF analysis if available
  comprehensive_path <- file.path(outdata_dir, "comprehensive_emf_analysis.rds")
  if (file.exists(comprehensive_path)) {
    comprehensive_results <- readRDS(comprehensive_path)
    comprehensive_stats <- comprehensive_results$unique_counts
  } else {
    # Try to calculate comprehensive analysis if not already done
    sample_points_path <- file.path(outdata_dir, "emf_samples_canada.gpkg")
    if (file.exists(sample_points_path)) {
      comprehensive_results <- calculate_comprehensive_emf_analysis(sample_points_path, outdata_dir)
      comprehensive_stats <- comprehensive_results$unique_counts
    } else {
      comprehensive_stats <- list()
    }
  }
  
  # Load spatial overlap statistics
  spatial_summary_path <- file.path(outdata_dir, "metrics_summary.rds")
  if (file.exists(spatial_summary_path)) {
    metrics_summary <- readRDS(spatial_summary_path)
    overlap_stats <- metrics_summary$ecoregion_overlap$stats
    coverage_stats <- metrics_summary$sampling_coverage$stats
  } else {
    overlap_stats <- list()
    coverage_stats <- list()
  }
  
  # Create comprehensive summary table
  summary_data <- data.frame(
    Metric = c(
      "Total EMF sequence records",
      "Total unique sampling locations (all data)",
      "Unique locations with EMF data",
      "Unique Other_ID values",
      "Unique UNITE_ID values", 
      "Unique EMF genera",
      "Unique EMF species",
      "Canadian EMF host tree species",
      "Host species with sequence data",
      "Percentage of host species with data",
      "Unique host genera with data",
      "Percentage of host genera with data",
      "Total ecoregions",
      "Ecoregions with potential EMF habitat",
      "Ecoregions sampled for EMF",
      "Percentage of habitat ecoregions unsampled"
    ),
    Value = c(
      ifelse("total_sequence_records" %in% names(basic_stats), 
             basic_stats$total_sequence_records, "N/A"),
      ifelse("total_sampling_locations_all_data" %in% names(basic_stats), 
             basic_stats$total_sampling_locations_all_data, "N/A"),
      ifelse("unique_locations_with_emf_data" %in% names(basic_stats), 
             basic_stats$unique_locations_with_emf_data, "N/A"),
      ifelse("unique_other_id" %in% names(comprehensive_stats), 
             comprehensive_stats$unique_other_id, "N/A"),
      ifelse("unique_unite_id" %in% names(comprehensive_stats), 
             comprehensive_stats$unique_unite_id, "N/A"),
      ifelse("unique_em_genera" %in% names(comprehensive_stats), 
             comprehensive_stats$unique_em_genera, "N/A"),
      ifelse("unique_em_species" %in% names(comprehensive_stats), 
             comprehensive_stats$unique_em_species, "N/A"),
      ifelse("emf_host_species_count" %in% names(basic_stats), 
             basic_stats$emf_host_species_count, "N/A"),
      ifelse("native_hosts_with_sequences" %in% names(basic_stats), 
             basic_stats$native_hosts_with_sequences, "N/A"),
      ifelse("percent_hosts_with_data" %in% names(basic_stats), 
             paste0(basic_stats$percent_hosts_with_data, "%"), "N/A"),
      ifelse("unique_host_genera" %in% names(basic_stats), 
             basic_stats$unique_host_genera, "N/A"),
      ifelse("percent_genera_with_data" %in% names(basic_stats), 
             paste0(basic_stats$percent_genera_with_data, "%"), "N/A"),
      ifelse("total_ecoregions" %in% names(overlap_stats), 
             overlap_stats$total_ecoregions, "N/A"),
      ifelse("ecoregions_with_habitat" %in% names(overlap_stats), 
             overlap_stats$ecoregions_with_habitat, "N/A"),
      ifelse("habitat_ecoregions_sampled" %in% names(coverage_stats), 
             coverage_stats$habitat_ecoregions_sampled, "N/A"),
      ifelse("percent_habitat_ecoregions_unsampled" %in% names(coverage_stats), 
             paste0(coverage_stats$percent_habitat_ecoregions_unsampled, "%"), "N/A")
    ),
    stringsAsFactors = FALSE
  )
  
  # Save summary table
  write.csv(summary_data, file.path(output_dir, "emf_summary_statistics.csv"), row.names = FALSE)
  
  # Create enhanced knowledge gaps table
  knowledge_gaps <- data.frame(
    Shortfall = c("Eltonian Shortfall", "Hutchinsonian Shortfall", "Wallacean Shortfall"),
    Description = c(
      "Gap in knowledge about species interactions",
      "Gap in knowledge about species' abiotic niches",
      "Gap in knowledge about species distributions"
    ),
    Key_Finding = c(
      paste0(basic_stats$percent_hosts_with_data, "% of Canadian EMF host tree species have associated fungal sequence data; ",
             ifelse("unique_em_species" %in% names(comprehensive_stats), 
                    paste(comprehensive_stats$unique_em_species, "unique EMF species identified"), "")),
      paste0(coverage_stats$percent_habitat_ecoregions_unsampled, "% of habitat ecoregions remain unsampled"),
      paste0(ifelse("total_sampling_locations_all_data" %in% names(basic_stats), 
                    paste(basic_stats$total_sampling_locations_all_data, "total unique sampling locations; "), ""),
             ifelse("unique_locations_with_emf_data" %in% names(basic_stats), 
                    paste(basic_stats$unique_locations_with_emf_data, "locations with EMF sequence data; "), ""),
             ifelse("unique_other_id" %in% names(comprehensive_stats),
                    paste(comprehensive_stats$unique_other_id, "unique Other_ID and", 
                          comprehensive_stats$unique_unite_id, "unique UNITE_ID sequence-based taxa"), ""))
    ),
    stringsAsFactors = FALSE
  )
  
  # Save knowledge gaps table
  write.csv(knowledge_gaps, file.path(output_dir, "knowledge_gaps_summary.csv"), row.names = FALSE)
  
  # Create proportion statistics table if available
  proportion_table <- NULL
  proportion_stats_path <- file.path(outdata_dir, "proportion_stats.rds")
  if (file.exists(proportion_stats_path)) {
    proportion_stats <- readRDS(proportion_stats_path)
    
    proportion_table <- data.frame(
      Statistic = c(
        "Cells with data",
        "Cells with zero EMF proportion",
        "Cells with non-zero EMF proportion",
        "Mean EMF proportion",
        "Maximum EMF proportion",
        "Minimum EMF proportion"
      ),
      Value = c(
        proportion_stats$cells_with_data,
        proportion_stats$cells_zero_proportion,
        proportion_stats$cells_nonzero_proportion,
        proportion_stats$mean_proportion,
        proportion_stats$max_proportion,
        proportion_stats$min_proportion
      ),
      stringsAsFactors = FALSE
    )
    
    write.csv(proportion_table, file.path(output_dir, "proportion_statistics.csv"), row.names = FALSE)
  }
  
  message("✓ All summary tables created")
  
  # Return the tables for use in RMD
  result <- list(
    summary_statistics = summary_data,
    knowledge_gaps = knowledge_gaps
  )
  
  if (!is.null(proportion_table)) {
    result$proportion_statistics <- proportion_table
  }
  
  if (length(comprehensive_stats) > 0) {
    result$comprehensive_stats = comprehensive_stats
  }
  
  return(result)
}

#' Clean up temporary files and objects
#' 
#' @param outdata_dir Character. Directory containing data files
#' @param keep_essential Logical. Whether to keep essential files
clean_temp_files <- function(outdata_dir = here::here("outdata"), keep_essential = TRUE) {
  
  message("Cleaning up temporary files...")
  
  # Define path to temporary directory
  temp_dir <- file.path(outdata_dir, "temp")
  
  files_removed <- 0
  
  # Remove entire temp directory if it exists
  if (dir.exists(temp_dir)) {
    # Count files before removal
    temp_files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
    files_removed <- length(temp_files)
    
    # Remove the entire temp directory
    unlink(temp_dir, recursive = TRUE)
    message(paste("✓ Removed temporary directory:", temp_dir))
  }
  
  # Also clean up any old temporary files that might exist in the main directory
  # (for backward compatibility)
  old_temp_patterns <- c(
    "*_raw.gpkg",    # Raw boundary files
    "*.zip",         # Downloaded zip files
    "temp_*",        # Any temporary files
    "*checkpoint*"   # Checkpoint files
  )
  
  for (pattern in old_temp_patterns) {
    old_temp_files <- list.files(outdata_dir, pattern = glob2rx(pattern), 
                                 full.names = TRUE, recursive = FALSE)
    
    if (length(old_temp_files) > 0) {
      file.remove(old_temp_files)
      files_removed <- files_removed + length(old_temp_files)
    }
  }
  
  # Clean up R environment
  # Remove large objects but keep essential ones
  if (keep_essential) {
    essential_objects <- c("crdref", "canada_albers")
    all_objects <- ls(envir = .GlobalEnv)
    objects_to_remove <- setdiff(all_objects, essential_objects)
    
    if (length(objects_to_remove) > 0) {
      rm(list = objects_to_remove, envir = .GlobalEnv)
    }
  }
  
  # Force garbage collection
  gc()
  
  message(paste("✓ Cleanup complete -", files_removed, "temporary files removed"))
  
  return(invisible(files_removed))
}

#' Generate analysis report summary
#' 
#' @param outdata_dir Character. Directory containing data files
#' @return Character. Summary text for report
generate_analysis_summary <- function(outdata_dir = here::here("outdata")) {
  
  # Load basic statistics
  basic_stats <- readRDS(file.path(outdata_dir, "basic_emf_stats.rds"))
  
  # Load comprehensive analysis if available
  comprehensive_path <- file.path(outdata_dir, "comprehensive_emf_analysis.rds")
  if (file.exists(comprehensive_path)) {
    comprehensive_results <- readRDS(comprehensive_path)
    comprehensive_stats <- comprehensive_results$unique_counts
  } else {
    comprehensive_stats <- list()
  }
  
  # Load metrics if available
  metrics_path <- file.path(outdata_dir, "metrics_summary.rds")
  if (file.exists(metrics_path)) {
    metrics_summary <- readRDS(metrics_path)
    overlap_stats <- metrics_summary$ecoregion_overlap$stats
    coverage_stats <- metrics_summary$sampling_coverage$stats
  } else {
    overlap_stats <- coverage_stats <- list()
  }
  
  # Generate summary text
  summary_text <- paste0(
    "## Analysis Summary\n\n",
    "This comprehensive analysis reveals key shortfalls:\n\n",
    "**Eltonian Shortfall (Species Interactions):**\n",
    "- ", basic_stats$percent_hosts_with_data, "% of Canadian EMF host tree species have associated fungal sequence data\n",
    "- ", ifelse("unique_other_id" %in% names(comprehensive_stats), 
                 comprehensive_stats$unique_other_id, "N/A"), " unique Other_ID sequence-based taxa\n",
    "- ", ifelse("unique_unite_id" %in% names(comprehensive_stats), 
                 comprehensive_stats$unique_unite_id, "N/A"), " unique UNITE_ID sequence-based taxa\n",
    "- ", ifelse("unique_em_genera" %in% names(comprehensive_stats), 
                 comprehensive_stats$unique_em_genera, "N/A"), " unique EMF genera\n",
    "- ", ifelse("unique_em_species" %in% names(comprehensive_stats), 
                 comprehensive_stats$unique_em_species, "N/A"), " unique EMF species\n",
    "- Host-taxon association matrices created for all taxonomic levels\n\n",
    
    "**Hutchinsonian Shortfall (Environmental Niches):**\n",
    if (length(coverage_stats) > 0) {
      paste0("- ", coverage_stats$percent_habitat_ecoregions_unsampled, "% of habitat ecoregions remain unsampled\n",
             "- ", coverage_stats$habitat_ecoregions_sampled, " of ", coverage_stats$habitat_ecoregions_total, " potential habitat ecoregions have been sampled\n\n")
    } else {
      "- Sampling coverage analysis pending\n\n"
    },
    
    "**Wallacean Shortfall (Species Distributions):**\n",
    "- ", ifelse("total_sampling_locations_all_data" %in% names(basic_stats), 
                 basic_stats$total_sampling_locations_all_data, "N/A"), " total unique sampling locations across Canada\n",
    "- ", ifelse("unique_locations_with_emf_data" %in% names(basic_stats), 
                 basic_stats$unique_locations_with_emf_data, "N/A"), " locations with detailed EMF sequence data\n",
    "- Sampling intensity calculated for all taxonomic levels (Other_ID, UNITE_ID, genera, species)\n",
    "- Location-per-taxon statistics provide insights into sampling gaps\n\n"
  )
  
  return(summary_text)
}

#' Save workspace state
#' 
#' @param outdata_dir Character. Directory to save workspace
#' @param workspace_name Character. Name for workspace file
save_workspace_state <- function(outdata_dir = here::here("outdata"), workspace_name = "emf_analysis_workspace.RData") {
  
  message("Saving workspace state...")
  
  # Define essential objects to keep
  essential_objects <- c(
    "crdref", "canada_albers"  # Keep coordinate systems
  )
  
  # Get current objects
  current_objects <- ls(envir = .GlobalEnv)
  objects_to_save <- intersect(current_objects, essential_objects)
  
  # Save workspace
  workspace_path <- file.path(outdata_dir, workspace_name)
  save(list = objects_to_save, file = workspace_path, envir = .GlobalEnv)
  
  message(paste("✓ Workspace saved:", workspace_path))
  
  return(workspace_path)
}

#' Load workspace state
#' 
#' @param outdata_dir Character. Directory containing workspace
#' @param workspace_name Character. Name of workspace file
load_workspace_state <- function(outdata_dir = here::here("outdata"), workspace_name = "emf_analysis_workspace.RData") {
  
  workspace_path <- file.path(outdata_dir, workspace_name)
  
  if (!file.exists(workspace_path)) {
    warning(paste("Workspace file not found:", workspace_path))
    return(FALSE)
  }
  
  message("Loading workspace state...")
  load(workspace_path, envir = .GlobalEnv)
  message("✓ Workspace loaded")
  
  return(TRUE)
}

#' Get file size in human readable format
#' 
#' @param filepath Character. Path to file
#' @return Character. File size with units
get_file_size <- function(filepath) {
  
  if (!file.exists(filepath)) {
    return("File not found")
  }
  
  size_bytes <- file.size(filepath)
  
  if (size_bytes < 1024) {
    return(paste(size_bytes, "B"))
  } else if (size_bytes < 1024^2) {
    return(paste(round(size_bytes / 1024, 1), "KB"))
  } else if (size_bytes < 1024^3) {
    return(paste(round(size_bytes / 1024^2, 1), "MB"))
  } else {
    return(paste(round(size_bytes / 1024^3, 1), "GB"))
  }
}

#' Create file inventory
#' 
#' @param outdata_dir Character. Directory to inventory
#' @return Data frame with file information
create_file_inventory <- function(outdata_dir = here::here("outdata")) {
  
  message("Creating file inventory...")
  
  # Get all files in outdata directory (excluding temp directory)
  all_files <- list.files(outdata_dir, recursive = TRUE, full.names = TRUE)
  
  # Exclude files in temp directory from inventory
  temp_dir <- file.path(outdata_dir, "temp")
  all_files <- all_files[!grepl(paste0("^", temp_dir), all_files)]
  
  # Create inventory
  inventory <- data.frame(
    File = basename(all_files),
    Path = dirname(all_files),
    Size = sapply(all_files, get_file_size),
    Modified = file.mtime(all_files),
    stringsAsFactors = FALSE
  )
  
  # Sort by modification time (newest first)
  inventory <- inventory[order(inventory$Modified, decreasing = TRUE), ]
  
  # Save inventory
  write.csv(inventory, file.path(outdata_dir, "file_inventory.csv"), row.names = FALSE)
  
  message("✓ File inventory created")
  
  return(inventory)
}

#' Format table for kable display
#' 
#' @param table_data Data frame to format
#' @param caption Character. Table caption
#' @param col_names Character vector. Column names
#' @return kable object
format_kable_table <- function(table_data, caption = "", col_names = NULL) {
  
  if (!require(kableExtra, quietly = TRUE)) {
    # Fall back to basic kable if kableExtra not available
    return(knitr::kable(table_data, caption = caption, col.names = col_names))
  }
  
  table_data %>%
    knitr::kable(caption = caption, col.names = col_names, format = "html") %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE,
      position = "left"
    ) %>%
    kableExtra::column_spec(1, bold = TRUE) %>%
    kableExtra::row_spec(0, bold = TRUE, background = "#f0f0f0")
}

#' Create interactive table with DT
#' 
#' @param table_data Data frame to display
#' @param caption Character. Table caption
#' @return DT datatable object
create_interactive_table <- function(table_data, caption = "") {
  
  if (!require(DT, quietly = TRUE)) {
    warning("DT package not available, returning basic table")
    return(table_data)
  }
  
  DT::datatable(
    table_data,
    caption = caption,
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'))
    ),
    class = 'cell-border stripe'
  )
}

#' Check and validate required analysis files
#' 
#' @param outdata_dir Character. Directory containing data files
#' @return Logical. TRUE if all files exist
validate_analysis_files <- function(outdata_dir = here::here("outdata")) {
  
  required_files <- c(
    # Processed spatial files
    file.path(outdata_dir, "canada_simple.gpkg"),
    file.path(outdata_dir, "tree_ranges_canada.gpkg"),
    file.path(outdata_dir, "ecoregions_processed.gpkg"),
    
    # Processed fungal files
    file.path(outdata_dir, "emf_host_tree_ranges.gpkg"),
    file.path(outdata_dir, "emf_samples_canada.gpkg"),
    file.path(outdata_dir, "emf_samples_aggregated_1000m.gpkg"),
    
    # Statistics files
    file.path(outdata_dir, "basic_emf_stats.rds")
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  
  if (length(missing_files) > 0) {
    warning("Missing required analysis files:")
    for (file in missing_files) {
      warning(paste(" -", file))
    }
    return(FALSE)
  }
  
  message("✓ All required analysis files found")
  return(TRUE)
}