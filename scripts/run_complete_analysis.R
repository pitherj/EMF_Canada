# =============================================================================
# EMF Canada Analysis - Master Script
# =============================================================================
# Purpose: Run complete analysis pipeline
# Author: Jason Pither, Clara Qin, Isaac Eckert
# Usage: source("run_complete_analysis.R")
# =============================================================================

# Clear environment
rm(list = ls())

# Load required libraries for the master script
if (!require(here, quietly = TRUE)) {
  stop("Package 'here' is required but not installed.")
}

# Set up paths
project_root <- here::here()
rawdata_dir <- file.path(project_root, "rawdata")
outdata_dir <- file.path(project_root, "outdata")
scripts_dir <- file.path(project_root, "scripts")
figures_dir <- file.path(project_root, "figures")

cat("=================================================================\n")
cat("EMF CANADA ANALYSIS - COMPLETE PIPELINE\n")
cat("=================================================================\n")
cat("Project directory:", project_root, "\n")
cat("Analysis started at:", as.character(Sys.time()), "\n")
cat("=================================================================\n\n")

# Track processing time
start_time <- Sys.time()

# =============================================================================
# STEP 1: SETUP AND INITIALIZATION
# =============================================================================

cat("STEP 1: SETUP AND INITIALIZATION\n")
cat("---------------------------------\n")

setup_success <- FALSE
tryCatch({
  source(file.path(scripts_dir, "01_setup.R"))
  
  # Load utility functions early for temp directory management
  source(file.path(scripts_dir, "utils.R"))
  
  # Create temporary directory for analysis
  temp_dir <- create_temp_dir(outdata_dir)
  cat("Temporary directory created at:", temp_dir, "\n")
  
  cat("✓ Environment Setup completed successfully\n\n")
  setup_success <- TRUE
}, error = function(e) {
  cat("✗ Error in Environment Setup:", e$message, "\n\n")
  setup_success <<- FALSE
})

if (!setup_success) {
  stop("Setup failed. Cannot continue with analysis.")
}

# =============================================================================
# STEP 2: DATA DOWNLOAD
# =============================================================================

cat("STEP 2: DATA DOWNLOAD\n")
cat("---------------------\n")

download_success <- FALSE
download_summary <- NULL
tryCatch({
  source(file.path(scripts_dir, "02_download_data.R"))
  download_summary <- download_all_data(rawdata_dir, outdata_dir)
  cat("✓ Data Download completed successfully\n\n")
  download_success <- TRUE
}, error = function(e) {
  cat("✗ Error in Data Download:", e$message, "\n\n")
  download_success <<- FALSE
})

if (!download_success) {
  warning("Data download failed. Check internet connection and try again.")
}

# =============================================================================
# STEP 3: SPATIAL DATA PROCESSING
# =============================================================================

cat("STEP 3: SPATIAL DATA PROCESSING\n")
cat("--------------------------------\n")

spatial_success <- FALSE
spatial_summary <- NULL
tryCatch({
  source(file.path(scripts_dir, "03_process_spatial.R"))
  spatial_summary <- process_all_spatial(rawdata_dir, outdata_dir, download_summary)
  cat("✓ Spatial Data Processing completed successfully\n\n")
  spatial_success <- TRUE
}, error = function(e) {
  cat("✗ Error in Spatial Data Processing:", e$message, "\n\n")
  spatial_success <<- FALSE
})

if (!spatial_success) {
  stop("Spatial processing failed. Cannot continue with analysis.")
}

# =============================================================================
# STEP 4: FUNGAL DATA PROCESSING
# =============================================================================

cat("STEP 4: FUNGAL DATA PROCESSING\n")
cat("-------------------------------\n")

fungal_success <- FALSE
tryCatch({
  source(file.path(scripts_dir, "04_process_fungal.R"))
  fungal_summary <- process_all_fungal(rawdata_dir, outdata_dir)
  cat("✓ Fungal Data Processing completed successfully\n\n")
  fungal_success <- TRUE
}, error = function(e) {
  cat("✗ Error in Fungal Data Processing:", e$message, "\n\n")
  fungal_success <<- FALSE
})

if (!fungal_success) {
  stop("Fungal data processing failed. Cannot continue with analysis.")
}

# =============================================================================
# STEP 5: METRICS CALCULATION
# =============================================================================

cat("STEP 5: METRICS CALCULATION\n")
cat("----------------------------\n")

metrics_success <- FALSE
tryCatch({
  source(file.path(scripts_dir, "05_calculate_metrics.R"))
  metrics_summary <- calculate_all_metrics(outdata_dir)
  cat("✓ Metrics Calculation completed successfully\n\n")
  metrics_success <- TRUE
}, error = function(e) {
  cat("✗ Error in Metrics Calculation:", e$message, "\n\n")
  metrics_success <<- FALSE
})

if (!metrics_success) {
  warning("Metrics calculation failed. Some analyses may be incomplete.")
}

# =============================================================================
# STEP 6: MAP CREATION
# =============================================================================

cat("STEP 6: MAP CREATION\n")
cat("--------------------\n")

maps_success <- FALSE
tryCatch({
  source(file.path(scripts_dir, "06_create_maps.R"))
  maps_summary <- create_all_maps(outdata_dir, figures_dir)
  cat("✓ Map Creation completed successfully\n\n")
  maps_success <- TRUE
}, error = function(e) {
  cat("✗ Error in Map Creation:", e$message, "\n\n")
  maps_success <<- FALSE
})

if (!maps_success) {
  warning("Map creation failed. Visualizations may be incomplete.")
}

# =============================================================================
# STEP 7: GENERATE SUMMARY TABLES
# =============================================================================

cat("STEP 7: GENERATE SUMMARY TABLES\n")
cat("--------------------------------\n")

cat("Creating summary tables...\n")

tables_success <- FALSE
tryCatch({
  tables_summary <- create_summary_tables(outdata_dir, outdata_dir)
  cat("✓ Summary tables created successfully\n\n")
  tables_success <- TRUE
}, error = function(e) {
  cat("✗ Error creating summary tables:", e$message, "\n\n")
  tables_success <<- FALSE
})

# =============================================================================
# STEP 8: CLEANUP TEMPORARY FILES
# =============================================================================

cat("STEP 8: CLEANUP TEMPORARY FILES\n")
cat("--------------------------------\n")

cleanup_success <- FALSE
tryCatch({
  files_removed <- clean_temp_files(outdata_dir = outdata_dir)
  cat("✓ Temporary files cleaned up successfully\n\n")
  cleanup_success <- TRUE
}, error = function(e) {
  cat("✗ Error during cleanup:", e$message, "\n")
  cat("Note: Temporary files may need manual cleanup\n\n")
  cleanup_success <<- FALSE
})

# =============================================================================
# ANALYSIS COMPLETION SUMMARY
# =============================================================================

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("=================================================================\n")
cat("ANALYSIS PIPELINE COMPLETED\n")
cat("=================================================================\n")
cat("Total processing time:", round(as.numeric(total_time), 2), "minutes\n")
cat("Completed at:", as.character(end_time), "\n")

# Check which components completed successfully
components <- c("Setup", "Download", "Spatial", "Fungal", "Metrics", "Maps", "Tables", "Cleanup")
success_flags <- c(setup_success, download_success, spatial_success, 
                   fungal_success, metrics_success, maps_success, tables_success, cleanup_success)

cat("\nComponent Status:\n")
for (i in seq_along(components)) {
  status <- if (success_flags[i]) "✓ SUCCESS" else "✗ FAILED"
  cat(sprintf("%-15s: %s\n", components[i], status))
}

# Generate file summary
cat("\nGenerated Files:\n")
cat("- Raw data files:", rawdata_dir, "\n")
cat("- Processed data and tables:", outdata_dir, "\n")
cat("- Figures:", figures_dir, "\n")

# Check for key output files
key_files <- c(
  file.path(outdata_dir, "basic_emf_stats.rds"),
  file.path(outdata_dir, "comprehensive_emf_analysis.rds"),
  file.path(outdata_dir, "comprehensive_emf_summary.csv"),
  file.path(outdata_dir, "locations_per_other_id.csv"),
  file.path(outdata_dir, "locations_per_unite_id.csv"),
  file.path(outdata_dir, "locations_per_genus.csv"),
  file.path(outdata_dir, "locations_per_species.csv"),
  file.path(outdata_dir, "host_taxon_matrix_other_id.csv"),
  file.path(outdata_dir, "host_taxon_matrix_unite_id.csv"),
  file.path(outdata_dir, "host_taxon_matrix_genus.csv"),
  file.path(outdata_dir, "host_taxon_matrix_species.csv"),
  file.path(figures_dir, "canada_emf_bivariate_map.png"),
  file.path(outdata_dir, "emf_summary_statistics.csv")
)

cat("\nKey Output Files:\n")
for (file in key_files) {
  status <- if (file.exists(file)) "✓ EXISTS" else "✗ MISSING"
  cat(sprintf("%-60s: %s\n", basename(file), status))
}

# =============================================================================
# NEXT STEPS
# =============================================================================

cat("\n=================================================================\n")
cat("NEXT STEPS\n")
cat("=================================================================\n")
cat("1. Review generated summary tables in 'outdata/' directory\n")
cat("2. Examine maps and figures in 'figures/' directory\n")
cat("3. Review comprehensive EMF taxon analysis outputs:\n")
cat("   - comprehensive_emf_summary.csv (summary across all taxonomic levels)\n")
cat("   - locations_per_*.csv (sampling locations per taxon for each level)\n")
cat("   - host_taxon_matrix_*.csv (host-taxon association matrices)\n")
cat("   - emf_taxa_per_host_*.csv (EMF taxa counts per host species)\n")
cat("   - hosts_per_emf_taxon_*.csv (host species counts per EMF taxon)\n")
cat("4. Knit the RMD file to generate the final report:\n")
cat("   rmarkdown::render('EMF_Analysis_Clean.Rmd')\n")
cat("5. Review any error messages above and re-run failed components\n")

if (all(success_flags)) {
  cat("\n🎉 ANALYSIS COMPLETED SUCCESSFULLY! 🎉\n")
  cat("All components processed without errors.\n")
} else {
  failed_components <- components[!success_flags]
  cat("\n⚠️  ANALYSIS COMPLETED WITH WARNINGS ⚠️\n")
  cat("Failed components:", paste(failed_components, collapse = ", "), "\n")
  cat("Check error messages above and re-run failed scripts individually.\n")
}

cat("=================================================================\n")

# Save final workspace
final_workspace_path <- file.path(outdata_dir, "final_analysis_workspace.RData")
tryCatch({
  save_workspace_state(outdata_dir, "final_analysis_workspace.RData")
  cat("Final workspace saved to:", final_workspace_path, "\n")
}, error = function(e) {
  cat("Note: Could not save workspace:", e$message, "\n")
})

cat("\nAnalysis pipeline complete!\n")

# Print final summary of comprehensive EMF analysis if available
comprehensive_path <- file.path(outdata_dir, "comprehensive_emf_analysis.rds")
if (file.exists(comprehensive_path)) {
  cat("\n=================================================================\n")
  cat("COMPREHENSIVE EMF TAXON ANALYSIS SUMMARY\n")
  cat("=================================================================\n")
  
  comprehensive_results <- readRDS(comprehensive_path)
  unique_counts <- comprehensive_results$unique_counts
  
  # Load basic stats for sampling coverage
  basic_stats_path <- file.path(outdata_dir, "basic_emf_stats.rds")
  if (file.exists(basic_stats_path)) {
    basic_stats <- readRDS(basic_stats_path)
  } else {
    basic_stats <- list()
  }
  
  cat("TAXONOMIC DIVERSITY:\n")
  cat("Unique Other_ID values:", unique_counts$unique_other_id, "\n")
  cat("Unique UNITE_ID values:", unique_counts$unique_unite_id, "\n")
  cat("Unique EMF genera:", unique_counts$unique_em_genera, "\n")
  cat("Unique EMF species:", unique_counts$unique_em_species, "\n")
  cat("Total EMF records:", unique_counts$total_records, "\n")
  cat("Unique sampling locations:", unique_counts$unique_locations, "\n")
  
  if (!is.null(comprehensive_results$host_associations)) {
    cat("\nHOST-TAXON ASSOCIATIONS:\n")
    cat("Host species with EMF data:", 
        comprehensive_results$host_associations$species$summary$total_host_species, "\n")
    cat("Mean EMF species per host:", 
        comprehensive_results$host_associations$species$summary$mean_taxa_per_host, "\n")
    cat("Mean hosts per EMF species:", 
        comprehensive_results$host_associations$species$summary$mean_hosts_per_taxon, "\n")
  }
  
  cat("\nSAMPLING INTENSITY (mean locations per taxon):\n")
  if (!is.null(comprehensive_results$sampling_intensity$other_id)) {
    cat("Other_ID level:", comprehensive_results$sampling_intensity$other_id$mean_locations, "\n")
  }
  if (!is.null(comprehensive_results$sampling_intensity$unite_id)) {
    cat("UNITE_ID level:", comprehensive_results$sampling_intensity$unite_id$mean_locations, "\n")
  }
  if (!is.null(comprehensive_results$sampling_intensity$genus)) {
    cat("Genus level:", comprehensive_results$sampling_intensity$genus$mean_locations, "\n")
  }
  if (!is.null(comprehensive_results$sampling_intensity$species)) {
    cat("Species level:", comprehensive_results$sampling_intensity$species$mean_locations, "\n")
  }
  
  cat("=================================================================\n")
}