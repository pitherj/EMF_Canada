# =============================================================================
# EMF Canada Analysis - Mapping Script
# =============================================================================
# Purpose: Create all maps and visualizations
# Author: Jason Pither, Clara Qin, Isaac Eckert
# =============================================================================

# Required packages
required_packages <- c("here", "terra", "sf", "cowplot", "ggplot2", "dplyr", "tidyterra", "patchwork")
options(warn = FALSE)
# Load packages with auto-installation for missing ones
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    install.packages(pkg)
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "could not be installed."))
    }
  }
}

# Load coordinate systems and spatial functions
if (!exists("canada_albers")) {
  source(here::here("scripts", "01_setup.R"))
}
if (!exists("safe_load_terra")) {
  source(here::here("scripts", "03_process_spatial.R"))
}

#' Create base map of Canada with ecoregions
#' 
#' @param outdata_dir Character. Directory containing data files
#' @return ggplot object of base map
create_base_map <- function(outdata_dir = here::here("outdata")) {
  
  message("Creating base map...")
  
  # Load spatial data
  canada_albers <- safe_load_terra(file.path(outdata_dir, "canada_ne_albers.gpkg"))
  ecoregions_processed <- safe_load_terra(file.path(outdata_dir, "ecoregions_processed.gpkg"))
  lakes_canada <- safe_load_terra(file.path(outdata_dir, "lakes_canada_albers.gpkg"))
  
  # Project ecoregions to Albers and clip to Canada
  ecoregions_albers <- terra::project(ecoregions_processed, canada_albers)
  ecoregions_clipped <- terra::intersect(ecoregions_albers, canada_albers)
  
  # Load ecozone names if available
  ecozone_names_path <- file.path(outdata_dir, "ecozone_names.csv")
  if (file.exists(ecozone_names_path)) {
    ecozone_names <- read.csv(ecozone_names_path)
    
    # Join with ecoregions
    ecoregions_with_names <- tidyterra::left_join(
      ecoregions_clipped,
      ecozone_names,
      by = "ECOZONE"
    )
  } else {
    ecoregions_with_names <- ecoregions_clipped
    ecoregions_with_names$NAME_EN <- paste("Ecozone", ecoregions_with_names$ECOZONE)
  }
  
  # Create distinct colors for ecozones
  unique_ecozones <- sort(unique(ecoregions_with_names$NAME_EN))
  n_ecozones <- length(unique_ecozones)
  
  # Create color palette
  set.seed(123)
  ecozone_colors <- rainbow(n_ecozones, start = 0, end = 0.85)
  names(ecozone_colors) <- unique_ecozones
  
  # Convert to sf for ggplot and ensure all data is clipped to Canada
  canada_sf <- sf::st_as_sf(canada_albers)
  ecoregions_sf <- sf::st_as_sf(ecoregions_with_names)
  lakes_sf <- sf::st_as_sf(lakes_canada)
  
  # Clip all layers to Canada boundary to ensure clean edges
  ecoregions_sf <- sf::st_intersection(ecoregions_sf, canada_sf)
  lakes_sf <- sf::st_intersection(lakes_sf, canada_sf)
  
  # Create base map
  base_map <- ggplot() +
    # Ecoregions (clipped to Canada)
    geom_sf(data = ecoregions_sf,
            aes(fill = NAME_EN),
            color = "grey30",
            linewidth = 0.2,
            alpha = 0.3) +
    
    # Lakes (clipped to Canada)
    geom_sf(data = lakes_sf,
            fill = "dodgerblue",
            color = "darkblue",
            linewidth = 0.2,
            alpha = 0.9) +
    
    # Canada boundary (defines map extent)
    geom_sf(data = canada_sf,
            fill = NA,
            color = "black",
            linewidth = 0.5) +
    
    # Color scale
    scale_fill_manual(values = ecozone_colors,
                      name = "Ecozones") +
    
    # Set map limits to Canada extent
    coord_sf(xlim = c(sf::st_bbox(canada_sf)[1], sf::st_bbox(canada_sf)[3]),
             ylim = c(sf::st_bbox(canada_sf)[2], sf::st_bbox(canada_sf)[4]),
             expand = FALSE) +
    
    # Theme
    theme_void() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.background = element_rect(fill = "white", color = NA),  # no border
      plot.background = element_rect(fill = "white", color = NA)    # no outer frame
    ) +
    
    
    # Labels
    labs(title = "Sampling locations and ecoregions")
  
  message("✓ Base map created")
  
  return(base_map)
}

#' Create map with sampling locations
#' 
#' @param data_dir Character. Directory containing data files
#' @param figures_dir Character. Directory to save figures
#' @return ggplot object of sampling map
create_sampling_map <- function(data_dir = here::here("data"), figures_dir = here::here("figures")) {
  
  message("Creating sampling locations map...")
  
  # Get base map (without default legend position)
  base_map <- create_base_map(data_dir) +
    theme(legend.position = "none")  # disable default legends
  
  # Load aggregated sample points
  aggregated_samples <- safe_load_terra(file.path(data_dir, "emf_samples_aggregated_1000m.gpkg"))
  
  # Project to Albers
  canada_albers_crs <- terra::crs(safe_load_terra(file.path(data_dir, "canada_ne_albers.gpkg")))
  aggregated_samples_albers <- terra::project(aggregated_samples, canada_albers_crs)
  
  # Convert to sf
  samples_sf <- sf::st_as_sf(aggregated_samples_albers)
  
  # Count samples by dataset
  dataset_counts <- table(samples_sf$dataset)
  
  # Create subtitle with counts
  #  subtitle_text <- paste0("Samples unique to Stephanie: ", dataset_counts["Stephanie"], 
  #                         "; Samples unique to van Galen: ", dataset_counts["vanGalen"], 
  #                          "; Samples common to both: ", dataset_counts["Both"])
  
  # Main map without legend
  map_body <- base_map +
    geom_sf(data = samples_sf,
            aes(color = dataset),
            size = 2,
            alpha = 0.7) +
    scale_color_manual(values = c("Both" = "black",
                                  "Stephanie" = "red",
                                  "vanGalen" = "darkgrey"),
                       name = "Sample source") +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(
        override.aes = list(size = 4),
        order = 2
      )
    ) # +
  # labs(subtitle = subtitle_text)
  
  # Extract legend using cowplot (only the legend, not the whole plot)
  legend_only <- get_legend(
    map_body +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8)
      )
  )
  
  # Combine map and legend using inset_element
  sampling_map <- map_body +
    patchwork::inset_element(
      legend_only,
      left = 0.70, bottom = 0.53,
      right = 0.97, top = 0.97
    )
  
  
  # Save map
  ggsave(file.path(figures_dir, "canada_sampling_locations_map.png"),
         plot = sampling_map, width = 12, height = 10, dpi = 300, bg = "white")
  
  message("✓ Sampling locations map created and saved")
  
  return(sampling_map)
}


#' Create EMF richness map
#' 
#' @param outdata_dir Character. Directory containing data files
#' @param figures_dir Character. Directory to save figures
#' @param resolution_deg Numeric. Resolution used for raster
#' @return ggplot object of EMF richness map
create_emf_richness_map <- function(outdata_dir = here::here("outdata"), 
                                    figures_dir = here::here("figures"), 
                                    resolution_deg = 1) {
  
  message("Creating EMF richness map...")
  
  # Load EMF richness raster
  emf_richness_path <- file.path(outdata_dir, paste0("emf_richness_", resolution_deg, "deg.tif"))
  if (!file.exists(emf_richness_path)) {
    stop("EMF richness raster not found. Calculate metrics first.")
  }
  
  emf_richness_raster <- terra::rast(emf_richness_path)
  canada_boundary <- safe_load_terra(file.path(outdata_dir, "canada_simple.gpkg"))
  
  # Load Albers-projected Canada boundary (same as sampling map)
  canada_albers <- safe_load_terra(file.path(outdata_dir, "canada_ne_albers.gpkg"))
  canada_albers_crs <- terra::crs(canada_albers)
  
  # Ensure raster is properly clipped to Canada
  emf_richness_raster <- terra::mask(emf_richness_raster, canada_boundary)
  
  # Reproject raster to Albers
  emf_richness_albers <- terra::project(emf_richness_raster, canada_albers_crs, method = "near")
  
  # Convert reprojected raster to dataframe for ggplot
  emf_df <- as.data.frame(emf_richness_albers, xy = TRUE)
  colnames(emf_df) <- c("x", "y", "emf_richness")
  emf_df <- emf_df[!is.na(emf_df$emf_richness) & emf_df$emf_richness > 0, ]
  
  # Convert Albers boundary to sf
  canada_albers_sf <- sf::st_as_sf(canada_albers)
  
  # Max value for scaling
  max_emf <- max(emf_df$emf_richness, na.rm = TRUE)
  emf_colors <- c("white", "lightblue", "blue", "darkblue", "purple", "red")
  
  # Main map without legend
  emf_richness_mp <- ggplot() +
    geom_sf(data = canada_albers_sf, fill = "whitesmoke", color = "darkgrey", linewidth = 0.3) +
    geom_tile(data = emf_df, aes(x = x, y = y, fill = emf_richness), 
              color = "lightgrey", linewidth = 0.05) +
    geom_sf(data = canada_albers_sf, fill = NA, color = "black", linewidth = 0.2) +
    scale_fill_gradientn(colors = emf_colors,
                         name = "EMF Host\nSpeciesRichness",
                         breaks = 0:max_emf,
                         labels = 0:max_emf,
                         limits = c(0, max_emf),
                         guide = "none") +  # hide legend for this plot
    coord_sf(xlim = c(st_bbox(canada_albers_sf)[1], st_bbox(canada_albers_sf)[3]),
             ylim = c(st_bbox(canada_albers_sf)[2], st_bbox(canada_albers_sf)[4]),
             expand = FALSE) +
    theme_void() +
    # labs(title = "EMF Host Species Richness") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.background = element_rect(fill = "white", color = NA)  # no border
    )
  
  # Dummy plot to extract legend only
  library(ggplot2)
  
  rich_legend_plot <- ggplot(data = emf_df, aes(x = x, y = y, fill = emf_richness)) +
    geom_raster(alpha = 0) +  # ← makes data invisible but keeps legend
    scale_fill_gradientn(
      colors = emf_colors,
      name = "EMF Host\nSpecies Richness",
      breaks = c(0, max_emf),
      labels = c("0", as.character(max_emf)),
      limits = c(0, max_emf),
      guide = guide_colorbar(
        title.position = "right",
        title.vjust = 0.5,
        title.theme = element_text(angle = 0, hjust = 0.5),
        frame.colour = "black",
        ticks.colour = "black",
        barwidth = unit(0.35, "cm"),
        barheight = unit(3.6, "cm")
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  
  # Extract legend using patchwork
  legend_grob <- patchwork::wrap_elements(ggplotGrob(rich_legend_plot))
  
  # Inset the legend into top-right area
  emf_richness_map <- emf_richness_mp +
    inset_element(legend_grob,
                  left = 0.80, bottom = 0.78,
                  right = 0.97, top = 0.85)
  # Save map
  ggsave(file.path(figures_dir, "canada_emf_richness_map.png"),
         plot = emf_richness_map, bg = "white", width = 10, height = 8, dpi = 300)
  
  message("✓ EMF richness map created and saved")
  
  return(emf_richness_map)
}

#' Create EMF proportion map
#' 
#' @param data_dir Character. Directory containing data files
#' @param figures_dir Character. Directory to save figures
#' @param resolution_deg Numeric. Resolution used for raster
#' @return ggplot object of EMF proportion map
create_emf_proportion_map <- function(data_dir = here::here("data"), 
                                      figures_dir = here::here("figures"), 
                                      resolution_deg = 1) {
  
  message("Creating EMF proportion map...")
  
  # Load proportion raster
  proportion_path <- file.path(data_dir, paste0("emf_proportion_", resolution_deg, "deg.tif"))
  if (!file.exists(proportion_path)) {
    stop("EMF proportion raster not found. Calculate metrics first.")
  }
  
  proportion_raster <- terra::rast(proportion_path)
  canada_boundary <- safe_load_terra(file.path(data_dir, "canada_simple.gpkg"))
  
  # Load Albers-projected Canada boundary (same as sampling map)
  canada_albers <- safe_load_terra(file.path(data_dir, "canada_ne_albers.gpkg"))
  canada_albers_crs <- terra::crs(canada_albers)
  
  # Ensure raster is properly clipped to Canada
  proportion_raster <- terra::mask(proportion_raster, canada_boundary)
  
  # Reproject raster to Albers
  proportion_albers <- terra::project(proportion_raster, canada_albers_crs, method = "bilinear")
  
  # Convert reprojected raster to dataframe for ggplot
  prop_df <- as.data.frame(proportion_albers, xy = TRUE)
  colnames(prop_df) <- c("x", "y", "proportion")
  prop_df <- prop_df[!is.na(prop_df$proportion), ]
  
  # Convert Albers boundary to sf
  canada_albers_sf <- sf::st_as_sf(canada_albers)
  
  # Create color palette
  prop_colors <- c("white", "lightblue", "blue", "darkblue", "red")
  
  # Main map without legend
  proportion_mp <- ggplot() +
    geom_sf(data = canada_albers_sf, fill = "whitesmoke", color = "darkgrey", linewidth = 0.3) +
    geom_tile(data = prop_df, aes(x = x, y = y, fill = proportion), 
              color = "lightgrey", linewidth = 0.05) +
    geom_sf(data = canada_albers_sf, fill = NA, color = "black", linewidth = 0.2) +
    scale_fill_gradientn(colors = prop_colors,
                         name = "EMF Data\nCoverage (Proportion)",
                         breaks = c(0, 1),
                         labels = c("0", "1"),
                         limits = c(0, 1),
                         guide = "none") +
    coord_sf(xlim = c(st_bbox(canada_albers_sf)[1], st_bbox(canada_albers_sf)[3]),
             ylim = c(st_bbox(canada_albers_sf)[2], st_bbox(canada_albers_sf)[4]),
             expand = FALSE) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Legend-only plot
  prop_legend_plot <- ggplot(data = prop_df, aes(x = x, y = y, fill = proportion)) +
    geom_raster(alpha = 0) +
    scale_fill_gradientn(
      colors = prop_colors,
      name = "EMF Data\nCoverage (Proportion)",
      breaks = c(0, 1),
      labels = c("0", "1"),
      limits = c(0, 1),
      guide = guide_colorbar(
        title.position = "right",
        title.vjust = 0.5,
        title.theme = element_text(angle = 0, hjust = 0.5),
        frame.colour = "black",
        ticks.colour = "black",
        barwidth = unit(0.35, "cm"),
        barheight = unit(3.6, "cm")
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  # Extract and inset legend
  legend_grob <- patchwork::wrap_elements(ggplotGrob(prop_legend_plot))
  proportion_map <- proportion_mp +
    inset_element(legend_grob,
                  left = 0.80, bottom = 0.78,
                  right = 0.97, top = 0.85)
  
  # Save final map
  ggsave(file.path(figures_dir, "canada_emf_proportion_map.png"),
         plot = proportion_map, bg = "white", width = 10, height = 8, dpi = 300)
  
  message("✓ EMF proportion map created and saved")
  
  return(proportion_map)
}


#' Create bivariate map combining richness and proportion
#' 
#' @param data_dir Character. Directory containing data files
#' @param figures_dir Character. Directory to save figures
#' @param resolution_deg Numeric. Resolution used for raster
#' @return ggplot object of bivariate map
create_bivariate_map <- function(data_dir = here::here("data"), 
                                 figures_dir = here::here("figures"), 
                                 resolution_deg = 1) {
  
  message("Creating bivariate map...")
  
  # Load rasters
  richness_path <- file.path(data_dir, paste0("species_richness_", resolution_deg, "deg.tif"))
  proportion_path <- file.path(data_dir, paste0("emf_proportion_", resolution_deg, "deg.tif"))
  
  if (!file.exists(richness_path) || !file.exists(proportion_path)) {
    stop("Required rasters not found. Calculate metrics first.")
  }
  
  richness_raster <- terra::rast(richness_path)
  proportion_raster <- terra::rast(proportion_path)
  canada_boundary <- safe_load_terra(file.path(data_dir, "canada_simple.gpkg"))
  
  # Load Albers-projected Canada boundary (same as sampling map)
  canada_albers <- safe_load_terra(file.path(data_dir, "canada_ne_albers.gpkg"))
  canada_albers_crs <- terra::crs(canada_albers)
  
  # Ensure rasters are properly clipped to Canada
  richness_raster <- terra::mask(richness_raster, canada_boundary)
  proportion_raster <- terra::mask(proportion_raster, canada_boundary)
  
  # Reproject rasters to Albers
  richness_albers <- terra::project(richness_raster, canada_albers_crs, method = "near")
  proportion_albers <- terra::project(proportion_raster, canada_albers_crs, method = "bilinear")
  
  # Get values and create classification data
  richness_vals <- terra::values(richness_albers)
  proportion_vals <- terra::values(proportion_albers)
  
  # Create dataframe for classification
  bivar_data <- data.frame(
    richness = richness_vals[, 1],
    proportion = ifelse(is.na(proportion_vals[, 1]) & richness_vals[, 1] > 0, 0, proportion_vals[, 1])
  )
  
  # Filter to valid cells
  bivar_data_clean <- bivar_data[!is.na(bivar_data$richness) & bivar_data$richness > 0, ]
  
  # Create breaks for 3x3 classification
  max_richness <- max(bivar_data_clean$richness, na.rm = TRUE)
  richness_breaks <- c(0, max_richness/3, 2*max_richness/3, max_richness)
  proportion_breaks <- c(0, 0.33, 0.66, 1)
  
  # Classify data
  bivar_data_clean$richness_class <- cut(bivar_data_clean$richness,
                                         breaks = richness_breaks,
                                         labels = 1:3,
                                         include.lowest = TRUE)
  bivar_data_clean$proportion_class <- cut(bivar_data_clean$proportion,
                                           breaks = proportion_breaks,
                                           labels = 1:3,
                                           include.lowest = TRUE)
  
  # Create bi_class
  bivar_data_clean$bi_class <- paste0(as.numeric(bivar_data_clean$richness_class),
                                      "-",
                                      as.numeric(bivar_data_clean$proportion_class))
  
  # Create coordinates for plotting (now in Albers projection)
  coords <- terra::xyFromCell(richness_albers, 1:terra::ncell(richness_albers))
  valid_cells <- which(!is.na(richness_vals[, 1]) & richness_vals[, 1] > 0)
  
  bivar_df <- data.frame(
    x = coords[valid_cells, 1],
    y = coords[valid_cells, 2],
    bi_class = bivar_data_clean$bi_class
  )
  
  # Define color palette
  bivar_colors <- c(
    "1-1" = "#e8e8e8", "2-1" = "#ace4e4", "3-1" = "#5ac8c8",
    "1-2" = "#dfb0d6", "2-2" = "#a5add3", "3-2" = "#5698b9",
    "1-3" = "#d272aa", "2-3" = "#ad6aad", "3-3" = "#7759a1"
  )
  
  # Convert Albers boundary to sf
  canada_albers_sf <- sf::st_as_sf(canada_albers)
  
  # Create main map
  main_map <- ggplot() +
    geom_sf(data = canada_albers_sf, fill = "whitesmoke", color = "darkgrey", linewidth = 0.3) +
    geom_tile(data = bivar_df, aes(x = x, y = y, fill = bi_class), 
              color = "grey70", linewidth = 0.1) +
    geom_sf(data = canada_albers_sf, fill = NA, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = bivar_colors, guide = "none") +
    coord_sf(xlim = c(st_bbox(canada_albers_sf)[1], st_bbox(canada_albers_sf)[3]),
             ylim = c(st_bbox(canada_albers_sf)[2], st_bbox(canada_albers_sf)[4]),
             expand = FALSE) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Create legend
  legend_df <- expand.grid(x = 1:3, y = 1:3)
  legend_df$bi_class <- paste0(legend_df$x, "-", legend_df$y)
  
  legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = bi_class)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_manual(values = bivar_colors, guide = "none") +
    scale_x_continuous(breaks = 1:3,
                       labels = paste0(round(richness_breaks[1:3], 0), "-",
                                       round(richness_breaks[2:4], 0)),
                       position = "bottom",
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:3,
                       labels = paste0(proportion_breaks[1:3], "-",
                                       proportion_breaks[2:4]),
                       position = "left",
                       expand = c(0, 0)) +
    labs(x = "Tree Species Richness",
         y = "EMF Coverage") +
    coord_fixed() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )
  
  # Combine maps using patchwork
  final_map <- main_map + 
    patchwork::inset_element(legend_plot,
                             left = 0.74, bottom = 0.68,
                             right = 0.99, top = 0.87)
  
  # Save map
  ggsave(file.path(figures_dir, "canada_emf_bivariate_map.png"),
         plot = final_map, width = 14, height = 10, dpi = 300)
  
  message("✓ Bivariate map created and saved")
  
  return(final_map)
}

#' Create all maps
#' 
#' @param data_dir Character. Directory containing data
#' @param figures_dir Character. Directory to save figures
#' @param resolution_deg Numeric. Resolution used for raster analysis
#' @return List of all created maps
create_all_maps <- function(data_dir = here::here("data"), 
                            figures_dir = here::here("figures"), 
                            resolution_deg = 1) {
  
  message("Creating all maps...")
  
  # Ensure figures directory exists
  if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
  }
  
  # Create all maps
  base_map <- create_base_map(data_dir)
  sampling_map <- create_sampling_map(data_dir, figures_dir)
  emf_richness_map <- create_emf_richness_map(data_dir, figures_dir, resolution_deg)
  emf_proportion_map <- create_emf_proportion_map(data_dir, figures_dir, resolution_deg)
  bivariate_map <- create_bivariate_map(data_dir, figures_dir, resolution_deg)
  
  # Create summary
  maps_summary <- list(
    base_map = base_map,
    sampling_map = sampling_map,
    emf_richness_map = emf_richness_map,
    emf_proportion_map = emf_proportion_map,
    bivariate_map = bivariate_map,
    figures_directory = figures_dir,
    timestamp = Sys.time()
  )
  
  # Save maps summary
  saveRDS(maps_summary, file.path(data_dir, "maps_summary.rds"))
  
  message("✓ All maps created and saved")
  
  return(maps_summary)
}