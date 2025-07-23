# EMF Canada Analysis

An assessment of knowledge gaps regarding ectomycorrhizal fungi (EMF) biodiversity across Canada's forest ecosystems.

## Overview

This project was funded by a working group grant from the Canadian Institute for Ecology & Evolution in 2025.

Contributors (members of the working group) include:
Jason Pither (UBC)
Justine Karst (U Alberta)
Bronte Shelton (UBC)
Diego Yusta (UBC)
Isaac Eckert (McGill U)
Monika Fischer (UBC)
Stephanie Kivlin (U Tennessee)
Clara Qin (SPUN)

This R project addresses three of the "shortfalls" discussed in the manuscript:
 
- **Eltonian Shortfall**: Gaps in species interactions knowledge
- **Hutchinsonian Shortfall**: Gaps in environmental niche understanding  
- **Wallacean Shortfall**: Gaps in species distribution data

## Authors

This R project was assembled by Jason Pither, with additional coding contributions from Isaac Eckert and Clara Qin

## Contact info

Jason Pither: jason<dot>pither[at]ubc<dot>ca

## Project Structure

```
EMF_Canada/
├── _README.md                    # This file
├── scripts/                     # R analysis scripts
│   ├── _README.md
│   ├── 01_setup.R              # Environment setup
│   ├── 02_download_data.R      # Data acquisition
│   ├── 03_process_spatial.R    # Spatial data processing
│   ├── 04_process_fungal.R     # EMF data processing
│   ├── 05_calculate_metrics.R  # Analysis metrics
│   ├── 06_create_maps.R        # Visualization
│   ├── run_complete_analysis.R # Master pipeline
│   └── utils.R                 # Helper functions
├── rawdata/                    # Input data (see rawdata/README.md)
│   ├── _README.md
│   ├── species_ranges.csv      # Tree species metadata
│   ├── vanGalen_Data/          # van Galen EMF dataset
│   ├── fungalroot/             # FungalRoot database
│   ├── Stephanie/              # Stephanie's EMF data
│   ├── ecoregions/             # Canadian ecoregions (downloaded)
│   ├── USTreeAtlas/            # US Tree Atlas (downloaded)
│   ├── admin_boundaries/       # Administrative boundaries (downloaded)
│   └── natural_earth/          # Natural Earth data (downloaded)
├── outdata/                    # Processed outputs
│   └── _README.md
└── figures/                    # Generated visualizations
    └── _README.md
```

## Quick Start

### Prerequisites

- R (≥ 4.0.0)
- Required packages (automatically installed):
  ```r
  c("here", "terra", "sf", "tidyverse", "geodata", "rnaturalearth", 
    "tidyterra", "ggplot2", "patchwork", "dplyr", "magrittr")
  ```

### Running the Complete Analysis

1. **Clone or download this repository**
2. **Prepare input data** (see `rawdata/_README.md` for requirements)
3. **Run the complete pipeline**:
   ```r
   source("run_complete_analysis.R")
   ```

The pipeline will automatically:
- Set up the analysis environment
- Download external datasets  
- Process spatial and fungal data
- Calculate biodiversity metrics
- Generate maps and visualizations
- Create summary tables

### Individual Script Usage

Scripts can also be run individually in order:

```r
# 1. Environment setup
source("scripts/01_setup.R")

# 2. Data download  
source("scripts/02_download_data.R")

# 3. Spatial processing
source("scripts/03_process_spatial.R")

# 4. Fungal data processing
source("scripts/04_process_fungal.R")

# 5. Calculate metrics
source("scripts/05_calculate_metrics.R")

# 6. Create visualizations
source("scripts/06_create_maps.R")
```

## Key Outputs

### Data Products
- **Processed EMF host tree ranges**: `outdata/emf_host_tree_ranges.gpkg`
- **EMF sampling locations**: `outdata/emf_samples_aggregated_1000m.gpkg`
- **Species richness rasters**: `outdata/species_richness_1deg.tif`
- **EMF coverage rasters**: `outdata/emf_proportion_1deg.tif`

### Visualizations
- **Sampling coverage map**: `figures/canada_sampling_locations_map.png`
- **EMF richness map**: `figures/canada_emf_richness_map.png`
- **Data coverage map**: `figures/canada_emf_proportion_map.png`  
- **Bivariate analysis map**: `figures/canada_emf_bivariate_map.png`

### Summary Statistics
- **Basic EMF statistics**: `outdata/emf_summary_statistics.csv`
- **Knowledge gaps summary**: `outdata/knowledge_gaps_summary.csv`
- **Spatial coverage metrics**: `outdata/proportion_statistics.csv`

## Data Sources

### Required Input Data
- **FungalRoot Database**: Mycorrhizal type assignments
- **EMF Sequence Data**: Stephanie Kivlin's Canadian EMF samples dataset
- **Tree Species Ranges**: Metadata for North American tree species

### External Data (Auto-downloaded)
- **GADM Administrative Boundaries**: Country/province boundaries
- **Canadian Ecoregions**: Ecological classification system
- **US Tree Atlas**: Species range maps
- **Natural Earth**: Base mapping data

## Analysis Workflow

1. **Data Integration**: Links tree species ranges with mycorrhizal type data
2. **Spatial Processing**: Creates unified spatial datasets for Canada
3. **EMF Mapping**: Identifies ectomycorrhizal host species distributions
4. **Gap Analysis**: Quantifies sampling coverage by ecoregion
5. **Metric Calculation**: Computes richness and coverage statistics
6. **Visualization**: Generates comprehensive maps and summary figures

## Key Findings

The analysis reveals significant knowledge gaps in Canadian EMF research:

- **Species Interactions**: Only a fraction of EMF host species have associated sequence data
- **Geographic Coverage**: Large areas of potential EMF habitat remain unsampled
- **Taxonomic Coverage**: Genus-level coverage varies substantially across regions

## Dependencies

### R Packages
- **Spatial**: `terra`, `sf`, `geodata`, `rnaturalearth`
- **Data**: `tidyverse`, `dplyr`, `magrittr`, `readr`
- **Visualization**: `ggplot2`, `tidyterra`, `patchwork`
- **Utilities**: `here`, `knitr`

### System Requirements
- Internet connection (for data download)
- Minimum 8GB RAM recommended
- ~5GB free disk space

## Citation

If you use this analysis framework, please cite:

```
TBD
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- FungalRoot database contributors
- Natural Resources Canada for ecoregions data
- US Forest Service for Tree Atlas data
- Natural Earth for base mapping data