## Prerequisites

- R (>= 4.x), Julia (>= 1.x), WhiteboxTools, GDAL/PROJ.
- R packages: sf, terra, dplyr, readxl, here, stringr, optparse, ini,
  rnaturalearth, RCurl, smoothr, stars, rmapshaper, whitebox,
  ggplot2, tidyterra, ggsci, ggpubr, patchwork, ggsankey, parallel, tidyr.

Clone the repo and set your working directory to the project root
(e.g., via an RStudio project). Scripts assume paths are relative to this root.

### Required data in `data/`

- `mamiferos/` – BioModelos SDM outputs, one folder per species.
- `IUCN/MAMMALS/MAMMALS.shp` – IUCN mammals shapefile.
- `dispersal/Generation Lenght for Mammals.xlsx` – trait data (body mass & generation length).
- `colombia.geojson` – Colombia boundary polygon.
- `config/omniscape_setting_template.ini` – Omniscape configuration template.

Other folders (`Env/`, `dispersal/species_dispersal_rate.csv`,
`results/`, etc.) are created by the scripts below.


## Select species and identify species traits

1. Run `1.1_scrutinize_species.R` to remote marine or freshwater only species from the species list of BioModelos.
2. Run `1.2_dispersal_rate.R` to calculate dispersal rate for each species.

## Model habitat quality

1. Run `2.1_prepare_variables.R` to prepare all environmental layers for modeling and future projections.
2. Run `2.2_reconstruct_models.R` to reconstruct BioModelos SDMs for selected species.
3. (Optional) Run `2.3_check_reconstruct.R` to check the reconstructed results and remove temporary files.

## Block size analysis

1. Run `3.1_block_size_sensitivity.R` to run Omniscape with 60 randomly selected species with different block sizes.
2. Run `3.2_pair_size_radius.R` to calculate the maximum block size for different dispersal distance.

## Assess connectivity

1. Run `4.1_summarize_species.R` to check how many time stamps all species will maintain connectivity to new habitats.
2. Run `4.2_bmd_connect_warp.R` to model omni-directional and multi-temporal habitat connectivity for all species. This script calls function `clim_connect` defined in file `clim_connect_omni_cmd.R`, which runs the model by calling the Julia script `run_omniscape.jl`. Alternatively, `clim_connect_omni.R` defines another version of `clim_connect` using R package `JuliaCall` to run Julia script.

## Typify connectivity

1. Run `5_ensemble_typify_connectivity.R` to ensemble connectivity layers and typify based on hydrological analysis.

## Other functions

- `utils.R` has function to clean small patches.
- `figures.R` make figures for the manuscript.