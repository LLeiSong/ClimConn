# ============================================================
# Script: 4.2_bmd_connect_warp.R
#
# Purpose:
#   Run Omniscape-based climate connectivity for selected BioModelos species:
#     - Prepare current and future habitat patches and suitability rasters
#     - Use block_size–radius relationship (allow_min_radius) to choose settings
#     - Call clim_connect to generate time-warped connectivity surfaces.
#
# Inputs:
#   Command-line:
#     -s / --sp        species name(s), comma-separated (folder names in data/mamiferos)
#     -o / --scenario  climate scenario (e.g., ssp126, ssp370, ssp585)
#
#   Files:
#     data/mamiferos/<sp>/reconstruct/current/MAXENT/       current binary/suitability
#     data/mamiferos/<sp>/reconstruct/future/...            future binary/suitability
#     data/mamiferos/<sp>/interest_areas/shape_M.shp        species AOI
#     data/colombia.geojson                                Colombia boundary
#     data/dispersal/species_dispersal_rate.csv            species dispersal rates
#     results/block_radius/allow_min_radius.csv            block_size–radius lookup
#     data/config/omniscape_setting_template.ini           Omniscape config template
#
# Outputs:
#   results/omni/<scenario>/<sp>/...  Omniscape connectivity outputs (created by clim_connect)
#
# Author: Lei Song <lei.song@rutgers.edu>
# Last updated: 2025-11-21
# ============================================================

# Load libraries
library(ini)
library(readxl)
library(dplyr)
library(sf)
library(terra)
library(optparse)
sf_use_s2(FALSE)

# Parse inline parameters
option_list <- list(
  make_option(c("-s", "--sp"), 
              action = "store", type = 'character',
              help = "The species to process."),
  make_option(c("-o", "--scenario"), 
              action = "store", type = 'character',
              help = "The scenario to process."))
opt <- parse_args(OptionParser(option_list = option_list))
sp <- opt$sp
sp_list <- strsplit(sp, ",")[[1]]
scenario <- opt$scenario

# Define paths and parameters
code_dir <- "/home/ls1686/ClimConn"
root_dir <- "/scratch/ls1686/ClimConn"
# Local: /Users/leisong/.juliaup/bin
julia_home <- "/home/ls1686/.conda/envs/sdm/bin"
crs_analysis <- "ESRI:54030"
years <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
ssps <- c("ssp126", "ssp370", "ssp585")
data_dir <- file.path(root_dir, "data/mamiferos")
cur_path <- "reconstruct/current/MAXENT/%s_10_MAXENT.tif"
cur_suit_path <- "reconstruct/current/MAXENT/%s_MAXENT.tif"
fut_paths <- sprintf(
  "%s/%s_10.tif", 
  "reconstruct/future/MAXENT/future_ensemble", 
  apply(expand.grid(years, ssps), 1, paste, collapse = "_"))
fut_suit_paths <- sprintf(
  "%s/%s_maxent.tif", 
  "reconstruct/future/MAXENT/future_ensemble", 
  apply(expand.grid(years, ssps), 1, paste, collapse = "_"))
aoi_path <- "interest_areas/shape_M.shp"
dispersal_path <- file.path(
  root_dir, "data/dispersal", "species_dispersal_rate.csv")

# Load radius_block relationship
allow_min_radius <- read.csv(
  file.path(root_dir, "results/block_radius", "allow_min_radius.csv"))

# Load function
source(file.path(code_dir, "R/clim_connect_omni_cmd.R"))

# Load files 
col <- st_read(file.path(root_dir, "data/colombia.geojson")) %>% 
  st_transform(crs_analysis)

for (sp in sp_list){
  # Define the time periods
  time_periods <- c("1970-2000", years)
  
  # Get the file names
  ## Current patch and suitability
  cur_fname <- file.path(
    data_dir, sp, sprintf(cur_path, gsub("\\.", "_", sp)))
  cur_suit_fn <- file.path(
    data_dir, sp, sprintf(cur_suit_path, gsub("\\.", "_", sp)))
  
  ## Future patches and suitabilities
  fut_fnames <- file.path(data_dir, sp, fut_paths[grep(scenario, fut_paths)])
  fut_suit_fns <- file.path(
    data_dir, sp, fut_suit_paths[grep(scenario, fut_suit_paths)])
  aoi_fname <- file.path(data_dir, sp, aoi_path)
  
  # Get potential dispersal distances for different time periods
  dispersal_rate <- read.csv(dispersal_path) %>% 
    filter(species == gsub("\\.", " ", sp)) %>% pull(dispersal_rate)
  if (length(dispersal_rate) == 0) {
    dispersal_rate <- read.csv(dispersal_path) %>% 
      pull(dispersal_rate) %>% median()}
  disersal_dists <- dispersal_rate * 1000 * c(40, 60, 80, 100) # in meter
  
  # AOIs with dispersal buffer for different time periods
  aoi <- st_read(aoi_fname) %>% st_transform(crs_analysis)
  aoi_spreads <- lapply(disersal_dists, function(disersal_dist){
    aoi_spread <- st_buffer(aoi, disersal_dist)
    col_spread <- st_buffer(col, disersal_dist * 2)
    st_intersection(aoi_spread, col_spread)
  })
  
  # Remove patches smaller than 100 square km
  cur_range <- rast(cur_fname) %>% project(crs_analysis) %>% 
    mask(aoi) %>% mask(st_buffer(col, disersal_dists[1] * 2))
  if (global(cur_range, "notNA")[[1]] != 0){
    cur_patches <- patches(cur_range, zeroAsNA = TRUE)
    rz <- zonal(cellSize(cur_patches, unit = "km"), 
                cur_patches, sum, as.raster = TRUE)
    cur_patches <- ifel(rz < 100, NA, cur_patches)
    
    if (global(cur_patches, "notNA")[[1]] != 0){
      fut_ranges <- lapply(1:length(fut_fnames), function(i){
        fut_fname <- fut_fnames[i]
        aoi_spread <- aoi_spreads[[i]]
        fut_range <- mask(rast(fut_fname) %>% project(crs_analysis), aoi_spread)
        if (global(fut_range, "notNA")[[1]] != 0){
          # Remove patches smaller than 100 square km
          fut_patches <- patches(fut_range, zeroAsNA = TRUE)
          rz <- zonal(cellSize(fut_patches, unit = "km"), 
                      fut_patches, sum, as.raster = TRUE)
          lyr <- ifel(rz < 100, NA, fut_patches)
          if (global(lyr, "notNA")[[1]] != 0){
            lyr
          } else NULL
        } else NULL
      })
      
      null_id <- which(sapply(fut_ranges, is.null))
      if (length(null_id) == 0){
        null_id <- length(time_periods)
      } else null_id <- min(null_id)
      
      if (null_id > 1){
        fut_ranges <- fut_ranges[1:(null_id - 1)]
        fut_ranges <- do.call(c, fut_ranges)
        time_periods <- time_periods[1:null_id]
        n_tp <- length(time_periods)
        
        hbts <- c(cur_patches, fut_ranges)
        names(hbts) <- time_periods
        
        # Extract suitability for patches
        cur_suit <- rast(cur_suit_fn) %>% project(crs_analysis) %>% 
          mask(aoi_spreads[[4]])
        fut_suits <- rast(fut_suit_fns[1:(null_id - 1)]) %>% 
          project(crs_analysis) %>% mask(aoi_spreads[[4]])
        suits <- c(cur_suit, fut_suits)
        names(suits) <- time_periods
        
        lyrs <- c(hbts, suits) %>% trim()
        
        hbts <- subset(lyrs, 1:n_tp)
        hbts[hbts == 0] <- NA
        names(hbts) <- time_periods
        
        suits <- subset(lyrs, (n_tp + 1):(2 * n_tp))
        names(suits) <- time_periods
        
        # Update dispersal distances to between each time points
        disersal_dists <- dispersal_rate * 1000 * c(40, 20, 20, 20)
        disersal_dists <- ceiling(disersal_dists / res(suits)[1])[1:(null_id - 1)]
        
        circ_dir <- file.path(root_dir, "results", "omni", scenario)
        if (!dir.exists(circ_dir)) dir.create(circ_dir, recursive = TRUE)
        
        config_template <- read.ini(file.path(
          root_dir, "data/config", "omniscape_setting_template.ini"))
        
        tryCatch({
          clim_connect(
            sp, disersal_dists, suits, hbts, time_periods, 
            circ_dir, config_template, allow_min_radius, julia_home, TRUE)
        }, error = function(e){print(e)})
      } else{
        message(sprintf("%s will lose all connection to future in Colombia.", sp))
      }
    } else {
      message(sprintf("%s has no valid big enough habitat in Colombia.", sp))
    }
  } else{
    message(sprintf("%s has no valid habitat in Colombia.", sp))
  }
}
