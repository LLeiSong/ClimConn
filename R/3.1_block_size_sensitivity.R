# ============================================================
# Script: 3.1_block_size_sensitivity.R
#
# Purpose:
#   Sensitivity analysis of Omniscape block size vs. dispersal radius:
#     - Select representative species with intermediate dispersal distances
#     - Run Omniscape repeatedly for each species with varying block_size
#     - Generate cumulative current maps to derive a rule of thumb
#       for choosing block size given a radius.
#
# Inputs:
#   data/smr_species.csv                    summary species table with time periods & SSPs
#   data/dispersal/species_dispersal_rate.csv  species-level dispersal rates
#   data/colombia.geojson                   Colombia boundary
#   data/mamiferos/<sp>/reconstruct/...     current & future habitat / suitability rasters
#   data/config/omniscape_setting_template.ini  Omniscape configuration template
#
# Outputs:
#   results/omni_block_size/<sp>/2021-2040_*  cum_currmap_*.tif for multiple block sizes
#
# Usage guidance:
# Objective: of course using block = =1 would give us the most accurate pixel
# level results, but it is very time consuming. If the future possible habitat
# is large, it will take forever to run. This script is to run a few example
# species to detect a rule of thumb to select block size based on radius.
#
# Author: Lei Song <lei.song@rutgers.edu>
# Last updated: 2025-11-21
# ============================================================

# Load libraries
library(dplyr)
library(sf)
library(terra)
library(ini)
library(classInt)
library(parallel)
sf_use_s2(FALSE)

# Setting directories and parameters
code_dir <- "/home/ls1686/ClimConn"
proj_dir <- "/scratch/ls1686/ClimConn"

julia_home <- "/home/ls1686/.conda/envs/sdm/bin"
crs_analysis <- "ESRI:54030"
data_dir <- file.path(proj_dir, "data/mamiferos")
cur_path <- "reconstruct/current/MAXENT/%s_10_MAXENT.tif"
cur_suit_path <- "reconstruct/current/MAXENT/%s_MAXENT.tif"
years <- "2021-2040"
ssps <- "ssp126"
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
  proj_dir, "data/dispersal", "species_dispersal_rate.csv")
scenario <- "ssp126"

# Step1: select some representative species
smr_species <- read.csv(file.path(proj_dir, "data/smr_species.csv")) %>% 
  filter(!is.na(time_periods)) %>% 
  filter(ssp == scenario & time_periods == "1970-2000")

smr_species <- read.csv(dispersal_path) %>% 
  mutate(species = gsub(" ", ".", species)) %>% 
  select(species, dispersal_rate) %>% 
  left_join(smr_species, ., by = c("sp" = "species")) %>% 
  mutate(dispersal_dist = floor(dispersal_rate * 1000 * 40 / 832.3368)) %>% 
  # 5 to have at least 3 runs, 83 is below 90 percentile
  filter(dispersal_dist >= 5 & dispersal_dist <= 83)

# Select species
set.seed(123)
smr_species <- smr_species %>%
  sample_n(size = 60) %>%
  select(sp, time_periods, dispersal_rate)

smr_species <- smr_species %>% 
  rbind(smr_species %>% mutate(time_periods = years))

message(sprintf("%s species are selected.", nrow(smr_species)))

# Load function
source(file.path(code_dir, "R/clim_connect_omni_cmd.R"))

# Load files 
col <- st_read(file.path(proj_dir, "data/colombia.geojson")) %>% 
  st_transform(crs_analysis)

round_down_to_odd <- function(x) {
  f <- floor(x)
  if (f %% 2 == 0) f - 1 else f
}

seq_int <- function(dispersal_dists, max_n = 10) {
  if (length(dispersal_dists) != 1L) stop("Use one value at a time.")
  dispersal_dists <- round_down_to_odd(dispersal_dists)
  
  s <- seq(1, dispersal_dists, by = 2)
  
  if (length(s) > max_n){
    s <- c(s[1:5], s[round(seq(6, length(s) - 1, length.out = max_n - 6))], 
           s[length(s)])
  }
  
  s
}

# Step2: for each species, use block size = 1 up to the radius size
# Step3: run models with all block sizes
# Step4: collect and save results
parallel::mclapply(unique(smr_species$sp), function(sp){
  # Define the time periods
  time_periods <- unique(smr_species$time_periods)
  
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
  disersal_dists <- dispersal_rate * 1000 * 40 # in meter
  
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
      if (length(null_id) != 0){
        stop(sprintf("%s will lose all connection to future in Colombia.", sp))
      }
      
      fut_ranges <- do.call(c, fut_ranges)
      n_tp <- length(time_periods)
      
      hbts <- c(cur_patches, fut_ranges)
      names(hbts) <- time_periods
      
      # Extract suitability for patches
      cur_suit <- rast(cur_suit_fn) %>% project(crs_analysis) %>% 
        mask(aoi_spreads[[1]])
      fut_suits <- rast(fut_suit_fns) %>% 
        project(crs_analysis) %>% mask(aoi_spreads[[1]])
      suits <- c(cur_suit, fut_suits)
      names(suits) <- time_periods
      
      lyrs <- c(hbts, suits) %>% trim()
      
      hbts <- subset(lyrs, 1:n_tp)
      hbts[hbts == 0] <- NA
      names(hbts) <- time_periods
      
      suits <- subset(lyrs, (n_tp + 1):(2 * n_tp))
      names(suits) <- time_periods
      
      disersal_dists <- ceiling(disersal_dists / res(suits)[1])
      
      circ_dir <- file.path(proj_dir, "results", "omni_block_size")
      if (!dir.exists(circ_dir)) dir.create(circ_dir, recursive = TRUE)
      
      config_template <- read.ini(file.path(
        proj_dir, "data/config", "omniscape_setting_template.ini"))
      
      for (i in seq_int(disersal_dists)){
        config_template$Options$block_size <- i
        
        tryCatch({
          clim_connect(
            sp, disersal_dists, suits, hbts, time_periods, 
            circ_dir, config_template, julia_home, TRUE)
        }, error = function(e){print(e)})
        
        # Rename the layer
        file.rename(
          file.path(circ_dir, sp, "2021-2040", "cum_currmap.tif"),
          file.path(circ_dir, sp, "2021-2040", sprintf("cum_currmap_%s.tif", i)))
        
        # Rename the folder
        old_folder <- file.path(circ_dir, sp, "2021-2040")
        new_folder <- file.path(circ_dir, sp, sprintf("2021-2040_%s", i))
        dir.create(new_folder)
        file.copy(
          list.files(old_folder, full.names = TRUE), 
          new_folder, recursive = TRUE)
        
        unlink(old_folder, recursive = TRUE)
      }
    } else {
      message(sprintf("%s has no valid big enough habitat in Colombia.", sp))
    }
  } else{
    message(sprintf("%s has no valid habitat in Colombia.", sp))
  }
}, mc.cores = 30)
