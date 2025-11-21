# ============================================================
# Script: clim_connect_omni_cmd.R
#
# Purpose:
#   Helper functions for Omniscape-based climate connectivity:
#     - get_suit_block(): choose block_size from radius–block lookup
#     - clim_connect(): climate-warped connectivity (chains time steps)
#     - clim_connect_single(): baseline connectivity (no chaining).
#
# Used by:
#   - 3.1_block_size_sensitivity.R
#   - 4.2_bmd_connect_warp.R
#   - (and other scripts that call clim_connect / clim_connect_single)
#
# Inputs (per call):
#   sp                species name (used as folder name)
#   disersal_dists    vector of dispersal radii between time steps (in cells)
#   suit_list         raster stack of suitability per time period
#   habitat_list      raster stack of habitat patches per time period
#   time_periods      character vector of time labels (same length as stacks)
#   work_dir          base output directory for Omniscape runs
#   config            Omniscape config list (from ini::read.ini)
#   allow_min_radius  data.frame with allow_min_radius × block_size
#   julia_home        path to Julia binary (for run_omniscape.jl)
#
# Outputs (side effects):
#   work_dir/<sp>/<time_period>/cum_currmap.tif and other Omniscape products
#   (created or skipped if already present).
#
# Author: Lei Song <lei.song@rutgers.edu>
# Last updated: 2025-11-21
# ============================================================

# Load libraries
library(ini)
library(terra)
library(sf)
library(dplyr)

# Function to get suitable block size
get_suit_block <- function(x, df){
  df <- df[order(df$allow_min_radius), ]
  matched_row <- tail(df[df$allow_min_radius <= x, ], 1)
  
  if (nrow(matched_row) == 0) return(1)
  return(matched_row$block_size)
}

# Define the function to do climate connectivity modeling
clim_connect <- function(sp,
                         disersal_dists,
                         suit_list,
                         habitat_list,
                         time_periods,
                         work_dir = ".",
                         config,
                         allow_min_radius,
                         julia_home = NULL,
                         ...){
  # Create working directory
  sp_dir <- file.path(work_dir, sp)
  if (!dir.exists(sp_dir)) dir.create(sp_dir)
  
  # Check time periods and layers
  if ((!nlyr(suit_list) == length(time_periods)) |
      (!nlyr(habitat_list) == length(time_periods))){
    stop("Different layer numbers for time periods.")
  } else {
    names(suit_list) <- time_periods
    names(habitat_list) <- time_periods
  }
  
  # Process the layers for omniscape
  for (i in 1:(length(time_periods) - 1)){
    # Directory path
    project_name <- file.path(sp_dir, time_periods[[i + 1]])
    
    # Do the job if not yet
    if (!file.exists(file.path(project_name, "cum_currmap.tif"))){
      disersal_dist <- disersal_dists[i]
      config$Options$radius <- disersal_dist
      config$Options$buffer <- ceiling(disersal_dist / 2)
      config$Options$block_size <- get_suit_block(disersal_dist, allow_min_radius)
      
      # Set names for runs
      nm <- time_periods[[i]]
      ground_name <- file.path(sp_dir, sprintf("ground_%s.asc", nm))
      source_name <- file.path(sp_dir, sprintf("source_%s.asc", nm))
      suit_name <- file.path(sp_dir, sprintf("suitability_%s.asc", nm))
      fname <- file.path(sp_dir, sprintf("config_%s.ini", nm))
      
      # Load layers
      suit_t1 <- suit_list[[time_periods[[i]]]] + 0.0001
      habitat_t1 <- habitat_list[[time_periods[[i]]]]
      suit_t2 <- suit_list[[time_periods[[i + 1]]]] + 0.0001
      habitat_t2 <- habitat_list[[time_periods[[i + 1]]]]
      grounds <- mask(suit_t1, habitat_t1)
      sources <- mask(suit_t2, habitat_t2)
      
      # Update suit_t1 if not the first run
      if (i > 1){
        cum_dir <- file.path(sp_dir, time_periods[[i]])
        connt <- rast(file.path(cum_dir, "cum_currmap.tif"))
        connt <- stretch(connt, minv = 0, maxv = 1) + 1
        suit_t1 <- suit_t1 * connt
      }
      
      # Write out layers and run 
      # Grounds
      grounds[is.na(grounds)] <- -9999
      writeRaster(grounds, ground_name, overwrite = TRUE,
                  wopt = list(NAflag = -9999))
      # Sources
      sources[is.na(sources)] <- -9999
      writeRaster(sources, source_name, overwrite = TRUE,
                  wopt = list(NAflag = -9999))
      
      # Suitability
      suit_t1[is.na(suit_t1)] <- -9999
      writeRaster(suit_t1, suit_name, overwrite = TRUE,
                  wopt = list(NAflag = -9999))
      
      # Update config
      config$`Input files`$resistance_file <- suit_name
      config$`Input files`$ground_file <- ground_name
      config$`Input files`$source_file <- source_name
      config$Options$project_name <- project_name
      write.ini(config, fname)
      
      # Run omniscape
      cmd <- sprintf(
        "%s/julia --threads %s R/run_omniscape.jl -i %s", 
        julia_home, config$`Output options`$parallel_batch_size, fname)
      result <- system(cmd)
      
      if (result != 0) {
        unlink(sp_dir, recursive = TRUE)
        stop(sprintf("Run failed for species %s.", sp))
      }
    }
  }
}

# Do not link prior time stamp
# This function provides a comparison to the proposed approach.
clim_connect_single <- function(sp,
                         disersal_dists,
                         suit_list,
                         habitat_list,
                         time_periods,
                         work_dir = ".",
                         config,
                         allow_min_radius,
                         julia_home = NULL,
                         ...){
  # Create working directory
  sp_dir <- file.path(work_dir, sp)
  if (!dir.exists(sp_dir)) dir.create(sp_dir)
  
  # Check time periods and layers
  if ((!nlyr(suit_list) == length(time_periods)) |
      (!nlyr(habitat_list) == length(time_periods))){
    stop("Different layer numbers for time periods.")
  } else {
    names(suit_list) <- time_periods
    names(habitat_list) <- time_periods
  }
  
  # Process the layers for omniscape
  for (i in 1:(length(time_periods) - 1)){
    # Directory path
    project_name <- file.path(sp_dir, time_periods[[i + 1]])
    
    # Do the job if not yet
    if (!file.exists(file.path(project_name, "cum_currmap.tif"))){
      disersal_dist <- disersal_dists[i]
      config$Options$radius <- disersal_dist
      config$Options$buffer <- ceiling(disersal_dist / 2)
      config$Options$block_size <- get_suit_block(disersal_dist, allow_min_radius)
      
      # Set names for runs
      nm <- time_periods[[i]]
      ground_name <- file.path(sp_dir, sprintf("ground_%s.asc", nm))
      source_name <- file.path(sp_dir, sprintf("source_%s.asc", nm))
      suit_name <- file.path(sp_dir, sprintf("suitability_%s.asc", nm))
      fname <- file.path(sp_dir, sprintf("config_%s.ini", nm))
      
      # Load layers
      suit_t1 <- suit_list[[time_periods[[1]]]] + 0.0001
      habitat_t1 <- habitat_list[[time_periods[[1]]]]
      suit_t2 <- suit_list[[time_periods[[i + 1]]]] + 0.0001
      habitat_t2 <- habitat_list[[time_periods[[i + 1]]]]
      grounds <- mask(suit_t1, habitat_t1)
      sources <- mask(suit_t2, habitat_t2)
      
      # Write out layers and run 
      # Grounds
      grounds[is.na(grounds)] <- -9999
      writeRaster(grounds, ground_name, overwrite = TRUE,
                  wopt = list(NAflag = -9999))
      # Sources
      sources[is.na(sources)] <- -9999
      writeRaster(sources, source_name, overwrite = TRUE,
                  wopt = list(NAflag = -9999))
      
      # Suitability
      suit_t1[is.na(suit_t1)] <- -9999
      writeRaster(suit_t1, suit_name, overwrite = TRUE,
                  wopt = list(NAflag = -9999))
      
      # Update config
      config$`Input files`$resistance_file <- suit_name
      config$`Input files`$ground_file <- ground_name
      config$`Input files`$source_file <- source_name
      config$Options$project_name <- project_name
      write.ini(config, fname)
      
      # Run omniscape
      cmd <- sprintf(
        "%s/julia --threads %s R/run_omniscape.jl -i %s", 
        julia_home, config$`Output options`$parallel_batch_size, fname)
      result <- system(cmd)
      
      if (result != 0) {
        unlink(sp_dir, recursive = TRUE)
        stop(sprintf("Run failed for species %s.", sp))
      }
    }
  }
}