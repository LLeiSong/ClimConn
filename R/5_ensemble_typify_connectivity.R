# ============================================================
# Script: 5_ensemble_typify_connectivity.R
#
# Purpose:
#   Ensemble and typify connectivity maps:
#     - Average Omniscape cumulative current maps across species
#       for each year × SSP combination
#     - Create a single-species comparison ensemble (ssp126)
#     - Derive a categorical connectivity typology (F1/F2/I1/I2)
#       using sinks/peaks and upper/lower quantiles.
#
# Inputs:
#   data/colombia.geojson                        Colombia boundary
#   data/Env/M_variable/vars.tif                 template raster (for grid)
#   data/smr_species.csv                         species summary (to select sps)
#   results/omni_single/<ssp>/<sp>/<year>/cum_currmap.tif  single-species runs
#   results/omni/<ssp>/<sp>/<year>/cum_currmap.tif         main Omniscape runs
#
# Outputs:
#   results/omni_single/conn_lyrs.tif            ensemble (ssp126) connectivity stack
#   results/omni/conn_lyrs.tif                   multi-SSP ensemble connectivity stack
#   results/typology/conn_types_<year>_<ssp>.tif connectivity typology rasters
#   results/typology/connectivity_typology.csv   class catalog (Background, F1, F2, I1, I2)
#
# Author: Lei Song <lei.song@rutgers.edu>
# Last updated: 2025-11-21
# ============================================================

# Load libraries and functions
library(here)
library(smoothr)
library(stars)
library(sf)
library(terra)
library(dplyr)
library(rmapshaper)
library(whitebox)
sf_use_s2(FALSE)
source("R/utils.R")

# Setting
root_dir <- here()
data_dir <- file.path(root_dir, "data")
result_dir <- file.path(root_dir, "results/omni")
result_cmp <- file.path(root_dir, "results/omni_single")
typology_dir <- file.path(root_dir, "results/typology")
if (!dir.exists(typology_dir)) dir.create(typology_dir)

# Parameters
years <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
ssps <- c("ssp126", "ssp370", "ssp585")
crs_analysis <- "ESRI:54030"
smr_species <- read.csv(file.path(data_dir, "smr_species.csv"))
sps <- smr_species %>% filter(big_habitat == 1) %>% 
  mutate(time_periods = ifelse(
    is.na(time_periods), "1970-2000", time_periods)) %>% 
  filter(time_periods == "1970-2000") %>% pull(sp) %>% unique()

#### Connectivity #####
## Make a template to burn the results in
col <- st_read(file.path(data_dir, "colombia.geojson"))
template <- rast(file.path(data_dir, "Env/M_variable/vars.tif"), lyrs = 1)
values(template) <- 0
template <- mask(crop(template, st_buffer(col, 0.5)), st_buffer(col, 0.5))
template <- template %>% project(crs_analysis, method = "near")

#### Connectivity comparison #####
ssps_cmp <- "ssp126"
conn_lyrs <- do.call(c, lapply(years, function(yr){
  lyrs <- do.call(c, lapply(ssps_cmp, function(ssp){
    message(sprintf("%s %s", yr, ssp))
    lyrs <- lapply(sps, function(sp){
      fname <- file.path(result_cmp, ssp, sp, yr, "cum_currmap.tif")
      if (file.exists(fname)){
        curmap <- rast(fname)
        curmap <- crop(curmap, template)
        if (global(curmap, sum, na.rm = TRUE)[[1]] != 0 |
            is.nan(global(curmap, sum, na.rm = TRUE)[[1]])){
          curmap <- stretch(curmap, minv = 0, maxv = 1)
          terra::extend(curmap, ext(template), fill = NA)
        }
      }
    })
    
    while(is(lyrs, "list")){
      lyrs <- do.call(c, lyrs)
    }
    
    mask(mean(lyrs, na.rm = TRUE), col %>% st_transform(crs_analysis))
  }))
  
  names(lyrs) <- sprintf("%s_%s", yr, ssps_cmp)
  lyrs
}))

# With a buffer
writeRaster(conn_lyrs, file.path(result_cmp, "conn_lyrs.tif"),
            overwrite = TRUE)

# Start the burning
conn_lyrs <- do.call(c, lapply(years, function(yr){
  lyrs <- do.call(c, lapply(ssps, function(ssp){
    message(sprintf("%s %s", yr, ssp))
    lyrs <- lapply(sps, function(sp){
      fname <- file.path(result_dir, ssp, sp, yr, "cum_currmap.tif")
      if (file.exists(fname)){
        curmap <- rast(fname)
        curmap <- crop(curmap, template)
        if (global(curmap, sum, na.rm = TRUE)[[1]] != 0 |
            is.nan(global(curmap, sum, na.rm = TRUE)[[1]])){
          curmap <- stretch(curmap, minv = 0, maxv = 1)
          terra::extend(curmap, ext(template), fill = NA)
        }
      }
    })
    
    while(is(lyrs, "list")){
      lyrs <- do.call(c, lyrs)
    }
    
    mask(mean(lyrs, na.rm = TRUE), col %>% st_transform(crs_analysis))
  }))
  
  names(lyrs) <- sprintf("%s_%s", yr, ssps)
  lyrs
}))

# With a buffer
writeRaster(conn_lyrs, file.path(result_dir, "conn_lyrs.tif"),
            overwrite = TRUE)

#### Reclassify connectivity #####
for (nm in names(conn_lyrs)){
  message(nm)
  
  # Get the layer, stretch it, does not really matter
  lyr <- conn_lyrs[[nm]]
  
  # STEP 0. Sinks and plateaus
  conn_file <- file.path(tempdir(), "conn.tif")
  conn_rev_file <- file.path(tempdir(), "conn_rev.tif")
  sink_file <- file.path(tempdir(), "sink.tif")
  peak_file <- file.path(tempdir(), "peak.tif")
  writeRaster(lyr, conn_file)
  writeRaster(1 - lyr, conn_rev_file)
  wbt_sink(conn_file, sink_file)
  wbt_sink(conn_rev_file, peak_file)
  
  pits <- rast(sink_file)
  pits <- patch_clean(pits)
  
  peaks <- rast(peak_file)
  peaks <- patch_clean(peaks)
  
  # Clean up
  file.remove(conn_file, conn_rev_file, sink_file, peak_file)
  
  # STEP 1. peaks and ridges
  top_all <- lyr >= quantile(values(lyr), 0.9, na.rm = TRUE)
  top_all <- patch_clean(top_all)
  top <- top_all * peaks * 10 + peaks
  
  # STEP 2. lowlands
  low_all <- lyr <= quantile(values(lyr), 0.1, na.rm = TRUE)
  low_all <- patch_clean(low_all)
  low <- low_all * pits * 20 + pits * 2
  
  # STEP 4. combine
  msk <- (top != 0) + (low != 0) > 1
  classes <- mask(top + low, col %>% st_transform(crs_analysis))
  classes <- mask(classes, msk, maskvalues = 1, updatevalue = 0)
  classes <- classes %>% trim()
  
  # Step 3: Set it as a categorical (factor) raster
  levels(classes) <- data.frame(
    ID = c(0, 1, 11, 2, 22),
    class = c("Background", "F1", "F2", "I1", "I2"))
  
  fname <- file.path(typology_dir, sprintf("conn_types_%s.tif", nm))
  writeRaster(classes, fname, overwrite = TRUE)
}

# Write out the class catalog
classes <- data.frame(
  ID = c(0, 1, 11, 2, 22),
  class = c("Background", "F1", "F2", "I1", "I2"),
  name = c("Background", "Facilitative area (Level 1)", 
           "Facilitative area (Level 2)", 
           "Impeded area (Level 1)",
           "Impeded area (Level 2)"))

fname <- file.path(typology_dir, "connectivity_typology.csv")
write.csv(classes, fname, row.names = FALSE)
