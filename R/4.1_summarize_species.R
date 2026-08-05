# ============================================================
# Script: 4.1_summarize_species.R
#
# Purpose:
#   Summarize species-level habitat status and connectivity potential:
#     - Quantify current and minimum future habitat area in Colombia
#     - Check if species have any / sufficiently large habitat patches
#     - Assess whether current and future habitats remain connected
#       within species-specific dispersal abilities.
#
# Inputs:
#   data/mamiferos/<sp>/reconstruct/current/MAXENT/   current binary/suitability rasters
#   data/mamiferos/<sp>/reconstruct/future/...        future binary rasters (all SSPs/periods)
#   data/mamiferos/<sp>/interest_areas/shape_M.shp    species AOI
#   data/colombia.geojson                             Colombia boundary
#   data/dispersal/species_dispersal_rate.csv         species-level dispersal rates
#
# Outputs:
#   data/smr_species.csv          running table of species × SSP summary (append mode)
#   data/smr_species_final.csv    full combined summary across all SSPs
#
# ============================================================

# Load libraries
library(readxl)
library(dplyr)
library(sf)
library(terra)
library(parallel)
sf_use_s2(FALSE)

# Define paths and parameters
code_dir <- "/home/ls1686/ClimConn"
root_dir <- "/scratch/ls1686/ClimConn"
result_dir <- file.path(root_dir, "results/omni")
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

# Load files 
col <- st_read(file.path(root_dir, "data/colombia.geojson")) %>% 
  st_transform(crs_analysis)

sps <- list.files(file.path(root_dir, "data/mamiferos"))

# Remove ongoing ones
sp_list <- lapply(ssps, function(ssp){
  data.frame(sp = sps, ssp = ssp)
}) %>% bind_rows()

# Calculated area of Colombia in square km.
col_area <- 1132385

df <- data.frame(
  sp = character(0), 
  ssp = character(0),
  habitat = numeric(0), 
  big_habitat = numeric(0),
  connect = numeric(0), 
  time_periods = character(0),
  area_current = numeric(0), 
  area_future = numeric(0),
  area_overlap = numeric(0))

write.csv(df, row.names = FALSE,
          file.path(code_dir, "data/smr_species.csv"))

smr_species <- mclapply(1:nrow(sp_list), function(i){
  record <- sp_list %>% slice(i)
  sp <- record$sp
  scenario <- record$ssp
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
    st_intersection(aoi_spread, col) %>% st_make_valid()
  })
  
  # Remove patches smaller than 100 square km
  cur_range <- rast(cur_fname) %>% project(crs_analysis) %>% 
    mask(aoi) %>% mask(col)
  
  # Define some values
  area_cur <- NA
  area_fut <- NA
  area_overlap <- NA
  
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
        connect <- 1
        habitat <- 1
        big_habitat <- 1
        time_periods <- time_periods[1:null_id]
        
        fut_ranges <- subset(fut_ranges, !sapply(fut_ranges, is.null))
        
        # Current area
        area_cur <- global(
          mask(cellSize(cur_patches, unit = "km"), cur_patches), 
          "sum", na.rm = TRUE)[[1]] / col_area * 100
        
        # Future areas
        area_fut <- sapply(fut_ranges, function(fut_range){
          global(
            mask(cellSize(fut_range, unit = "km"), fut_range), 
            "sum", na.rm = TRUE) / col_area * 100
        }) %>% unlist() %>% min()
        
        # Overlap
        lyrs <- c(cur_patches, do.call(c, fut_ranges))
        lyr <- app(lyrs, fun = function(x) all(!is.na(x)))
        area_overlap <- global(
          cellSize(lyr, unit = "km") * lyr, 
          "sum", na.rm = TRUE)[[1]] / col_area * 100
        
      } else{
        message(sprintf("%s will lose all connection to future in Colombia.", sp))
        connect <- 0
        habitat <- 1
        big_habitat <- 1
        time_periods <- NA
      }
    } else {
      message(sprintf("%s has no valid big enough habitat in Colombia.", sp))
      connect <- 0
      habitat <- 1
      big_habitat <- 0
      time_periods <- NA
    }
  } else{
    message(sprintf("%s has no valid habitat in Colombia.", sp))
    connect <- 0
    habitat <- 0
    big_habitat <- 0
    time_periods <- NA
  }
  
  rst <- data.frame(sp = sp, ssp = scenario,
             habitat = habitat, big_habitat = big_habitat,
             connect = connect, time_periods = time_periods,
             area_current = area_cur, area_future = area_fut,
             area_overlap = area_overlap)
  
  write.table(rst, file.path(code_dir, "data/smr_species.csv"), 
              sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  
  rst
}, mc.cores = 20) %>% bind_rows()

write.csv(smr_species, row.names = FALSE,
          file.path(code_dir, "data/smr_species_final.csv"))
