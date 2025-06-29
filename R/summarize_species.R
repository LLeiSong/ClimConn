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
results <- lapply(ssps, function(ssp){
  lapply(sps, function(sp){
    fnames <- list.files(
      file.path(result_dir, ssp, sp), 
      pattern = "cum_currmap.tif", recursive = TRUE)
    
    fnum <- sum(dir.exists(
      file.path(result_dir, ssp, sp, years)))
    
    fnum_s <- sum(dir.exists(
      file.path(result_dir, ssp, sp, sprintf("%s_1", years))))
    
    fnum_good <- sum(
      file.exists(file.path(result_dir, ssp, sp, years, 
                            "cum_currmap.tif")))
    
    fnum_abrupt <- sum(
      file.exists(file.path(result_dir, ssp, sp, sprintf("%s_1", years), 
                            "cum_currmap.tif")))
    
    data.frame(
      sp = sp,
      fnum = length(fnames),
      abrupt = ifelse(fnum_s != 0 & fnum_abrupt != 0, 1, 0),
      all = ifelse(fnum == fnum_good, 1, 0),
      done = ifelse(fnum == length(fnames) & fnum == 4 & fnum_s == 0, 1, 0))
  }) %>% bind_rows() %>% mutate(ssp = ssp)
}) %>% bind_rows()

sp_list <- results %>% filter(done == 0)

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
        connect <- 1
        habitat <- 1
        big_habitat <- 1
        time_periods <- time_periods[1:null_id]
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
  
  data.frame(sp = sp, ssp = scenario,
             habitat = habitat, big_habitat = big_habitat,
             connect = connect, time_periods = time_periods)
}, mc.cores = 3) %>% bind_rows()

write.csv(smr_species, row.names = FALSE,
          file.path(code_dir, "data/smr_species.csv"))
