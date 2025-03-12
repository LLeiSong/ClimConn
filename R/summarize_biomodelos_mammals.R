# Load libraries
library(sf)
library(terra)
library(dplyr)
library(parallel)
library(pbmcapply)
library(ggplot2)
library(tidyterra)
library(stringr)
sf_use_s2(FALSE)

# Set directories and paths
root_dir <- here::here()
crs_analysis <- "ESRI:54030"
data_dir <- file.path(root_dir, "data/mamiferos")
dispersal_path <- file.path(
  root_dir, "data/dispersal", "species_dispersal_rate.csv")
species_list <- list.files(data_dir)
cur_path <- "ensembles/current/MAXENT/%s_10_MAXENT.tif"
fut_paths <- sprintf(
  "%s/2021-2040_ssp%s_10.tif", 
  "ensembles/future/MAXENT/future_ensemble", c("126", "585"))
aoi_path <- "interest_areas/shape_M.shp"

# Template for csv
fname <- file.path(root_dir, "data/sp_numbers.csv")
columns <- c("species", "ssp", "aoi_area", "cur_area", "fut_area", "IoU")
nums <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(nums) <- columns
write.csv(nums, fname, row.names = FALSE)

# Calculate the numbers for filtering
sp_numbers <- pbmclapply(species_list, function(sp){
  cur_fname <- file.path(data_dir, sp, sprintf(cur_path, gsub("\\.", "_", sp)))
  fut_fnames <- file.path(data_dir, sp, fut_paths)
  aoi_fname <- file.path(data_dir, sp, aoi_path)
  
  if (all(file.exists(cur_fname, fut_fnames, aoi_fname))){
    cur_range <- rast(cur_fname) %>% project(crs_analysis)
    
    # Remove patches smaller than 100 square km
    cur_patches <- patches(cur_range, zeroAsNA = TRUE)
    rz <- zonal(cellSize(cur_patches, unit = "km"), 
                cur_patches, sum, as.raster = TRUE)
    cur_patches <- ifel(rz < 100, NA, cur_patches)
    
    # Get potential dispersal distance
    dispersal_rate <- read.csv(dispersal_path) %>% 
      filter(species == gsub("\\.", " ", sp)) %>% pull(dispersal_rate)
    if (length(dispersal_rate) == 0) {
      dispersal_rate <- read.csv(dispersal_path) %>% 
        pull(dispersal_rate) %>% median()
    }
    dis_dist <- dispersal_rate * 40 * 1000 # in meter
    
    # Read AOI
    aoi <- st_read(aoi_fname) %>% st_transform(crs_analysis)
    aoi_spread <- st_buffer(aoi, dis_dist)
    
    for (fut_fname in fut_fnames){
      ssp <- str_extract(basename(fut_fname), "ssp[0-9]{3}")
      
      fut_range <- mask(rast(fut_fname) %>% project(crs_analysis), aoi_spread)
      
      # Remove patches smaller than 100 square km
      fut_patches <- patches(fut_range, zeroAsNA = TRUE)
      rz <- zonal(cellSize(fut_patches, unit = "km"), 
                  fut_patches, sum, as.raster = TRUE)
      fut_patches <- ifel(rz < 100, NA, fut_patches)
      
      aoi_area <- st_area(aoi) %>% units::set_units("km2") %>% sum()
      cur_area <- global(mask(cellSize(cur_patches, unit = "km"), cur_patches), 
                         "sum", na.rm = TRUE)[[1]]
      fut_area <- global(mask(cellSize(fut_patches, unit = "km"), fut_patches), 
                         "sum", na.rm = TRUE)[[1]]
      
      # Calculate IoU to indicate the overlap between current and future
      a_overlap <- mask(cur_patches, fut_patches)
      a_union <- sum(cur_patches, fut_patches, na.rm = TRUE)
      
      IoU <- (global(mask(cellSize(a_overlap, unit = "km"), a_overlap), 
                     "sum", na.rm = TRUE)[[1]] / 
                global(mask(cellSize(a_union, unit = "km"), a_union), 
                       "sum", na.rm = TRUE)[[1]]) * 100
      
      num_row <- data.frame(species = sp, ssp = ssp, aoi_area = aoi_area,
                            cur_area = cur_area, fut_area = fut_area,
                            IoU = IoU)
      
      write.table(num_row, fname, sep = ",", row.names = FALSE, 
                  col.names = FALSE, append = TRUE)
    }
  }
}, mc.cores = 6)
