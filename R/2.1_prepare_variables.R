# ============================================================
# Script: 2.1_prepare_variables.R
#
# Purpose:
#   Prepare environmental predictor variables for BioModelos
#   MaxEnt/enmeval/kuenm workflow:
#     - Parse MaxEnt logs to get variable sets and GCMs
#     - Download WorldClim v2.1 current and CMIP6 future layers
#     - Clip to Colombia / study region
#     - Organize files into MaxEnt-ready M_variable and G_variable folders.
#
# Inputs:
#   data/mamiferos/                 species folders with final_models_enmeval and interest_areas/shape_M.shp
#   data/colombia.geojson           study area boundary
#   WorldClim v2.1 + CMIP6 (online) climate and elevation layers
#
# Outputs:
#   data/Env/                       clipped and stacked baseline/future rasters
#   data/Env/M_variable/vars.tif    training stack
#   data/Env/G_variable/...         scenario-specific bio_*.asc layers
#   data/mamiferos/*/M_variables/   per-species M-variable ASCII grids
#
# Author: Lei Song <lei.song@rutgers.edu>
# Last updated: 2025-11-21
# ============================================================

# Load libraries
library(sf)
library(dplyr)
library(terra)
library(here)
library(stringr)
library(tools)
library(rnaturalearth)
library(RCurl)

dat_dir <- here("data")

# Get variable set and names for modeling and projection
bmls <- file.path(dat_dir, "mamiferos")
fnames <- list.files(bmls, full.names = TRUE)

# Check model package
use_kuenm <- lapply(fnames, function(fname){
  log <- readLines(file.path(fname, "log_file.txt"))
  if (length(log[str_detect(log, "kuenm")]) > 0) fname
}) %>% unlist()
if (is.null(use_kuenm)) message("No kuenm used for tuning these models.")

# Get the model catalog
models_info <- lapply(fnames, function(fname){
  mods <- list.files(file.path(fname, "final_models_enmeval"),
                     full.names = TRUE)
  
  # just one model has everything.
  for (mod in mods){
    log <- readLines(file.path(mod, "maxent.log"))
    
    # Get var set and vars for modeling
    m_vars <- log[str_detect(log, "Checking header of ")]
    m_vars <- str_extract(m_vars, "[A-Z]:.*M_variables.*.asc")
    m_vars <- gsub("\\\\", .Platform$file.sep, m_vars)
    m_set <- unique(basename(dirname(m_vars)))
    m_vars <- file_path_sans_ext(basename(m_vars))
    
    # Get var set and models for projection
    p_vars <- log[str_detect(log, "Converting file .*G_variables.*")]
    if (length(p_vars) == 0){
      next
    } else{
      p_vars <- str_extract(p_vars, "[A-Z]:.*G_variables.*.asc")
      p_vars <- gsub("\\\\", .Platform$file.sep, p_vars)
      p_set <- unique(basename(dirname(dirname(p_vars))))
      p_mods <- unique(basename(dirname(p_vars)))
      p_mods <- p_mods[p_mods != "M"]
      
      # Return
      statistics <- list(vars = m_vars, m_set = m_set, 
                         p_set = p_set, p_mods = p_mods)
    }
  }
  
  statistics
})

# Get the list of variables
var_list <- lapply(models_info, function(sts) sts$vars)
if (length(unique(var_list)) == 1) var_list <- var_list[[1]]

# Get the var set
m_set <- lapply(models_info, function(sts) sts$m_set)
if (length(unique(m_set)) == 1) message(sprintf("M set is %s", m_set[[1]]))

p_set <- lapply(models_info, function(sts) sts$p_set)
if (length(unique(p_set)) == 1) message(sprintf("P set is %s", p_set[[1]]))

# Get the Worldclim projection models
p_mods <- lapply(models_info, function(sts) sts$p_mods)
if (length(unique(p_mods)) == 1) p_mods <- p_mods[[1]]
p_mods <- lapply(p_mods, function(x) strsplit(x, "_")[[1]][1]) %>% 
  unlist() %>% unique()
p_mods[p_mods == "Earth3-Veg"] = "EC-Earth3-Veg"

# Start to download the Worldclim layers
## Create folder structure
env_dir <- file.path(dat_dir, "Env")
if (!dir.exists(env_dir)) dir.create(env_dir)

M_dir <- file.path(env_dir, "M_variable")
if (!dir.exists(M_dir)) dir.create(M_dir)

G_dir <- file.path(env_dir, "G_variable", p_set[[1]])
if (!dir.exists(G_dir)) dir.create(G_dir, recursive = TRUE)

## Define a shape for the layers for projection
sf_use_s2(FALSE)
col <- st_read(file.path(dat_dir, "colombia.geojson")) %>% 
  st_bbox() %>% st_as_sfc() %>% 
  st_buffer(11, endCapStyle = "SQUARE") %>% 
  st_bbox() %>% st_as_sfc()

sa <- ne_countries(scale = 10) %>% 
  filter(continent %in% c("South America", "North America")) %>% 
  st_union() %>% st_cast("POLYGON") %>% st_as_sf() %>% 
  mutate(area = st_area(.)) %>% filter(area == max(area)) %>% select(-area)
col <- st_intersection(col, sa) %>% vect()

## Reformat the items to download
years <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
ssps <- c("ssp126", "ssp370", "ssp585")
url_cur <- "https://geodata.ucdavis.edu/climate/worldclim/2_1/base"
url_proj <- "https://geodata.ucdavis.edu/cmip6/30s"

## Start
options_orig <- options()
options(timeout = 60 * 60)

## - current DEM and bio
items <- c("wc2.1_30s_elev", "wc2.1_30s_bio")
for (item in items){
  # Download
  destfile <- file.path(env_dir, sprintf("%s.zip", item))
  tryCatch({
    download.file(file.path(url_cur, sprintf("%s.zip", item)), 
                  destfile = destfile)
  }, error = function(msg){
    Sys.sleep(10)
    download.file(file.path(url_cur, sprintf("%s.zip", item)), 
                  destfile = destfile)
  })
  
  # Extract the files
  temp_dir <- file.path(env_dir, "temp")
  unzip(destfile, exdir = temp_dir)
  file.remove(destfile)
  
  # Clip
  lyrs <- rast(list.files(temp_dir, recursive = TRUE, full.names = TRUE))
  writeRaster(lyrs, file.path(env_dir, sprintf("%s.tif", item)))
  
  lyrs <- crop(lyrs, col)
  
  fname <- file.path(
    env_dir, sprintf("%s.tif", gsub("wc2.1_30s_", "", item)))
  writeRaster(lyrs, fname)
  unlink(temp_dir, recursive = TRUE)
}

for(mod in p_mods){
  for (yr in years){
    for (ssp in ssps){
      fname <- sprintf("wc2.1_30s_bioc_%s_%s_%s.tif", mod, ssp, yr)
      dl_url <- file.path(url_proj, mod, ssp, fname)
      fname_to <- file.path(env_dir, gsub("wc2.1_30s_bioc_", "", fname))
      
      if(url.exists(dl_url) & !file.exists(fname_to)){
        destfile <- file.path(env_dir, fname)
        tryCatch({
          download.file(dl_url, destfile = destfile)
        }, error = function(msg){
          Sys.sleep(10)
          download.file(dl_url, destfile = destfile)
        })
        
        # Clip
        lyrs <- rast(destfile)
        lyrs <- crop(lyrs, col)
        
        writeRaster(lyrs, fname_to)
        file.remove(destfile)
      }
    }
  }
}

# Reverse back to default options
options(options_orig)

# Reorganize files and save as MaxEnt ready format
## Define functions

## Get the number of decimals
## https://github.com/PEM-Humboldt/biomodelos-sdm/modelling/
## R/process_env_variables.R#L438
peel_vars <- function(env_Fstack, dst_dir, dnum = 3){
  # writing future layers
  for (d in 1:nlyr(env_Fstack)) {
    env_Fd <- env_Fstack[[d]]
    
    vars_int <- c(sprintf("bio_%s", c(12:14, 16, 17)), "wc21elev_s")
    
    if (names(env_Fd) %in% vars_int) {
      datTyp <- "INT2S"
      decinum <- 0
    } else {
      datTyp <- "FLT4S"
      decinum <- dnum
    }
    
    writeRaster(
      x = round(env_Fd, digits = decinum),
      filename = file.path(dst_dir, sprintf("%s.asc", names(env_Fd))),
      overwrite = T,
      NAflag = -9999,
      datatype = datTyp)
  }
}

## Save the current layer with larger extent for training
## This layer will be crop to interest_area for each species later
aoi <- do.call(rbind, lapply(fnames, function(fname){
  st_read(file.path(fname, "interest_areas/shape_M.shp"), quiet = TRUE)
})) %>% st_union() %>% vect() %>% buffer(1000)

ele <- rast(file.path(env_dir, "wc2.1_30s_elev.tif"))
ele <- mask(crop(ele, aoi), aoi)
ele <- c(ele, terra::terrain(ele, v = "slope"))
names(ele) <- var_list[16:17]
vars <- rast(file.path(env_dir, "wc2.1_30s_bio.tif"))
names(vars) <- str_extract(names(vars), "bio_[0-9]+")
vars <- mask(crop(vars, aoi), aoi)
vars <- subset(c(vars, ele), var_list)

writeRaster(vars, file.path(M_dir, "vars.tif"))
rm(aoi, ele, vars)

## Do G_variable first
## file structure should be like: Set_1\Earth3-Veg_2021-2040_ssp585\bio_1.asc
## 1. Prepare terrain vars
ele <- rast(file.path(env_dir, "elev.tif"))
ele <- c(ele, terra::terrain(ele, v = "slope"))
names(ele) <- var_list[16:17]

p_mods <- c("M", p_mods)
for(mod in p_mods){
  ## 2. current vars
  if (mod == "M"){
    # Load and subset the layers
    fname <- file.path(env_dir, "bio.tif")
    lyrs <- rast(fname); names(lyrs) <- str_extract(names(lyrs), "bio_[0-9]+")
    lyrs <- subset(c(lyrs, ele), var_list)
    
    # Peel variables into separate asc files
    dst_dir <- file.path(G_dir, "M"); dir.create(dst_dir, showWarnings = FALSE)
    peel_vars(lyrs, dst_dir)
  ## 3. future vars
  } else{
    for (yr in years){
      for (ssp in ssps){
        # Load and subset the layers
        fname <- file.path(env_dir, sprintf("%s_%s_%s.tif", mod, ssp, yr))
        if (file.exists(fname)){
          lyrs <- rast(fname)
          names(lyrs) <- paste("bio", str_extract(names(lyrs), "[0-9]+$"), sep = "_")
          lyrs <- subset(c(lyrs, ele), var_list)
          
          # Peel variables into separate asc files
          dst_dir <- file.path(G_dir, sprintf("%s_%s_%s", mod, yr, ssp))
          dir.create(dst_dir, showWarnings = FALSE)
          peel_vars(lyrs, dst_dir)
        }
      }
    }
  }
}

## Do M_variable for modeling
vars <- rast(file.path(M_dir, "vars.tif"))
fnames <- list.files(bmls, full.names = TRUE)
for (fname in fnames){
  # Load data
  msk <- rast(file.path(
    fname, "ensembles/current/MAXENT/",
    sprintf("%s_MAXENT.tif", gsub("\\.", "_", basename(fname)))))
  vars_s <- crop(vars, msk, mask = TRUE) %>% trim()
  
  # Peel variables into separate asc files
  dst_dir <- file.path(fname, "M_variables/Set_1")
  if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE)
  peel_vars(vars_s, dst_dir)
}

# Clean up
rm(list = ls())

# So far, these variables are ready for projection.
