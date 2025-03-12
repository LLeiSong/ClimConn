# Load libraries
library(readxl)
library(dplyr)
library(sf)
library(terra)
library(optparse)

# Parse inline parameters
option_list <- list(
  make_option(c("-s", "--sp"), 
              action = "store", type = 'character',
              help = "The species to process."))
opt <- parse_args(OptionParser(option_list = option_list))
sp <- opt$sp

# Define paths and parameters
root_dir <- "/home/lsong/ClimConn"
julia_home <- "/Users/leisong/.juliaup/bin"
crs_analysis <- "ESRI:54030"
data_dir <- file.path(root_dir, "data/mamiferos")
cur_path <- "ensembles/current/MAXENT/%s_10_MAXENT.tif"
cur_suit_path <- "ensembles/current/MAXENT/%s_MAXENT.tif"
fut_paths <- sprintf(
  "%s/2021-2040_ssp%s_10.tif", 
  "ensembles/future/MAXENT/future_ensemble", c("126", "585"))
fut_suit_path <- "ensembles/future/MAXENT"
aoi_path <- "interest_areas/shape_M.shp"
dispersal_path <- file.path(
  root_dir, "data/dispersal", "species_dispersal_rate.csv")

# Load function
source(file.path(root_dir, "R/biomodelos_connect_omni.R"))

# Load files 
col <- st_read(file.path(root_dir, "data/colombia.geojson")) %>% 
  st_transform(crs_analysis)
col <- vect(col) %>% ext()

time_periods <- c("1970-2000", "2021-2040")
cur_fname <- file.path(data_dir, sp, sprintf(cur_path, gsub("\\.", "_", sp)))
cur_suit_fn <- file.path(data_dir, sp, sprintf(cur_suit_path, gsub("\\.", "_", sp)))

fut_fname <- file.path(data_dir, sp, fut_paths[grep(scenario, fut_paths)])
fut_suit_fn <- file.path(data_dir, sp, fut_suit_path)
aoi_fname <- file.path(data_dir, sp, aoi_path)

# Get potential dispersal distance
dispersal_rate <- read.csv(dispersal_path) %>% 
  filter(species == gsub("\\.", " ", sp)) %>% pull(dispersal_rate)
if (length(dispersal_rate) == 0) {
  dispersal_rate <- read.csv(dispersal_path) %>% 
    pull(dispersal_rate) %>% median()
}
dis_dist <- dispersal_rate * 40 * 1000 # in meter

aoi <- st_read(aoi_fname) %>% st_transform(crs_analysis)
aoi_spread <- st_buffer(aoi, dis_dist)

# Remove patches smaller than 100 square km
cur_range <- rast(cur_fname) %>% project(crs_analysis)
cur_patches <- patches(cur_range, zeroAsNA = TRUE)
rz <- zonal(cellSize(cur_patches, unit = "km"), 
            cur_patches, sum, as.raster = TRUE)
cur_patches <- ifel(rz < 100, NA, cur_patches)

fut_range <- mask(rast(fut_fname) %>% project(crs_analysis), aoi_spread)
# Remove patches smaller than 100 square km
fut_patches <- patches(fut_range, zeroAsNA = TRUE)
rz <- zonal(cellSize(fut_patches, unit = "km"), 
            fut_patches, sum, as.raster = TRUE)
fut_patches <- ifel(rz < 100, NA, fut_patches)

hbts <- c(cur_range, fut_range)

# NOTE: cannot believe they only save the AOI area.
cur_suit <- rast(cur_suit_fn) %>% project(crs_analysis)
fut_suit <- rast(list.files(
  fut_suit_fn, pattern = sprintf("2021-2040_%s_MAXENT.tif", scenario), 
  full.name = TRUE)) %>% median(na.rm = TRUE) %>% project(crs_analysis)
suits <- c(cur_suit, fut_suit)

lyrs <- c(hbts, suits) %>% trim()

hbts <- subset(lyrs, 1:2)
hbts[hbts == 0] <- NA
names(hbts) <- time_periods

suits <- subset(lyrs, 3:4)
names(suits) <- time_periods

dis_dist <- ceiling(dis_dist / res(suits)[1])

circ_dir <- file.path(root_dir, "results", "omni")
if (!dir.exists(circ_dir)) dir.create(circ_dir)

config_template <- file.path(
  root_dir, "data/config", "omniscape_setting_template.ini")

clim_connect(
  sp, dis_dist, suits, hbts, time_periods, 
  config_template, circ_dir, julia_home, TRUE)
