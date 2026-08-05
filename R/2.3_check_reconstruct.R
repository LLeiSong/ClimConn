# ============================================================
# Script: 2.3_check_reconstruct.R
#
# Purpose:
#   Check completeness and consistency of reconstructed MAXENT results:
#     - Verify M_variables, kuenm/MaxEnt projections, and ensembles
#     - Summarize missing outputs per species
#     - Compare original vs reconstructed suitability surfaces
#     - Clean temporary variables and intermediate prediction files.
#
# Inputs:
#   data/mamiferos/<species>/      reconstructed outputs from 2.2_reconstruct_models.R
#
# Outputs:
#   data/species_miss_var.csv              species missing M_variables
#   data/species_miss_maxent.csv           species missing MaxEnt projections
#   data/species_miss_current.csv          species missing current ensembles
#   data/species_miss_future.csv           species missing future ensembles
#   data/species_miss_future_ensemble.csv  species missing future_ensemble
#   data/species_miss.csv                  union of all missing species
#   data/suit_cors.csv                     correlation of original vs reconstructed suitability
#
# ============================================================

# Load libraries
library(dplyr)
library(here)
library(stringr)
library(terra)

# Directories
dat_dir <- "data"
bmls <- file.path(dat_dir, "mamiferos")

# Parameters to reference
var_list <- c(paste("bio", c(1:7, 10:17), sep = "_"), 
              "wc21elev_s", "wc21slope")
mods <- c("EC-Earth3-Veg", "GFDL-ESM4", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR")
years <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
ssps <- c("ssp126", "ssp370", "ssp585")

## Check the variables in M_variables folder
fnames <- list.files(bmls, full.names = TRUE)
sp_ms_var <- lapply(fnames, function(fname){
  var_dir <- file.path(fname, "M_variables/Set_1")
  files <- file.path(var_dir, sprintf("%s.asc", var_list))
  if (!all(file.exists(files))) basename(fname)
}) %>% unlist()

if (length(sp_ms_var) > 0){
  sp_ms_var <- data.frame(sp = sp_ms_var)
  write.csv(sp_ms_var, file.path(dat_dir, "species_miss_var.csv"),
            row.names = FALSE)
} else{
  message("All species have good variables.")
}

## Check the MAXENT results in folder final_models_enmeval_r
sp_ms_maxent <- lapply(fnames, function(fname){
  mxt_dir <- file.path(fname, "final_models_enmeval_r")
  mod_dir <- list.files(mxt_dir, full.names = TRUE)
  
  if (length(mod_dir) > 0){
    sp_ms_mxt <- lapply(mod_dir, function(mod){
      # All possible files
      files <- tidyr::crossing(mods, years, ssps) %>% 
        mutate(fname = paste(mods, years, ssps, sep = "_")) %>% 
        pull(fname)
      
      # Remove files should be missing
      to_move <- sprintf("GFDL-ESM4_%s_ssp585", years)
      files <- files[!files %in% to_move]
      
      # Mosaic to the full path
      files <- paste(basename(fname), files, sep = "_")
      files <- file.path(
        mod, sprintf("%s.asc", c(files, paste(basename(fname), "M", sep = "_"))))
      
      if (!all(file.exists(files))) basename(fname)
    }) %>% unlist()
  } else sp_ms_mxt <- basename(fname)
  
  if (!is.null(sp_ms_mxt)) sp_ms_mxt
}) %>% unlist()

if (length(sp_ms_maxent) > 0){
  sp_ms_maxent <- data.frame(sp = sp_ms_maxent)
  write.csv(sp_ms_maxent, file.path(dat_dir, "species_miss_maxent.csv"),
            row.names = FALSE)
} else{
  message("All species have good MAXENT results.")
}

## Check the ensemble results in folder reconstruct
## both current, future and future/future_ensemble
## Current
sp_ms_cur <- lapply(fnames, function(fname){
  cur_dir <- file.path(fname, "reconstruct/current/MAXENT")
  sp <- gsub("\\.", "_", basename(fname))
  
  files <- file.path(
    cur_dir, sprintf("%s_%s_MAXENT.tif", sp, c(0, 10, 20, 30, "E")))
  files <- c(files, file.path(cur_dir, sprintf("%s_MAXENT.tif", sp)))
  
  if (!all(file.exists(files))) basename(fname)
}) %>% unlist()

if (length(sp_ms_cur) > 0){
  sp_ms_cur <- data.frame(sp = sp_ms_cur)
  write.csv(sp_ms_cur, file.path(dat_dir, "species_miss_current.csv"),
            row.names = FALSE)
} else{
  message("All species have good current results.")
}

## Future
sp_ms_fut <- lapply(fnames, function(fname){
  fut_dir <- file.path(fname, "reconstruct/future/MAXENT")
  sp <- gsub("\\.", "_", basename(fname))
  
  # All possible files
  files <- tidyr::crossing(
    mods, years, ssps, thred = c(0, 10, 20, 30, "E")) %>% 
    filter(!(mods == "GFDL-ESM4" & ssps == "ssp585")) %>% 
    mutate(fname = paste(mods, years, ssps, thred, sep = "_")) %>% 
    pull(fname)
  
  # Mosaic to the full path
  files <- file.path(
    fut_dir, sprintf("%s_%s_MAXENT.tif", sp, files))
  
  if (!all(file.exists(files))) basename(fname)
}) %>% unlist()

if (length(sp_ms_fut) > 0){
  sp_ms_fut <- data.frame(sp = sp_ms_fut)
  write.csv(sp_ms_fut, file.path(dat_dir, "species_miss_future.csv"),
            row.names = FALSE)
} else{
  message("All species have good future results.")
}

## Future ensemble
sp_ms_fut_esb <- lapply(fnames, function(fname){
  fut_dir <- file.path(fname, "reconstruct/future/MAXENT/future_ensemble")
  
  # All possible files
  files <- tidyr::crossing(years, ssps, item = c("10", "maxent")) %>% 
    mutate(fname = paste(years, ssps, item, sep = "_")) %>% 
    pull(fname)
  
  # Mosaic to the full path
  files <- file.path(fut_dir, sprintf("%s.tif", files))
  
  if (!all(file.exists(files))) basename(fname)
}) %>% unlist()

if (length(sp_ms_fut_esb) > 0){
  sp_ms_fut_esb <- data.frame(sp = sp_ms_fut_esb)
  write.csv(
    sp_ms_fut_esb, file.path(dat_dir, "species_miss_future_ensemble.csv"),
    row.names = FALSE)
} else{
  message("All species have good future ensemble results.")
}

# Put everything together
sps_miss <- unique(c(sp_ms_maxent$sp, sp_ms_cur$sp, 
                     sp_ms_fut$sp, sp_ms_fut_esb$sp))
if (length(sps_miss) > 0){
  message(sprintf("%s species are missing.", length(sps_miss)))
  sps_miss <- data.frame(sp = sps_miss)
  write.csv(
    sps_miss, file.path(dat_dir, "species_miss.csv"),
    row.names = FALSE)
} else{
  message("All species have good results.")
}

# Check the similarity of the results
fnames <- list.files(bmls, full.names = TRUE)
cors <- lapply(fnames, function(fname){
  sp <- basename(fname)
  cur_org <- rast(
    file.path(fname, "ensembles/current/MAXENT",
              sprintf("%s_MAXENT.tif", gsub("\\.", "_", basename(fname)))))
  cur_rct <- rast(
    file.path(fname, "reconstruct/current/MAXENT",
              sprintf("%s_MAXENT.tif", gsub("\\.", "_", basename(fname)))))
  cur_rct <- mask(crop(cur_rct, cur_org), cur_org)
  
  r <- min(layerCor(c(cur_org, cur_rct), "cor", use = "complete.obs")$correlation)
  data.frame(sp = sp, cor = r)
}) %>% bind_rows()

write.csv(cors, file.path(dat_dir, "suit_cors.csv"), row.names = FALSE)

# Remove temporary files for reconstruct
for (fname in fnames){
  # Remove variables
  unlink(file.path(fname, "M_variables"), recursive = TRUE)
  
  # Remove independent future predictions
  preds <- list.files(file.path(fname, "reconstruct/future/MAXENT"),
                      pattern = ".tif", full.names = TRUE)
  file.remove(preds)
  
  # Remove original independent future predictions
  preds <- list.files(file.path(fname, "final_models_enmeval_r"),
                      pattern = sprintf("%s_.+.asc", basename(fname)), 
                      full.names = TRUE, recursive = TRUE)
  file.remove(preds)
  
  # Remove the ensemble folder
  unlink(file.path(fname, "ensemble"), recursive = TRUE)
}
