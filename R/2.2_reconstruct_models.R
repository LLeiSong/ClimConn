# ============================================================
# Script: 2.2_reconstruct_models.R
#
# Purpose:
#   Reconstruct BioModelos MAXENT models for selected species using kuenm:
#     - Reformat occurrences for kuenm
#     - Refit / project models with prepared M/G variables
#     - Build current and future ensembles and binary maps.
#
# Inputs:
#   data/mamiferos/<species>/occurrences/jointID_occ.csv
#   data/mamiferos/<species>/M_variables/Set_1/         M-variable grids
#   data/mamiferos/<species>/ensembles/current/MAXENT/binValues_MAXENT.csv
#   data/Env/G_variable/                                G-variable grids for projections
#
# Outputs (per species):
#   data/mamiferos/<species>/occurrences/occ_joint_kuenm.csv
#   data/mamiferos/<species>/final_models_enmeval_r/    kuenm MAXENT runs and projections
#   data/mamiferos/<species>/reconstruct/current/MAXENT/  current ensembles (cont. + binary)
#   data/mamiferos/<species>/reconstruct/future/MAXENT/   future ensembles by scenario
#   data/mamiferos/<species>/reconstruct/future/MAXENT/future_ensemble/  final future ensembles
#
# Usage:
#   Rscript 2.2_reconstruct_models.R -f "Species_one,Species_two"
#
# Reference:
# This script is constructed highly relying on 
# https://github.com/PEM-Humboldt/biomodelos-sdm.git
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
# devtools::install_github("marlonecobos/kuenm")
library(kuenm)
library(optparse)

# Directories
dat_dir <- "data"
bmls <- file.path(dat_dir, "mamiferos")
env_dir <- file.path(dat_dir, "Env")
G_dir <- file.path(env_dir, "G_variable")

# Project the models
## Define the functions
kuenm.occ <- function(data., spname, foldersp, occname) {
  occ_kuenm <- data.
  occ_kuenm$species <- rep(spname, nrow(occ_kuenm))
  occ_kuenm <- occ_kuenm [, c(3, 1, 2)]
  
  write.csv(
    x = occ_kuenm, file = paste0(foldersp, "/occurrences/", occname, ".csv"), 
    row.names = F)
}

## Ensemble layers from different models
currentEns_byAlg <- function(rasF.Stack, algorithm, foldersp,
                             tim, esc.nm, crs.proyect, bins) {
  # folder to save ensembles
  dir.create(paste0(foldersp, "/reconstruct/", tim, "/", algorithm), 
             showWarnings = F, recursive = T)
  
  # species name with hyphen not dot
  sp_hyphen <- gsub(basename(foldersp), pattern = "\\.", replacement = "_")
  
  esc.nm <- ifelse(tim == "current", "_", paste0("_", esc.nm, "_"))
  
  #----------------------
  # Binary and ensemble
  biomodelos.thresh <- bins
  if (terra::nlyr(rasF.Stack) > 1) {
    Ras.med <- terra::median(rasF.Stack)
  } else {
    Ras.med <- rasF.Stack
  }
  
  # apply thresholds to G layer
  BinsRas <- lapply(X = names(biomodelos.thresh), function(X) {
    Ras.med >= bins[1, X]
  })
  BinsRas <- rast(BinsRas)
  
  if (nlyr(rasF.Stack) > 1) {
    Ras.devstd <- terra::app(x = rasF.Stack, sd)
    Ras.cv <- Ras.devstd / Ras.med
    Ras.sum <- terra::app(x = rasF.Stack, sum)
    
    # stacking results of ensembles
    Resensembles <- c(Ras.med, Ras.devstd, Ras.cv, Ras.sum)
    names(Resensembles) <- c(
      paste0(sp_hyphen, esc.nm, algorithm), 
      paste0(sp_hyphen, "_devstd", esc.nm, algorithm),
      paste0(sp_hyphen, "_CV", esc.nm, algorithm), 
      paste0(sp_hyphen, "_sum", esc.nm, algorithm)
    )
    terra::crs(Resensembles) <- crs.proyect
  } else {
    Resensembles <- rasF.Stack
    names(Resensembles) <- c(
      paste0(sp_hyphen, esc.nm, algorithm)
    )
  }
  
  # writing binaries
  names(BinsRas) <- paste0(sp_hyphen, esc.nm, names(biomodelos.thresh), "_", algorithm)
  
  for (i in 1:nlyr(BinsRas)) {
    BinsRas.i <- BinsRas[[i]]
    terra::crs(BinsRas.i) <- crs.proyect
    terra::writeRaster(BinsRas.i, paste0(
      foldersp, "/reconstruct/", tim, "/", algorithm, "/",
      names(BinsRas[[i]]), ".tif"
    ),
    filetype = "GTiff",
    overwrite = T,
    datatype = "INT2S",
    NAflag = -9999,
    gdal = c("COMPRESS = LZW")
    )
  }
  
  # writing continuos
  
  for (i in 1:terra::nlyr(Resensembles)) {
    RasResen.i <- Resensembles[[i]]
    terra::crs(RasResen.i) <- crs.proyect
    writeRaster(RasResen.i, paste0(
      foldersp, "/reconstruct/", tim, "/", algorithm, "/",
      names(Resensembles[[i]]), ".tif"
    ),
    filetype = "GTiff",
    overwrite = T,
    NAflag = -9999,
    datatype = "FLT4S",
    gdal = c("COMPRESS = LZW")
    )
  }
}

#------------------------
# variation coeficient function
CV <- function(x, na.rm = TRUE) {
  x1 <- sd(x, na.rm = na.rm) / median(x, na.rm = na.rm)
  x2 <- x1 * 100
}

#------------------------------------------
## Slip future layers based on path and year
futAuxiliar <- function(fut.list.ras) {
  listras.lenght <- length(fut.list.ras)
  layers.lenght <- nlyr(fut.list.ras[[1]]) # example of layers length
  
  list.layersFEsc <- list()
  nm_vector <- as.numeric()
  
  for (c in 1:layers.lenght) {
    layersbyEscenario <- terra::rast()
    
    for (d in 1:listras.lenght) {
      stackc <- fut.list.ras[[d]]
      rasc <- stackc[[c]]
      nms <- unlist(regmatches(
        x = names(rasc), regexpr("_", names(rasc)), invert = TRUE))[2]
      nms <- gsub(".asc", "", nms)
      nm_vector[c] <- nms
      layersbyEscenario <- c(layersbyEscenario, rasc)
    }
    
    list.layersFEsc[[c]] <- layersbyEscenario
  }
  names(list.layersFEsc) <- nm_vector
  return(list.ras = list.layersFEsc)
}

# Ensemble future layers of all models
future_ensemble <- function(
    folder_in, algo = "MAXENT", threshold = 10, erase_old = TRUE, 
    ensemble_method = "binary", uncertainty = TRUE, continuos_median = TRUE,
    continuos_binary = TRUE) {
  
  # Define paths and create new folders if necessary
  path_fut <- "reconstruct/future"
  new_folder <- file.path(folder_in, path_fut, algo, "future_ensemble")
  
  if (dir.exists(new_folder) && erase_old) {
    unlink(x = new_folder, recursive = TRUE, force = TRUE)
  }
  
  dir.create(new_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Read and filter raster files
  rasters <- list.files(file.path(folder_in, path_fut, algo), pattern = "*.tif")
  
  if (sum(grepl(pattern = sprintf("_%s_", threshold), x = rasters)) < 2) {
    stop("Not enough rasters to process.")
  }
  
  variability <- grepl(pattern = "_sum_|_devstd_|_CV_", x = rasters)
  rasters <- rasters[!variability]
  
  rasters2 <- gsub(pattern = "_Cor", replacement = "", x = rasters)
  rasters_str <- stringr::str_split(rasters2, "_", simplify = TRUE)
  
  mdl <- unique(rasters_str[, 3])
  yr <- unique(rasters_str[, 4])
  path <- unique(rasters_str[, 5])
  bin <- unique(rasters_str[, 6])
  bin <- sapply(bin, function(x) gsub(pattern = ".tif", replacement = "", x))
  
  index_bin <- which(bin == as.character(threshold) | bin == paste0("X", threshold))
  bin <- bin[index_bin]
  
  # Read threshold values
  bin_table <- read.csv(paste0(folder_in, "/ensembles/current/", 
                               algo, "/binValues_", algo, ".csv"))
  if(threshold != 0){
    bin_value <- bin_table[1, grepl(pattern = threshold, x = names(bin_table))]  
  } else {
    bin_value <- bin_table[1, c(names(bin_table) == 0 | names(bin_table) == "X0")]
  }
  
  # Main Processing Loop
  for (i in seq_along(yr)) {
    rasters_yr_i <- rasters[grepl(pattern = yr[i], x = rasters)]
    
    for (a in seq_along(path)) {
      rasters_yr_path_a <- rasters_yr_i[grepl(pattern = path[a], x = rasters_yr_i)]
      
      if(ensemble_method == "continuos"){
        pattern_con <- paste0(path[a], "_", algo)
        rasters_yr_path_a <- rasters_yr_path_a[
          grepl(pattern = pattern_con, x = rasters_yr_path_a)]
        raster_data <- file.path(folder_in, path_fut, 
                                 algo, rasters_yr_path_a) %>% terra::rast()
        
        # Calculate standard deviation
        if (uncertainty) {
          ens_devStd <- round(terra::app(raster_data, sd), 3)
          writeRaster(ens_devStd, 
                      paste0(new_folder, "/", yr[i], "_", 
                             path[a], "_devstd_maxent.tif"), 
                      NAflag = -9999, overwrite = TRUE, datatype = "FLT4S")
        }
        
        # Calculate median
        if (continuos_median) {
          ens_median <- round(median(raster_data), 3)
          names(ens_median) <- paste0(folder_in, "_", yr[i], "_", path[a], "_maxent.tif")
          writeRaster(ens_median, paste0(new_folder, "/", yr[i], "_", path[a], "_maxent.tif"), 
                      NAflag = -9999, overwrite = TRUE, datatype = "FLT4S")
        }
        
        if (continuos_binary) {
          binary <- (ens_median >= bin_value) * 1
          binary <- terra::classify(binary, cbind(0, NA), cbind(1, 1))
          writeRaster(binary, paste0(new_folder, "/", yr[i], "_", path[a], "_", 
                                     threshold, ".tif"), 
                      datatype = "INT1U", NAflag = 255, overwrite = TRUE)
        }
      }
      
      if(ensemble_method == "binary"){
        for (b in seq_along(bin)) {
          pattern_bin <- paste0("_", bin[b], ".tif|_", bin[b], "_")
          rasters_yr_path_bin_b <- rasters_yr_path_a[
            grepl(pattern = pattern_bin, x = rasters_yr_path_a)]
          raster_data <- terra::rast(
            file.path(folder_in, path_fut, algo, rasters_yr_path_bin_b))
          
          ens <- sum(raster_data)
          ens_values <- unique(ens)
          index_max <- which(ens_values == max(ens_values))
          m <- matrix(data = c(t(ens_values), rep(0, nrow(ens_values))), 
                      nrow = nrow(ens_values),
                      ncol = 2,
                      byrow = FALSE)
          m[index_max, 2] <- 1
          binary <- terra::classify(ens, m)
          
          writeRaster(binary, paste0(
            new_folder, "/", yr[i], "_", path[a], "_", threshold, ".tif"), 
            datatype = "INT1U", NAflag = 255, overwrite = TRUE)
          
          if (uncertainty) {
            ens_sum <- sum(raster_data)
            writeRaster(
              ens_sum, paste0(new_folder, "/", yr[i], "_", 
                              path[a], "_sum_", bin[b], ".tif"), 
                        datatype = "INT1U", NAflag = 255, overwrite = TRUE)
          }
        }
      }
    }
  }
}

##### Start the process ########
# Parse inline parameters
option_list <- list(
  make_option(c("-f", "--folder_sps"), 
              action = "store", type = 'character',
              help = "The folder of species to process."))
opt <- parse_args(OptionParser(option_list = option_list))
folder_sps <- opt$folder_sps
folder_sps <- strsplit(folder_sps, ",")[[1]]

for (folder_sp in folder_sps){
  folder_sp <- file.path(bmls, folder_sp)
  occ <- read.csv(file.path(folder_sp, "/occurrences/jointID_occ.csv"))
  kuenm.occ(occ[, -c(1, 2)], spname = basename(folder_sp), 
            foldersp = folder_sp, occname = "occ_joint_kuenm")
  
  env_Mdir <- file.path(folder_sp, "M_variables")
  
  kuenm_mod(
    occ.joint = file.path(folder_sp, "occurrences/occ_joint_kuenm.csv"),
    M.var.dir = env_Mdir, out.eval = file.path(folder_sp, "eval_results_enmeval"),
    batch = file.path(folder_sp, "final_models"), rep.n = 1, rep.type = "Bootstrap",
    jackknife = FALSE, max.memory = 2000, out.format = "cloglog", project = TRUE, 
    out.dir = file.path(folder_sp, "final_models_enmeval_r"),
    G.var.dir = G_dir, ext.type = "no_ext", write.mess = FALSE,
    write.clamp = FALSE, maxent.path = here("java"), 
    args = c(NULL, paste0("maximumbackground=", 10000)), 
    wait = TRUE, run = TRUE)
  
  # Calculate ensembles
  ## Current
  current_M_files <- list.files(
    path = file.path(folder_sp, "final_models_enmeval_r"), 
    pattern = paste0(basename(folder_sp), "_M.asc"), 
    full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
  current_M_proj <- rast(current_M_files)
  names(current_M_proj) <- str_extract(current_M_files, "M_.*Set_1")
  
  bins <- read.csv(
    file.path(folder_sp, "ensembles/current/MAXENT/binValues_MAXENT.csv"))
  names(bins) <- c("E", "0", "10", "20", "30")
  
  enscurr <- currentEns_byAlg(
    rasF.Stack = current_M_proj, 
    algorithm = "MAXENT", foldersp = folder_sp, 
    tim = "current", esc.nm = "",
    crs.proyect = "+proj=longlat +datum=WGS84 +no_defs +type=crs", 
    bins = bins)
  
  ## Future
  env.folder <- list.dirs(file.path(folder_sp, "final_models_enmeval_r"), 
                          full.names = T, recursive = F)
  env.folderNm <- list.dirs(file.path(folder_sp, "final_models_enmeval_r"), 
                            full.names = F, recursive = F)
  
  fut_proj_list <- list()
  for (i in 1:length(env.folder)) {
    # folder of layers of each model
    env.listFolder <- list.files(
      env.folder[i], pattern = ".asc$", 
      full.names = T, recursive = T, include.dirs = F)
    env.listnms <- list.files(
      env.folder[i], pattern = ".asc$", 
      full.names = F, recursive = T, include.dirs = F)
    
    # no future asc files
    noFRas <- c(grep("*_M.asc$", env.listFolder), 
                grep("*_G.asc$", env.listFolder), 
                grep(paste0(basename(folder_sp), ".asc$"), 
                     env.listFolder))
    
    # removing no future files
    env.FlistRas <- env.listFolder[-noFRas]
    fut_proj <- terra::rast(env.FlistRas)
    names(fut_proj) <- env.listnms[-noFRas]
    
    fut_proj_list[[i]] <- fut_proj
  }
  names(fut_proj_list) <- env.folderNm
  
  layersF <- futAuxiliar(fut.list.ras = fut_proj_list)
  
  for (f in 1:length(layersF)) {
    enscurr <- currentEns_byAlg(
      rasF.Stack = layersF[[f]], 
      algorithm = "MAXENT", foldersp = folder_sp, 
      tim = "future", esc.nm = names(layersF[f]),
      crs.proyect = "+proj=longlat +datum=WGS84 +no_defs +type=crs", 
      bins = bins)
  }
  
  future_ensemble(folder_in = folder_sp, ensemble_method = "continuos")
}
