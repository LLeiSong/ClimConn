library(ini)
library(JuliaCall)
library(terra)
library(sf)
library(dplyr)

clim_connect <- function(sp,
                         dispersal_dist,
                         suit_list,
                         habitat_list,
                         time_periods,
                         config_template,
                         work_dir = ".",
                         julia_home = NULL,
                         verbose = FALSE){
  # Create working directory
  sp_dir <- file.path(circ_dir, sp)
  if (!dir.exists(sp_dir)) dir.create(sp_dir)
  
  # Check time periods and layers
  if ((!nlyr(suit_list) == length(time_periods)) |
      (!nlyr(habitat_list) == length(time_periods))){
    stop("Different layer numbers for time periods.")
  } else {
    names(suit_list) <- time_periods
    names(habitat_list) <- time_periods
  }
  
  # set up Julia
  if (is.null(julia_home)){
    stop("No speficied Julia. Suggest to install Julia and Omniscape with the right version.")
  } else {
    # julia_home <- "/Users/leisong/.juliaup/bin"
    julia_setup(JULIA_HOME = julia_home)
  }
  
  # Load Omniscape
  julia_library("Omniscape")
  
  # Let Omniscape shut up if request
  if (!verbose){
    julia_library("Logging")
    julia_command("Logging.disable_logging(Logging.Info)", show_value = FALSE)
  }
  
  config <- read.ini(config_template)
  config$Options$radius <- dispersal_dist
  config$Options$buffer <- ceiling(dispersal_dist / 2)
  # config$Options$block_size <- 5
  
  # Process the layers for omniscape
  i <- 1
  # Set names for runs
  nm <- time_periods[[i]]
  ground_name <- file.path(sp_dir, sprintf("ground_%s.asc", nm))
  source_name <- file.path(sp_dir, sprintf("source_%s.asc", nm))
  suit_name <- file.path(sp_dir, sprintf("suitability_%s.asc", nm))
  project_name <- file.path(sp_dir, time_periods[[i + 1]])
  fname <- file.path(sp_dir, sprintf("config_%s.ini", nm))
  
  # Load layers
  suit_t1 <- suit_list[[time_periods[[1]]]] + 0.0001 # non-zero
  habitat_t1 <- habitat_list[[time_periods[[1]]]]
  suit_t2 <- suit_list[[time_periods[[2]]]] + 0.0001
  habitat_t2 <- habitat_list[[time_periods[[2]]]]
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
  julia_command(sprintf('run_omniscape("%s")', fname), show_value = FALSE)
}
