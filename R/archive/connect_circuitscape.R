library(ini)
library(JuliaCall)
library(terra)
library(sf)
library(dplyr)
library(readr)
library(spThin)

clim_connect <- function(dis_dist, # in meters
                         suit_t1, 
                         habitat_t1,
                         habitat_t2,
                         circ_dir, # full working directory
                         n_iteration = 1000) {
  sp_dir <- file.path(circ_dir, sp)
  if (!dir.exists(sp_dir)) dir.create(sp_dir)
  
  # Make a reasonable iteration number if not set
  if (is.null(n_iteration)){
    n_iteration <- floor(
      global(cellSize(habitat_t1), "sum", na.rm = TRUE)[[1]] / 
        (pi * dis_dist^2) * 2)}
  
  # set up Julia
  # julia_setup(installJulia = TRUE)
  julia_install_package_if_needed("Circuitscape")
  julia_library("Circuitscape")
  julia_library("Logging")
  julia_command("Logging.disable_logging(Logging.Info)", show_value = FALSE)
  config <- read.ini(file.path("data", "circuitscape_setting_template.ini"))
  
  suit_name <- file.path(sp_dir, "suit.asc")
  node_name <- file.path(sp_dir, "nodes.asc")
  
  # Update config for circuitscape
  config$`Output options`$write_cur_maps <- 0
  config$`Output options`$write_volt_maps <- 0
  config$`Output options`$write_cum_cur_map_only <- "True"
  config$`Output options`$write_max_cur_maps <- "False"
  config$`Logging Options`$log_file <-
    file.path(sp_dir, "log.log")
  config$`Options for pairwise and one-to-all and all-to-one modes`$point_file <- node_name
  config$`Habitat raster or graph`$habitat_file <- suit_name
  config$`Calculation options`$parallelize <- 'True'
  config$`Calculation options`$print_timings <- 'False'
  config$`Calculation options`$print_rusages <- 'False'
  config$`Calculation options`$preemptive_memory_release <- 'True'
  
  cum_rasters <- lapply(1:n_iteration, function(n){
    nm <- sprintf('iter_%s', n)
    config$`Output options`$output_file <- file.path(sp_dir, nm)
    fname <- file.path(sp_dir, "config.ini")
    write.ini(config, fname)
    
    # Make a list of source nodes
    set.seed(n)
    src_nodes <- spatSample(
      habitat_t1, size = min(100, global(habitat_t1, "sum", na.rm = TRUE)[[1]]),
      na.rm = TRUE, values = FALSE, xy = TRUE) %>% data.frame() %>% 
      st_as_sf(coords = c("x", "y"), crs = crs(habitat_t1)) %>% 
      st_transform(4326) %>% st_coordinates() %>% data.frame() %>% 
      mutate(SPEC = sp)
    set.seed(n)
    src_nodes <- thin(src_nodes, lat.col = "Y", long.col = "X", 
                      thin.par = dis_dist * 2.1 / 1000, reps = 1,
                      write.files = FALSE, verbose = FALSE,
                      locs.thinned.list.return = TRUE)[[1]]
    src_nodes <- src_nodes %>% 
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
      st_transform(crs(habitat_t1))
    
    pseudo_window <- src_nodes %>% st_buffer(dist = dis_dist) %>% 
      st_union()
    
    suit_t1 <- mask(suit_t1, vect(pseudo_window)) %>% trim()
    suit_t1[is.na(suit_t1)] <- -9999
    writeRaster(suit_t1, suit_name, overwrite = TRUE,
                wopt = list(NAflag = -9999))
    
    habitat_t1 <- mask(habitat_t1, vect(pseudo_window)) %>% trim()
    habitat_t2 <- mask(habitat_t2, vect(pseudo_window)) %>% trim()
    
    # Pair the source nodes and random five sink nodes
    sink_nodes <- lapply(1:nrow(src_nodes), function(sid){
      src_node <- src_nodes %>% slice(sid)
      msk <- src_node %>% st_buffer(dis_dist)
      
      dst_nodes <- habitat_t2 %>% crop(msk) %>% mask(msk)
      dst_nodes <- mask(
        dst_nodes, rasterize(src_node, dst_nodes), inverse = TRUE)
      
      set.seed(n)
      spatSample(
        dst_nodes, size = min(20, global(dst_nodes, "sum", na.rm = TRUE)[[1]]), 
        na.rm = TRUE, values = FALSE, xy = TRUE) %>% data.frame() %>% 
        st_as_sf(coords = c("x", "y"), crs = crs(dst_nodes))
    }) %>% bind_rows()
    
    # Prepare nodes
    nodes <- rbind(src_nodes %>% mutate(id = 1), 
                   sink_nodes %>% mutate(id = 2) %>% select(id))
    nodes <- rasterize(nodes, suit_t1, field = "id")
    nodes[is.na(nodes)] <- -9999
    writeRaster(nodes, node_name, overwrite = TRUE,
                wopt = list(NAflag = -9999))
    
    # Run the experiment
    if (file.exists(fname)) {
      julia_command(sprintf('compute("%s")', fname), show_value = FALSE)
    } else {
      warning("Failed to write out config file.")
    }
    
    cum_raster <- rast(
      file.path(sp_dir, sprintf("%s_cum_curmap.asc", nm)))
    mask(cum_raster, nodes, maskvalues = -9999, 
         inverse = TRUE, updatevalue = NA)
  })
  
  cum_rasters <- do.call(
    c, lapply(cum_rasters, function(lyr) extend(lyr, suit_t1)))
  cum_rasters <- mean(cum_rasters, na.rm = TRUE)
  
  writeRaster(cum_rasters, 
              file.path(sp_dir, "cum_curmap.tif"))
}
