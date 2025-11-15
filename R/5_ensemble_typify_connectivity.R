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

#### Reclassify connectivity #####
for (nm in names(conn_lyrs)){
  message(nm)
  
  # Get the layer, stretch it, does not really matter
  lyr <- stretch(conn_lyrs[[nm]], 0, 8000)
  
  # Define a circle smoothing windows
  w7 <- focalMat(
    lyr, d = res(lyr)[1] * 3.5, 
    type = "circle", fillNA = TRUE)
  w7[!is.na(w7)] <- 1
  
  # STEP 0. Sinks and plateaus
  conn_file <- file.path(tempdir(), "conn.tif")
  conn_rev_file <- file.path(tempdir(), "conn_rev.tif")
  sink_file <- file.path(tempdir(), "sink.tif")
  peak_file <- file.path(tempdir(), "peak.tif")
  writeRaster(lyr, conn_file)
  writeRaster(8000 - lyr, conn_rev_file)
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
  top <- top_all * peaks * 100 + top_all * 10 + peaks
  
  # STEP 2. lowlands
  low_all <- lyr <= quantile(values(lyr), 0.1, na.rm = TRUE)
  low_all <- patch_clean(low_all)
  low <- low_all * pits * 200 + low_all * 20 + pits * 2
  
  # STEP 3. scarp slope
  lyr_slope <- terrain(lyr, v = "slope", unit = "degrees", neighbors = 8)
  scarp <- lyr_slope >= quantile(values(lyr_slope), 0.9, na.rm = TRUE)
  
  conn_high <- focal(lyr, w = w7, fun = "max", na.rm = TRUE) >= 
    quantile(values(lyr), 0.5, na.rm = TRUE)
  
  conn_low <- focal(lyr, w = w7, fun = "min", na.rm = TRUE) < 
    quantile(values(lyr), 0.5, na.rm = TRUE)
  
  scarp <- scarp * conn_high * conn_low
  scarp <- focal(scarp, w = w7, fun = "modal", na.rm = TRUE)
  
  # STEP 4. combine
  overlap <- (top != 0) & (low != 0)
  classes <- cover(scarp * 300, ifel(overlap, 0, top + low), 0)
  classes <- focal(classes, w = w7, fun = "modal", na.rm = TRUE)
  
  # Step 3: Set it as a categorical (factor) raster
  levels(classes) <- data.frame(
    ID = c(0, 1, 10, 111, 2, 20, 222, 300),
    class = c("Background", "Through Lane", "Freeway mainline", "Express Lane",
              "Detour Zone", "Avoidance Zone", "Bottleneck", "Merge Lane"))
  
  fname <- file.path(typology_dir, sprintf("conn_types_%s.tif", nm))
  writeRaster(classes, fname, overwrite = TRUE)
}

classes <- data.frame(
  ID = c(0, 1, 10, 111, 2, 20, 222, 300),
  class = c("Background", "Through Lane", "Freeway mainline", "Express Lane",
            "Detour Zone", "Avoidance Zone", "Bottleneck", "Merge Lane"))

fname <- file.path(typology_dir, "connectivity_typology.csv")
write.csv(classes, fname, row.names = FALSE) 

# #### Extract contour lines #####
# # Parameters
# smoothness <- 1
# contour_interval <- 0.02
# contour_min <- round(min(global(conn_lyrs, "min", na.rm = TRUE)[["min"]]), 1)
# contour_max <- round(max(global(conn_lyrs, "max", na.rm = TRUE)[["max"]]), 1)
# 
# for (nm in names(conn_lyrs)){
#   message(nm)
#   
#   # Get the layer
#   lyr <- conn_lyrs[[nm]]
#   
#   # Get the raw contours
#   contours <- st_as_stars(lyr) %>% 
#     st_contour(breaks = seq(
#       from = contour_min, to = contour_max, by = contour_interval))
#   
#   # Drop small pieces gradually until 100 square km
#   areas <- seq(10, 100, 10)
#   bry <- contours %>% st_union() %>% 
#     st_buffer(units::set_units(-0.01, "m")) %>% st_cast("MULTILINESTRING")
#   
#   for (area in areas){
#     # Target small pieces along the border
#     tiny_plys <- st_cast(contours, "MULTIPOLYGON") %>% 
#       st_cast("POLYGON", warn = FALSE) %>% 
#       mutate(ar = st_area(., "km2")) %>% 
#       filter(ar <= units::set_units(area, "km2")) %>% 
#       slice(unique(unlist(st_intersects(bry, .)))) %>% 
#       select(-ar)
#     
#     # Drop crumbs, add back the border ones, then fill_holes
#     if (nrow(tiny_plys) > 0){
#       contours <- contours %>% 
#         drop_crumbs(units::set_units(area, "km2")) %>% 
#         vect() %>% combineGeoms(vect(tiny_plys), distance = TRUE) %>% 
#         st_as_sf() %>% fill_holes(units::set_units(area * 10, "km2"))
#     } else{
#       contours <- contours %>% 
#         drop_crumbs(units::set_units(area, "km2")) %>% 
#         fill_holes(units::set_units(area * 10, "km2"))
#     }
#   }
#   
#   # Smooth it
#   contours <- contours %>% smooth(method = "ksmooth", smoothness = smoothness)
#   contours <- st_cast(contours, "MULTIPOLYGON") %>% 
#     st_cast("POLYGON", warn = FALSE)
#   rm(areas, bry); names(contours)[1] <- "Range"
#   contours <- contours %>% mutate(ID = 1:nrow(.)) %>% 
#     select(ID, Range, Min, Max)
#   
#   # Save out
#   fname <- file.path(
#     typology_dir, sprintf('contours_%s_%s_%s.shp', 
#                           contour_interval, smoothness, nm))
#   if (file.exists(fname)) st_delete(fname)
#   st_write(contours, fname)
# }
