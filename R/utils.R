# Clean tiny patches from a raster
patch_clean <- function(rst){
  rst <- patches(rst, zeroAsNA = TRUE)
  rz <- zonal(cellSize(rst, unit = "km"), 
              rst, sum, as.raster = TRUE)
  rst <- ifel(rz < 1000, NA, rst)
  rst[!is.na(rst)] <- 1
  rst[is.na(rst)] <- 0
  
  rst
}
