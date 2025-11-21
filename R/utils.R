# ============================================================
# Script: utils.R
#
# Purpose:
#   Utility helpers for raster processing.
#   Currently includes:
#     - patch_clean(): remove very small patches from a raster
#       (keeps only patches ≥ 1,000 km² and returns a 0/1 mask).
#
# Inputs (for patch_clean):
#   rst    SpatRaster (integer / categorical) with patch IDs or presence values
#
# Output:
#   SpatRaster mask with:
#     1 = retained (large) patches
#     0 = background / removed small patches
#
# Author: Lei Song <lei.song@rutgers.edu>
# Last updated: 2025-11-21
# ============================================================

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
