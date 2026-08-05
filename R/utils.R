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

# The function to calculate geometric mean
geom_m <- function(v){
  v[ v <= 0] <- 1e-6
  exp(mean(log(v), na.rm = TRUE))
}

# The function to get effective mask
f_mask <- function(lyr){
  vals <- values(lyr)
  vals_pos <- vals[vals > 0]  # ignore exact zeros
  
  # if everything is zero, bail out
  if (length(vals_pos) > 0) {
    vals_sorted <- sort(vals_pos)
    cum_sum <- cumsum(vals_sorted) / sum(vals_sorted)
    
    # choose threshold so that values below 
    # contribute at most 1% of total connectivity
    idx <- which(cum_sum > 0.01)[1] 
    Tmin <- vals_sorted[idx]
    
    # optional safety: don't let Tmin be absurdly small
    Tmin <- max(Tmin, 1e-6)
    
    # make an "effective" connectivity raster that ignores near-zero junk
    lyr_eff <- lyr
    lyr_eff[lyr_eff < Tmin] <- NA
  } else {
    lyr_eff <- lyr
  }
  
  lyr_eff
}
