# ============================================================
# Script: 3.2_pair_size_radius.R
#
# Purpose:
#   Derive a rule of thumb linking Omniscape block size to dispersal radius:
#     - Compute correlations between block_size runs and block_size = 1 baseline
#     - Fit logistic curves of cor ~ radius_size for each block size
#     - Identify minimum radius allowing a given block size at target correlation.
#
# Inputs:
#   results/omni_block_size/<sp>/2021-2040_*   cum_currmap_*.tif + config.ini
#
# Outputs (in results/block_radius/):
#   correlation_block_size.csv   correlations by species × block_size
#   fit_curves.rda               nls objects for each block_size
#   allow_min_radius.csv         min radius allowed per block_size at cor_target
#   fitted_curves.csv            smoothed fitted curves over radius_size
#
# ============================================================

# Loading parameters
library(here)
library(dplyr)
library(terra)
library(stringr)

# Set directories
root_dir <- here()
bs_sens_dir <- file.path(root_dir, "results/omni_block_size")
results_dir <- file.path(root_dir, "results/block_radius")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Other settings
# as this is the biggest regular one, we don't want to use a too large block
max_block_size <- 11
cor_target <- 0.95

# Start the process
sps <- list.files(bs_sens_dir)

sens_cors <- lapply(sps, function(sp){
  dirs <- list.dirs(file.path(bs_sens_dir, sp), recursive = FALSE)
  
  baseline_dir <- dirs[str_detect(dirs, "2021-2040_1$")]
  baseline <- rast(file.path(baseline_dir, "cum_currmap_1.tif"))
  
  lapply(dirs[!str_detect(dirs, "2021-2040_1$")], function(fld){
    i <- readLines(file.path(fld, "config.ini"), warn = FALSE)
    i <- i[str_detect(i, "block_size = ")]
    i <- as.integer(str_extract(i, "[0-9]+"))
    
    lyr <- rast(file.path(fld, sprintf("cum_currmap_%s.tif", i)))
    vals <- values(c(lyr, baseline))
    vals <- vals[complete.cases(vals), ]
    r <- cor(vals[, 1], vals[, 2])
    
    data.frame(block_size = i, cor = r)
  }) %>% bind_rows() %>% 
    rbind(data.frame(block_size = 1, cor = 1)) %>% 
    mutate(sp = gsub("\\.", " ", sp)) %>% arrange(block_size)
}) %>% bind_rows()

# Get the max size for each species
radius_size <- sens_cors %>% 
  group_by(sp) %>% summarise(radius_size = max(block_size))

# Attach the max size back to the data.frame
sens_cors <- sens_cors %>% 
  left_join(radius_size, by = "sp")
rm(radius_size)

# Filter out the data to be used
sens_cors <- sens_cors %>% 
  filter(block_size <= max_block_size & block_size > 1)

write.csv(sens_cors, file.path(results_dir, "correlation_block_size.csv"),
          row.names = FALSE)

# fit logistic growth model
mods <- lapply(seq(3, max_block_size, 2), function(i){
  dat <- sens_cors %>% filter(block_size == i) %>% select(cor, radius_size)
  nls(cor ~ a / (1 + b * exp(-k * radius_size)),
      data = dat,
      start = list(a = max(dat$cor), b = 1, k = 0.1))
}); names(mods) <- seq(3, max_block_size, 2)

save(mods, file = file.path(results_dir, "fit_curves.rda"))

# Calculate the intersect with the target cor
allow_min_radius <- lapply(seq(3, max_block_size, 2), function(i){
  mod <- mods[[as.character(i)]]
  a <- mod$m$getPars()[[1]]
  b <- mod$m$getPars()[[2]]
  k <- mod$m$getPars()[[3]]
  
  -log((a / (cor_target * b)) - (1 / b)) / k
}) %>% unlist() %>% ceiling() %>% 
  data.frame(allow_min_radius = .,
             block_size = seq(3, max_block_size, 2))

write.csv(allow_min_radius, file.path(results_dir, "allow_min_radius.csv"),
          row.names = FALSE)

# Generate smooth curve
fitted_curves <- data.frame(
  radius_size = seq(min(sens_cors$radius_size), max(sens_cors$radius_size), 
                 length.out = 200))
fitted_curves <- lapply(seq(3, max_block_size, 2), function(i){
  predict(mods[[as.character(i)]], newdat) %>% 
    data.frame(pred = .) %>% mutate(block_size = i) %>% cbind(fitted_curves)
}) %>% bind_rows()

write.csv(fitted_curves, file.path(results_dir, "fitted_curves.csv"),
          row.names = FALSE)
