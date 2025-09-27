## Check and clean the Biomodelos species
## Because some species might not be suitable to model connectivity
## we need to check and remove these species
## - Marine species
## - Aquatic species
## - Domestic species
## - Extinct species
## Author: Lei Song (lei.song@rutgers.edu)

# Load libraries
library(sf)
library(dplyr)
library(terra)
library(here)

dat_dir <- here("data")
bmls <- file.path(dat_dir, "mamiferos")

# 1. Remove empty models, which are not included in BioModelos (e.g. marine).
fnames <- list.files(bmls, full.names = TRUE)

for (fname in fnames){
  n_files <- length(list.files(file.path(fname, "final_models_enmeval")))
  
  if (n_files == 0) unlink(fname, recursive = TRUE)
}

# 2. Remove the remaining freshwater or other non-suitable species
fnames <- list.files(bmls, full.names = TRUE)
species_info <- st_read(file.path(dat_dir, "IUCN/MAMMALS/MAMMALS.shp"))
species_info <- species_info %>% st_drop_geometry() %>% 
  select(sci_name, marine, terrestria, freshwater) %>% unique()

sp_names <- gsub("\\.", " ", basename(fnames)) %>% 
  data.frame("biomodelos" = .)

species_info <- right_join(
  species_info, sp_names, by = c("sci_name" = "biomodelos"))

# Split to species to manually check and automatic filtering 
species_info_mc <- species_info %>% filter(is.na(marine))

# Save out to edit, because X11 is not capable 
# and I cannot reinstall R right now.
write.csv(species_info_mc, "data/species_info_mc.csv", row.names = FALSE)

species_info_mc <- read.csv("data/species_info_mc.csv") %>% 
  filter(Wild == "Yes")
species_info_af <- species_info %>% 
  filter(marine != "true" & terrestria != "false")

# 3. Save out species list.
species <- c(species_info_mc$sci_name, species_info_af$sci_name)
save(species, file = file.path(dat_dir, "mammal_species_conn.rda"))

# 4. Clean unwanted species
for (fname in fnames){
  if (!basename(fname) %in% gsub(" ", ".", species)){
    message(sprintf("Remove species folder: %s", basename(fname)))
    unlink(fname, recursive = TRUE)
  }
}

# Clean up
rm(list = ls())

# So far, only the relevant species left.
