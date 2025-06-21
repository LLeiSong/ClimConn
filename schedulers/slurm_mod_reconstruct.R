# Directories
dat_dir <- "data"
bmls <- file.path(dat_dir, "mamiferos")

fname <- file.path(dat_dir, "species_miss.csv")
if (file.exists(fname)){
  folder_sps <- read.csv(fname)[["sp"]]
} else{
  folder_sps <- list.files(bmls)
}

message(sprintf("%s species to process.", length(folder_sps)))

chunk_length <- 10
folder_sps_list <- split(
  folder_sps, ceiling(seq_along(folder_sps) / chunk_length))

for (folder_sps in folder_sps_list){
  folder_sp <- paste(folder_sps, collapse = ",")
  system(sprintf("sbatch schedulers/mod_reconstruct.sh %s", folder_sp))
}
