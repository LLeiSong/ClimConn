library(dplyr)

root_dir <- "/scratch/ls1686/ClimConn"
sp_list <- list.files(file.path(root_dir, "data/mamiferos"))

chunk_length <- 5
sp_list <- split(
  sp_list, ceiling(seq_along(sp_list) / chunk_length))

for (scenario in c("ssp126", "ssp370", "ssp585")){
  for (sp in sp_list){
    sp <- paste(sp, collapse = ",")
    system(sprintf("sbatch schedulers/bmd_connect.sh %s %s", sp, scenario))
  }
}
