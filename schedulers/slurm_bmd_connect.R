library(dplyr)

root_dir <- "/home/lsong/ClimConn"
fname <- file.path(root_dir, "data/sp_numbers.csv")
sp_numbers <- read.csv(fname)

for (scenario in c("ssp126", "ssp585")){
  sp_selected <- sp_numbers %>% filter(ssp == scenario) %>% 
    filter(IoU <= 50 & aoi_area >= 1e6) %>% 
    arrange(cur_area)
  
  sp_list <- sp_selected$species
  
  for (sp in sp_list){
    system(sprintf("sbatch schedulers/bmd_connect.sh %s %s", sp, scenario))
  }
}
