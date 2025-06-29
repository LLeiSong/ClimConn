## Select species and identify species traits

1. Run `scrutinize_species.R` to remote marine or freshwater only species from the species list of BioModelos.
2. Run `dispersal_rate.R` to calculate dispersal rate for each species.

## Model habitat quality

1. Run `prepare_variables.R` to prepare all environmental layers for modeling and future projections.
2. Run `reconstruct_models.R` to reconstruct BioModelos SDMs for selected species.
3. (Optional) Run `check_reconstruct.R` to check the reconstructed results and remove temporary files.

## Assess connectivity

1. Run `bmd_connect_warp.R` to model omni-directional and multi-temporal habitat connectivity for all species. This script calls function `clim_connect` defined in file `clim_connect_omni_cmd.R`, which runs the model by calling the Julia script `run_omniscape.jl`. Alternatively, `clim_connect_omni.R` defines another version of `clim_connect` using R package `JuliaCall` to run Julia script.
2. Run `summarize_species.R` to check how many time stamps all species will maintain connectivity to new habitats.
