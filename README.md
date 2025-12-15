# ClimConn
Climate adaptive habitat connectivity of mammal speicies in Colombia.

## Overview

This study uses SHapley Additive exPlanations (SHAP) to assess how climate and land cover drivers influence future mammal distributions across different global change scenarios. We quantify both directional changes (favoring/disfavoring turnovers) and gradual changes (favorability shifts) to evaluate species-specific spatial risks and opportunities under global change.

The repository includes:

- Scripts to project habitat suitability using BioModelos SDMs.
- Code for simulating climate adaptive connectivity using a multi-step framework
- Reproducible figures and tables

## Repository structure

```text
├── docs/                    # Technical reports and environmental setting.
├── R/                       # All analysis scripts
├── Schedulers/              # Supportive HPC schedule scripts for R scripts.
├── data/
│   └── config/              # Omniscape configs
│   └── dispersal/           # Species dispersal information
│   └── Env/                 # Environmental variables used for modeling
│   └── IUCN/                # IUCN redlist expert range maps
│   └── mamiferos/           # BioModelos SDMs for mammals
├── java/                    # MaxEnt files for BioModelos
├── results/
│   ├── block_radius/        # Summrized results for experiments on block size
│   └── omni/                # The Omniscape results for connectivity
│   └── omni_single/         # Comparison experiments with single step modeling
│   └── omni_block_size/     # The experiments results for block size
│   └── classification/      # Connectivity classification maps
│   └── figures/             # All visual figures and tables
└── README.md
```

## Start

Please check documents: [docs/biomodelos_connec.Rmd](docs/biomodelos_connec.Rmd) to set up the environment.

*NOTE:* due to some updates in Omniscape, users need to install our customized version.

## Session info

All analyses were conducted in R. The full R session information, including package versions and system details, is saved in [session_info.txt](session_info.txt) for reproducibility.

## Acknowledgement

This work is supported by project “DISES: Decision Making for Land Use Planning under Future Climate Scenarios through Engaged Research via Co-Design” funded by a grant to A.E.F. from the U.S. National Science Foundation (award number: 2401273).
