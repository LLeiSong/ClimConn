#!/bin/bash
#SBATCH -c 2
#SBATCH --mem 16G
#SBATCH -t 30-00:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong@ucsb.edu

SP=$1

conda activate conn
cd ~/ClimConn

srun Rscript R/reconstruct_models.R -f $SP