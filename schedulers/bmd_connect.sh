#!/bin/bash
#SBATCH -c 10
#SBATCH --mem 40G
#SBATCH -t 30-00:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lsong@ucsb.edu

SP=$1
SCENARIO=$2

conda activate conn
cd ~/ClimConn

srun Rscript R/bmd_connect_warp.R -s $SP -o $SCENARIO