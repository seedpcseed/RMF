#!/bin/sh
#SBATCH -A b1051
#SBATCH -p buyin
#SBATCH -J RMetaflow
#SBATCH -t 7-0:0
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --error=./log/R-%x.%j.err 
#SBATCH --output=./log/R-%x.%j.out

bash

cd $SLURM_SUBMIT_DIR

module purge all
module load singularity
module load R/4.0.3

Rscript --vanilla R/startup.R