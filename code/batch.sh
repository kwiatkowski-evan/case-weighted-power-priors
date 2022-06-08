#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 47:59:59
#SBATCH --mem=1g
#SBATCH --output=../slurm/%a.slurm
#SBATCH --error=../error/%a.err
#SBATCH --array=1-2400

## add R module
module add r/4.1.0

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" 02-full-sims.R ../logs/log$SLURM_ARRAY_TASK_ID.Rout
