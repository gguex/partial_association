#!/bin/bash -l

##SBATCH --account 1
#SBATCH --mail-type ALL 
#SBATCH --mail-user guillaume.guex@unil.ch

#SBATCH --chdir ./
#SBATCH --job-name cluster_3kernels
#SBATCH --output cluster_outputs/cluster_3kernels_%a.out

#SBATCH --partition cpu
#SBATCH --ntasks 1 

#SBATCH --cpus-per-task 16
#SBATCH --mem 16G 
#SBATCH --time 12:00:00 
#SBATCH --export NONE

#SBATCH --array=0-10

module load r-light

Rscript 2.1_cluster_3kernels.R $SLURM_ARRAY_TASK_ID
