#!/bin/bash
#SBATCH --job-name find_versions
#SBATCH --cpus-per-task 2
#SBATCH -c 2
#SBATCH --mem 16g
#SBATCH --partition allnodes 
#SBATCH --output scripts/find_versions.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.0.5 R CMD BATCH scripts/find_versions.R scripts/find_versions.Rout
