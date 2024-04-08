#!/bin/bash
#SBATCH --job-name CreateSeuratObjectWithCCA
#SBATCH --cpus-per-task 16
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Wu_etal_2021_BRCA_scRNA-CreateSeuratObjectWithCCA.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.0.5 R CMD BATCH scripts/Wu_etal_2021_BRCA_scRNA-CreateSeuratObjectWithCCA.R scripts/Wu_etal_2021_BRCA_scRNA-CreateSeuratObjectWithCCA.Rout
