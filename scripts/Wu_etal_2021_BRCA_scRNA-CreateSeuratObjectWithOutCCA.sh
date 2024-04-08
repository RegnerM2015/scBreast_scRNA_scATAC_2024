#!/bin/bash
#SBATCH --job-name CreateSeuratObjectWithOutCCA
#SBATCH --cpus-per-task 16
#SBATCH -c 16
#SBATCH --mem 64g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Wu_etal_2021_BRCA_scRNA-CreateSeuratObjectWithOutCCA.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.0.5 R CMD BATCH scripts/Wu_etal_2021_BRCA_scRNA-CreateSeuratObjectWithOutCCA.R scripts/Wu_etal_2021_BRCA_scRNA-CreateSeuratObjectWithOutCCA.Rout
