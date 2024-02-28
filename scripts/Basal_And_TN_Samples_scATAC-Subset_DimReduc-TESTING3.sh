#!/bin/bash
#SBATCH --job-name Basal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3
#SBATCH -c 16
#SBATCH --mem 64g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Basal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.out 
#SBATCH --error scripts/Basal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.err

# Activate conda env that has singularity installed
# source /home/regnerm/anaconda3/etc/profile.d/conda.sh
# conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.1.2 Rscript scripts/Basal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.R
