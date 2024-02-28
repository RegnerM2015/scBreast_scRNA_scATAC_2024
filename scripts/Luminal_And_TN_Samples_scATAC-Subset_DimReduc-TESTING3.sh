#!/bin/bash
#SBATCH --job-name Luminal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3
#SBATCH -c 8
#SBATCH --mem 128g
#SBATCH --partition allnodes 
#SBATCH --error scripts/Luminal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.err
#SBATCH --output scripts/Luminal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.out 

# Activate conda env that has singularity installed
# source /home/regnerm/anaconda3/etc/profile.d/conda.sh
# conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.1.2 Rscript scripts/Luminal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.R
