#!/bin/bash
#SBATCH --job-name Basal_And_TN_Samples_scATAC-Transfer_Labels_from_scRNA_Call_Peaks-TESTING3
#SBATCH -c 32
#SBATCH --mem 128g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Basal_And_TN_Samples_scATAC-Transfer_Labels_from_scRNA_Call_Peaks-TESTING3.out 

# Activate conda env that has singularity installed
# source /home/regnerm/anaconda3/etc/profile.d/conda.sh
# conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.1.2 Rscript scripts/Basal_And_TN_Samples_scATAC-Transfer_Labels_from_scRNA_Call_Peaks-TESTING3.R
