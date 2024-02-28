#!/bin/bash
#SBATCH --job-name All_Samples_scATAC-QC_DoubletRemoval_Preprocessing
#SBATCH --cpus-per-task 16
#SBATCH -c 16
#SBATCH --mem 256g
#SBATCH --partition allnodes 
#SBATCH --output scripts/All_Samples_scATAC-QC_DoubletRemoval_Preprocessing.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.0.5 R CMD BATCH scripts/All_Samples_scATAC-QC_DoubletRemoval_Preprocessing.R scripts/All_Samples_scATAC-QC_DoubletRemoval_Preprocessing.Rout
