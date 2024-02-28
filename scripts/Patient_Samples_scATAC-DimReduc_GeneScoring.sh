#!/bin/bash
#SBATCH --job-name Patient_Samples_scATAC-DimReduc_GeneScoring
#SBATCH -c 8
#SBATCH --mem 256g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Patient_Samples_scATAC-DimReduc_GeneScoring.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.0.5 R CMD BATCH scripts/Patient_Samples_scATAC-DimReduc_GeneScoring.R scripts/Patient_Samples_scATAC-DimReduc_GeneScoring.Rout
