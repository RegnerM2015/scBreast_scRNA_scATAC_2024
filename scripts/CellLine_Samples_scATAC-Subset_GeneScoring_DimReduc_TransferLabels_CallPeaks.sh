#!/bin/bash
#SBATCH --job-name CellLine_Samples_scATAC-Subset_GeneScoring_DimReduc_TransferLabels_CallPeaks
#SBATCH -c 16
#SBATCH --mem 64g
#SBATCH --partition allnodes 
#SBATCH --output scripts/CellLine_Samples_scATAC-Subset_GeneScoring_DimReduc_TransferLabels_CallPeaks.out 

# # Activate conda env that has singularity installed
# source /home/regnerm/anaconda3/etc/profile.d/conda.sh
# conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.0.5 Rscript scripts/CellLine_Samples_scATAC-Subset_GeneScoring_DimReduc_TransferLabels_CallPeaks.R 