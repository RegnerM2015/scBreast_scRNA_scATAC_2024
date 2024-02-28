#!/usr/bin/env bash
#SBATCH --job-name Individual_Samples_scRNA-RemoveLowMappingRatePopulation_Reprocess_3FCDEL
#SBATCH -c 32
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output scripts/Individual_Samples_scRNA-RemoveLowMappingRatePopulation_Reprocess_3FCDEL.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore \
                 --home /datastore/nextgenout5/share/labs/francolab/scBreast_scRNA_scATAC_2023 \
                 docker://regnerm/scbreast_2023:1.0.5 \
                 R CMD BATCH '--args 3FCDEL' scripts/Individual_Samples_scRNA-RemoveLowMappingRatePopulation_Reprocess_3FCDEL.R scripts/Individual_Samples_scRNA-RemoveLowMappingRatePopulation_Reprocess_3FCDEL.Rout
