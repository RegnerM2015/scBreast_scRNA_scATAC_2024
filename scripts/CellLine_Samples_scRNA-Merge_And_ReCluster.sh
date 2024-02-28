#!/usr/bin/env bash
#SBATCH --job-name CellLine_Samples_scRNA-Merge_And_ReCluster
#SBATCH -c 32
#SBATCH --mem 96g
#SBATCH --partition allnodes
#SBATCH --output scripts/CellLine_Samples_scRNA-Merge_And_ReCluster.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore \
                 --home /datastore/nextgenout5/share/labs/francolab/scBreast_scRNA_scATAC_2023 \
                 docker://regnerm/scbreast_2023:1.0.5 \
                 R CMD BATCH scripts/CellLine_Samples_scRNA-Merge_And_ReCluster.R scripts/CellLine_Samples_scRNA-Merge_And_ReCluster.Rout
