#!/usr/bin/env bash
#SBATCH --job-name Patient_Samples_scRNA-Merge_And_ReCluster
#SBATCH -c 32
#SBATCH --mem 128g
#SBATCH --partition allnodes
#SBATCH --output scripts/Patient_Samples_scRNA-Merge_And_ReCluster.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore \
                 --home /datastore/nextgenout5/share/labs/francolab/scBreast_scRNA_scATAC_2023 \
                 docker://regnerm/scbreast_2023:1.0.5 \
                 R CMD BATCH scripts/Patient_Samples_scRNA-Merge_And_ReCluster.R scripts/Patient_Samples_scRNA-Merge_And_ReCluster.Rout
