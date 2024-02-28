#!/usr/bin/env bash
#SBATCH --job-name Basal_And_TN_Samples_scRNA-Subset_And_ReCluster-TESTING
#SBATCH -c 8
#SBATCH --mem 16g
#SBATCH --partition allnodes
#SBATCH --output scripts/Basal_And_TN_Samples_scRNA-Subset_And_ReCluster-TESTING.out 

# Activate conda env that has singularity installed
# source /home/regnerm/anaconda3/etc/profile.d/conda.sh
# conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.1.2 Rscript scripts/Basal_And_TN_Samples_scRNA-Subset_And_ReCluster-TESTING.R
