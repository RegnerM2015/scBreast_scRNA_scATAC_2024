#!/bin/bash
#SBATCH --job-name scLME_update-metacells-Basal-SingFits_OLS
#SBATCH -c 96
#SBATCH --mem 256g
#SBATCH --partition allnodes 
#SBATCH --output scripts/scLME_update-metacells-Basal-SingFits_OLS.out 
#SBATCH --error scripts/scLME_update-metacells-Basal-SingFits_OLS.err

# Activate conda env that has singularity installed
# source /home/regnerm/anaconda3/etc/profile.d/conda.sh
# conda activate singularity

# Run R script in docker container 
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.1.2 Rscript scripts/scLME_update-metacells-Basal-SingFits_OLS.R

