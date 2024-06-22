#!/bin/bash
#SBATCH --job-name DepMap_Analysis_Basal_Luminal-updates
#SBATCH -c 4
#SBATCH --mem 16g
#SBATCH --partition allnodes 
#SBATCH --output scripts/DepMap_Analysis_Basal_Luminal-updates.out 

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:2.0.0 Rscript scripts/DepMap_Analysis_Basal_Luminal-updates.R 
