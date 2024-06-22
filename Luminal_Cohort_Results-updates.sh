#!/bin/bash
#SBATCH --job-name Luminal_Cohort_Results-updates
#SBATCH -c 32
#SBATCH --mem 128g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Luminal_Cohort_Results-updates.out 

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.2.0 Rscript scripts/Luminal_Cohort_Results-updates.R 