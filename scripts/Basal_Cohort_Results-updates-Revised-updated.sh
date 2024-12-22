#!/bin/bash
#SBATCH --job-name Basal_Cohort_Results-updates-Revised-updated
#SBATCH -c 32
#SBATCH --mem 128g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Basal_Cohort_Results-updates-Revised-updated.out 

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.2.0 Rscript scripts/Basal_Cohort_Results-updates-Revised-updated.R 