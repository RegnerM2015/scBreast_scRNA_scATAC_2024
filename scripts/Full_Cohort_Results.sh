#!/bin/bash
#SBATCH --job-name Full_Cohort_Results
#SBATCH -c 32
#SBATCH --mem 128g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Full_Cohort_Results.out 

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.4.0 Rscript scripts/Full_Cohort_Results.R 
