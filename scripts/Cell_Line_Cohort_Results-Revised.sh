#!/bin/bash
#SBATCH --job-name Cell_Line_Cohort_Results-Revised
#SBATCH -c 32
#SBATCH --mem 128g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Cell_Line_Cohort_Results-Revised.out 

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.6.0 Rscript scripts/Cell_Line_Cohort_Results-Revised.R 
