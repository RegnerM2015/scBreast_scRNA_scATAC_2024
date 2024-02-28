#!/bin/bash
#SBATCH --job-name Table_1_S1
#SBATCH -c 4
#SBATCH --mem 16g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Table_1_S1.out 

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.4.0 Rscript scripts/Table_1_S1.R 
