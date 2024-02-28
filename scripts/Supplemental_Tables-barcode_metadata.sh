#!/bin/bash
#SBATCH --job-name Supplemental_Tables-barcode_metadata
#SBATCH -c 4
#SBATCH --mem 16g
#SBATCH --partition allnodes 
#SBATCH --output scripts/Supplemental_Tables-barcode_metadata.out 

# Run R script in docker container 
singularity exec --bind /datastore --home $PWD docker://regnerm/scbreast_2023:1.8.0 Rscript scripts/Supplemental_Tables-barcode_metadata.R 
