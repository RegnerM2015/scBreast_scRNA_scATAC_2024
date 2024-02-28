#!/bin/bash
#SBATCH --job-name make_cellranger_directories
#SBATCH --cpus-per-task 1
#SBATCH -c 1
#SBATCH --mem 4g
#SBATCH --partition allnodes
#SBATCH --output scripts/make_cellranger_directories.out 

# Make cellranger output directories
mkdir ./cellranger-atac_outputs
mkdir ./cellranger_outputs

# Store sample/container IDs
declare -a StringArray=("35A4AL 35EE8L 3821AL 3B3E9L 3C7D1L 3D388L 3FCDEL 43E7BL 43E7CL 44F0AL 45CB0L 49758L 49CFCL 4AF75L 4B146L 4C2E5L 4D0D2L HCC1143 MCF7 SUM149PT T47D")

# Iterate the string array using for loop, make subdirectories for each sample's filtered feature barcode matrix
for i in ${StringArray[@]}; do
   mkdir ./cellranger_outputs/filtered_feature_bc_matrix_${i}/
done
