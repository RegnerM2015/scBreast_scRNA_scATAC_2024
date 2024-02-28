#!/bin/bash
#SBATCH --job-name get_SCSubtype_training_data
#SBATCH --cpus-per-task 1
#SBATCH -c 1
#SBATCH --mem 4g
#SBATCH --partition allnodes
#SBATCH --output scripts/get_SCSubtype_training_data.out 

# Make directory to store
mkdir SCSubtype_Training_Data

# Copy data from Franco Lab directory
cp /datastore/nextgenout5/share/labs/francolab/SCSubtype_Training_Data/* ./SCSubtype_Training_Data/