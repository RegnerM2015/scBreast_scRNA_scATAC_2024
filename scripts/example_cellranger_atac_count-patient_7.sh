#!/usr/bin/env bash

#SBATCH --job-name HT-35EE8L-ATAC_MAY2019_A6
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output HT-35EE8L-ATAC_MAY2019_A6.cellranger-count.job.update.out 
#SBATCH --error HT-35EE8L-ATAC_MAY2019_A6.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger-atac count --id=HT-35EE8L-ATAC_MAY2019_A6_Update \
                      --fastqs=./fastq_path/H2NFLBGXB/HT-35EE8L-ATAC_MAY2019_A6 \
                      --reference=${DATA}/refdata-cellranger-atac-GRCh38-1.2.0 \
                      --sample=HT-35EE8L-ATAC_MAY2019_A6 \
                      --localcores=16
