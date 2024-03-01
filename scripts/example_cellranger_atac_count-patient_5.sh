#!/usr/bin/env bash

#SBATCH --job-name HT-35A4AL-ATAC_MAY2019_A5
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output HT-35A4AL-ATAC_MAY2019_A5.cellranger-count.job.update.out 
#SBATCH --error HT-35A4AL-ATAC_MAY2019_A5.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger-atac count --id=HT-35A4AL-ATAC_MAY2019_A5_Update \
                      --fastqs=./fastq_path/H2NFLBGXB/HT-35A4AL-ATAC_MAY2019_A5 \
                      --reference=${DATA}/refdata-cellranger-atac-GRCh38-1.2.0 \
                      --sample=HT-35A4AL-ATAC_MAY2019_A5 \
                      --localcores=16
