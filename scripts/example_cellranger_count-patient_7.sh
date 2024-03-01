#!/usr/bin/env bash

#SBATCH --job-name HT-35EE8L-RNA_MAY2019_G6
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output HT-35EE8L-RNA_MAY2019_G6.cellranger-count.job.update.out 
#SBATCH --error HT-35EE8L-RNA_MAY2019_G6.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger count  --id=HT-35EE8L-RNA_MAY2019_G6_Update \
                      --fastqs=./fastq_path/H2YT7BGXB/HT-35EE8L-RNA_MAY2019_G6 \
                      --transcriptome=${DATA}/refdata-cellranger-GRCh38-3.0.0 \
                      --sample=HT-35EE8L-RNA_MAY2019_G6 \
                      --localcores=16 \
                      --localmem=80
