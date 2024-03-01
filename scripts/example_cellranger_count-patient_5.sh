#!/usr/bin/env bash

#SBATCH --job-name HT-35A4AL-RNA_MAY2019_G5
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output HT-35A4AL-RNA_MAY2019_G5.cellranger-count.job.update.out 
#SBATCH --error HT-35A4AL-RNA_MAY2019_G5.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger count  --id=HT-35A4AL-RNA_MAY2019_G5_Update \
                      --fastqs=./fastq_path/H2YT7BGXB/HT-35A4AL-RNA_MAY2019_G5 \
                      --transcriptome=${DATA}/refdata-cellranger-GRCh38-3.0.0 \
                      --sample=HT-35A4AL-RNA_MAY2019_G5 \
                      --localcores=16 \
                      --localmem=80
