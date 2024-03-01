#!/usr/bin/env bash

#SBATCH --job-name H2YT7BGXB
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output H2YT7BGXB_demultiplex.job.out
#SBATCH --error H2YT7BGXB_demultiplex.job.err

DATA=/datastore/nextgenout5/share/labs/bioinformatics/seqware/francolab_10x_copy
OUT=/datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.05.28.2019

cellranger mkfastq --id=H2YT7BGXB \
                   --run=${DATA}/190524_NS500270_0300_AH2YT7BGXB \
                   --csv=${OUT}/Samplesheet.csv \
                   --qc \
                   --localcores=16
