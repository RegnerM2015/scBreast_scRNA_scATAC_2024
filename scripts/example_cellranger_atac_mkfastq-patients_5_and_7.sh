#!/usr/bin/env bash

#SBATCH --job-name H2NFLBGXB
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output H2NFLBGXB_demultiplex.job.out
#SBATCH --error H2NFLBGXB_demultiplex.job.err

DATA=/datastore/nextgenout5/share/labs/bioinformatics/seqware/francolab_10x_copy
OUT=/datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.05.29.2019

cellranger-atac mkfastq --id=H2NFLBGXB \
                        --run=${DATA}/190528_NS500270_0301_AH2NFLBGXB \
                        --csv=${OUT}/Samplesheet.csv \
                        --qc \
                        --localcores=16
