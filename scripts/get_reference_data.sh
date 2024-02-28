#!/bin/bash
#SBATCH --job-name get_reference_data
#SBATCH --cpus-per-task 1
#SBATCH -c 1
#SBATCH --mem 4g
#SBATCH --partition allnodes
#SBATCH --output scripts/get_reference_data.out 

# Get data from GEO Swarbrick paper
wget -qO- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz | tar xfz -

# Rename files to conform to Seurat's Read10x() function
mv Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx Wu_etal_2021_BRCA_scRNASeq/matrix.mtx
mv Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv Wu_etal_2021_BRCA_scRNASeq/features.tsv
mv Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv Wu_etal_2021_BRCA_scRNASeq/barcodes.tsv

# Gzip to conform to Seurat's Read10x() function
gzip Wu_etal_2021_BRCA_scRNASeq/matrix.mtx
gzip Wu_etal_2021_BRCA_scRNASeq/features.tsv
gzip Wu_etal_2021_BRCA_scRNASeq/barcodes.tsv
