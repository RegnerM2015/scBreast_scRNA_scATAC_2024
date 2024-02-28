#!/bin/bash
#SBATCH --job-name copy_cellranger_outputs
#SBATCH --cpus-per-task 2
#SBATCH -c 2
#SBATCH --mem 16g
#SBATCH --partition allnodes
#SBATCH --output scripts/copy_cellranger_outputs.out 

# Copy cellranger-atac fragments files to location

# Patient samples
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.02.02.2021/49758L-ATAC/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_49758L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.02.02.2021/49758L-ATAC/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_49758L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.25.2021/49CFCL-ATAC_update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_49CFCL.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.25.2021/49CFCL-ATAC_update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_49CFCL.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.25.2021/4AF75L-ATAC_update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_4AF75L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.25.2021/4AF75L-ATAC_update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_4AF75L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.29.2021/HN4B146L_Mar2021_ATAC_update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_4B146L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.29.2021/HN4B146L_Mar2021_ATAC_update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_4B146L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.05.29.2019/HT-35A4AL-ATAC_MAY2019_A5_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_35A4AL.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.05.29.2019/HT-35A4AL-ATAC_MAY2019_A5_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_35A4AL.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.05.29.2019/HT-35EE8L-ATAC_MAY2019_A6_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_35EE8L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.05.29.2019/HT-35EE8L-ATAC_MAY2019_A6_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_35EE8L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.10.11.2019/3821AL-ATAC_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_3821AL.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.10.11.2019/3821AL-ATAC_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_3821AL.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.01.15.2020/3B3E9L-ATAC_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_3B3E9L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.01.15.2020/3B3E9L-ATAC_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_3B3E9L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.02.10.2020/3C7D1L-ATAC_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_3C7D1L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.02.10.2020/3C7D1L-ATAC_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_3C7D1L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.19.2020/3D388L-ATAC/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_3D388L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.19.2020/3D388L-ATAC/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_3D388L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.19.2020/3FCDEL-ATAC/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_3FCDEL.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.19.2020/3FCDEL-ATAC/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_3FCDEL.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.08.07.2020/43E7BL-ATAC/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_43E7BL.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.08.07.2020/43E7BL-ATAC/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_43E7BL.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.08.07.2020/43E7CL-ATAC/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_43E7CL.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.08.07.2020/43E7CL-ATAC/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_43E7CL.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.09.04.2020/44F0AL-ATAC/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_44F0AL.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.09.04.2020/44F0AL-ATAC/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_44F0AL.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.02.02.2021/45CB0L-ATAC/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_45CB0L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.02.02.2021/45CB0L-ATAC/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_45CB0L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.29.2021/42CE5L-ATAC_update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_4C2E5L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.03.29.2021/42CE5L-ATAC_update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_4C2E5L.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.05.03.2021/4D0D2L-ATAC/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_4D0D2L.tsv.gz
cp /datastore/nextgenout5/share/labs/francolab/scATAC-seq_Breast.05.03.2021/4D0D2L-ATAC/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_4D0D2L.tsv.gz.tbi

# Cell lines
cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-8-04_200803_NS500270_0354_AHHJKYBGXF/MCF7_ATAC_July2020_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_MCF7.tsv.gz
cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-8-04_200803_NS500270_0354_AHHJKYBGXF/MCF7_ATAC_July2020_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_MCF7.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-8-04_200803_NS500270_0354_AHHJKYBGXF/T47D_ATAC_July2020_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_T47D.tsv.gz
cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-8-04_200803_NS500270_0354_AHHJKYBGXF/T47D_ATAC_July2020_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_T47D.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-3-18_200314_NS500270_0347_AH5H22BGXF/HCC1143_Mar2020_ATAC_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_HCC1143.tsv.gz
cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-3-18_200314_NS500270_0347_AH5H22BGXF/HCC1143_Mar2020_ATAC_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_HCC1143.tsv.gz.tbi

cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-3-18_200314_NS500270_0347_AH5H22BGXF/SUM149PT_Mar2020_ATAC_Update/outs/fragments.tsv.gz ./cellranger-atac_outputs/fragments_SUM149PT.tsv.gz
cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-3-18_200314_NS500270_0347_AH5H22BGXF/SUM149PT_Mar2020_ATAC_Update/outs/fragments.tsv.gz.tbi ./cellranger-atac_outputs/fragments_SUM149PT.tsv.gz.tbi

# Copy cellranger filtered_feature_bc_matrices to location

# Patient samples
cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.02.02.2021/49758L-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_49758L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.03.30.2021/49CFCL-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_49CFCL/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.03.30.2021/4AF75L-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_4AF75L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.03.31.2021/HN4B146L_Mar2021_RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_4B146L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.05.28.2019/HT-35A4AL-RNA_MAY2019_G5_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_35A4AL/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.05.28.2019/HT-35EE8L-RNA_MAY2019_G6_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_35EE8L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.09.11.2019/3821AL-RNA_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_3821AL/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.01.15.2020/3B3E9L-RNA_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_3B3E9L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.02.07.2020/3C7D1L-RNA_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_3C7D1L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.03.19.2020/3D388L-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_3D388L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.03.19.2020/3FCDEL-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_3FCDEL/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.08.10.2020/43E7BL-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_43E7BL/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.08.10.2020/43E7CL-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_43E7CL/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.09.04.2020/44F0AL-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_44F0AL/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.02.02.2021/45CB0L-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_45CB0L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.03.31.2021/4C2E5L-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_4C2E5L/

cp /datastore/nextgenout5/share/labs/francolab/scRNA-seq_Breast.05.03.2021/4D0D2L-RNA/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_4D0D2L/

# Cell lines
cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-8-04_200803_NS500276_0188_AHHGGWBGXF/MCF7_RNA_July2020_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_MCF7/

cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-8-04_200803_NS500276_0188_AHHGGWBGXF/T47D_RNA_July2020_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_T47D/

cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-3-18_200314_NS500270_0346_AHTLGYBGXC/HCC1143_Mar2020_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_HCC1143/

cp /datastore/nextgenout5/share/labs/bioinformatics/seqware/hu_10x_nextseq_copy/hu_10x_2020-3-18_200314_NS500270_0346_AHTLGYBGXC/SUM149PT_Mar2020_Update/outs/filtered_feature_bc_matrix/* ./cellranger_outputs/filtered_feature_bc_matrix_SUM149PT/
