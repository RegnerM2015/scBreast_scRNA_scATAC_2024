# scBreast_scRNA_scATAC_2024
Repository to store code used for the analyses presented in "Defining the Regulatory Logic of Breast Cancer Using Single-Cell Epigenetic And Transcriptome Profiling"

### Data availability:
Processed data will be made publicly available at GEO (https://www.ncbi.nlm.nih.gov/geo/) after publication. 

Raw data will be available with controlled access at dbGaP (https://www.ncbi.nlm.nih.gov/gap/).

### Repository format:
The following bash scripts were run in this order to execute specific tasks under the each heading. More specifically, these are slurm scripts that either sourced external data or submitted an R script to run in a Docker container executed with Singularity in an HPC environment. 

Please visit the [wiki](https://github.com/RegnerM2015/scBreast_scRNA_scATAC_2024/wiki) for an in-depth walkthrough of these scripts and how they were used in the making of each figure. 

#### Source external data and process scRNA-seq data:
1. scripts/make_cellranger_directories.sh
2. scripts/copy_cellranger_outputs.sh
3. scripts/find_versions.sh
4. scripts/get_reference_data.sh
5. scripts/Sumbit-Individual_Samples_scRNA-QC_DoubletRemoval_Preprocessing.sh
6. scripts/Submit-Individual_Samples_scRNA-MultiKClustering.sh,
   
   scripts/Submit-Individual_Samples_scRNA-MultiKClustering_AlternateSeed.sh,

   scripts/Submit-Individual_Samples_scRNA-MultiKClustering_AlternateSeed_SecondAttempt.sh,

   scripts/Submit-Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt.sh
7. scripts/Submit-Individual_Samples_scRNA-FindClusterMarkerGenes.sh
8. scripts/Individual_Samples_scRNA-RemoveLowMappingRatePopulation_Reprocess_3FCDEL.sh
9. scripts/Wu_etal_2021_BRCA_scRNA-CreateSeuratObjectWithCCA.sh
10. scripts/Wu_etal_2021_BRCA_scRNA-CreateSeuratObjectWithOutCCA.sh
11. scripts/Submit-Individual_Samples_scRNA-CellTypeAnnotation.sh
12. scripts/Submit-Individual_Samples_scRNA-inferCNV_CancerCellDetection.sh
13. scripts/get_SCSubtype_training_data.sh
14. scripts/Submit-Individual_Samples_scRNA-SCSubtype_Classification.sh
15. scripts/Patient_Samples_scRNA-Merge_And_ReCluster.sh

#### Process scATAC-seq data:
16. scripts/All_Samples_scATAC-QC_DoubletRemoval_Preprocessing.sh
17. scripts/Patient_Samples_scATAC-Subset.sh
18. scripts/Patient_Samples_scATAC-DimReduc_GeneScoring.sh
19. scripts/Patient_Samples_scATAC-Transfer_Labels_from_scRNA.sh

#### Subset to Basal-like BC in scRNA-seq and scATAC-seq:
20. scripts/Basal_And_TN_Samples_scRNA-Subset_And_ReCluster-TESTING.sh
21. scripts/Basal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.sh
22. scripts/Basal_And_TN_Samples_scATAC-Transfer_Labels_from_scRNA_Call_Peaks-TESTING3.sh

#### Subset to Luminal BC in scRNA-seq and scATAC-seq:
23. scripts/Luminal_And_TN_Samples_scRNA-Subset_And_ReCluster-TESTING.sh
24. scripts/Luminal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.sh
25. scripts/Luminal_And_TN_Samples_scATAC-Transfer_Labels_from_scRNA_Call_Peaks-TESTING3.sh
    
#### Subset to BC cell lines in scRNA-seq and scATAC-seq:
26. scripts/CellLine_Samples_scRNA-Merge_And_ReCluster.sh
27. scripts/CellLine_Samples_scATAC-Subset_GeneScoring_DimReduc_TransferLabels_CallPeaks.sh

#### Perform peak-to-gene association analyses:
28. scripts/scLME_update-metacells-Basal-SingFits_OLS.sh
29. scripts/scLME_update-metacells-Luminal-SingFits_OLS.sh
30. scripts/scLME_update-metacells-SingFits_OLS-CellLines.sh

#### Generate visualizations and tabular data:
31. scripts/Full_Cohort_Results.sh
32. scripts/Basal_Cohort_Results.sh
33. scripts/Luminal_Cohort_Results.sh
34. scripts/Cell_Line_Cohort_Results.sh
35. scripts/Table_1_S1.sh
36. scripts/Supplemental_Tables-P2Gs.sh
37. scripts/Supplemental_Tables-barcode_metadata.sh

### Docker container to help replicate the computational environment:
R scripts were run in a [Docker](https://www.docker.com/resources/what-container/) container executed with [Singularity](https://github.com/sylabs/singularity) in an HPC environment.

To help users replicate our computational environment, we have uploaded the Docker image to DockerHub for public access: https://hub.docker.com/r/regnerm/scbreast_2023.

To pull the latest version of the Docker image, you may run the following command: 
```
docker pull regnerm/scbreast_2023:1.8.0
```
