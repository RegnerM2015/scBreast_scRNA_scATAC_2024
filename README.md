# scBreast_scRNA_scATAC_2024
Repository to store code used for the analyses presented in "Single-Cell ATAC and RNA Sequencing of Human Breast Cancer Reveals Salient Cancer-Specific Enhancers"

### Run these scripts in order while in the main directory:

#### Get data and process/QC individual scRNA datasets
1. scripts/make_cellranger_directories.sh
2. scripts/copy_cellranger_outputs.sh
3. scripts/find_versions.sh
4. scripts/get_reference_data.sh
5. scripts/Sumbit-Individual_Samples_scRNA-QC_DoubletRemoval_Preprocessing.sh
6. scripts/Submit-Individual_Samples_scRNA-MultiKClustering.sh, scripts/Submit-Individual_Samples_scRNA-MultiKClustering_AlternateSeed.sh, scripts/Submit-Individual_Samples_scRNA-MultiKClustering_AlternateSeed_SecondAttempt.sh, and scripts/Submit-Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt.sh
7. scripts/Submit-Individual_Samples_scRNA-FindClusterMarkerGenes.sh
8. scripts/Individual_Samples_scRNA-RemoveLowMappingRatePopulation_Reprocess_3FCDEL.sh
9. scripts/Submit-Individual_Samples_scRNA-CellTypeAnnotation.sh
10. scripts/Submit-Individual_Samples_scRNA-inferCNV_CancerCellDetection.sh
11. scripts/get_SCSubtype_training_data.sh
12. scripts/Submit-Individual_Samples_scRNA-SCSubtype_Classification.sh,
scripts/Submit-Individual_CellLine_Samples_scRNA-SCSubtype_Classification.sh

#### QC all scATAC data
13. scripts/All_Samples_scATAC-QC_DoubletRemoval_Preprocessing.sh

#### Subset out patient scATAC data and process 
14. scripts/Patient_Samples_scATAC-Subset.sh
15. scripts/Patient_Samples_scATAC-DimReduc_GeneScoring.sh

#### Create Patient scRNA dataset
16. scripts/Patient_Samples_scRNA-Merge_And_ReCluster.sh

#### Transfer labels between patient scRNA and patient scATAC
17. scripts/Patient_Samples_scATAC-Transfer_Labels_from_scRNA.sh

#### Create cell line scRNA dataset
18. scripts/CellLine_Samples_scRNA-Merge_And_ReCluster.sh

#### Subset out cell line scATAC data, process, and transfer labels from cell line scRNA dataset
19. scripts/CellLine_Samples_scATAC-Subset_GeneScoring_DimReduc_TransferLabels_CallPeaks.sh

#### Create Basal subtype scRNA dataset 
20. scripts/Basal_And_TN_Samples_scRNA-Subset_And_ReCluster-TESTING.sh

#### Subset out Basal scATAC data, process, and transfer labels from Basal scRNA dataset
21. scripts/Basal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.sh
22. scripts/Basal_And_TN_Samples_scATAC-Transfer_Labels_from_scRNA_Call_Peaks-TESTING3.sh

#### Create Luminal subtype scRNA dataset 
23. scripts/Luminal_And_TN_Samples_scRNA-Subset_And_ReCluster-TESTING.sh

#### Subset out Luminal scATAC data, process, and transfer labels from Luminal scRNA dataset
24. scripts/Luminal_And_TN_Samples_scATAC-Subset_DimReduc-TESTING3.sh
25. scripts/Luminal_And_TN_Samples_scATAC-Transfer_Labels_from_scRNA_Call_Peaks-TESTING3.sh

#### Perform P2G association analyses in Basal, Luminal, and cell lines
26. scripts/scLME_update-metacells-Basal-SingFits_OLS.sh
27. scripts/scLME_update-metacells-Basal_basalBackground-SingFits_OLS.sh
28. scripts/scLME_update-metacells-Luminal-SingFits_OLS.sh
29. scripts/scLME_update-metacells-SingFits_OLS-CellLines.sh

#### Generate visuals and tabular data 
30. scripts/Table_1_S1.sh
31. scripts/Full_Cohort_Results.sh
32. scripts/Basal_Cohort_Results.sh
33. scripts/Luminal_Cohort_Results.sh
34. scripts/Cell_Line_Cohort_Results.sh
35. scripts/Supplemental_Tables-P2Gs.sh
36. scripts/Supplemental_Tables-barcode_metadata.sh
