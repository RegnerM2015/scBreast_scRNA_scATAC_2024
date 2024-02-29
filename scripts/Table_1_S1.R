###################################################
# Matt Regner
# Franco Lab 
###################################################
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(ArchR)
library(stringr)
library(stringi)
library(RColorBrewer)
library(dplyr)
library(parallel)
library(ArchR)
library(tidyr)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(ComplexHeatmap)
library(msigdbr)
library(EnhancedVolcano)
library(readxl)

# Make figure output folder
if(dir.exists("Table_1_S1_outputs")==TRUE){
  print("Directory Table_1_S1_outputs already exists!")
}else{
  dir.create("Table_1_S1_outputs")
}

# Read in clinical data
clinical <- read_excel("miscellaneous/Franco-Perou-SingleCellBreastCancerDataset-July2022_MR-MJR_KW.xlsx")
clinical$Patient_Age_at_Surgery <- as.numeric(clinical$Patient_Age_at_Surgery)
clinical$BMI <- as.numeric(clinical$BMI)
colnames(clinical)[c(10,11,16)] <- c("ER_IHC","PR_IHC","Menopause_Status")
clinical$Menopause_Status <- gsub("\\ .*","",clinical$Menopause_Status)
clinical <- clinical[,-grep("Gender",colnames(clinical))]

clinical[,grep("IHC",colnames(clinical))[1:3]][ clinical[,grep("IHC",colnames(clinical))[1:3]] == "NA" ] <- "Not Performed"
clinical[,grep("IHC",colnames(clinical))[4]][ clinical[,grep("IHC",colnames(clinical))[4]] == "NA" ] <- "Not Determined"

clinical[,grep("FISH",colnames(clinical))][ clinical[,grep("FISH",colnames(clinical))] == "NA" ] <- "Not Performed"

clinical$Type <- c(rep("Primary Breast Tumor",11),rep("Normal Breast Tissue",4),"Primary Breast Tumor")

clinical$Patient <- c(5,7,8,9,10,13,11,14,14,12,15,1,2,3,4,6)
clinical <- dplyr::arrange(clinical,Patient)
clinical <- clinical[,c("Sample_ID","Patient","Type",colnames(clinical)[-c(1,ncol(clinical),ncol(clinical)-1)])]

fwrite(clinical,"./Table_1_S1_outputs/Supplemental_Table_1.csv",sep = ",",col.names = TRUE)

# Make Main Table 1
main <- clinical[,-c(4:10,14:ncol(clinical))]
main$Subject <- paste0("Patient ",main$Patient)
main <- main[,c("Sample_ID","Subject","Type","ER_IHC","PR_IHC","Her2_IHC")]

cell_lines <- data.frame(Sample_ID = c("MCF7","T47D","HCC1143","SUM149PT"),
                         Subject = c("MCF7","T47D","HCC1143","SUM149PT"),
                         Type = rep("Cell Line",4),
                         ER_IHC=rep("Not Performed",4),
                         PR_IHC=rep("Not Performed",4),
                         Her2_IHC=rep("Not Performed",4))

main <- rbind(main,cell_lines)

main$ER_IHC <- ifelse(main$ER_IHC == "Positive (>10% or stated as positive)", "+", main$ER_IHC)
main$ER_IHC <- ifelse(main$ER_IHC == "Negative (<1% or stated as negative)", "-", main$ER_IHC)

main$PR_IHC <- ifelse(main$PR_IHC == "Positive (>10% or stated as positive)", "+", main$PR_IHC)
main$PR_IHC <- ifelse(main$PR_IHC == "Negative (<1% or stated as negative)", "-", main$PR_IHC)

main$Her2_IHC <- ifelse(main$Her2_IHC == "Negative", "-", main$Her2_IHC)

# report number of scRNA cells
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
patients <- table(rna$orig.ident)[order(match(names(table(rna$orig.ident)),main$Sample_ID[1:16]))]

rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_CellLine_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
cellLines <- table(rna$orig.ident)[order(match(names(table(rna$orig.ident)),main$Sample_ID[17:20]))]

all.equal(names(c(patients,cellLines)),main$Sample_ID)
main$scRNA <- c(as.numeric(patients),as.numeric(cellLines))

# report number of scATAC cells
atac <- loadArchRProject(path =  "./Patient_Samples_scATAC")
patients <- table(atac$Sample)[order(match(names(table(atac$Sample)),main$Sample_ID[1:16]))]

atac <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")
cellLines <- table(atac$Sample)[order(match(names(table(atac$Sample)),main$Sample_ID[17:20]))]

all.equal(names(c(patients,cellLines)),main$Sample_ID)
main$scATAC <- c(as.numeric(patients),as.numeric(cellLines))

fwrite(main,"./Table_1_S1_outputs/Table_1_main.csv",sep = ",",col.names = TRUE)
