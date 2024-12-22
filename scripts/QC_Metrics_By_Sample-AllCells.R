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
library(patchwork)
source("scripts/Matrix.utils.R")

set.seed(1234)

# Make output folder
if(dir.exists("QC_Metrics_By_Sample-AllCells")){
  print("Directory QC_Metrics_By_Sample-AllCells already exists!")
}else{
  dir.create("QC_Metrics_By_Sample-AllCells")
}

# Read in first patient sample
IDs <- c("35A4AL",
         "4C2E5L",
         "35EE8L",
         "3821AL",
         "3B3E9L",
         "3C7D1L",
         "3D388L",
         "3FCDEL",
         "43E7BL",
         "43E7CL",
         "44F0AL",
         "45CB0L",
         "49758L",
         "49CFCL",
         "4AF75L",
         "4B146L")

counts <- Read10X(data.dir = paste0("cellranger_outputs/filtered_feature_bc_matrix_",IDs[1]))
rna.aggr <- CreateSeuratObject(counts = counts, 
                          min.cells = 0, 
                          min.features = 0,
                          project = IDs[1])

for(i in 2:length(IDs)){
  counts <- Read10X(data.dir = paste0("cellranger_outputs/filtered_feature_bc_matrix_",IDs[i]))
  rna.new <- CreateSeuratObject(counts = counts, 
                            min.cells = 0, 
                            min.features = 0,
                            project = IDs[i])
  
  rna.aggr <- merge(x=rna.aggr,
                    y=rna.new)
  
}
# Rename aggregated Seurat object
rna <- rna.aggr

# Plot QC metrics stratified by sample

df <- as.data.frame(rna@meta.data)

df$Sample <- factor(df$orig.ident,levels=rev(c("49758L",
                                                "49CFCL",
                                                "4AF75L",
                                                "4B146L",
                                                "35A4AL",
                                                "4C2E5L",
                                                "35EE8L",
                                                "3821AL",
                                                "3B3E9L",
                                                "3C7D1L",
                                                "3FCDEL",
                                                "44F0AL",
                                                "3D388L",
                                                "43E7BL",
                                                "43E7CL",
                                                "45CB0L"
                                                )))

library(RColorBrewer)
cols <- RColorBrewer::brewer.pal(12,"Paired")
cols <- rev(c("#B4B4B4","#828282","#505050","#323232",cols[1:6],cols[8],cols[10],cols[7],"#a97ac2","#ccb4d9",cols[12]))

qcRNAnCountByPatientSample <- ggplot(df, aes(y=Sample,
                              x=log10(nCount_RNA),
                              fill=Sample)) + 
  scale_fill_manual(values = cols) + 
  geom_boxplot()+
  theme_bw()+
  NoLegend()+
  ggtitle(paste("PRE-QC log10(RNA counts) by sample","\n","n=",nrow(df)," cells"))

qcRNAnFeatureByPatientSample <- ggplot(df, aes(y=Sample,
                                      x=nFeature_RNA,
                                      fill=Sample)) + 
  scale_fill_manual(values = cols) + 
  geom_boxplot()+
  theme_bw()+
  ggtitle(paste("PRE-QC RNA features by sample","\n","n=",nrow(df)," cells"))


# Read in first cell line sample
IDs <- c("HCC1143",
         "SUM149PT",
         "MCF7",
         "T47D")

counts <- Read10X(data.dir = paste0("cellranger_outputs/filtered_feature_bc_matrix_",IDs[1]))
rna.aggr <- CreateSeuratObject(counts = counts, 
                               min.cells = 0, 
                               min.features = 0,
                               project = IDs[1])

for(i in 2:length(IDs)){
  counts <- Read10X(data.dir = paste0("cellranger_outputs/filtered_feature_bc_matrix_",IDs[i]))
  rna.new <- CreateSeuratObject(counts = counts, 
                                min.cells = 0, 
                                min.features = 0,
                                project = IDs[i])
  
  rna.aggr <- merge(x=rna.aggr,
                    y=rna.new)
  
}
# Rename aggregated Seurat object
rna <- rna.aggr

# Plot QC metrics stratified by sample

df <- as.data.frame(rna@meta.data)

cols <- rev(c("#662377","#A82973","#EF5064","#FC867D"))

df$Sample <- factor(df$orig.ident,levels=rev(c("HCC1143",
                                               "SUM149PT",
                                               "MCF7",
                                               "T47D")))

qcRNAnCountByCellLine <- ggplot(df, aes(y=Sample,
                                      x=log10(nCount_RNA),
                                      fill=Sample)) + 
  scale_fill_manual(values = cols) + 
  geom_boxplot()+
  theme_bw()+
  NoLegend()+
  ggtitle(paste("PRE-QC log10(RNA counts) by sample","\n","n=",nrow(df)," cells"))

qcRNAnFeatureByCellLine <- ggplot(df, aes(y=Sample,
                                        x=nFeature_RNA,
                                        fill=Sample)) + 
  scale_fill_manual(values = cols) + 
  geom_boxplot()+
  theme_bw()+
  ggtitle(paste("PRE-QC RNA features by sample","\n","n=",nrow(df)," cells"))

# Plot boxplots from patients and cell lines (All Cells called by cellranger)
qcRNAnCountByPatientSample+qcRNAnFeatureByPatientSample+qcRNAnCountByCellLine+qcRNAnFeatureByCellLine+plot_layout(ncol=2)
ggsave("./QC_Metrics_By_Sample-AllCells/QC_Boxplots_RNA_Patients_and_CellLines.pdf",width = 9,height = 8)
