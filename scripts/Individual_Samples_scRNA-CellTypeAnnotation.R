################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: This script performs the following tasks  
#         1) Annotate MultiK clusters
################################################################################
library(Seurat)
library(tidyverse)
library(dplyr)
library(sigclust)
library(MultiK)
library(DoubletFinder)
library(S4Vectors)
library(dendextend)
library(ggdendro)
library(stringr)
library(stringi)
library(parallel)
library(future)
library(grDevices)
options(future.globals.maxSize = 6000 * 1024^2)
plan("multicore")
set.seed(1)

args=(commandArgs(TRUE))
ID <- args[[1]]
print(ID)

rna <- readRDS(paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering.rds"))
Idents(object = rna) <- "MultiK_clusters"

markers <- readRDS(paste0("./DifferentialExpression_ClusterAnnotation_scRNA/WilcoxDEGS_",ID,".rds"))

# Check for Mast cell clusters
idx.tpsb2 <- grep("TPSB2",markers$gene)
idx.tpsab1 <- grep("TPSAB1",markers$gene)

if (length(idx.tpsb2) >0 | length(idx.tpsab1) >0){
  markers.tpsb2 <- markers[idx.tpsb2,]
  cluster.tpsb2 <- as.character(markers.tpsb2$cluster)
  
  markers.tpsab1 <- markers[idx.tpsab1,]
  cluster.tpsab1 <- as.character(markers.tpsab1$cluster)
  
  mast.clusters <- c(cluster.tpsab1,cluster.tpsb2)
  
  rna@meta.data$in_MastCell_cluster <- ifelse(rna@meta.data$MultiK_clusters %in% mast.clusters,"Yes","No")
  
}else{
  markers <- markers
  rna$in_MastCell_cluster <- "No"
}

# Calculate Enrichment of cell type gene signatures per cluster from panglaodb
tsv=gzfile("./miscellaneous/PanglaoDB_markers_27_Mar_2020.tsv")  
panglaodb <- read.csv(tsv,header=T,sep = "\t") 
panglaodb <- dplyr::filter(panglaodb,species == "Hs" | species == "Mm Hs")# Human subset 
panglaodb <- split(as.character(panglaodb$official.gene.symbol), panglaodb$cell.type)

# Epithelial 
rna <- AddModuleScore(rna,features = list(panglaodb$`Epithelial cells`,
                                          panglaodb$`Mammary epithelial cells`,
                                          panglaodb$`Myoepithelial cells`),
                      name = c("Panglaodb_Epithelial_cells.",
                                "Panglaodb_Mammary_epithelial_cells.",
                                "Panglaodb_Myoepithelial_cells."),search = T)
# B cells
rna <- AddModuleScore(rna,features = list(panglaodb$`B cells`,
                                          panglaodb$`B cells memory`,
                                          panglaodb$`B cells naive`),
                      name = c("Panglaodb_B_cells.",
                                "Panglaodb_B_cells_memory.",
                                "Panglaodb_B_cells_naive."),search = T)

# T/NK cells
rna <- AddModuleScore(rna,features = list(panglaodb$`T cells`,
                                          panglaodb$`T cells naive`,
                                          panglaodb$`T cytotoxic cells`,
                                          panglaodb$`T helper cells`,
                                          panglaodb$`T memory cells`,
                                          panglaodb$`T regulatory cells`,
                                          panglaodb$`NK cells`,
                                          panglaodb$`Natural killer T cells`),
                      name = c("Panglaodb_T_cells.",
                                "Panglaodb_T_cells_naive.",
                                "Panglaodb_T_cytotoxic_cells.",
                                "Panglaodb_T_helper_cells.",
                                "Panglaodb_T_memory_cells.",
                                "Panglaodb_T_regulatory_cells.",
                                "Panglaodb_NK_cells.",
                                "Panglaodb_Natural_killer_T_cells."),search = T)

# Myeloid
rna <- AddModuleScore(rna,features = list(panglaodb$Macrophages,
                                          panglaodb$Monocytes,
                                          panglaodb$`Mast cells`,
                                          panglaodb$`Dendritic cells`,
                                          panglaodb$Neutrophils),
                      name = c("Panglaodb_Macrophages.",
                                "Panglaodb_Monocytes.",
                                "Panglaodb_Mast_cells.",
                                "Panglaodb_Dendritic_cells.",
                                "Panglaodb_Neutrophils."),search = T)

# Vascular
rna <- AddModuleScore(rna,features = list(panglaodb$`Endothelial cells`,
                                          panglaodb$Pericytes),
                      name = c("Panglaodb_Endothelial_cells.",
                                "Panglaodb_Pericytes."),search = T)

# Stromal/Mesenchymal
rna <- AddModuleScore(rna,features = list(panglaodb$Fibroblasts,
                                          panglaodb$Myofibroblasts,
                                          panglaodb$`Smooth muscle cells`,
                                          panglaodb$Adipocytes),
                      name = c("Panglaodb_Fibroblasts.",
                                "Panglaodb_Myofibroblasts.",
                                "Panglaodb_Smooth_muscle_cells.",
                                "Panglaodb_Adipocytes."),search = T)

# Plot stacked violins
source("./scripts/stacked_violin.R")

StackedVlnPlot(rna,features = c("Panglaodb_Epithelial_cells.1",
                                "Panglaodb_Mammary_epithelial_cells.2",
                                "Panglaodb_Myoepithelial_cells.3",
                                "Panglaodb_B_cells.1",
                                "Panglaodb_B_cells_memory.2",
                                "Panglaodb_B_cells_naive.3",
                                "Panglaodb_T_cells.1",
                                "Panglaodb_T_cells_naive.2",
                                "Panglaodb_T_cytotoxic_cells.3",
                                "Panglaodb_T_helper_cells.4",
                                "Panglaodb_T_memory_cells.5",
                                "Panglaodb_T_regulatory_cells.6",
                                "Panglaodb_NK_cells.7",
                                "Panglaodb_Natural_killer_T_cells.8",
                                "Panglaodb_Macrophages.1",
                                "Panglaodb_Monocytes.2",
                                "Panglaodb_Mast_cells.3",
                                "Panglaodb_Dendritic_cells.4",
                                "Panglaodb_Neutrophils.5",
                                "Panglaodb_Endothelial_cells.1",
                                "Panglaodb_Pericytes.2",
                                "Panglaodb_Fibroblasts.1",
                                "Panglaodb_Myofibroblasts.2",
                                "Panglaodb_Smooth_muscle_cells.3",
                                "Panglaodb_Adipocytes.4"))
ggsave(paste0("./DifferentialExpression_ClusterAnnotation_scRNA/SignatureEnrichment_Panglaodb_",
              ID,".pdf"),width =12,height = 48)

# Add reference-based annotations

# CCA corrected reference (celltype_minor)
ref.cca <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Wuetal2021_with_CCA.rds")
anchors <- FindTransferAnchors(reference = ref.cca, query = rna,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref.cca$celltype_minor,
                            dims = 1:30)
rna <- AddMetaData(rna, metadata = predictions[,c("predicted.id","prediction.score.max")])
rna$GSE176078_ref_with_CCA_predicted.id <- rna$predicted.id
rna$GSE176078_ref_with_CCA_prediction.score.max <- rna$prediction.score.max

Idents(object = rna) <- "MultiK_clusters"
# Assign celltype label to each cluster based on majority label
cells <- rna@meta.data %>% dplyr::group_by(MultiK_clusters) %>% dplyr::count(GSE176078_ref_with_CCA_predicted.id) 

cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in factor(cells$MultiK_clusters)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(MultiK_clusters ==i) %>% arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$GSE176078_ref_with_CCA_predicted.id[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$celltypeLabel <- Idents(rna)

rna$Cluster_Majority_CCA_Reference_Cell_Type_Minor <- paste0(rna$MultiK_clusters,"-",rna$celltypeLabel)

# Edit Mast cell annotations
rna$Cluster_Majority_CCA_Reference_Cell_Type_Minor <- ifelse(rna$in_MastCell_cluster == "Yes",
                                                       paste0(rna$MultiK_clusters,"-Mast cell"),
                                                       rna$Cluster_Majority_CCA_Reference_Cell_Type_Minor)

# CCA corrected reference (celltype_major)
anchors <- FindTransferAnchors(reference = ref.cca, query = rna,
                               dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref.cca$celltype_major,
                            dims = 1:30)
rna <- AddMetaData(rna, metadata = predictions[,c("predicted.id","prediction.score.max")])
rna$GSE176078_ref_with_CCA_predicted.id <- rna$predicted.id
rna$GSE176078_ref_with_CCA_prediction.score.max <- rna$prediction.score.max

Idents(object = rna) <- "MultiK_clusters"
# Assign celltype label to each cluster based on majority label
cells <- rna@meta.data %>% dplyr::group_by(MultiK_clusters) %>% dplyr::count(GSE176078_ref_with_CCA_predicted.id) 

cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in factor(cells$MultiK_clusters)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(MultiK_clusters ==i) %>% arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$GSE176078_ref_with_CCA_predicted.id[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$celltypeLabel <- Idents(rna)

rna$Cluster_Majority_CCA_Reference_Cell_Type_Major <- paste0(rna$MultiK_clusters,"-",rna$celltypeLabel)

# Edit Mast cell annotations
rna$Cluster_Majority_CCA_Reference_Cell_Type_Major <- ifelse(rna$in_MastCell_cluster == "Yes",
                                                             paste0(rna$MultiK_clusters,"-Mast cell"),
                                                             rna$Cluster_Majority_CCA_Reference_Cell_Type_Major)

# No batch correction reference (celltype_minor)
ref <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Wuetal2021_without_CCA.rds")
anchors <- FindTransferAnchors(reference = ref, query = rna,
                               dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$celltype_minor,
                            dims = 1:30)
rna <- AddMetaData(rna, metadata = predictions[,c("predicted.id","prediction.score.max")])
rna$GSE176078_ref_without_CCA_predicted.id <- rna$predicted.id
rna$GSE176078_ref_without_CCA_prediction.score.max <- rna$prediction.score.max

Idents(object = rna) <- "MultiK_clusters"
# Assign celltype label to each cluster based on majority label
cells <- rna@meta.data %>% dplyr::group_by(MultiK_clusters) %>% dplyr::count(GSE176078_ref_without_CCA_predicted.id) 

cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in factor(cells$MultiK_clusters)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(MultiK_clusters ==i) %>% arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$GSE176078_ref_without_CCA_predicted.id[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$celltypeLabel <- Idents(rna)

rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor <- paste0(rna$MultiK_clusters,"-",rna$celltypeLabel)

# Edit Mast cell annotations
rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor <- ifelse(rna$in_MastCell_cluster == "Yes",
                                                             paste0(rna$MultiK_clusters,"-Mast cell"),
                                                             rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor)

# No batch correction reference (celltype_major)
anchors <- FindTransferAnchors(reference = ref, query = rna,
                               dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$celltype_major,
                            dims = 1:30)
rna <- AddMetaData(rna, metadata = predictions[,c("predicted.id","prediction.score.max")])
rna$GSE176078_ref_without_CCA_predicted.id <- rna$predicted.id
rna$GSE176078_ref_without_CCA_prediction.score.max <- rna$prediction.score.max

Idents(object = rna) <- "MultiK_clusters"
# Assign celltype label to each cluster based on majority label
cells <- rna@meta.data %>% dplyr::group_by(MultiK_clusters) %>% dplyr::count(GSE176078_ref_without_CCA_predicted.id) 

cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in factor(cells$MultiK_clusters)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(MultiK_clusters ==i) %>% arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$GSE176078_ref_without_CCA_predicted.id[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$celltypeLabel <- Idents(rna)

rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Major <- paste0(rna$MultiK_clusters,"-",rna$celltypeLabel)

# Edit Mast cell annotations
rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Major <- ifelse(rna$in_MastCell_cluster == "Yes",
                                                                paste0(rna$MultiK_clusters,"-Mast cell"),
                                                                rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Major)

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering-CellTypeAnnotations.rds"))
