################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: Merge CellLine Samples and Re-cluster
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
library(harmony)
library(grDevices)
# options(future.globals.maxSize = 6000 * 1024^2)
# plan("multicore")
set.seed(321)

ID <- "CellLine_Cohort"

SUM149PT <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_SUM149PT-MultiKClustering-SCSubtype.rds")
T47D <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_T47D-MultiKClustering-SCSubtype.rds")
MCF7 <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_MCF7-MultiKClustering-SCSubtype.rds")
HCC1143 <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_HCC1143-MultiKClustering-SCSubtype.rds")

# Merge Seurat objects
rna <- merge(x = SUM149PT,
             y = list(T47D,
                      MCF7,
                      HCC1143))

#normalizing the data
rna <- NormalizeData(object = rna, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)

# Find variable features
rna <- FindVariableFeatures(object = rna, selection.method = "vst", nfeatures = 2000,
                            loess.span = 0.3, clip.max = "auto",
                            num.bin = 20, binning.method = "equal_width", verbose = F)

# Scaling data
rna <- ScaleData(object = rna, verbose=F)

# Run PCA to reduce dimensions
rna <- RunPCA(object = rna, features = VariableFeatures(object = rna), npcs = 50, verbose=F)

# Run UMAP
rna <- RunUMAP(rna,reduction = "pca", dims = 1:30,seed.use = 1,umap.method = "uwot",
               n.neighbors = 30L,
               metric = "cosine",
               learning.rate = 1,
               min.dist = 0.3)
# Run clustering
rna <- FindNeighbors(rna,reduction = "pca",dims = 1:30)
rna <- FindClusters(rna,resolution = 0.4)

if(dir.exists("CellLine_Cohort_scRNA")==TRUE){
  print("Directory CellLine_Cohort_scRNA already exists!")
}else{
  dir.create("CellLine_Cohort_scRNA")
}

DimPlot(rna,group.by = "SCSubtype")
ggsave(paste0("./CellLine_Cohort_scRNA/UMAP_SCSubtype_",
              ID,".pdf"),
       width = 8,
       height = 8)
DimPlot(rna,group.by = "RNA_snn_res.0.4")
ggsave(paste0("./CellLine_Cohort_scRNA/UMAP_res04_clusters_",
              ID,".pdf"),
       width = 8,
       height = 8)
# Begin plotting
meta <- rna@meta.data
df <- meta %>% dplyr::group_by_at("RNA_snn_res.0.4") %>% dplyr::count(SCSubtype)
colnames(df) <- c("Cluster","SCSubtype","Cells")
df %>%
  ggplot(aes(fill=SCSubtype, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()
ggsave(paste0("./CellLine_Cohort_scRNA/SCSubtype_Proportion_BarChart_",
              ID,".pdf"),
       width = 6,
       height = 8)

meta <- rna@meta.data
df <- meta %>% dplyr::group_by_at("RNA_snn_res.0.4") %>% dplyr::count(orig.ident)
colnames(df) <- c("Cluster","Sample","Cells")
df %>%
  ggplot(aes(fill=Sample, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()
ggsave(paste0("./CellLine_Cohort_scRNA/Sample_Proportion_BarChart_",
              ID,".pdf"),
       width = 6,
       height = 8)

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds"))

