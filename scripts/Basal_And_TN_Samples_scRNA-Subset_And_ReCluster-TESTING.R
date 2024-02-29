################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: Merge Basal_TN Samples and Re-cluster
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
library(ArchR)
# options(future.globals.maxSize = 6000 * 1024^2)
# plan("multicore")
set.seed(321)

ID <- "Basal_TN_Subset-TESTING"

# Read in scRNA
rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

# Subset
rna <- rna[,rna$RNA_snn_res.0.4 %in% c("11",
                                       "14",
                                       "4",
                                       "7")]
normals <- c("49758L","4AF75L","4B146L","49CFCL")
rna$From_Normal_Sample <- ifelse(rna$orig.ident %in% normals,"Yes","No")

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 %in% c("4","7") & rna$From_Normal_Sample == "No","Drop","Keep" )

rna <- rna[,rna$Keep.Cells != "Drop"]

basals <- c("35A4AL","4C2E5L")
rna$From_Basal_Sample <- ifelse(rna$orig.ident %in% basals,"Yes","No")

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 %in% c("11","14") & rna$From_Basal_Sample == "No","Drop","Keep" )

rna <- rna[,rna$Keep.Cells != "Drop"]

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 %in% c("11","14") & rna$SCSubtype != "Basal_SC","Drop","Keep" )
rna <- rna[,rna$Keep.Cells != "Drop"]

sub <- rna[,rna$RNA_snn_res.0.4 %in% c("4","7")]
table(sub$orig.ident)

sub <- rna[,rna$RNA_snn_res.0.4 %in% c("11","14")]
table(sub$orig.ident)

table(rna$SCSubtype)

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

# Perform clustering
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna, resolution = 0.015)

if(dir.exists("Basal_TN_Cohort_scRNA-TESTING")==TRUE){
  print("Directory Basal_TN_Cohort_scRNA-TESTING already exists!")
}else{
  dir.create("Basal_TN_Cohort_scRNA-TESTING")
}

DimPlot(rna,group.by = "SCSubtype")
ggsave(paste0("./Basal_TN_Cohort_scRNA-TESTING/UMAP_SCSubtype_",
              ID,".pdf"),
       width = 8,
       height = 8)
DimPlot(rna,group.by = "orig.ident")
ggsave(paste0("./Basal_TN_Cohort_scRNA-TESTING/UMAP_Sample_",
              ID,".pdf"),
       width = 8,
       height = 8)
DimPlot(rna,group.by = "RNA_snn_res.0.015")
ggsave(paste0("./Basal_TN_Cohort_scRNA-TESTING/UMAP_Cluster_",
              ID,".pdf"),
       width = 8,
       height = 8)
# Begin plotting
meta <- rna@meta.data
df <- meta %>% dplyr::group_by_at("orig.ident") %>% dplyr::count(SCSubtype)
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
ggsave(paste0("./Basal_TN_Cohort_scRNA-TESTING/SCSubtype_Proportion_BarChart_",
              ID,".pdf"),
       width = 6,
       height = 8)

meta <- rna@meta.data
df <- meta %>% dplyr::group_by_at("RNA_snn_res.0.015") %>% dplyr::count(orig.ident)
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
ggsave(paste0("./Basal_TN_Cohort_scRNA-TESTING/Cluster_Proportion_BarChart_",
              ID,"-Sample.pdf"),
       width = 6,
       height = 8)

meta <- rna@meta.data
df <- meta %>% dplyr::group_by_at("RNA_snn_res.0.015") %>% dplyr::count(SCSubtype)
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
ggsave(paste0("./Basal_TN_Cohort_scRNA-TESTING/Cluster_Proportion_BarChart_",
              ID,"-SCSubtype.pdf"),
       width = 6,
       height = 8)

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,".rds"))

