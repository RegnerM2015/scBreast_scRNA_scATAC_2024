################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: Merge Patient Samples and Re-cluster
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

ID <- "Patient_Cohort"

HT35A4AL <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_35A4AL-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT35EE8L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_35EE8L-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT3821AL <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_3821AL-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT3B3E9L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_3B3E9L-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT3C7D1L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_3C7D1L-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT3D388L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_3D388L-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT3FCDEL <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_3FCDEL-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT43E7BL <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_43E7BL-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT43E7CL <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_43E7CL-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT44F0AL <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_44F0AL-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT45CB0L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_45CB0L-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HT4C2E5L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_4C2E5L-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
#HT4D0D2L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_4D0D2L-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
HN49758L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_49758L-MultiKClustering-CellTypeAnnotations.rds")
HN49CFCL <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_49CFCL-MultiKClustering-CellTypeAnnotations.rds")
HN4AF75L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_4AF75L-MultiKClustering-CellTypeAnnotations.rds")
HN4B146L <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_4B146L-MultiKClustering-CellTypeAnnotations.rds")

# Merge Seurat objects
rna <- merge(x = HT35A4AL,
             y = list(HT35EE8L,
                      HT3821AL,
                      HT3B3E9L,
                      HT3C7D1L,
                      HT3D388L,
                      HT3FCDEL,
                      HT43E7BL,
                      HT43E7CL,
                      HT44F0AL,
                      HT45CB0L,
                      HT4C2E5L,
                      HN49758L,
                      HN49CFCL,
                      HN4AF75L,
                      HN4B146L))

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

if(dir.exists("Patient_Cohort_scRNA")==TRUE){
  print("Directory Patient_Cohort_scRNA already exists!")
}else{
  dir.create("Patient_Cohort_scRNA")
}

DimPlot(rna,group.by = "SCSubtype")
ggsave(paste0("./Patient_Cohort_scRNA/UMAP_SCSubtype_",
              ID,".pdf"),
       width = 8,
       height = 8)
DimPlot(rna,group.by = "RNA_snn_res.0.4")
ggsave(paste0("./Patient_Cohort_scRNA/UMAP_res04_clusters_",
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
ggsave(paste0("./Patient_Cohort_scRNA/SCSubtype_Proportion_BarChart_",
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
ggsave(paste0("./Patient_Cohort_scRNA/Sample_Proportion_BarChart_",
              ID,".pdf"),
       width = 6,
       height = 8)

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds"))

