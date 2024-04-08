################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: This script performs the following tasks  
#         1) Read in Swarbrick breast cancer dataset
#         2) Process aggregated Swarbrick dataset without CCA
################################################################################
library(Seurat)
library(tidyverse)
library(parallel)
library(future)
options(future.globals.maxSize = 6000 * 1024^2)
plan("multicore")
set.seed(1)

# Read in data
brca.ref <- Read10X(data.dir = "Wu_etal_2021_BRCA_scRNASeq",gene.column=1)
meta <- read.csv("Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
rownames(meta) <- meta$X
length(which(rownames(meta)==colnames(brca.ref)))
brca.ref <- CreateSeuratObject(brca.ref,meta.data = meta)

# Process data
# normalizing the data
brca.ref <- NormalizeData(object = brca.ref, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)

# Find variable features
brca.ref <- FindVariableFeatures(object = brca.ref, selection.method = "vst", nfeatures = 2000,
                            loess.span = 0.3, clip.max = "auto",
                            num.bin = 20, binning.method = "equal_width", verbose = F)

# Scaling data
all.genes <- rownames(x = brca.ref)
brca.ref <- ScaleData(object = brca.ref, features = all.genes, verbose=F)

# Run PCA to reduce dimensions
brca.ref <- RunPCA(object = brca.ref, features = VariableFeatures(object = brca.ref), npcs = 50, verbose=F)

# Run UMAP
brca.ref <- RunUMAP(brca.ref,reduction = "pca", dims = 1:30,seed.use = 1,umap.method = "uwot",
               n.neighbors = 30L,
               metric = "cosine",
               learning.rate = 1,
               min.dist = 0.3)

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(brca.ref,"Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Wuetal2021_without_CCA.rds")
