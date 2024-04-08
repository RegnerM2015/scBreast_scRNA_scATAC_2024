################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: This script performs the following tasks  
#         1) Read in Swarbrick breast cancer dataset
#         2) Integrate via Seurat's CCA
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

list <- SplitObject(brca.ref, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(object = x, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)

  x <- FindVariableFeatures(object = x, selection.method = "vst", nfeatures = 2000,
                              loess.span = 0.3, clip.max = "auto",
                              num.bin = 20, binning.method = "equal_width", verbose = F)
})

anchors <- FindIntegrationAnchors(object.list = list, dims =1:30)

# this command creates an 'integrated' data assay
integrated <- IntegrateData(anchorset = anchors,dims=1:30)


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
# Scaling data
all.genes <- rownames(x = integrated)
integrated <- ScaleData(object = integrated, features = all.genes, verbose=F)

# Run PCA to reduce dimensions
integrated <- RunPCA(object = integrated, features = VariableFeatures(object = integrated), npcs = 50, verbose=F)

# Run UMAP
integrated <- RunUMAP(integrated,reduction = "pca", dims = 1:30,seed.use = 1,umap.method = "uwot",
               n.neighbors = 30L,
               metric = "cosine",
               learning.rate = 1,
               min.dist = 0.3)

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(integrated,"Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Wuetal2021_with_CCA.rds")
