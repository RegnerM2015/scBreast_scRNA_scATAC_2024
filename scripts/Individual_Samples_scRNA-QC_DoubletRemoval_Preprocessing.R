################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: This script performs the following tasks  
#         1) scRNA-seq quality control and doublet removal
#         2) Standard Seurat pre-processing and dim reduction
################################################################################
library(Seurat)
library(tidyverse)
library(dplyr)
library(DoubletFinder)
library(S4Vectors)
library(stringr)
library(stringi)
library(parallel)
library(future)
options(future.globals.maxSize = 6000 * 1024^2)
plan("multicore")
set.seed(1)

args=(commandArgs(TRUE))
ID <- args[[1]]
print(ID)

counts <- Read10X(data.dir = paste0("cellranger_outputs/filtered_feature_bc_matrix_",ID))
rna <- CreateSeuratObject(counts = counts, 
                          min.cells = 3, 
                          min.features = 0,
                          project = ID)

# Add mitochondrial percentage
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

if(dir.exists("QualityControl_scRNA") == TRUE){
  print("Directory QualityControl_scRNA already exists!")
}else{
  dir.create("QualityControl_scRNA")
}

# Visualize QC metrics as a violin plot
plot1 <- VlnPlot(rna, features = "nFeature_RNA",pt.size = 0)
plot1+geom_hline(yintercept = 500)
ggsave(paste0("./QualityControl_scRNA/VlnPlot_nFeature_RNA_",ID,".pdf"),
       width = 6,
       height = 8)
plot2 <- VlnPlot(rna, features = "nCount_RNA",pt.size = 0)
plot2+geom_hline(yintercept = 1000)
ggsave(paste0("./QualityControl_scRNA/VlnPlot_nCount_RNA_",ID,".pdf"),
       width = 6,
       height = 8)
plot3 <- VlnPlot(rna, features = "percent.mt",pt.size = 0)
plot3+geom_hline(yintercept = 20)
ggsave(paste0("./QualityControl_scRNA/VlnPlot_percent.mt_",ID,".pdf"),
       width = 6,
       height = 8)

# Visualize QC metrics as scatter plots
plot1 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1+geom_hline(yintercept = 20)+geom_vline(xintercept = 1000)
ggsave(paste0("./QualityControl_scRNA/FeatureScatter_nCount_RNA_percent.mt_",
              ID,".pdf"),
       width = 12,
       height = 8)
plot2 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2+geom_hline(yintercept = 500)+geom_vline(xintercept = 1000)
ggsave(paste0("./QualityControl_scRNA/FeatureScatter_nCount_RNA_nFeature_RNA_",
              ID,".pdf"),
       width = 12,
       height = 8)

# Apply constant QC filters described in Slyper et al. (Title: A single-cell and 
# single-nucleus RNA-Seq toolbox for fresh and frozen human tumors)
rna <- subset(rna, 
              subset = nFeature_RNA >=500 & nCount_RNA >=1000 & percent.mt <=20)

# normalizing the data
rna <- NormalizeData(object = rna, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)

# Find variable features
rna <- FindVariableFeatures(object = rna, selection.method = "vst", nfeatures = 2000,
                             loess.span = 0.3, clip.max = "auto",
                             num.bin = 20, binning.method = "equal_width", verbose = F)

# Scaling data
all.genes <- rownames(x = rna)
rna <- ScaleData(object = rna, features = all.genes, verbose=F)

# Run PCA to reduce dimensions
rna <- RunPCA(object = rna, features = VariableFeatures(object = rna), npcs = 50, verbose=F)

# Run UMAP
rna <- RunUMAP(rna,reduction = "pca", dims = 1:30,seed.use = 1,umap.method = "uwot",
               n.neighbors = 30L,
               metric = "cosine",
               learning.rate = 1,
               min.dist = 0.3)

# Doublet detection with DoubletFinder
## pK Identification (no ground-truth) -----------------------------------------
sweep.res.list <- paramSweep_v3(rna, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- as.numeric(unfactor(dplyr::arrange(bcmvn,desc(BCmetric))$pK[1]))
print(pK)
doublet.table <- read.table("./miscellaneous/DoubletRates_10xGenomicsUserGuideV3.csv",
                            sep = ",",header = T)
doublet.table$Rate <- str_remove(string=doublet.table$Rate,pattern = "~")
doublet.table$Rate <- str_remove(string=doublet.table$Rate,pattern = "%")
doublet.table$Rate <- as.numeric(doublet.table$Rate)/100
doublet.table$Loaded <- str_remove(string=doublet.table$Loaded,pattern = "~")
doublet.table$Loaded <- str_remove(string=doublet.table$Loaded,pattern = ",")
doublet.table$Loaded <- as.numeric(doublet.table$Loaded)
doublet.table$Recovered <- str_remove(string=doublet.table$Recovered,pattern = "~")
doublet.table$Recovered <- str_remove(string=doublet.table$Recovered,pattern = ",")
doublet.table$Recovered <- as.numeric(doublet.table$Recovered)

idx <- which.min(abs(doublet.table$Recovered-ncol(counts)))
rate <- doublet.table$Rate[idx]
print(rate)
nExp_poi <- round(rate*ncol(counts)) 
print(nExp_poi)

# Run DoubletFinder
rna <- doubletFinder_v3(rna, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
rna$Doublet.Call <- rna[[paste0("DF.classifications_0.25_",pK,"_",nExp_poi)]]
print(table(rna$Doublet.Call))

# Generate visuals showing predicted doublets
DimPlot(rna,group.by = "Doublet.Call")
ggsave(paste0("./QualityControl_scRNA/UMAP_predicted_doublets_",
              ID,".pdf"),
       width = 8,
       height = 8)
# Visualize QC metrics as scatter plots
plot1 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        group.by = "Doublet.Call")
plot1+geom_hline(yintercept = 20)+geom_vline(xintercept = 1000)
ggsave(paste0("./QualityControl_scRNA/FeatureScatter_nCount_RNA_percent.mt_predicted_doublets_",
              ID,".pdf"),
       width = 12,
       height = 8)
plot2 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        group.by = "Doublet.Call")
plot2+geom_hline(yintercept = 500)+geom_vline(xintercept = 1000)
ggsave(paste0("./QualityControl_scRNA/FeatureScatter_nCount_RNA_nFeature_RNA_predicted_doublets_",
              ID,".pdf"),
       width = 12,
       height = 8)
saveRDS(rna@meta.data,paste0("./QualityControl_scRNA/Doublet_Annotations_",ID,".rds"))
# Remove predicted doublets
rna <- subset(rna, 
              subset = Doublet.Call == "Singlet")

# normalizing the data
rna <- NormalizeData(object = rna, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)

# Find variable features
rna <- FindVariableFeatures(object = rna, selection.method = "vst", nfeatures = 2000,
                            loess.span = 0.3, clip.max = "auto",
                            num.bin = 20, binning.method = "equal_width", verbose = F)

# Scaling data
all.genes <- rownames(x = rna)
rna <- ScaleData(object = rna, features = all.genes, verbose=F)

# Run PCA to reduce dimensions
rna <- RunPCA(object = rna, features = VariableFeatures(object = rna), npcs = 50, verbose=F)

# Run UMAP
rna <- RunUMAP(rna,reduction = "pca", dims = 1:30,seed.use = 1,umap.method = "uwot",
               n.neighbors = 30L,
               metric = "cosine",
               learning.rate = 1,
               min.dist = 0.3)

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,".rds"))
