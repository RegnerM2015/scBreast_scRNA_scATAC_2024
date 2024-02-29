################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: The remove population with suspected cDNA degradation and low 
# mapping rate in reads
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
set.seed(321)

ID <- "3FCDEL"
print(ID)

rna <- readRDS(paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering.rds"))

Idents(rna) <- "MultiK_clusters"
idents <- unique(rna$MultiK_clusters)
idents <- idents[-c(1:3,8)]
rna <- subset(x=rna,idents=idents)

#normalizing the data
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

if(dir.exists("MultiK_scRNA_clustering")==TRUE){
  print("Directory MultiK_scRNA_clustering already exists!")
}else{
  dir.create("MultiK_scRNA_clustering")
}

# Run MultiK to estimate robust groups 
multik <- MultiK(rna, reps=100,seed = 321,resolution = seq(0.5,2,0.05))
saveRDS(multik,paste0("./MultiK_scRNA_clustering/multik_object_",ID,".rds"))
pdf(paste0("./MultiK_scRNA_clustering/InitialDiagnosticPlot_",ID,".pdf"),
    width = 16,
    height = 8)
DiagMultiKPlot(multik$k, multik$consensus)
dev.off()

# get optimal k
tog <- as.data.frame(table(multik$k)[table(multik$k) > 1])
colnames(tog)[1] <- "ks"
pacobj <- CalcPAC(x1=0.1, x2=0.9, xvec = tog$ks, ml = multik$consensus)
tog$rpac <- pacobj$rPAC
tog$one_minus_rpac  <- 1-tog$rpac
optK <- findOptK(tog)
print(optK)

if(length(optK) > 1){
  num <- as.numeric(optK)
  optK <- as.character(min(num))
}else{
  optK <- optK
}

# get clusters
clusters <- getClusters(rna, optK)

# Assign cluster resolution to seurat object
if(length(which(rownames(rna@meta.data)== rownames(clusters$clusters))) == nrow(rna@meta.data)){
  df <- as.data.frame(clusters$clusters)
  rna$MultiK_clusters <- as.character(df[,1])
}else{
  print("BARCODES DO NOT MATCH!")
}
Idents(object = rna) <- "MultiK_clusters"

if(dir.exists("DifferentialExpression_ClusterAnnotation_scRNA")==TRUE){
  print("Directory DifferentialExpression_ClusterAnnotation_scRNA already exists!")
}else{
  dir.create("DifferentialExpression_ClusterAnnotation_scRNA")
}

markers <- FindAllMarkers(rna, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25,
                          test.use = "wilcox")
markers <- markers[markers$p_val_adj <= 0.01,]
saveRDS(markers,paste0("./DifferentialExpression_ClusterAnnotation_scRNA/WilcoxDEGS_",ID,".rds"))

markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(rna, features = top10$gene) + NoLegend()
ggsave(paste0("./DifferentialExpression_ClusterAnnotation_scRNA/Top10DEGsHeatmp_",
              ID,".pdf"),
       width = 12,
       height = 8)

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering.rds"))

