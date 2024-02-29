################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: This script performs the following tasks  
#         1) scRNA-seq clustering with MultiK
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

rna <- readRDS(paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,".rds"))

if(dir.exists("MultiK_scRNA_clustering")==TRUE){
  print("Directory MultiK_scRNA_clustering already exists!")
}else{
  dir.create("MultiK_scRNA_clustering")
}

# Run MultiK to estimate robust groups 
multik <- MultiK(rna, reps=100,seed = 1,resolution = seq(0.5,2,0.05))
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

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering.rds"))
