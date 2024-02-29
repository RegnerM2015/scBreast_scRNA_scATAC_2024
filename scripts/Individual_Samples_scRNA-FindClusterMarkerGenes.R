################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: This script performs the following tasks  
#         1) Find DEGs with Seurat's FindAllMarkers()
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

args=(commandArgs(TRUE))
ID <- args[[1]]
print(ID)

rna <- readRDS(paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering.rds"))
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