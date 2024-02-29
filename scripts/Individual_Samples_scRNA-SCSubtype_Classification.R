################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
################################################################################
library(Seurat)
library(tidyverse)
library(dplyr)
library(scales)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(stringi)
library(infercnv)
library(EnsDb.Hsapiens.v86)
library(parallel)
library(future)
options(future.globals.maxSize = 6000 * 1024^2)
plan("multicore")
set.seed(1)

args=(commandArgs(TRUE))
ID <- args[[1]]
print(ID)

rna <- readRDS(paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering-CellTypeAnnotations-inferCNV.rds"))

#Reading in the training datasets
Her2combinedNov2019 <- readRDS("./SCSubtype_Training_Data/Her2combinedNov2019.rds")
LumAcombinedNov2020 <- readRDS("./SCSubtype_Training_Data/LumAcombinedNov2020.rds")
LumBcombinedNov2019 <- readRDS("./SCSubtype_Training_Data/LumBcombinedNov2019.rds")
BasalTrainingdataset2020 <- readRDS("./SCSubtype_Training_Data/BasalTrainingdataset2020.rds")

#Merge the testing with the training datasets(RDS objects)
Tumors.combined <- merge(Her2combinedNov2019,y=c(LumAcombinedNov2020,LumBcombinedNov2019,BasalTrainingdataset2020,rna[,rna$normal_cell_call=="cancer"]))
Idents(object = Tumors.combined) <- "orig.ident"

#Performing normalizing and scaling the whole dataset
DefaultAssay(Tumors.combined) <- "RNA"
# Tumors.combined[["percent.mt"]] <- PercentageFeatureSet(object = Tumors.combined, pattern = "^MT.")
# plot1 <- FeatureScatter(object = Tumors.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(object = Tumors.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
Tumors.combined <- NormalizeData(object = Tumors.combined)
Tumors.combined <- FindVariableFeatures(object = Tumors.combined, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(x = Tumors.combined)
Tumors.combined <- ScaleData(object = Tumors.combined, features = all.genes)

#SCTyper Signature Calculations
sigdat <- read.csv("./miscellaneous/NatGen_Supplementary_table_S4.csv",sep=',',na.strings=c("NA","NaN", ""))
sigdat <- t(sigdat)
tocalc<-as.data.frame(Tumors.combined@assays$RNA@scale.data)
outdat <- matrix(0,nrow=nrow(sigdat),
                 ncol=ncol(tocalc),
                 dimnames=list(rownames(sigdat),
                               colnames(tocalc)))
for(i in 1:nrow(sigdat)){
  sigdat[i,!is.na(sigdat[i,])]->module
  row <- as.character(unlist(module))
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc) %in% row)
  temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
  outdat[i,]<-as.numeric(temp)
}
final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
final<-as.data.frame(final)
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
finalm<-as.matrix(final)
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
finalm.sweep.t <- finalm.sweep.t[rownames(finalm.sweep.t) %in% colnames(rna[,rna$normal_cell_call == "cancer"]),]
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtype <- Finalnames
finalm.sweep.t$barcode <- rownames(finalm.sweep.t)

# Merge SCSubtype predictions to Seurat object
rna@meta.data <- merge(rna@meta.data,finalm.sweep.t,by="barcode",all.x=TRUE)
rownames(rna@meta.data) <- rna$barcode

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds"))