################################################################################
# Matt Regner
# Franco Lab
# Description: This script performs the following tasks  
# subset to patient samples
################################################################################
library(scater)
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(msigdbr)
library(fgsea)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(ArchR)
library(SingleR)
library(viridis)
library(parallel)
set.seed(1)
addArchRThreads(threads =16) 
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject(path =  "./All_Samples_scATAC")

# Subset ArchR project to patient samples
samples <- c("35A4AL",
             "35EE8L",
             "3821AL",
             "3B3E9L",
             "3C7D1L",
             "3D388L",
             "3FCDEL",
             "43E7BL",
             "43E7CL",
             "44F0AL",
             "45CB0L",
             "49758L",
             "49CFCL",
             "4AF75L",
             "4B146L",
             "4C2E5L")

idxSample <- BiocGenerics::which(proj$Sample %in% samples)
cellsSample <- proj$cellNames[idxSample]
  
subsetArchRProject(
  ArchRProj = proj,
  cells = cellsSample,
  outputDirectory = "Patient_Samples_scATAC",
  dropCells = FALSE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)
