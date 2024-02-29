################################################################################
# Matt Regner
# Franco Lab
# Description: This script performs the following tasks  
# Transfer labels from scRNA patient cohort to scATAC patient cohort
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
addArchRThreads(threads =8)
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject(path =  "./Patient_Samples_scATAC")

# Read in scRNA
rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

groupList <- SimpleList()
for (i in levels(factor(proj$Sample))){
  
  rna.sub <- rna[,rna$orig.ident == i]
  RNA.cells <- colnames(rna.sub)
  
  proj.meta <- as.data.frame(proj@cellColData)
  proj.meta.sub <- proj.meta[proj.meta$Sample == i,]
  ATAC.cells <- rownames(proj.meta.sub)
  
  groupList[[i]] <- SimpleList(
    ATAC = ATAC.cells,
    RNA = RNA.cells
  )
}

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = rna,
  groupATAC = NULL,
  groupRNA = "RNA_snn_res.0.4",
  groupList = groupList,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  embeddingATAC = NULL,
  embeddingRNA = NULL,
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  plotUMAP = TRUE,
  UMAPParams = list(n_neighbors = 30, min_dist = 0.3, metric = "cosine", verbose =
                      FALSE),
  nGenes = 2000,
  useImputation = FALSE,
  reduction = "cca",
  addToArrow = TRUE,
  scaleTo = 10000,
  genesUse = NULL,
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  transferParams = list(),
  threads = getArchRThreads(),
  verbose = TRUE,
  force = TRUE,
  logFile = createLogFile("addGeneIntegrationMatrix")
)


p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-TSS.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "log10(nFrags)", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-logFrags.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-predictedGroup.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedScore", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-predictedScore.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-Sample.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = getOutputDirectory(proj),
  overwrite = TRUE,
  load = FALSE,
  dropCells = T,
  logFile = createLogFile("saveArchRProject"),
  threads = getArchRThreads()
)