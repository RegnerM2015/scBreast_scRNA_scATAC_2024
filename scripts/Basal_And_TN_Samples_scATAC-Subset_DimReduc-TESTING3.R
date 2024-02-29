################################################################################
# Matt Regner
# Franco Lab
# Description: This script performs the following tasks  
# Subset, DimReduc 
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

proj <- loadArchRProject(path =  "./Patient_Samples_scATAC")

# Subset ArchR project to patient samples
idxSample <- BiocGenerics::which(proj$predictedGroup %in% c("11",
                                                            "14",
                                                            "4",
                                                            "7"))

normals <- c("49758L","4AF75L","4B146L","49CFCL")
proj$From_Normal_Sample <- ifelse(proj$Sample %in% normals,"Yes","No")

proj$Keep.Cells <- ifelse(proj$predictedGroup %in% c("4","7") & proj$From_Normal_Sample == "No","Drop","Keep" )

idxSample.remove.1 <- BiocGenerics::which(proj$Keep.Cells == "Drop")

basals <- c("35A4AL","4C2E5L")
proj$From_Basal_Sample <- ifelse(proj$Sample %in% basals,"Yes","No")

proj$Keep.Cells <- ifelse(proj$predictedGroup %in% c("11","14") & proj$From_Basal_Sample == "No","Drop","Keep" )

idxSample.remove.2 <- BiocGenerics::which(proj$Keep.Cells == "Drop")

# Read in scRNA
rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
rna$barcode <- rownames(rna@meta.data)

meta <- as.data.frame(proj@cellColData)
all.equal(rownames(meta),proj$cellNames)
meta$cellNames <- rownames(meta)
all.equal(meta$cellNames,proj$cellNames)
length(intersect(meta$predictedCell,colnames(rna)))
length(intersect(meta$predictedCell,rownames(rna@meta.data)))
length(intersect(meta$predictedCell,rna$barcode))
meta$barcode <- meta$predictedCell
meta$id  <- 1:nrow(meta)

meta <- merge(meta,rna@meta.data,by="barcode")
meta <- meta[order(meta$id), ]

all.equal(meta$cellNames,proj$cellNames)

proj$Keep.Cells <- ifelse(meta$predictedGroup %in% c("11",
                                                     "14") & meta$SCSubtype %ni% c("Basal_SC"),"Drop","Keep" )

idxSample.remove.3 <- BiocGenerics::which(proj$Keep.Cells == "Drop")

idxSample <- idxSample[idxSample %ni% c(idxSample.remove.1,idxSample.remove.2,idxSample.remove.3)]
cellsSample <- proj$cellNames[idxSample]

df <- as.data.frame(proj@cellColData)
df$barcode <- rownames(df)
df <- df[df$barcode %in% cellsSample,]

sub <- df[df$predictedGroup %in% c("4","7"),]
table(sub$Sample)

sub <- df[df$predictedGroup %in% c("11","14"),]
table(sub$Sample)

subsetArchRProject(
  ArchRProj = proj,
  cells = cellsSample,
  outputDirectory = "Basal_TN_Samples_scATAC-TESTING3",
  dropCells = FALSE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

proj <- loadArchRProject(path =  "./Basal_TN_Samples_scATAC-TESTING3")

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start
                       = 10),
  firstSelection = "top",
  depthCol = "nFrags",
  varFeatures = 5000,
  dimsToUse = 1:30,
  LSIMethod = 2,
  scaleDims = TRUE,
  corCutOff = 0.75,
  binarize = TRUE,
  outlierQuantiles = c(0.02, 0.98),
  filterBias = TRUE,
  sampleCellsPre = 10000,
  projectCellsPre = FALSE,
  sampleCellsFinal = NULL,
  selectionMethod = "var",
  scaleTo = 10000,
  totalFeatures = 5e+05,
  filterQuantile = 0.995,
  excludeChr = c(),
  saveIterations = TRUE,
  UMAPParams = list(n_neighbors = 30, min_dist = 0.3, metric = "cosine", verbose =
                      FALSE),
  nPlot = 10000,
  outDir = getOutputDirectory(proj),
  threads = getArchRThreads(),
  seed = 2,
  verbose = TRUE,
  force = TRUE,
  logFile = createLogFile("addIterativeLSI")
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.3,
  metric = "cosine",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  sampleCells = NULL,
  outlierQuantile = 0.9,
  saveModel = TRUE,
  verbose = TRUE,
  seed = 1,
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-TSS.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "log10(nFrags)", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-logFrags.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

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
