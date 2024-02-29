################################################################################
# Matt Regner
# Franco Lab
# Description: This script performs the following tasks  
# DimReduc, and Gene Scoring
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
set.seed(3)
addArchRThreads(threads =16) 
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject(path =  "./Patient_Samples_scATAC")

proj <- addTileMatrix(
  input = proj,
  chromSizes = if (inherits(proj, "ArchRProject")) getChromSizes(proj) else NULL,
  blacklist = if (inherits(proj, "ArchRProject")) getBlacklist(proj) else NULL,
  tileSize = 500,
  binarize = TRUE,
  excludeChr = c("chrM", "chrY"),
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addTileMatrix")
)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start
                       = 10),
  firstSelection = "top",
  depthCol = "nFrags",
  varFeatures = 10000,
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
  seed = 3,
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

# Add Gene Score matrix
proj <- addGeneScoreMatrix(
  input = proj,
  genes = getGenes(proj),
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  matrixName = "GeneScoreMatrix",
  extendUpstream = c(1000, 1e+05),
  extendDownstream = c(1000, 1e+05),
  geneUpstream = 5000,
  geneDownstream = 0,
  useGeneBoundaries = TRUE,
  useTSS = FALSE,
  extendTSS = FALSE,
  tileSize = 500,
  ceiling = 4,
  geneScaleFactor = 5,
  scaleTo = 10000,
  excludeChr = c("chrY", "chrM"),
  blacklist = getBlacklist(proj),
  threads = getArchRThreads(),
  parallelParam = NULL,
  subThreading = TRUE,
  force = TRUE,
  logFile = createLogFile("addGeneScoreMatrix")
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
