################################################################################
# Matt Regner
# Franco Lab
# Description: This script performs the following tasks  
# Subset
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
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject(path =  "./All_Samples_scATAC")

# Subset ArchR project to patient samples
samples <- c("MCF7",
             "T47D",
             "HCC1143",
             "SUM149PT")

idxSample <- BiocGenerics::which(proj$Sample %in% samples)
cellsSample <- proj$cellNames[idxSample]
  
subsetArchRProject(
  ArchRProj = proj,
  cells = cellsSample,
  outputDirectory = "CellLine_Samples_scATAC-update",
  dropCells = FALSE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

proj <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")

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

# Read in scRNA
rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_CellLine_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = rna,
  groupATAC = NULL,
  groupRNA = "orig.ident",
  groupList = NULL,
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

print(all.equal(proj$Sample,proj$predictedGroup))

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

# Find path to Macs2
pathToMacs2 <- findMacs2()

# Add peakset according to predicted labels
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "predictedGroup",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 500,
  maxFragments = 25 * 10^6,
  minReplicates = 2,
  maxReplicates = 5,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),
  returnGroups = FALSE,
  parallelParam = NULL,
  force = FALSE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "predictedGroup",
  peakMethod = "Macs2",
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 25,
  excludeChr = c("chrM", "chrY"),
  pathToMacs2 = pathToMacs2,
  genomeSize = NULL,
  shift = -75,
  extsize = 150,
  method = "q",
  cutOff = 0.1,
  additionalParams = "--nomodel --nolambda",
  extendSummits = 250,
  promoterRegion = c(2000, 100),
  genomeAnnotation = getGenomeAnnotation(proj),
  geneAnnotation = getGeneAnnotation(proj),
  plot = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = FALSE,
  verbose = TRUE,
  logFile = createLogFile("addReproduciblePeakSet")
)

proj <- addPeakMatrix(
  ArchRProj = proj,
  ceiling = 4,
  binarize = FALSE,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addPeakMatrix")
)

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = getOutputDirectory(proj),
  overwrite = TRUE,
  load = FALSE,
  dropCells = T,
  logFile = createLogFile("saveArchRProject"),
  threads = getArchRThreads()
)