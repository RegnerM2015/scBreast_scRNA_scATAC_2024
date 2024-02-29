################################################################################
# Matt Regner
# Franco Lab
# Description: This script performs the following tasks  
#         1) scATAC-seq quality control
#         2) Creates Arrow files for the ArchR Software
#         3) Create ArchR Project for all Arrow Files (All_Samples_scATAC)
#         4) Plot QC metrics
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
addArchRThreads(threads = 32) 
addArchRGenome("hg38")
h5disableFileLocking()

# Set up inputFiles
inputFiles <- list.files(path = "./cellranger-atac_outputs",
                         pattern = ".tsv.gz")
inputFiles <- str_remove(string = inputFiles, 
                         pattern = ".tbi")
inputFiles <- unique(inputFiles)

# Set up sampleNames/outputNames
sampleNames <- str_remove(string = inputFiles,
                          pattern = ".tsv.gz")
sampleNames <- str_remove(string = sampleNames,
                          pattern = "fragments_")
outputNames <- sampleNames

# Create Arrow and ArchR project
################################################################################
ArrowFiles <- createArrowFiles(
  inputFiles = paste0("./cellranger-atac_outputs/",inputFiles),
  sampleNames = sampleNames,
  outputNames = outputNames,
  validBarcodes = NULL,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  minTSS = 8,
  minFrags = 1000,
  maxFrags = 1e+05,
  QCDir = "QualityControl",
  nucLength = 147,
  promoterRegion = c(2000, 100),
  TSSParams = list(),
  excludeChr = c("chrM", "chrY"),
  nChunk = 5,
  bcTag = "qname",
  gsubExpression = NULL,
  bamFlag = NULL,
  offsetPlus = 4,
  offsetMinus = -5,
  addTileMat = TRUE,
  TileMatParams = list(),
  addGeneScoreMat = F,
  GeneScoreMatParams = list(),
  force = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  subThreading = TRUE,
  verbose = TRUE,
  cleanTmp = TRUE,
  logFile = createLogFile("createArrows"),
  filterFrags = NULL,
  filterTSS = NULL
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  useMatrix = "TileMatrix",
  k = 10,
  nTrials = 5,
  dimsToUse = 1:30,
  LSIMethod = 1,
  scaleDims = FALSE,
  corCutOff = 0.75,
  knnMethod = "UMAP",
  UMAPParams = list(n_neighbors = 30, min_dist = 0.3, metric = "cosine", verbose = FALSE),
  LSIParams = list(outlierQuantiles = NULL, filterBias = FALSE),
  threads = getArchRThreads(),
  force = TRUE,
  parallelParam = NULL,
  verbose = TRUE,
  logFile = createLogFile("addDoubletScores")
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "All_Samples_scATAC",
  copyArrows = T 
)

saveRDS(proj@cellColData,"./All_Samples_scATAC/predoublet_barcode_metadata.rds")

proj <- filterDoublets(proj,filterRatio = 1,cutEnrich = 1,cutScore = -Inf)

saveRDS(proj@cellColData,"./All_Samples_scATAC/postdoublet_barcode_metadata.rds")

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = getOutputDirectory(proj),
  overwrite = TRUE,
  load = FALSE,
  dropCells = T,
  logFile = createLogFile("saveArchRProject"),
  threads = getArchRThreads()
)

proj <- loadArchRProject(path =  "./All_Samples_scATAC")

# Plot QC measures
p1 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p2 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p4 <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", 
        ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)
