findOptimaloverlapCutoff <- function(ArchRProj,
                                     reducedDims,
                                     useMatrix,
                                     dimsToUse,
                                     scaleDims,
                                     corCutOff,
                                     k,
                                     knnIteration,
                                     restrictKNN,
                                     predictionCutoff,
                                     seed,
                                     threads,
                                     verbose,
                                     groupBy,
                                     groupByScore,
                                     cellsToUse
                                     ){
  collect <- data.frame(overlapCutoff = 0,
                        groupBy = 0,
                        restrictKNN = 0,
                        NumberOfMetacells = 0,
                        NumberOfSingleCells = 0)
  df <- as.data.frame(ArchRProj@cellColData)
  df$cellNames <- rownames(df)
  df <- df[df$predictedScore >  predictionCutoff,]
  df <- df[df$cellNames %in% cellsToUse,]
  for ( i in seq(0,1,.05)){
    
    for ( j in levels(factor(df[[groupBy]]))){
      
      for ( l in levels(factor(df[[restrictKNN]][df[[groupBy]] == j]))){
        
        meta <- as.data.frame(ArchRProj@cellColData)
        meta$cellNames <- rownames(meta)
        meta <- meta[meta$cellNames %in% df$cellNames,]
        print(dim(meta))
        
        # Get cell names of specific cell group identity
        print(dim(meta))
        meta <- meta[meta[[groupBy]] == j & meta[[restrictKNN]] == l,]
        cellNames <- rownames(meta)
        
        print(dim(meta))
        print(l)
        print(j)
        print(length(cellNames))
        
        # Get Reduced Dims
        rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
        rD <- rD[cellNames, ,drop=FALSE]
        
        print(dim(rD))
        
        #Subsample 
        set.seed(seed)
        idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
        
        #KNN Matrix 
        knnObj.temp <- .computeKNN(data = rD, query = rD[idx,], k = k)
        
        #Determine Overlap 
        print(dim(knnObj.temp))
        keepKnn <- determineOverlapCpp(knnObj.temp, floor(i * k))
        
        #Keep Above Cutoff 
        knnObj.temp <- as.data.frame(knnObj.temp)
        knnObj.temp <- knnObj.temp[keepKnn==0,]
        
        collect <- rbind(collect,data.frame(overlapCutoff = i,
                                            groupBy = j,
                                            restrictKNN = l,
                                            NumberOfMetacells = nrow(knnObj.temp),
                                            NumberOfSingleCells = length(cellNames)))
        print(collect)
      }
      
    }
    print(paste0("Finished overlapCutoff: ",i))
  }  
  collect <- collect[-1,]
  collect.sub <- collect[collect$NumberOfMetacells >=50,]
  res <- sort(as.numeric(names(table(collect.sub$overlapCutoff)[table(collect.sub$overlapCutoff) == length(unique(df[[restrictKNN]]))])))[1]
  p <- ggplot(collect,aes(overlapCutoff,NumberOfMetacells,color=restrictKNN))+geom_point()+geom_line()
  return(list(res,collect,p))
}

getMetaCells <- function(ArchRProj,
                         reducedDims,
                         useMatrix,
                         dimsToUse,
                         scaleDimsL,
                         corCutOff,
                         k,
                         knnIteration,
                         overlapCutoff,
                         maxDist,
                         scaleTo,
                         log2Norm,
                         predictionCutoff,
                         seed,
                         threads,
                         verbose,
                         cellsToUse,
                         groupBy,
                         restrictKNN,
                         linkageMethod,
                         compareGroup1,
                         compareGroup2
                         ){
  
  # Get Arrow files
  ArrowFiles <- getArrowFiles(ArchRProj)
  
  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  
  #Gene Info
  geneSet <- .getFeatureDF(ArrowFiles, useMatrix, threads = threads)
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
  
  # Get metacells
  knnObj <- SimpleList()
  knnObj.groupBy <- SimpleList()
  knnObj.restrictKNN <- SimpleList()
  print("Finding metacell groups...")
  for ( i in levels(factor(ArchRProj@cellColData[[groupBy]][rownames(ArchRProj@cellColData) %in% cellsToUse]))){
    
    for ( j in levels(factor(ArchRProj@cellColData[[restrictKNN]][rownames(ArchRProj@cellColData) %in% cellsToUse &
                                                                  ArchRProj@cellColData[[groupBy]] == i]))){
      # Get cell names of specific cell group identity
      meta <- as.data.frame(ArchRProj@cellColData)
      meta$cellNames <- rownames(meta)
      meta <- meta[meta$cellNames %in% cellsToUse,]
      meta <- meta[meta[groupBy] == i & meta[restrictKNN] == j,]
      cellNames <- rownames(meta)
      
      # Get Reduced Dims
      rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
      rD <- rD[cellNames, ,drop=FALSE]
      
      #Subsample 
      set.seed(seed)
      idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
      
      #KNN Matrix 
      knnObj.temp <- .computeKNN(data = rD, query = rD[idx,], k = k)
      
      #Determine Overlap 
      keepKnn <- determineOverlapCpp(knnObj.temp, floor(overlapCutoff * k))
      
      #Keep Above Cutoff 
      knnObj.temp <- knnObj.temp[keepKnn==0,]
      
      #Convert To Names List 
      
      # Store cellNames
      knnObj.append <- lapply(seq_len(nrow(knnObj.temp)), function(x){
        rownames(rD)[knnObj.temp[x, ]]
      }) %>% SimpleList
      knnObj <- append(knnObj,knnObj.append)
      
      # Store groupBy labels
      knnObj.append <- lapply(seq_len(nrow(knnObj.temp)), function(x){
        meta[[groupBy]][knnObj.temp[x, ]]
      }) %>% SimpleList
      knnObj.groupBy <- append(knnObj.groupBy,knnObj.append)
      
      # Store restrictKNN labels
      knnObj.append <- lapply(seq_len(nrow(knnObj.temp)), function(x){
        meta[[restrictKNN]][knnObj.temp[x, ]]
      }) %>% SimpleList
      knnObj.restrictKNN <- append(knnObj.restrictKNN,knnObj.append)
    }
    
  }
  
  #Check Chromosomes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)
  
  #Features
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- seqnames(geneStart)
  peakDF$seqnames <- seqnames(peakSet)
  
  print("Making RNA metacell group matrix...")
  #Group Matrix RNA
  groupMatRNA <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
    useMatrix = useMatrix,
    threads = threads,
    verbose = FALSE
  )
  #groupMatRNA <- groupMatRNA[rowSums(groupMatRNA) != 0,]
  rawMatRNA <- groupMatRNA
  
  print("Making ATAC metacell group matrix...")
  #Group Matrix ATAC
  groupMatATAC <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = FALSE
  )
  #groupMatATAC <- groupMatATAC[rowSums(groupMatATAC) != 0,]
  rawMatATAC <- groupMatATAC
  
  groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo
  groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo
  
  if(log2Norm){
    groupMatRNA  <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)    
  }
  
  names(geneStart) <- NULL
  
  print("Making RNA SummarizedExperiment...")
  seRNA <- SummarizedExperiment(
    assays = SimpleList(RNA = groupMatRNA, RawRNA = rawMatRNA), 
    rowRanges = geneStart[as.numeric(sub('.', '', rownames(groupMatRNA)))]
  )
  metadata(seRNA)$KNNList <- knnObj
  
  names(peakSet) <- NULL
  
  print("Making ATAC SummarizedExperiment...")
  seATAC <- SummarizedExperiment(
    assays = SimpleList(ATAC = groupMatATAC, RawATAC = rawMatATAC), 
    rowRanges = peakSet[as.numeric(sub('.', '', rownames(groupMatATAC)))]
  )
  metadata(seATAC)$KNNList <- knnObj
  
  rm(groupMatRNA, groupMatATAC)
  gc()
  
  print("Storing metacell information...")
  # Find order of metacells
  store.majorities <- c("fill")
  store.proportions <- c(0)
  for(i in 1:length(knnObj.groupBy)){
    
    store.proportions <- append(store.proportions,max(table(knnObj.groupBy[[i]]))/sum(table(knnObj.groupBy[[i]])))
    store.majorities <- append(store.majorities,names(which.max(table(knnObj.groupBy[[i]]))))
    
  }
  store.majorities.groupBy <- store.majorities[-1]# Remove placeholder
  store.proportions.groupBy <- store.proportions[-1]# Remove placeholder
  
  print(summary(store.proportions.groupBy))
  
  # Find order of metacells
  store.majorities <- c("fill")
  store.proportions <- c(0)
  for(i in 1:length(knnObj.restrictKNN)){
    
    store.proportions <- append(store.proportions,max(table(knnObj.restrictKNN[[i]]))/sum(table(knnObj.restrictKNN[[i]])))
    store.majorities <- append(store.majorities,names(which.max(table(knnObj.restrictKNN[[i]]))))
    
  }
  store.majorities.restrictKNN <- store.majorities[-1]# Remove placeholder
  store.proportions.restrictKNN <- store.proportions[-1]# Remove placeholder
  
  print(summary(store.proportions.restrictKNN))
  
  return(list(seRNA,
              seATAC,
              store.majorities.groupBy,
              store.proportions.groupBy,
              store.majorities.restrictKNN,
              store.proportions.restrictKNN))
}
source("scripts/HiddenUtils-ArchR.R")
source("scripts/ValidationUtils-ArchR.R")
source("scripts/RcppExports-ArchR.R")
source("scripts/ArrowRead-ArchR.R")
source("scripts/ArrowUtils-ArchR.R")
source("scripts/LoggerUtils-ArchR.R")
source("scripts/runDiffLME_alt.R")

library(ArchR) 
library(parallel)
addArchRThreads(threads = 64) 
addArchRGenome("hg38")
proj <- loadArchRProject(path =  "./Basal_TN_Samples_scATAC-TESTING3")
################################################################################
# Remove samples with low # of cells
samplesToUse <- names(table(proj$predictedGroup)[table(proj$predictedGroup) >150])

# Make cellsToUse
idxSample <- BiocGenerics::which(proj$predictedGroup %in% c("1","2","3"))
cellsToUse <- proj$cellNames[idxSample]

# Make cancer normal status
proj$cancer_normal_status <- ifelse(proj$Sample %in% c("49758L", "49CFCL", "4AF75L", "4B146L"),"Normal","Cancer" )

# Make cancer normal status
proj$patient <- ifelse(proj$Sample %in% c("43E7BL","43E7CL"),"43E7BL_43E7CL",proj$Sample )

# Determine optimal overlapCutoff param
res <- findOptimaloverlapCutoff(ArchRProj = proj,
                                reducedDims = "IterativeLSI",
                                useMatrix = "GeneIntegrationMatrix",
                                dimsToUse = 1:30,
                                scaleDims = NULL,
                                corCutOff = 0.75,
                                k = 100,
                                knnIteration = 100,
                                restrictKNN = "patient",
                                predictionCutoff = 0.4,
                                seed = 1,
                                threads = max(floor(getArchRThreads()/2), 1),
                                verbose = TRUE,
                                groupBy = "cancer_normal_status",
                                groupByScore = "predictedScore",
                                cellsToUse = cellsToUse)
print(res[[1]])
meta_cell_data <- getMetaCells( ArchRProj = proj ,
                                reducedDims = "IterativeLSI",
                                useMatrix = "GeneIntegrationMatrix",
                                dimsToUse = 1:30,
                                scaleDims = NULL,
                                corCutOff = 0.75,
                                k = 100,
                                knnIteration = 100,
                                overlapCutoff = 0.8,
                                maxDist = 500000,
                                scaleTo = 10^4,
                                log2Norm = TRUE,
                                predictionCutoff = 0.5,
                                seed = 5,
                                threads = max(floor(getArchRThreads()/2), 1),
                                verbose = TRUE,
                                cellsToUse = cellsToUse,
                                groupBy = "cancer_normal_status",
                                restrictKNN = "patient",
                                linkageMethod = "regression",
                                compareGroup1 = "cancer",
                                compareGroup2 = "normal"
                                )
saveRDS(meta_cell_data,"./LME_Basal_out-SingFits_OLS/meta_cell_data.rds")

# order of meta_cell_data list
# seRNA,
# seATAC,
# store.majorities.groupBy,
# store.proportions.groupBy,
# store.majorities.restrictKNN,
# store.proportions.restrictKNN


#Overlaps
o <- DataFrame(
  findOverlaps(
    resize(meta_cell_data[[1]], 2 * 500000 + 1, "center"), 
    resize(rowRanges(meta_cell_data[[2]]), 1, "center"), 
    ignore.strand = TRUE
  )
)

#Get Distance from Fixed point A B 
colnames(o) <- c("idxRNA", "idxATAC")

# Add gene and peak names
o$geneName <- meta_cell_data[[1]]@rowRanges$name[o$idxRNA]
o$peakName <- paste0(meta_cell_data[[2]]@rowRanges@seqnames[o$idxATAC],":",meta_cell_data[[2]]@rowRanges@ranges[o$idxATAC])
o$P2G <- paste0(o$peakName,"|",o$geneName)
o <- o[,c(4,3)]

p2gs_to_test <- as.data.frame(o)
saveRDS(p2gs_to_test,"./LME_Basal_out-SingFits_OLS/p2gs_to_test-metacells.rds")

# Create RNA and ATAC matrices
ATAC <- assay(meta_cell_data[[2]])
rownames(ATAC) <- paste0(meta_cell_data[[2]]@rowRanges@seqnames,":",meta_cell_data[[2]]@rowRanges@ranges)
colnames(ATAC) <- paste0("metacell_",1:ncol(ATAC))

RNA <- assay(meta_cell_data[[1]])
rownames(RNA) <- meta_cell_data[[1]]@rowRanges$name
colnames(RNA) <- paste0("metacell_",1:ncol(RNA))

# Make metadata
meta <- data.frame(cell_type = meta_cell_data[[3]],
                   patient = meta_cell_data[[5]])
rownames(meta) <- paste0("metacell_",1:ncol(ATAC))# identical metacells between RNA and ATAC

# Fit P2G models
start <- Sys.time()
runDiffLME_alt(peak_mat = ATAC,
               gene_mat = RNA,
               p2gs_to_test = p2gs_to_test,
               meta_data = meta,
               covariates = NULL,
               random_effect = "patient",
               cell_type = "cell_type",
               cell_type_1 = "Cancer",
               cell_type_2 = "Normal",
               meta_cells = TRUE,
               cores = 96,
               out_path = "LME_Basal_out-SingFits_OLS",
               id = "metacells")
end <- Sys.time()

end-start

# # Confirmatory testing
# df <- data.frame(gene = RNA[grep("TPSAB1",rownames(RNA)),],
#                  peak = ATAC[grep("chr16:964638-965138",rownames(ATAC)),])
# df <- cbind(df,meta[,colnames(meta) %in% c("cell_type","patient")])
# 
# df[["cell_type"]] <- ifelse(df[["cell_type"]] == "Cancer",1,0)
# df$gene <- base::scale(df$gene)
# df$peak <- base::scale(df$peak)
# 
# library(lme4)
# library(lmerTest)
# mod <- lmerTest::lmer(gene~peak+cell_type+(1|patient),data=df,REML = FALSE)
# mod.int <- lmerTest::lmer(gene~peak+cell_type+peak:cell_type+(1|patient),data=df,REML = FALSE)
# 
# 
# peak_mat = ATAC
# gene_mat = RNA
# p2gs_to_test = p2gs_to_test
# meta_data = meta
# covariates = NULL
# random_effect = "patient"
# cell_type = "cell_type"
# cell_type_1 = "Cancer"
# cell_type_2 = "Normal"
# meta_cells = TRUE
# cores = 96
# out_path = "scLME_out_repro_check_v2"
# id = "metacells"
