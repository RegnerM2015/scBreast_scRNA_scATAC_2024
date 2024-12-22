###################################################
# Matt Regner
# Franco Lab 
###################################################

library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(ArchR)
library(stringr)
library(stringi)
library(RColorBrewer)
library(dplyr)
library(parallel)
library(ArchR)
source("scripts/Matrix.utils.R")
library(tidyr)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(ComplexHeatmap)
library(msigdbr)
library(EnhancedVolcano)
source("scripts/plotBrowserTrack-modified-2_LoopTracks-PseudobulkExprBoxPlots-TwoBulkTracks-CustomHighlightDistance.R")
source("scripts/HiddenUtils-ArchR.R")
source("scripts/ValidationUtils-ArchR.R")
source("scripts/RcppExports-ArchR.R")
source("scripts/ArrowRead-ArchR.R")
source("scripts/ArrowUtils-ArchR.R")
source("scripts/LoggerUtils-ArchR.R")
source("scripts/IntegrativeAnalysis-ArchR.R")
source("scripts/runDiffLME_alt.R")

set.seed(1234)
addArchRThreads(threads = 32)
addArchRGenome("hg38")

# Make output folder
if(dir.exists("Basal_Cohort_Results-updates-Revised-updated")){
  print("Directory Basal_Cohort_Results-updates-Revised-updated already exists!")
}else{
  dir.create("Basal_Cohort_Results-updates-Revised-updated")
}

################################################################################
# Plot UMAPs
################################################################################

# Plot seurat UMAP
seurat <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Basal_TN_Subset-TESTING.rds")
seurat$cancer <- ifelse(seurat$orig.ident %in% c("49758L","49CFCL", "4AF75L", "4B146L"),"normal","cancer")

# Plot UMAP colored by cell type cluster
cols <- c("#828282","#323232","#D26262","#B51515")

df <- as.data.frame(seurat@reductions$umap@cell.embeddings)
all.equal(rownames(df),rownames(seurat@meta.data))
df$cellType <- as.character(seurat$RNA_snn_res.0.015)
df$cellType <- ifelse(df$cellType == "0","Normal_Basal",df$cellType)
df$cellType <- ifelse(df$cellType == "1","Normal_LP",df$cellType)
df$cellType <- ifelse(df$cellType == "2","Basal_BC_35A4AL",df$cellType)
df$cellType <- ifelse(df$cellType == "3","Basal_BC_4C2E5L",df$cellType)
df$cellType <- factor(df$cellType,levels=rev(levels(factor(df$cellType))))

p1 <- ggplot(df,aes(x =UMAP_1,y=UMAP_2,color = cellType))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  scale_color_manual(values=cols)+
  ggtitle(paste0("scRNA-seq by cell type cluster\nn=",nrow(df)," cells"))

# Plot UMAP by sample
cols <- c("#A6CEE3", "#1F78B4","#B4B4B4", "#828282", "#505050", "#323232" )

df <- as.data.frame(seurat@reductions$umap@cell.embeddings)
all.equal(rownames(df),rownames(seurat@meta.data))
df$sample <- factor(seurat$orig.ident,levels=c("35A4AL","4C2E5L",
                                               "49758L",
                                               "49CFCL",
                                               "4AF75L",
                                               "4B146L"))

p2 <- ggplot(df,aes(x =UMAP_1,y=UMAP_2,color = sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  scale_color_manual(values=cols)+
  ggtitle(paste0("scRNA-seq by patient sample \nn=",nrow(df)," cells"))

# Plot UMAP colored by cell type cluster
cols <- c("#828282","#323232","#D26262","#B51515")

proj <- loadArchRProject(path =  "./Basal_TN_Samples_scATAC-TESTING3")
df <- plotEmbedding(proj,colorBy = "cellColData",name = "predictedGroup",embedding = "UMAP")
df <- as.data.frame(df$data)
df$color <- sub(".*?-","",df$color)
df$cellType <- as.character(df$color)
df$cellType <- ifelse(df$cellType == "0","Normal_Basal",df$cellType)
df$cellType <- ifelse(df$cellType == "1","Normal_LP",df$cellType)
df$cellType <- ifelse(df$cellType == "2","Basal_BC_35A4AL",df$cellType)
df$cellType <- ifelse(df$cellType == "3","Basal_BC_4C2E5L",df$cellType)
df$cellType <- factor(df$cellType,levels=rev(levels(factor(df$cellType))))

p3 <- ggplot(df,aes(x =x,y=y,color = cellType))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  scale_color_manual(values=cols)+
  ggtitle(paste0("scATAC-seq by cell type cluster\nn=",nrow(df)," cells"))


# Plot UMAP by sample
cols <- c("#A6CEE3", "#1F78B4","#B4B4B4", "#828282", "#505050", "#323232" )

df <- plotEmbedding(proj,colorBy = "cellColData",name = "Sample",embedding = "UMAP")
df <- as.data.frame(df$data)
df$color <- sub(".*?-","",df$color)
df$sample <- as.character(df$color)
df$sample <- factor(df$sample,levels=c("35A4AL","4C2E5L",
                                       "49758L",
                                       "49CFCL",
                                       "4AF75L",
                                       "4B146L"))

p4 <- ggplot(df,aes(x =x,y=y,color = sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  scale_color_manual(values=cols)+
  ggtitle(paste0("scATAC-seq by patient sample \nn=",nrow(df)," cells"))

p1 + p3 

ggsave(paste0("./Basal_Cohort_Results-updates-Revised-updated/UMAPs_",
              "Basal_Cell_Type_Cluster",".pdf"),
       width = 14,
       height = 6 )

p2 + p4 

ggsave(paste0("./Basal_Cohort_Results-updates-Revised-updated/UMAPs_",
              "Basal_Patient_Sample",".pdf"),
       width = 14,
       height = 6 )

################################################################################
# Differential gene expression
# Using code adapted from https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
################################################################################
# Plot seurat UMAP
seurat <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Basal_TN_Subset-TESTING.rds")
seurat$cancer <- ifelse(seurat$orig.ident %in% c("49758L","49CFCL", "4AF75L", "4B146L"),"normal","cancer")
seurat$patient <- ifelse(seurat$orig.ident %in% c("43E7CL","43E7BL"),"43E7BL",seurat$orig.ident)
seurat$full <- "full"
seurat <- seurat[,seurat$RNA_snn_res.0.015 != "0"]
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts

metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

cluster_names <- levels(factor(colData(sce)$full))
cluster_names

sample_names <- levels(factor(colData(sce)$patient))
sample_names

groups <- colData(sce)[, c("full", "patient")]
head(groups)


# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)),
                                groupings = groups, fun = "sum")

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
# aggr_counts[1:6, 1:6]

aggr_counts <- t(aggr_counts)
# aggr_counts[1:6, 1:6]


## Exploring structure of function output (list)
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)


# Using which() to look up tstrsplit() output
cell_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "full")
cell_idx

colnames(aggr_counts)[cell_idx]
aggr_counts[1:10, cell_idx]

cluster_names


# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)


head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::select(cancer, full, patient)

dim(metadata)
head(metadata)

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]

dim(metadata)
head(metadata)

rownames(metadata) <- metadata$patient
head(metadata)

t <- table(colData(sce)$patient,
           colData(sce)$full)
# t[1:6, 1:6]


## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata,
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)

# Select cell type of interest
cluster_names

# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))

idx <- which(names(counts_ls) == "full")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]


# Check contents of extracted objects
# cluster_counts[1:6, 1:6]
head(cluster_metadata)

# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))

# Read in clincal info

library(readxl)
# Read in clinical data
clinical <- read_excel("miscellaneous/Franco-Perou-SingleCellBreastCancerDataset-July2022_MR-MJR_KW.xlsx")
clinical$Sample <- clinical$Sample_ID
clinical$Patient_Age_at_Surgery <- as.numeric(clinical$Patient_Age_at_Surgery)
clinical$BMI <- as.numeric(clinical$BMI)
colnames(clinical)[16] <- "Menopause_Status"
clinical$Menopause_Status <- gsub("\\ .*","",clinical$Menopause_Status)
clinical$sample_id <- ifelse(clinical$Sample %in% c("43E7CL","43E7BL"),"43E7BL_43E7CL",clinical$Sample)
clinical <- clinical[,-c(grep("Sample_ID",colnames(clinical)),grep("Sample",colnames(clinical)))]
clinical <- clinical[!duplicated(clinical),]

cluster_metadata <- merge(cluster_metadata,clinical,by="sample_id")

# Create DESeq2 object
cluster_metadata$group_id <- ifelse(cluster_metadata$sample_id %in% c("49758L","49CFCL", "4AF75L", "4B146L"),"normal","cancer")
cluster_metadata$group_id <- factor(cluster_metadata$group_id,levels=c("cancer","normal"))

dds <- DESeqDataSetFromMatrix(cluster_counts[rowSums(cluster_counts) >0,],
                              colData = cluster_metadata,
                              design = ~ group_id)

# Transform counts for data visualization
rlog <- rlog(dds, blind=TRUE)

topProp <- 0.10
n_features = round(topProp*length(rlog))

# Plot PCA
p1 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "sample_id")+  ggtitle(paste0("PCA colored by sample\nTop 10% variable genes: ",n_features))
p2 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "group_id")+  ggtitle(paste0("PCA colored by group\nTop 10% variable genes: ",n_features))
p3 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "cell_count")+  ggtitle(paste0("PCA colored by cell count\nTop 10% variable genes: ",n_features))
p4 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Primary_Diagnosis")+  ggtitle(paste0("PCA colored by primary diagnosis\nTop 10% variable genes: ",n_features))
p5 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Patient_Age_at_Surgery")+  ggtitle(paste0("PCA colored by age\nTop 10% variable genes: ",n_features))
p6 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Patient_Race")+  ggtitle(paste0("PCA colored by race\nTop 10% variable genes: ",n_features))
p7 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "BMI")+  ggtitle(paste0("PCA colored by BMI\nTop 10% variable genes: ",n_features))
p8 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "T_value")+  ggtitle(paste0("PCA colored by T value\nTop 10% variable genes: ",n_features))
p9 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "N_value")+  ggtitle(paste0("PCA colored by N value\nTop 10% variable genes: ",n_features))
p10 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Grade")+  ggtitle(paste0("PCA colored by Grade\nTop 10% variable genes: ",n_features))
p11 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "ER")+  ggtitle(paste0("PCA colored by ER\nTop 10% variable genes: ",n_features))
p12 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "PR")+  ggtitle(paste0("PCA colored by PR\nTop 10% variable genes: ",n_features))
p13 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Her2_IHC")+  ggtitle(paste0("PCA colored by Her2\nTop 10% variable genes: ",n_features))
p14 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "BRCA")+  ggtitle(paste0("PCA colored by BRCA\nTop 10% variable genes: ",n_features))
p15 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Menopause_Status")+  ggtitle(paste0("PCA colored by Menopause status\nTop 10% variable genes: ",n_features))


plot_list <- list(p1,
                  p2,
                  p3,
                  p4,
                  p5,
                  p6,
                  p7,
                  p8,
                  p9,
                  p10,
                  p11,
                  p12,
                  p13,
                  p14,
                  p15)

ggsave(
  filename = paste0("./Basal_Cohort_Results-updates-Revised-updated/",
                    "Pseudobulk_RNA_PCA",".pdf"),
  plot = marrangeGrob(plot_list, nrow=1, ncol=1),
  width = 10, height = 10
)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

plotDispEsts(dds)

# Write out dds
saveRDS(dds,"./Basal_Cohort_Results-updates-Revised-updated/pseudobulk_expression_dds.rds")

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds,
               alpha = 0.05,
               contrast = c("group_id","cancer","normal"))

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res <- lfcShrink(dds,
                 res=res,
                 type="normal",
                 contrast = c("group_id","cancer","normal"))

pdf("Basal_Cohort_Results-updates-Revised-updated/EnhancedVolcano_DEGs.pdf",width=10,height=10)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Basal-like cancer versus normal LP',
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 6.0)
dev.off()


# Turn the DESeq2 results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl

padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

# Set thresholds
log2fc_cutoff <- 0.58

# Count significantly up/down genes above threshold
n_sig_up <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff) %>%
  nrow()
n_sig_dn <- dplyr::filter(sig_res, log2FoldChange <= -log2fc_cutoff) %>%
  nrow()

# Write out DEGs
degs <- dplyr::filter(sig_res, abs(log2FoldChange) >= log2fc_cutoff)
saveRDS(degs,paste0("./Basal_Cohort_Results-updates-Revised-updated/DEGs_padj_",padj_cutoff,"_log2FC_",log2fc_cutoff,".rds"))

## Extract normalized counts from dds object
normalized_counts <- counts(dds, normalized = TRUE)
saveRDS(normalized_counts,"./Basal_Cohort_Results-updates-Revised-updated/pseudobulk_expression_normalized_counts.rds")

## Extract top 20 DEG from resLFC (make sure to order by padj)
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  dplyr::pull(gene) %>%
  head(n = 4)

## Extract matching normalized count values from matrix
top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
top20_sig_counts

## Convert wide matrix to long data frame for ggplot2
top20_sig_df <- data.frame(top20_sig_counts)
top20_sig_df$gene <- rownames(top20_sig_counts)

top20_sig_df <- melt(setDT(top20_sig_df),
                     id.vars = c("gene"),
                     variable.name = "cluster_sample_id") %>%
  data.frame()

## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
top20_sig_df

## Join counts data frame with metadata
top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                           by = "cluster_sample_id")
top20_sig_df

## Generate plot
ggplot(top20_sig_df, aes(y = value, x = group_id, col = group_id)) +
  geom_jitter(height = 0, width = 0.15) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("condition") +
  ggtitle("Top 4 Significant DE Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  facet_wrap(~ gene)

ggsave(paste0("./Basal_Cohort_Results-updates-Revised-updated/DEGs_",
              "Basal",".pdf"),
       width = 5,
       height = 5 )


# Get upregulated genes in cancer
up.genes <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff)
saveRDS(up.genes,"./Basal_Cohort_Results-updates-Revised-updated/up_genes.rds")

library(fgsea)
library(msigdbr)
gset = msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::distinct(gs_name, gene_symbol)
gset <- split(gset$gene_symbol,gset$gs_name)

sig_res <- dplyr::filter(res_tbl, padj < 0.05) %>%
  dplyr::arrange(log2FoldChange)

ranks <- as.numeric(sig_res$log2FoldChange)
names(ranks) <- sig_res$gene
head(ranks)

fgseaRes <- fgsea::fgsea(gset, ranks)
fgseaRes.sig <- fgseaRes[fgseaRes$padj < 0.01,]
fgseaRes.sig <- dplyr::arrange(fgseaRes.sig,desc(NES))
head(fgseaRes.sig)

pdf("Basal_Cohort_Results-updates-Revised-updated/GSEA_Basal_DEGs.pdf",width=10,height=4)
plotGseaTable(gset[fgseaRes.sig$pathway], ranks, fgseaRes.sig,
              gseaParam=0.5)
dev.off()

# Plot heatmap of degs using rlog values

resSig = res[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 0.58), ]
saveRDS(resSig,"./Basal_Cohort_Results-updates-Revised-updated/resSig_genes.rds")
allSig_genes <- rownames(resSig)

# topVarGenes <- head(order(rowVars(assay(rlog)[allSig_genes,]), decreasing = TRUE), 100)

set.seed(1234)
mat <- assay(rlog)
colnames(mat) <- gsub(".*_","",colnames(mat))

pdf("Basal_Cohort_Results-updates-Revised-updated/Heatmap_of_DEGs_Basal.pdf",width=6,height=4)
ComplexHeatmap::Heatmap(t(scale(t(mat[allSig_genes,]))),cluster_rows = T,cluster_columns = TRUE)
dev.off()

################################################################################
# Differential peak accessibility
# Using code adapted from https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
################################################################################

proj <- loadArchRProject(path =  "./Basal_TN_Samples_scATAC-TESTING3")
################################################################################
# Remove samples with low # of cells
samplesToUse <- names(table(proj$predictedGroup)[table(proj$predictedGroup) >150])

# Make cellsToUse
idxSample <- BiocGenerics::which(proj$predictedGroup %in% c("1","2","3"))
cellsToUse <- proj$cellNames[idxSample]
proj <- proj[cellsToUse, ]

metadata <- as.data.frame(proj@cellColData)
all.equal(rownames(metadata),proj$cellNames)
metadata$cellNames <- proj$cellNames
metadata$patient <- metadata$Sample
metadata$cancer <- ifelse(metadata$Sample %in% c("49758L","49CFCL", "4AF75L", "4B146L"),"normal","cancer")

peaks <- getMatrixFromProject(ArchRProj = proj,useMatrix = "PeakMatrix",binarize = T)
peakData <- peaks@assays@data$PeakMatrix
rownames(peakData) <- paste0(peaks@rowRanges@seqnames,":",peaks@rowRanges@ranges)

all.equal(colnames(peakData),metadata$cellNames)
all.equal(colnames(peakData),rownames(metadata))

metadata <- metadata[order(match(rownames(metadata),colnames(peakData))),]

all.equal(colnames(peakData),metadata$cellNames)
all.equal(colnames(peakData),rownames(metadata))

metadata$full <- "full"

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = peakData),
                            colData = metadata)

cluster_names <- levels(factor(colData(sce)$full))
cluster_names

sample_names <- levels(factor(colData(sce)$patient))
sample_names

groups <- colData(sce)[, c("full", "patient")]
head(groups)


# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)),
                                groupings = groups, fun = "sum")

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
# aggr_counts[1:6, 1:6]

aggr_counts <- t(aggr_counts)
# aggr_counts[1:6, 1:6]


## Exploring structure of function output (list)
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)


# Using which() to look up tstrsplit() output
cell_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "full")
cell_idx

colnames(aggr_counts)[cell_idx]
aggr_counts[1:10, cell_idx]

cluster_names


# Loop over all cell types to extract corresponding counts, and store information in a list

## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)


head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::select(cancer, full, patient)

dim(metadata)
head(metadata)

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]

dim(metadata)
head(metadata)

rownames(metadata) <- metadata$patient
head(metadata)

t <- table(colData(sce)$patient,
           colData(sce)$full)
# t[1:6, 1:6]


## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata,
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)

# Select cell type of interest
cluster_names

# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))

idx <- which(names(counts_ls) == "full")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]


# Check contents of extracted objects
# cluster_counts[1:6, 1:6]
head(cluster_metadata)

# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))

# Read in clincal info

library(readxl)
# Read in clinical data
clinical <- read_excel("miscellaneous/Franco-Perou-SingleCellBreastCancerDataset-July2022_MR-MJR_KW.xlsx")
clinical$Sample <- clinical$Sample_ID
clinical$Patient_Age_at_Surgery <- as.numeric(clinical$Patient_Age_at_Surgery)
clinical$BMI <- as.numeric(clinical$BMI)
colnames(clinical)[16] <- "Menopause_Status"
clinical$Menopause_Status <- gsub("\\ .*","",clinical$Menopause_Status)
clinical$sample_id <- ifelse(clinical$Sample %in% c("43E7CL","43E7BL"),"43E7BL_43E7CL",clinical$Sample)
clinical <- clinical[,-c(grep("Sample_ID",colnames(clinical)),grep("Sample",colnames(clinical)))]
clinical <- clinical[!duplicated(clinical),]

cluster_metadata <- merge(cluster_metadata,clinical,by="sample_id")

# Create DESeq2 object
cluster_metadata$group_id <- ifelse(cluster_metadata$sample_id %in% c("49758L","49CFCL", "4AF75L", "4B146L"),"normal","cancer")
cluster_metadata$group_id <- factor(cluster_metadata$group_id,levels=c("cancer","normal"))

dds <- DESeqDataSetFromMatrix(cluster_counts[rowSums(cluster_counts) >0,],
                              colData = cluster_metadata,
                              design = ~ group_id)

# Transform counts for data visualization
rlog <- rlog(dds, blind=TRUE)

topProp <- 0.10
n_features = round(topProp*length(rlog))

# Plot PCA
p1 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "sample_id")+  ggtitle(paste0("PCA colored by sample\nTop 10% variable genes: ",n_features))
p2 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "group_id")+  ggtitle(paste0("PCA colored by group\nTop 10% variable genes: ",n_features))
p3 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "cell_count")+  ggtitle(paste0("PCA colored by cell count\nTop 10% variable genes: ",n_features))
p4 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Primary_Diagnosis")+  ggtitle(paste0("PCA colored by primary diagnosis\nTop 10% variable genes: ",n_features))
p5 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Patient_Age_at_Surgery")+  ggtitle(paste0("PCA colored by age\nTop 10% variable genes: ",n_features))
p6 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Patient_Race")+  ggtitle(paste0("PCA colored by race\nTop 10% variable genes: ",n_features))
p7 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "BMI")+  ggtitle(paste0("PCA colored by BMI\nTop 10% variable genes: ",n_features))
p8 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "T_value")+  ggtitle(paste0("PCA colored by T value\nTop 10% variable genes: ",n_features))
p9 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "N_value")+  ggtitle(paste0("PCA colored by N value\nTop 10% variable genes: ",n_features))
p10 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Grade")+  ggtitle(paste0("PCA colored by Grade\nTop 10% variable genes: ",n_features))
p11 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "ER")+  ggtitle(paste0("PCA colored by ER\nTop 10% variable genes: ",n_features))
p12 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "PR")+  ggtitle(paste0("PCA colored by PR\nTop 10% variable genes: ",n_features))
p13 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Her2_IHC")+  ggtitle(paste0("PCA colored by Her2\nTop 10% variable genes: ",n_features))
p14 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "BRCA")+  ggtitle(paste0("PCA colored by BRCA\nTop 10% variable genes: ",n_features))
p15 <- DESeq2::plotPCA(rlog, ntop = n_features, intgroup = "Menopause_Status")+  ggtitle(paste0("PCA colored by Menopause status\nTop 10% variable genes: ",n_features))

plot_list <- list(p1,
                  p2,
                  p3,
                  p4,
                  p5,
                  p6,
                  p7,
                  p8,
                  p9,
                  p10,
                  p11,
                  p12,
                  p13,
                  p14,
                  p15)

ggsave(
  filename = paste0("./Basal_Cohort_Results-updates-Revised-updated/",
                    "Pseudobulk_ATAC_PCA",".pdf"),
  plot = marrangeGrob(plot_list, nrow=1, ncol=1),
  width = 10, height = 10
)

# rlog_mat <- assay(rlog)
# rlog_cor <- cor(rlog_mat)
#
# # Plot heatmap
# pheatmap(rlog_cor, annotation = cluster_metadata[, c("group_id"), drop=F])

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

plotDispEsts(dds)

# Write out dds
saveRDS(dds,"./Basal_Cohort_Results-updates-Revised-updated/pseudobulk_accessibility_dds.rds")

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds,
               alpha = 0.05,
               contrast = c("group_id","cancer","normal"))

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res <- lfcShrink(dds,
                 res=res,
                 type="normal",
                 contrast = c("group_id","cancer","normal"))


pdf("Basal_Cohort_Results-updates-Revised-updated/EnhancedVolcano_DAPs.pdf",width=10,height=10)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Basal-like cancer versus normal LP',
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 6.0)
dev.off()

# Turn the DESeq2 results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl

padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

# Set thresholds
log2fc_cutoff <- 0.58

# Count significantly up/down genes above threshold
n_sig_up <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff) %>%
  nrow()
n_sig_dn <- dplyr::filter(sig_res, log2FoldChange <= -log2fc_cutoff) %>%
  nrow()

# Write out DAPs
daps <- dplyr::filter(sig_res, abs(log2FoldChange) >= log2fc_cutoff)
saveRDS(daps,paste0("./Basal_Cohort_Results-updates-Revised-updated/DAPs_padj_",padj_cutoff,"_log2FC_",log2fc_cutoff,".rds"))

## Extract normalized counts from dds object
normalized_counts <- counts(dds, normalized = TRUE)
saveRDS(normalized_counts,"./Basal_Cohort_Results-updates-Revised-updated/pseudobulk_accessibility_normalized_counts.rds")

## Extract top 20 DEG from resLFC (make sure to order by padj)
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  dplyr::pull(gene) %>%
  head(n = 4)

## Extract matching normalized count values from matrix
top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
top20_sig_counts

## Convert wide matrix to long data frame for ggplot2
top20_sig_df <- data.frame(top20_sig_counts)
top20_sig_df$gene <- rownames(top20_sig_counts)

top20_sig_df <- melt(setDT(top20_sig_df),
                     id.vars = c("gene"),
                     variable.name = "cluster_sample_id") %>%
  data.frame()

## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
top20_sig_df

## Join counts data frame with metadata
top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                           by = "cluster_sample_id")
top20_sig_df

## Generate plot
ggplot(top20_sig_df, aes(y = value, x = group_id, col = group_id)) +
  geom_jitter(height = 0, width = 0.15) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized accessibility level") +
  xlab("condition") +
  ggtitle("Top 4 Significant DE Peaks") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  facet_wrap(~ gene)

ggsave(paste0("./Basal_Cohort_Results-updates-Revised-updated/DAPs_",
              "Basal",".pdf"),
       width = 5,
       height = 5 )



# Get upregulated genes in cancer
up.peaks <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff)
saveRDS(up.peaks,"./Basal_Cohort_Results-updates-Revised-updated/up_peaks.rds")

# Plot heatmap of degs using rlog values

resSig = res[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 0.58), ]
saveRDS(resSig,"./Basal_Cohort_Results-updates-Revised-updated/resSig_peaks.rds")
allSig_genes <- rownames(resSig)

# topVarGenes <- head(order(rowVars(assay(rlog)[allSig_genes,]), decreasing = TRUE), 100)

set.seed(1234)
mat <- assay(rlog)
colnames(mat) <- gsub(".*_","",colnames(mat))

pdf("Basal_Cohort_Results-updates-Revised-updated/Heatmap_of_DAPs_Basal.pdf",width=6,height=4)
ComplexHeatmap::Heatmap(t(scale(t(mat[allSig_genes,]))),cluster_rows = T,cluster_columns = TRUE)
dev.off()


################################################################################
# Read and process peak-to-gene associations
################################################################################

up.genes <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/up_genes.rds")
up.peaks <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/up_peaks.rds")

# Find P2Gs 
peakSet <- getPeakSet(proj)
names(peakSet) <- NULL
peakSet <- as.data.frame(peakSet)
peakSet$peakName <- paste0(peakSet$seqnames,":",peakSet$start,"-",peakSet$end)

# Unbalanced
meta <- fread("LME_Basal_out-SingFits_OLS/interactionLMM_univariateLMM_and_OLS_results-metacells.csv")
meta <- as.data.frame(meta)
colnames(meta) <- meta[2,]
meta <- meta[-2,]
names <- grep("Name",colnames(meta))
singFits <- grep("singular_fit",colnames(meta))
names <- unique(c(names,singFits))
meta[,-names] <- sapply(meta[,-names], as.numeric)
meta <- meta[complete.cases(meta),]
meta$P2G <- paste0(meta$LME_cell_type_1_model_peakName,"|",meta$LME_cell_type_1_model_geneName)

# meta.cancer <- meta
# meta.cancer$FDR <- p.adjust(meta.cancer$LME_cell_type_1_model_peak_pval,method = "fdr")
# meta.cancer <- meta.cancer[meta.cancer$FDR < 1e-04,]
# meta.cancer$peakName <- meta.cancer$LME_cell_type_1_model_peakName
# meta.cancer <- merge(meta.cancer,peakSet,by="peakName")
# meta.cancer <- meta.cancer[meta.cancer$peakType %in% c("Intronic","Distal"),]

# meta.normal <- meta
# meta.normal$FDR <- p.adjust(meta.normal$LME_cell_type_2_model_peak_pval,method = "fdr")
# meta.normal <- meta.normal[meta.normal$FDR < 1e-04,]
# meta.normal$peakName <- meta.normal$LME_cell_type_2_model_peakName
# meta.normal <- merge(meta.normal,peakSet,by="peakName")
# meta.normal <- meta.normal[meta.normal$peakType %in% c("Intronic","Distal"),]

meta.int <- meta
meta.int$peakName <- meta.int$LME_cell_type_2_model_peakName
meta.int <- merge(meta.int,peakSet,by="peakName")
# meta.int <- meta.int[meta.int$peakType %in% c("Intronic","Distal"),]

meta.int$Cancer_FDR <- p.adjust(meta.int$LME_cell_type_1_model_peak_pval,method = "fdr")
meta.int$Normal_FDR <- p.adjust(meta.int$LME_cell_type_2_model_peak_pval,method = "fdr")

p2g <- meta.int

classes = rep(0, nrow(p2g))

pvA = (p2g$Cancer_FDR < 1e-04)
pvB = (p2g$Normal_FDR < 1e-04)
cAup = (p2g$LME_cell_type_1_model_peak_effect_size > 0)
cBup = (p2g$LME_cell_type_2_model_peak_effect_size > 0)

#UpUp
classes[which(pvA & pvB & cAup & cBup)] = 1
#UpNon
classes[which(pvA & !pvB & cAup)] = 2
#UpDown
classes[which(pvA & pvB & cAup & !cBup)] = 3
#NonUp
classes[which(!pvA & pvB & cBup)] = 4
#NonNon
classes[which(!pvA & !pvB)] = 5
#NonDown
classes[which(!pvA & pvB & !cBup)] = 6
#DownUp
classes[which(pvA & pvB & !cAup & cBup)] = 7
#DownNon
classes[which(pvA & !pvB & !cAup)] = 8
#DownDown
classes[which(pvA & pvB & !cAup & !cBup)] = 9


classes = factor(classes, levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                 labels = c("NonSig", "+/+", "+/0", "+/-",
                            "0/+", "0/0", "0/-", "-/+", "-/0", "-/-"))

p2g$DiffClass <- classes

p2g$Both <- ifelse(p2g$DiffClass %in% c("+/+",
                                        "+/-",
                                        "-/+",
                                        "-/-"),"Shared",p2g$DiffClass)
p2g$Both <- ifelse(p2g$DiffClass %in% c("+/0",
                                        "-/0"),"Cancer",p2g$Both)
p2g$Both <- ifelse(p2g$DiffClass %in% c("0/+",
                                        "0/-"),"Normal",p2g$Both)
p2g$Both <- ifelse(p2g$DiffClass %in% c("0/0"),"Neither",p2g$Both)

table(p2g$DiffClass)
table(p2g$Both)

saveRDS(p2g,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_univariate-prefilter.rds")

p2g <- p2g[p2g$DiffClass != "0/0",]
p2g$Both <- factor(p2g$Both,levels=c("Normal","Cancer","Shared"))
p2g$DiffClass <- factor(p2g$DiffClass,levels=c("0/+","0/-",
                                               "+/0","-/0",
                                               "+/+","-/-",
                                               "+/-","-/+"
                                               ))

saveRDS(p2g,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_univariate-postfilter.rds")

# Annotate ENCODE and TCGA regions

# Downloaded TCGA BRCA peaks from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG 
brca <- read.delim("./miscellaneous/BRCA_peakCalls.txt",sep = "\t")
brca.gr <- GRanges(paste0(brca$seqnames,":",brca$start,"-",brca$end))

# Downloaded ENCODE cCREs from https://screen.encodeproject.org 
encode <- read.delim("./miscellaneous/GRCh38-cCREs.bed",header = F)
encode.gr <- GRanges(paste0(encode$V1,":",encode$V2,"-",encode$V3))

p2g.gr <- GRanges(unique(p2g$peakName))

# Overlap with references
overlaps <- subsetByOverlaps(x = p2g.gr,
                             ranges = brca.gr)

brcaPct <- length(overlaps) / length(p2g.gr)

p2g$TCGA_BRCA_Peak_Overlap <- ifelse(p2g$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")

overlaps <- subsetByOverlaps(x = p2g.gr,
                             ranges = encode.gr)

encodePct <- length(overlaps) / length(p2g.gr)

p2g$ENCODE_Peak_Overlap <- ifelse(p2g$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")

# peakType proportion barchart
df <- p2g %>% dplyr::group_by_at("Both") %>% dplyr::count(peakType) %>% mutate(pct= prop.table(n) * 100)
colnames(df) <- c("conditionSignif","peakType","Cells","Pct")
df$peakType <- factor(df$peakType,levels=c("Intronic","Distal","Promoter","Exonic"))
df %>%
  ggplot(aes(fill=peakType, y=Pct, x= conditionSignif, label = Pct))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=rev(RColorBrewer::brewer.pal(5,"Set1")[-1]))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nNormal: ",table(p2g$Both)[1],"  |  Cancer: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

ggsave("./Basal_Cohort_Results-updates-Revised-updated/peakType_proportionBarchart-univariate.pdf",width = 4, height = 4)

peakTypePropUni <- df %>%
  ggplot(aes(fill=peakType, y=Pct, x= conditionSignif, label = Pct))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=rev(RColorBrewer::brewer.pal(5,"Set1")[-1]))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nNormal: ",table(p2g$Both)[1],"  |  Cancer: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

# ENCDOE proportion barchart
df <- p2g %>% dplyr::group_by_at("Both") %>% dplyr::count(ENCODE_Peak_Overlap) %>% mutate(pct= prop.table(n) * 100)
colnames(df) <- c("conditionSignif","ENCODE","Cells","Pct")
df$peakType <- factor(df$ENCODE)
df %>%
  ggplot(aes(fill=ENCODE, y=Pct, x= conditionSignif))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=c("gray10","gray65"))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nNormal: ",table(p2g$Both)[1],"  |  Cancer: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

ggsave("./Basal_Cohort_Results-updates-Revised-updated/encode_proportionBarchart-univariate.pdf",width = 4,height = 4)

encodePropUni <- df %>%
  ggplot(aes(fill=ENCODE, y=Pct, x= conditionSignif))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=c("gray10","gray65"))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nNormal: ",table(p2g$Both)[1],"  |  Cancer: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

pdf("Basal_Cohort_Results-updates-Revised-updated/peakType_and_ENCODE_overlap_proportion_bar_charts.pdf",width=7,height=4)
cowplot::plot_grid(peakTypePropUni,encodePropUni, ncol = 2, nrow = 1)
dev.off()

pdf("Basal_Cohort_Results-updates-Revised-updated/peakType_and_ENCODE_overlap_proportion_bar_charts-NoLegend.pdf",width=5,height=4)
cowplot::plot_grid(peakTypePropUni+NoLegend(),
                   encodePropUni+NoLegend(), ncol = 2, nrow = 1)
dev.off()

p2g <- p2g[p2g$peakType %in% c("Intronic","Distal"),]
p2g$FDR <- p.adjust(p2g$LME_interaction_model_peak_cell_type_interaction_pval,method = "fdr")
p2g <- p2g[p2g$FDR < 1e-04,]

saveRDS(p2g,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction.rds")

table(p2g$DiffClass)
table(p2g$Both)

# peakType proportion barchart
df <- p2g %>% dplyr::group_by_at("Both") %>% dplyr::count(peakType) %>% mutate(pct= prop.table(n) * 100)
colnames(df) <- c("conditionSignif","peakType","Cells","Pct")
df$peakType <- factor(df$peakType,levels=c("Intronic","Distal","Promoter","Exonic"))
df %>%
  ggplot(aes(fill=peakType, y=Pct, x= conditionSignif, label = Pct))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=rev(RColorBrewer::brewer.pal(5,"Set1")[-1]))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nNormal: ",table(p2g$Both)[1],"  |  Cancer: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

ggsave("./Basal_Cohort_Results-updates-Revised-updated/peakType_proportionBarchart-interaction.pdf",width = 4, height = 4)

peakTypePropInt <- df %>%
  ggplot(aes(fill=peakType, y=Pct, x= conditionSignif, label = Pct))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=rev(RColorBrewer::brewer.pal(5,"Set1")[-1]))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nNormal: ",table(p2g$Both)[1],"  |  Cancer: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

# ENCDOE proportion barchart
df <- p2g %>% dplyr::group_by_at("Both") %>% dplyr::count(ENCODE_Peak_Overlap) %>% mutate(pct= prop.table(n) * 100)
colnames(df) <- c("conditionSignif","ENCODE","Cells","Pct")
df$peakType <- factor(df$ENCODE)
df %>%
  ggplot(aes(fill=ENCODE, y=Pct, x= conditionSignif))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=c("gray10","gray65"))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nNormal: ",table(p2g$Both)[1],"  |  Cancer: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

ggsave("./Basal_Cohort_Results-updates-Revised-updated/encode_proportionBarchart-interaction.pdf",width = 4,height = 4)

encodePropInt <- df %>%
  ggplot(aes(fill=ENCODE, y=Pct, x= conditionSignif))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=c("gray10","gray65"))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nNormal: ",table(p2g$Both)[1],"  |  Cancer: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

# Diff P2G scatter plots and bar chart

# Custom function to get default ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(3)

p2gScatterBoth <- ggplot(p2g,aes(LME_cell_type_1_model_peak_effect_size,
                                 LME_cell_type_2_model_peak_effect_size,
                                 color=Both))+geom_point(size=0.25,alpha=0.5)+
  theme_bw()+
  xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
  ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values=c(cols[2],cols[1],cols[3]))+
  ggtitle(paste(nrow(p2g)," significant differential P2Gs"))

cols <- brewer.pal(n = 8, name = "Dark2")
cols[6] <- "gray40"
cols[8] <- "gray20"

colsToUse <- c(cols[3],cols[4],
               cols[1],cols[2],
               cols[7],cols[5],
               cols[8],cols[6])

df <- p2g %>% count(DiffClass) %>% arrange(desc(n))
df$DiffClass <- factor(df$DiffClass,levels=df$DiffClass)

DiffClassBar <- ggplot(df,aes(x=DiffClass,y=n))+
  geom_bar(stat='identity',
           fill=c(colsToUse[1],colsToUse[4],
                  colsToUse[3],colsToUse[2],
                  colsToUse[5],colsToUse[7],
                  colsToUse[6],colsToUse[8]))+
  geom_text(aes(label= n),
            position=position_stack(vjust=1.05)) +
  theme_bw()+
  ggtitle(paste(nrow(p2g)," significant differential P2Gs"))

# Relevel DiffClass order
p2g$DiffClass <- factor(p2g$DiffClass,levels=df$DiffClass)

p2gScatterDiffClass <- ggplot(p2g,aes(LME_cell_type_1_model_peak_effect_size,
                                      LME_cell_type_2_model_peak_effect_size,
                                      color=DiffClass))+
  geom_point(size=0.25,alpha=0.5)+
  theme_bw()+
  xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
  ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values=c(colsToUse[1],colsToUse[4],
                              colsToUse[3],colsToUse[2],
                              colsToUse[5],colsToUse[7],
                              colsToUse[6],colsToUse[8]))+
  ggtitle(paste(nrow(p2g)," significant differential P2Gs"))

pdf("Basal_Cohort_Results-updates-Revised-updated/p2g_scatter_plots_and_DiffClass_bar_chart.pdf",width=18,height=6)
cowplot::plot_grid(p2gScatterBoth, p2gScatterDiffClass, DiffClassBar, ncol = 3, nrow = 1)
dev.off()

pdf("Basal_Cohort_Results-updates-Revised-updated/p2g_scatter_plots_and_DiffClass_bar_chart-NoLegend.pdf",width=14,height=6)
cowplot::plot_grid( 
  p2gScatterBoth+NoLegend(), 
  p2gScatterDiffClass+NoLegend(), 
  DiffClassBar+NoLegend(), ncol = 3, nrow = 1)
dev.off()

# Split DiffClass Scatter into 8 different plots

p2gScatterDiffClassFacet <-  p2g %>% 
                                    group_by(DiffClass) %>%
                                    mutate(median_es_cancer = median(LME_cell_type_1_model_peak_effect_size)) %>%
                                    mutate(median_es_normal = median(LME_cell_type_2_model_peak_effect_size)) %>%
                                    ungroup() %>%
                                    ggplot(aes(LME_cell_type_1_model_peak_effect_size,
                                                                             LME_cell_type_2_model_peak_effect_size,
                                                                             color=DiffClass))+
                                    geom_point(size=0.25,alpha=0.5)+
                                    theme_bw()+facet_wrap(vars(DiffClass),scales="free")+
                                    geom_vline(aes(xintercept = median_es_cancer),linetype="dotted",linewidth=1)+
                                    geom_hline(aes(yintercept = median_es_normal),linetype="dotted",linewidth=1)+
                                    xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
                                    ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
                                    guides(colour = guide_legend(override.aes = list(size=3)))+
                                    scale_color_manual(values=c(colsToUse[1],colsToUse[4],
                                                                colsToUse[3],colsToUse[2],
                                                                colsToUse[5],colsToUse[7],
                                                                colsToUse[6],colsToUse[8]))+
                                    ggtitle(paste(nrow(p2g)," significant differential P2Gs"))
ggsave("Basal_Cohort_Results-updates-Revised-updated/p2g_scatter_plots-FacetByDiffClass.pdf",width=8,height=8)


# mat <- as.matrix(data.frame(cancer=p2g$LME_cell_type_1_model_peak_effect_size,
#                   normal=p2g$LME_cell_type_2_model_peak_effect_size))

# idx <- sample(1:nrow(mat),2000)

# split by a vector
# fh = function(x) fastcluster::hclust(dist(x))
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs-Heatmap.pdf",width=8,height=16)
# ComplexHeatmap::Heatmap(mat, name = "beta effect size", 
#         row_split = p2g$DiffClass,
#         cluster_rows = fh,
#         cluster_columns = fh)
# dev.off()

# Function to plot heatmap, scatter plots, gene set enrichment

plotDiffClassPlots <- function(p2g,
                               DiffClass,
                               up.genes,
                               up.peaks,
                               geneOfInterest,
                               brca.gr,
                               encode.gr,
                               path,
                               seed){
  
  # Screen for DiffClass P2Gs
  table(p2g$DiffClass)
  table(p2g$Both)
  
  p2g <- p2g[p2g$DiffClass == DiffClass,]
  
  saveRDS(p2g,paste0(path,"/p2gs_interaction-",gsub(pattern="/",replacement = "_",x=DiffClass),".rds"))
  
  if(!is.null(up.genes)){
    
    if(is.null(up.peaks)){
      
      # Screen for up genes
      p2g.up.gene <- p2g[p2g$LME_cell_type_1_model_geneName %in% up.genes$gene,]
      
      saveRDS(p2g.up.gene,paste0(path,"/p2gs_interaction-",gsub(pattern="/",replacement = "_",x=DiffClass),"-upgenes.rds"))
      
      p2g.use <- p2g.up.gene
    
    }else{
      
      # Screen for up genes and up peaks
      p2g.up.gene.up.peak <- p2g[p2g$LME_cell_type_1_model_geneName %in% up.genes$gene &
                                   p2g$LME_cell_type_1_model_peakName %in% up.peaks$gene,]
      
      saveRDS(p2g.up.gene.up.peak,paste0(path,"/p2gs_interaction-",gsub(pattern="/",replacement = "_",x=DiffClass),"-upgenes-uppeaks.rds"))
      
      p2g.use <- p2g.up.gene.up.peak
      
    }
    
  }else{
    p2g.use <- p2g
  }
  
  # Annotate ENCODE and TCGA regions
  
  p2g.use.gr <- GRanges(unique(p2g.use$peakName))
  
  # Overlap with references
  overlaps <- subsetByOverlaps(x = p2g.use.gr,
                               ranges = brca.gr)
  
  brcaPct <- length(overlaps) / length(p2g.use.gr)
  
  p2g.use$TCGA_BRCA_Peak_Overlap <- ifelse(p2g.use$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
  
  overlaps <- subsetByOverlaps(x = p2g.use.gr,
                               ranges = encode.gr)
  
  encodePct <- length(overlaps) / length(p2g.use.gr)
  
  p2g.use$ENCODE_Peak_Overlap <- ifelse(p2g.use$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
  
  if(nrow(p2g.use) > 10000){
    set.seed(seed)
    p2g.use <- p2g.use[sample(1:nrow(p2g.use),10000),]
  }
  
  # Plot Heatmap of effect sizes
  mat <- as.matrix(data.frame(normal=p2g.use$LME_cell_type_2_model_peak_effect_size,
                              cancer=p2g.use$LME_cell_type_1_model_peak_effect_size))
  
  ha = HeatmapAnnotation(
    ENCODE = p2g.use$ENCODE_Peak_Overlap, 
    TCGA_BRCA = p2g.use$TCGA_BRCA_Peak_Overlap,
    col = list(
      ENCODE=c("TRUE"="grey70","FALSE"="gray10"),
      TCGA_BRCA=c("TRUE"="gray70","FALSE"="gray10")
    ),
    which="row",
    annotation_name_side = "top",
    annotation_name_gp = grid::gpar(fontsize = 8)
  )
  
  library(circlize)
  col_fun = colorRamp2(c(quantile(mat, probs = seq(.01, .99, by = .01))[1],
                         0,
                         quantile(mat, probs = seq(.01, .99, by = .01))[99]),
                       c("#CCCCFF", "white", "red"))
  
  if(!is.null(up.genes)){
    if(is.null(up.peaks)){
      pdf(paste0(path,"/Significant_Differential_P2Gs-Heatmap_",gsub(pattern="/",replacement = "_",x=DiffClass),"_upgenes.pdf"),width=4.15,height=5.3)
     
    }else{
      pdf(paste0(path,"/Significant_Differential_P2Gs-Heatmap_",gsub(pattern="/",replacement = "_",x=DiffClass),"_upgenes_uppeaks.pdf"),width=4.15,height=5.3)
    }
  }else{
    pdf(paste0(path,"/Significant_Differential_P2Gs-Heatmap_",gsub(pattern="/",replacement = "_",x=DiffClass),".pdf"),width=4.15,height=5.3)
  }

  ht <- ComplexHeatmap::Heatmap(mat, 
                                name = "beta",
                                col = col_fun, 
                                right_annotation = ha,row_split=ifelse(p2g.use$ENCODE_Peak_Overlap == TRUE,"Annotated","Unannotated"),
                                top_annotation = HeatmapAnnotation(Condition = c("Normal","Cancer"),
                                                                   col = list(Condition=c("Normal"="#828282","Cancer"="#a84343")),
                                                                   which = "column",
                                                                   annotation_name_side = "left",
                                                                   annotation_name_gp = grid::gpar(fontsize = 8)),
                                cluster_columns = FALSE,
                                heatmap_legend_param = list(legend_direction = "vertical"),
                                column_title_gp = grid::gpar(fontsize = 8),
                                column_names_gp = grid::gpar(fontsize = 8),
                                column_title = paste0( round(encodePct*100,2),"% of ",
                                                       length(unique(p2g.use$peakName)),
                                                       " unique peaks and\n",
                                                       as.numeric(round((table(p2g.use$ENCODE_Peak_Overlap)[2]/nrow(mat))*100,2)),
                                                       "% of ", nrow(mat), " total P2Gs\noverlap with ENCODE annotations |\n",
                                                       round(brcaPct*100,2),"% of ",
                                                       length(unique(p2g.use$peakName)),
                                                       " unique peaks and\n",
                                                       as.numeric(round((table(p2g.use$TCGA_BRCA_Peak_Overlap)[2]/nrow(mat))*100,2)),
                                                       "% of ", nrow(mat), " total P2Gs\noverlap with TCGA BRCA bulk ATAC peaks"))
  draw(ht,heatmap_legend_side = "right",annotation_legend_side = "right")
  
  dev.off()
  
  if(!is.null(up.genes)){
    # Plot Heatmap of effect sizes ranked by differential gene FC 
    df.int <- data.frame(normal=p2g.use$LME_cell_type_2_model_peak_effect_size,
                         cancer=p2g.use$LME_cell_type_1_model_peak_effect_size,
                         gene=p2g.use$LME_cell_type_1_model_geneName)
    df.int <- merge(df.int,up.genes,by="gene")
    
    if(sign(df.int$log2FoldChange) < 0){
      df.int$log2FoldChange <- df.int$log2FoldChange*-1
    }else{
      df.int$log2FoldChange <- df.int$log2FoldChange
    }
    
    df.int <- dplyr::arrange(df.int,desc(log2FoldChange))
    mat <- as.matrix(data.frame(normal=df.int$normal,
                                cancer=df.int$cancer))
    
    ha = HeatmapAnnotation(
      log2FC = df.int$log2FoldChange, 
      col = list(
        log2FC=colorRamp2(c(min(df.int$log2FoldChange), max(df.int$log2FoldChange)), 
                          c("#fee4d2", "#c65102"))
      ),
      which="row",
      annotation_name_side = "top",
      annotation_name_gp = grid::gpar(fontsize = 8)
    )
    
    library(circlize)
    col_fun = colorRamp2(c(quantile(mat, probs = seq(.01, .99, by = .01))[1],
                           0,
                           quantile(mat, probs = seq(.01, .99, by = .01))[99]),
                         c("#CCCCFF", "white", "red"))
    
    if(!is.null(up.genes)){
      if(is.null(up.peaks)){
        pdf(paste0(path,"/Significant_Differential_P2Gs-Heatmap_",gsub(pattern="/",replacement = "_",x=DiffClass),"_upgenes-rankedByLog2FC.pdf"),width=4.15,height=5.3)
        
      }else{
        pdf(paste0(path,"/Significant_Differential_P2Gs-Heatmap_",gsub(pattern="/",replacement = "_",x=DiffClass),"_upgenes_uppeaks-rankedByLog2FC.pdf"),width=4.15,height=5.3)
      }
    }else{
      pdf(paste0(path,"/Significant_Differential_P2Gs-Heatmap_",gsub(pattern="/",replacement = "_",x=DiffClass),"-rankedByLog2FC.pdf"),width=4.15,height=5.3)
    }
    
    if(is.null(geneOfInterest)){
      ht <- ComplexHeatmap::Heatmap(mat, 
                                    name = "beta",
                                    col = col_fun, 
                                    right_annotation = ha, 
                                    top_annotation = HeatmapAnnotation(Condition = c("Normal","Cancer"),
                                                                       col = list(Condition=c("Normal"="#828282","Cancer"="#a84343")),
                                                                       which = "column",
                                                                       annotation_name_side = "left",
                                                                       annotation_name_gp = grid::gpar(fontsize = 8)),
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE,
                                    heatmap_legend_param = list(legend_direction = "vertical"),
                                    column_title_gp = grid::gpar(fontsize = 8),
                                    column_names_gp = grid::gpar(fontsize = 8),
                                    column_title = paste0( round(encodePct*100,2),"% of ",
                                                           length(unique(p2g.use$peakName)),
                                                           " unique peaks and\n",
                                                           as.numeric(round((table(p2g.use$ENCODE_Peak_Overlap)[2]/nrow(mat))*100,2)),
                                                           "% of ", nrow(mat), " total P2Gs\noverlap with ENCODE annotations |\n",
                                                           round(brcaPct*100,2),"% of ",
                                                           length(unique(p2g.use$peakName)),
                                                           " unique peaks and\n",
                                                           as.numeric(round((table(p2g.use$TCGA_BRCA_Peak_Overlap)[2]/nrow(mat))*100,2)),
                                                           "% of ", nrow(mat), " total P2Gs\noverlap with TCGA BRCA bulk ATAC peaks"))
    }else{
      ht <- ComplexHeatmap::Heatmap(mat, 
                                    name = "beta",
                                    col = col_fun, 
                                    right_annotation = c(ha,rowAnnotation(geneName = anno_mark(at = grep(paste0("\\b",geneOfInterest,"\\b"),df.int$gene), 
                                                                                                                              labels = rep(geneOfInterest,length(grep(paste0("\\b",geneOfInterest,"\\b"),df.int$gene)))))), 
                                    top_annotation = HeatmapAnnotation(Condition = c("Normal","Cancer"),
                                                                       col = list(Condition=c("Normal"="#828282","Cancer"="#a84343")),
                                                                       which = "column",
                                                                       annotation_name_side = "left",
                                                                       annotation_name_gp = grid::gpar(fontsize = 8)),
                                    cluster_columns = FALSE,
                                    cluster_rows = FALSE,
                                    heatmap_legend_param = list(legend_direction = "vertical"),
                                    column_title_gp = grid::gpar(fontsize = 8),
                                    column_names_gp = grid::gpar(fontsize = 8),
                                    column_title = paste0( round(encodePct*100,2),"% of ",
                                                           length(unique(p2g.use$peakName)),
                                                           " unique peaks and\n",
                                                           as.numeric(round((table(p2g.use$ENCODE_Peak_Overlap)[2]/nrow(mat))*100,2)),
                                                           "% of ", nrow(mat), " total P2Gs\noverlap with ENCODE annotations |\n",
                                                           round(brcaPct*100,2),"% of ",
                                                           length(unique(p2g.use$peakName)),
                                                           " unique peaks and\n",
                                                           as.numeric(round((table(p2g.use$TCGA_BRCA_Peak_Overlap)[2]/nrow(mat))*100,2)),
                                                           "% of ", nrow(mat), " total P2Gs\noverlap with TCGA BRCA bulk ATAC peaks"))
    }
    
    draw(ht,heatmap_legend_side = "right",annotation_legend_side = "right")
    
    dev.off()
    
  }
 
  # Scatter plots of subsetted P2Gs
  p1 <- ggplot(p2g.use,aes(LME_cell_type_1_model_peak_effect_size,
                               LME_cell_type_2_model_peak_effect_size))+geom_hex(bins = 200) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+
    xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
    ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
    ggtitle(paste(nrow(p2g.use)," significant differential P2Gs"))
  
  
  p2 <- ggplot(p2g.use,aes(LME_cell_type_1_model_peak_effect_size,
                               LME_cell_type_2_model_peak_effect_size,
                               color=Both))+geom_point(size=0.5,alpha=0.5)+
    theme_bw()+
    xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
    ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    ggtitle(paste(nrow(p2g.use)," significant differential P2Gs"))
  
  p3 <- ggplot(p2g.use,aes(LME_cell_type_1_model_peak_effect_size,
                               LME_cell_type_2_model_peak_effect_size,
                               color=DiffClass))+geom_point(size=0.5,alpha=0.5)+
    theme_bw()+
    xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
    ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    ggtitle(paste(nrow(p2g.use)," significant differential P2Gs"))
  
  
  if(!is.null(up.genes)){
    if(is.null(up.peaks)){
      cowplot::plot_grid(p1,
                         p2,
                         p3,
                         ncol = 3, nrow = 1)
      ggsave(paste0(path,"/Significant_Differential_P2Gs_",gsub(pattern="/",replacement = "_",x=DiffClass),"_upgenes.pdf"),width=18,height=4)
    }else{
      cowplot::plot_grid(p1,
                         p2,
                         p3,
                         ncol = 3, nrow = 1)
      ggsave(paste0(path,"/Significant_Differential_P2Gs_",gsub(pattern="/",replacement = "_",x=DiffClass),"_upgenes_uppeaks.pdf"),width=18,height=4)
    }
  }else{
    cowplot::plot_grid(p1,
                       p2,
                       p3,
                       ncol = 3, nrow = 1)
    ggsave(paste0(path,"/Significant_Differential_P2Gs_",gsub(pattern="/",replacement = "_",x=DiffClass),".pdf"),width=18,height=4)
  }

  # Gene set ORA 
  gset = msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::distinct(gs_name, gene_symbol)
  # Do +/0 first
  em.p2g.use <- clusterProfiler::enricher(gene=unique(p2g.use$LME_cell_type_1_model_geneName), 
                                              TERM2GENE=gset,
                                              pvalueCutoff=1,
                                              qvalueCutoff=1)
  if(!is.null(up.genes)){
    if(is.null(up.peaks)){
      saveRDS(em.p2g.use,paste0(path,"/em_p2g_",gsub(pattern="/",replacement = "_",x=DiffClass),"_up_gene.rds"))
    }else{
      saveRDS(em.p2g.use,paste0(path,"/em_p2g_",gsub(pattern="/",replacement = "_",x=DiffClass),"_up_gene_up_peak.rds"))
    }
  }else{
    saveRDS(em.p2g.use,paste0(path,"/em_p2g_",gsub(pattern="/",replacement = "_",x=DiffClass),".rds"))
  }
  
  if(!is.null(up.genes)){
    if(is.null(up.peaks)){
      enrichplot::dotplot(em.p2g.use,showCategory=3)+
        ggtitle(paste0(length(unique(p2g.use$LME_cell_type_1_model_geneName))," unique genes"))
      ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-upgenes-DotPlot.pdf"), 
             width = 8,
             height = 6)
    }else{
      enrichplot::dotplot(em.p2g.use,showCategory=3)+
        ggtitle(paste0(length(unique(p2g.use$LME_cell_type_1_model_geneName))," unique genes"))
      ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-upgenes-uppeaks-DotPlot.pdf"), 
                    width = 8,
                    height = 6)
    }
  }else{
    enrichplot::dotplot(em.p2g.use,showCategory=3)+
      ggtitle(paste0(length(unique(p2g.use$LME_cell_type_1_model_geneName))," unique genes"))
    ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-DotPlot.pdf"), 
                  width = 8,
                  height = 6)
  }

  if(!is.null(up.genes)){
    if(is.null(up.peaks)){
      
      barplot(em.p2g.use,showCategory=3)+
        ggtitle(paste0(length(unique(p2g.use$LME_cell_type_1_model_geneName))," unique genes"))
      
      ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-upgenes-BarPlot_Count.pdf"), 
                    width = 8,
                    height = 6)
    }else{
      
      barplot(em.p2g.use,showCategory=3)+
        ggtitle(paste0(length(unique(p2g.use$LME_cell_type_1_model_geneName))," unique genes"))
      
      ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-upgenes-uppeaks-BarPlot_Count.pdf"), 
                    width = 8,
                    height = 6)
    }
  }else{
    
    barplot(em.p2g.use,showCategory=3)+
      ggtitle(paste0(length(unique(p2g.use$LME_cell_type_1_model_geneName))," unique genes"))
    
    ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-BarPlot_Count.pdf"), 
                  width = 8,
                  height = 6)
  }

  if(!is.null(up.genes)){
    if(is.null(up.peaks)){
      mutate(em.p2g.use, qscore = -log(p.adjust, base=10)) %>% 
        barplot(x="qscore",showCategory=3)
      ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-upgenes-BarPlot_QScore.pdf"), 
                    width = 8,
                    height = 6)
    }else{
      mutate(em.p2g.use, qscore = -log(p.adjust, base=10)) %>% 
        barplot(x="qscore",showCategory=3)
      ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-upgenes-uppeaks-BarPlot_QScore.pdf"), 
                    width = 8,
                    height = 6)
    }
  }else{
    mutate(em.p2g.use, qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore",showCategory=3)
    ggsave(paste0(path,"/GeneSet_ORA_Hallmark_MSigDB_",gsub(pattern="/",replacement = "_",x=DiffClass),"-BarPlot_QScore.pdf"), 
                  width = 8,
                  height = 6)
  }

}

# Downloaded TCGA BRCA peaks from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG 
brca <- read.delim("./miscellaneous/BRCA_peakCalls.txt",sep = "\t")
brca.gr <- GRanges(paste0(brca$seqnames,":",brca$start,"-",brca$end))

# Downloaded ENCODE cCREs from https://screen.encodeproject.org 
encode <- read.delim("./miscellaneous/GRCh38-cCREs.bed",header = F)
encode.gr <- GRanges(paste0(encode$V1,":",encode$V2,"-",encode$V3))

# Plot visuals for cancer gained enhancers for up genes in cancer v. normal

# Read in genes upregulated in cancer v. normal 
up.genes <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/up_genes.rds")
# up.peaks <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/up_peaks.rds")

plotDiffClassPlots(p2g=p2g,
                   DiffClass="+/0",
                   up.genes=up.genes,
                   up.peaks=NULL,
                   geneOfInterest = "HEY1",
                   brca.gr=brca.gr,
                   encode.gr=encode.gr,
                   seed=123,
                   path="Basal_Cohort_Results-updates-Revised-updated")

# Plot visuals for switched enhancer (silencer in normal to enhancer in cancer)
plotDiffClassPlots(p2g=p2g,
                   DiffClass="+/-",
                   up.genes=up.genes,
                   up.peaks=NULL,
                   geneOfInterest = NULL,
                   brca.gr=brca.gr,
                   encode.gr=encode.gr,
                   seed=123,
                   path="Basal_Cohort_Results-updates-Revised-updated")

# Plot visuals for normal gained enhancers for genes upregulated in normal v. cancer 
degs <- readRDS(paste0("./Basal_Cohort_Results-updates-Revised-updated/DEGs_padj_",padj_cutoff,"_log2FC_",log2fc_cutoff,".rds"))
up.genes <- dplyr::filter(degs, log2FoldChange <= -log2fc_cutoff)

plotDiffClassPlots(p2g=p2g,
                   DiffClass="0/+",
                   up.genes=up.genes,
                   up.peaks=NULL,
                   geneOfInterest = NULL,
                   brca.gr=brca.gr,
                   encode.gr=encode.gr,
                   seed=123,
                   path="Basal_Cohort_Results-updates-Revised-updated")

# Plot visuals for all P2Gs in each DiffClass (if greater than 10,000, downsample to 10,000 for visualization purposes)
for( i in c( "+/+", "+/0", "+/-",
             "0/+", "0/-", "-/+", "-/0", "-/-")){
  plotDiffClassPlots(p2g=p2g,
                     DiffClass=i,
                     up.genes=NULL,
                     up.peaks=NULL,
                     geneOfInterest = NULL,
                     brca.gr=brca.gr,
                     encode.gr=encode.gr,
                     seed=123,
                     path="Basal_Cohort_Results-updates-Revised-updated")
  
  print(paste0("Completed DiffClass: ",i," !"))
  
}


# 
# 
# # Screen for +/0 P2Gs
# table(p2g$DiffClass)
# table(p2g$Both)
# 
# p2g <- p2g[p2g$DiffClass == "+/0",]
# 
# saveRDS(p2g,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction-+_0.rds")
# 
# # Screen for up genes
# p2g.up.gene <- p2g[p2g$LME_cell_type_1_model_geneName %in% up.genes$gene,]
# 
# saveRDS(p2g.up.gene,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction-+_0-upgenes.rds")
# 
# # Annotate ENCODE and TCGA regions
# 
# # Downloaded TCGA BRCA peaks from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG 
# brca <- read.delim("./miscellaneous/BRCA_peakCalls.txt",sep = "\t")
# brca.gr <- GRanges(paste0(brca$seqnames,":",brca$start,"-",brca$end))
# 
# # Downloaded ENCODE cCREs from https://screen.encodeproject.org 
# encode <- read.delim("./miscellaneous/GRCh38-cCREs.bed",header = F)
# encode.gr <- GRanges(paste0(encode$V1,":",encode$V2,"-",encode$V3))
# 
# p2g.up.gene.gr <- GRanges(unique(p2g.up.gene$peakName))
# 
# # Overlap with references
# overlaps <- subsetByOverlaps(x = p2g.up.gene.gr,
#                              ranges = brca.gr)
# 
# brcaPct <- length(overlaps) / length(p2g.up.gene.gr)
# 
# p2g.up.gene$TCGA_BRCA_Peak_Overlap <- ifelse(p2g.up.gene$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
# 
# overlaps <- subsetByOverlaps(x = p2g.up.gene.gr,
#                              ranges = encode.gr)
# 
# encodePct <- length(overlaps) / length(p2g.up.gene.gr)
# 
# p2g.up.gene$ENCODE_Peak_Overlap <- ifelse(p2g.up.gene$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
# 
# # Plot Heatmap of effect sizes
# mat <- as.matrix(data.frame(normal=p2g.up.gene$LME_cell_type_2_model_peak_effect_size,
#                             cancer=p2g.up.gene$LME_cell_type_1_model_peak_effect_size))
# 
# ha = HeatmapAnnotation(
#   ENCODE = p2g.up.gene$ENCODE_Peak_Overlap, 
#   TCGA_BRCA = p2g.up.gene$TCGA_BRCA_Peak_Overlap,
#   col = list(
#              ENCODE=c("TRUE"="grey70","FALSE"="gray10"),
#              TCGA_BRCA=c("TRUE"="gray70","FALSE"="gray10")
#   ),
#   which="row",
#   annotation_name_side = "top",
#   annotation_name_gp = grid::gpar(fontsize = 8)
# )
# 
# library(circlize)
# col_fun = colorRamp2(c(quantile(mat, probs = seq(.01, .99, by = .01))[1],
#                        0,
#                        quantile(mat, probs = seq(.01, .99, by = .01))[99]),
#                     c("#CCCCFF", "white", "red"))
# 
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs-Heatmap_+_0_upgenes.pdf",width=4.15,height=5.3)
# ht <- ComplexHeatmap::Heatmap(mat, 
#                               name = "beta",
#                               col = col_fun, 
#                               right_annotation = ha, 
#                               top_annotation = HeatmapAnnotation(Condition = c("Normal","Cancer"),
#                                                                  col = list(Condition=c("Normal"="#828282","Cancer"="#a84343")),
#                                                                  which = "column",
#                                                                  annotation_name_side = "left",
#                                                                  annotation_name_gp = grid::gpar(fontsize = 8)),
#                               cluster_columns = FALSE,
#                               heatmap_legend_param = list(legend_direction = "vertical"),
#                               column_title_gp = grid::gpar(fontsize = 8),
#                               column_names_gp = grid::gpar(fontsize = 8),
#                               column_title = paste0( round(encodePct*100,2),"% of ",
#                                                      length(unique(p2g.up.gene$peakName)),
#                                                      " unique peaks and\n",
#                                                      as.numeric(round((table(p2g.up.gene$ENCODE_Peak_Overlap)[2]/nrow(mat))*100,2)),
#                                                      "% of ", nrow(mat), " total P2Gs\noverlap with ENCODE annotations |\n",
#                                                      round(brcaPct*100,2),"% of ",
#                                                      length(unique(p2g.up.gene$peakName)),
#                                                      " unique peaks and\n",
#                                                      as.numeric(round((table(p2g.up.gene$TCGA_BRCA_Peak_Overlap)[2]/nrow(mat))*100,2)),
#                                                      "% of ", nrow(mat), " total P2Gs\noverlap with TCGA BRCA bulk ATAC peaks"))
# draw(ht,heatmap_legend_side = "right",annotation_legend_side = "right")
# 
# dev.off()
# 
# # Scatter plots of subsetted P2Gs
# p1 <- ggplot(p2g.up.gene,aes(LME_cell_type_1_model_peak_effect_size,
#                              LME_cell_type_2_model_peak_effect_size))+geom_hex(bins = 200) +
#   scale_fill_continuous(type = "viridis") +
#   theme_bw()+
#   xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
#   ggtitle(paste(nrow(p2g.up.gene)," significant differential P2Gs"))
# 
# 
# p2 <- ggplot(p2g.up.gene,aes(LME_cell_type_1_model_peak_effect_size,
#                              LME_cell_type_2_model_peak_effect_size,
#                              color=Both))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw()+
#   xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   ggtitle(paste(nrow(p2g.up.gene)," significant differential P2Gs"))
# 
# p3 <- ggplot(p2g.up.gene,aes(LME_cell_type_1_model_peak_effect_size,
#                              LME_cell_type_2_model_peak_effect_size,
#                              color=DiffClass))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw()+
#   xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   ggtitle(paste(nrow(p2g.up.gene)," significant differential P2Gs"))
# 
# 
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs_+_0_upgenes.pdf",width=18,height=4)
# p1+p2+p3
# dev.off()
# 
# # Gene set ORA 
# gset = msigdbr(species = "Homo sapiens", category = "H") %>% 
#   dplyr::distinct(gs_name, gene_symbol)
# # Do +/0 first
# em.p2g.up.gene <- clusterProfiler::enricher(gene=unique(p2g.up.gene$LME_cell_type_1_model_geneName), 
#                                             TERM2GENE=gset,
#                                             pvalueCutoff=1,
#                                             qvalueCutoff=1)
# 
# saveRDS(em.p2g.up.gene,"./Basal_Cohort_Results-updates-Revised-updated/em_p2g_+_0_up_gene.rds")
# 
# enrichplot::dotplot(em.p2g.up.gene,showCategory=3)+
#   ggtitle(paste0(length(unique(p2g.up.gene$LME_cell_type_1_model_geneName))," unique genes"))
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_+_0_upgenes-DotPlot.pdf", 
#        width = 8,
#        height = 6)
# 
# barplot(em.p2g.up.gene,showCategory=3)+
#   ggtitle(paste0(length(unique(p2g.up.gene$LME_cell_type_1_model_geneName))," unique genes"))
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_+_0_upgenes-BarPlot_Count.pdf", 
#        width = 8,
#        height = 6)
# 
# mutate(em.p2g.up.gene, qscore = -log(p.adjust, base=10)) %>% 
#   barplot(x="qscore",showCategory=3)
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_+_0_upgenes-BarPlot_QScore.pdf", 
#        width = 8,
#        height = 6)
# 
# # Screen for up genes and peaks
# p2g.up.gene.up.peak <- p2g[p2g$LME_cell_type_1_model_geneName %in% up.genes$gene &
#                              p2g$LME_cell_type_1_model_peakName %in% up.peaks$gene,]
# 
# saveRDS(p2g.up.gene.up.peak,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction-+_0-upgenes_uppeaks.rds")
# 
# # Annotate ENCODE and TCGA regions
# 
# # Downloaded TCGA BRCA peaks from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG 
# brca <- read.delim("./miscellaneous/BRCA_peakCalls.txt",sep = "\t")
# brca.gr <- GRanges(paste0(brca$seqnames,":",brca$start,"-",brca$end))
# 
# # Downloaded ENCODE cCREs from https://screen.encodeproject.org 
# encode <- read.delim("./miscellaneous/GRCh38-cCREs.bed",header = F)
# encode.gr <- GRanges(paste0(encode$V1,":",encode$V2,"-",encode$V3))
# 
# p2g.up.gene.up.peak.gr <- GRanges(unique(p2g.up.gene.up.peak$peakName))
# 
# # Overlap with references
# overlaps <- subsetByOverlaps(x = p2g.up.gene.up.peak.gr,
#                              ranges = brca.gr)
# 
# brcaPct <- length(overlaps) / length(p2g.up.gene.up.peak.gr)
# 
# p2g.up.gene.up.peak$TCGA_BRCA_Peak_Overlap <- ifelse(p2g.up.gene.up.peak$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
# 
# overlaps <- subsetByOverlaps(x = p2g.up.gene.up.peak.gr,
#                              ranges = encode.gr)
# 
# encodePct <- length(overlaps) / length(p2g.up.gene.up.peak.gr)
# 
# p2g.up.gene.up.peak$ENCODE_Peak_Overlap <- ifelse(p2g.up.gene.up.peak$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
# 
# # Plot Heatmap of effect sizes
# mat <- as.matrix(data.frame(normal=p2g.up.gene.up.peak$LME_cell_type_2_model_peak_effect_size,
#                             cancer=p2g.up.gene.up.peak$LME_cell_type_1_model_peak_effect_size))
# 
# ha = HeatmapAnnotation(
#   ENCODE = p2g.up.gene.up.peak$ENCODE_Peak_Overlap, 
#   TCGA_BRCA = p2g.up.gene.up.peak$TCGA_BRCA_Peak_Overlap,
#   col = list(
#     ENCODE=c("TRUE"="grey70","FALSE"="gray10"),
#     TCGA_BRCA=c("TRUE"="gray70","FALSE"="gray10")
#   ),
#   which="row",
#   annotation_name_side = "top",
#   annotation_name_gp = grid::gpar(fontsize = 8)
# )
# 
# library(circlize)
# col_fun = colorRamp2(c(quantile(mat, probs = seq(.01, .99, by = .01))[1],
#                        0,
#                        quantile(mat, probs = seq(.01, .99, by = .01))[99]),
#                      c("#CCCCFF", "white", "red"))
# 
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs-Heatmap_+_0_upgenes_uppeaks.pdf",width=4.15,height=5.3)
# ht <- ComplexHeatmap::Heatmap(mat, 
#                               name = "beta",
#                               col = col_fun, 
#                               right_annotation = ha, 
#                               top_annotation = HeatmapAnnotation(Condition = c("Normal","Cancer"),
#                                                                  col = list(Condition=c("Normal"="#828282","Cancer"="#a84343")),
#                                                                  which = "column",
#                                                                  annotation_name_side = "left",
#                                                                  annotation_name_gp = grid::gpar(fontsize = 8)),
#                               cluster_columns = FALSE,
#                               heatmap_legend_param = list(legend_direction = "vertical"),
#                               column_title_gp = grid::gpar(fontsize = 8),
#                               column_names_gp = grid::gpar(fontsize = 8),
#                               column_title = paste0( round(encodePct*100,2),"% of ",
#                                                      length(unique(p2g.up.gene.up.peak$peakName)),
#                                                      " unique peaks and\n",
#                                                      as.numeric(round((table(p2g.up.gene.up.peak$ENCODE_Peak_Overlap)[2]/nrow(mat))*100,2)),
#                                                      "% of ", nrow(mat), " total P2Gs\noverlap with ENCODE annotations |\n",
#                                                      round(brcaPct*100,2),"% of ",
#                                                      length(unique(p2g.up.gene.up.peak$peakName)),
#                                                      " unique peaks and\n",
#                                                      as.numeric(round((table(p2g.up.gene.up.peak$TCGA_BRCA_Peak_Overlap)[2]/nrow(mat))*100,2)),
#                                                      "% of ", nrow(mat), " total P2Gs\noverlap with TCGA BRCA bulk ATAC peaks"))
# draw(ht,heatmap_legend_side = "right",annotation_legend_side = "right")
# 
# dev.off()
# 
# # Scatter plots of subsetted P2Gs
# p1 <- ggplot(p2g.up.gene.up.peak,aes(LME_cell_type_1_model_peak_effect_size,
#                                      LME_cell_type_2_model_peak_effect_size))+geom_hex(bins = 200) +
#   scale_fill_continuous(type = "viridis") +
#   theme_bw()+
#   xlab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in normal univariate model"))+
#   ggtitle(paste(nrow(p2g.up.gene.up.peak)," significant differential P2Gs"))
# 
# 
# p2 <- ggplot(p2g.up.gene.up.peak,aes(LME_cell_type_1_model_peak_effect_size,
#                                      LME_cell_type_2_model_peak_effect_size,
#                                      color=Both))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw()+
#   xlab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in normal univariate model"))+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   ggtitle(paste(nrow(p2g.up.gene.up.peak)," significant differential P2Gs"))
# 
# p3 <- ggplot(p2g.up.gene.up.peak,aes(LME_cell_type_1_model_peak_effect_size,
#                                      LME_cell_type_2_model_peak_effect_size,
#                                      color=DiffClass))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw()+
#   xlab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in normal univariate model"))+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   ggtitle(paste(nrow(p2g.up.gene.up.peak)," significant differential P2Gs"))
# 
# 
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs_+_0_upgenes_uppeaks.pdf",width=18,height=6)
# p1+p2+p3
# dev.off()
# 
# # Gene set ORA 
# gset = msigdbr(species = "Homo sapiens", category = "H") %>% 
#   dplyr::distinct(gs_name, gene_symbol)
# # Do +/0 first
# em.p2g.up.gene.up.peak <- clusterProfiler::enricher(gene=unique(p2g.up.gene.up.peak$LME_cell_type_1_model_geneName), 
#                                                     TERM2GENE=gset,
#                                                     pvalueCutoff=1,
#                                                     qvalueCutoff=1)
# 
# saveRDS(em.p2g.up.gene.up.peak,"./Basal_Cohort_Results-updates-Revised-updated/em_p2g_+_0_up_gene_up_peak.rds")
# 
# enrichplot::dotplot(em.p2g.up.gene.up.peak,showCategory=3)+
#   ggtitle(paste0(length(unique(p2g.up.gene.up.peak$LME_cell_type_1_model_geneName))," unique genes"))
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_+_0_upgenes_uppeaks-DotPlot.pdf")
# 
# barplot(em.p2g.up.gene.up.peak,showCategory=3)+
#   ggtitle(paste0(length(unique(p2g.up.gene.up.peak$LME_cell_type_1_model_geneName))," unique genes"))
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_+_0_upgenes_uppeaks-BarPlot_Count.pdf")
# 
# mutate(em.p2g.up.gene.up.peak, qscore = -log(p.adjust, base=10)) %>% 
#   barplot(x="qscore",showCategory=3)
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_+_0_upgenes_uppeaks-BarPlot_QScore.pdf")
# 
# 
# # Screen for 0/+ P2Gs
# p2g <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction.rds")
# 
# table(p2g$DiffClass)
# table(p2g$Both)
# 
# p2g <- p2g[p2g$DiffClass == "0/+",]
# 
# saveRDS(p2g,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction-0_+.rds")
# 
# # Screen for up genes
# # Get up genes in normal v. cancer
# degs <- readRDS(paste0("./Basal_Cohort_Results-updates-Revised-updated/DEGs_padj_",padj_cutoff,"_log2FC_",log2fc_cutoff,".rds"))
# up.genes <- dplyr::filter(degs, log2FoldChange <= -log2fc_cutoff)
# 
# p2g.up.gene <- p2g[p2g$LME_cell_type_1_model_geneName %in% up.genes$gene,]
# 
# saveRDS(p2g.up.gene,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction-0_+-upgenes.rds")
# 
# # Annotate ENCODE and TCGA regions
# 
# # Downloaded TCGA BRCA peaks from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG 
# brca <- read.delim("./miscellaneous/BRCA_peakCalls.txt",sep = "\t")
# brca.gr <- GRanges(paste0(brca$seqnames,":",brca$start,"-",brca$end))
# 
# # Downloaded ENCODE cCREs from https://screen.encodeproject.org 
# encode <- read.delim("./miscellaneous/GRCh38-cCREs.bed",header = F)
# encode.gr <- GRanges(paste0(encode$V1,":",encode$V2,"-",encode$V3))
# 
# p2g.up.gene.gr <- GRanges(unique(p2g.up.gene$peakName))
# 
# # Overlap with references
# overlaps <- subsetByOverlaps(x = p2g.up.gene.gr,
#                              ranges = brca.gr)
# 
# brcaPct <- length(overlaps) / length(p2g.up.gene.gr)
# 
# p2g.up.gene$TCGA_BRCA_Peak_Overlap <- ifelse(p2g.up.gene$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
# 
# overlaps <- subsetByOverlaps(x = p2g.up.gene.gr,
#                              ranges = encode.gr)
# 
# encodePct <- length(overlaps) / length(p2g.up.gene.gr)
# 
# p2g.up.gene$ENCODE_Peak_Overlap <- ifelse(p2g.up.gene$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
# 
# # Plot Heatmap of effect sizes
# mat <- as.matrix(data.frame(normal=p2g.up.gene$LME_cell_type_2_model_peak_effect_size,
#                             cancer=p2g.up.gene$LME_cell_type_1_model_peak_effect_size))
# 
# ha = HeatmapAnnotation(
#   ENCODE = p2g.up.gene$ENCODE_Peak_Overlap, 
#   TCGA_BRCA = p2g.up.gene$TCGA_BRCA_Peak_Overlap,
#   col = list(
#     ENCODE=c("TRUE"="grey70","FALSE"="gray10"),
#     TCGA_BRCA=c("TRUE"="gray70","FALSE"="gray10")
#   ),
#   which="row",
#   annotation_name_side = "top",
#   annotation_name_gp = grid::gpar(fontsize = 8)
# )
# 
# library(circlize)
# col_fun = colorRamp2(c(quantile(mat, probs = seq(.01, .99, by = .01))[1],
#                        0,
#                        quantile(mat, probs = seq(.01, .99, by = .01))[99]),
#                      c("#CCCCFF", "white", "red"))
# 
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs-Heatmap_0_+_upgenes.pdf",width=4.15,height=5.3)
# ht <- ComplexHeatmap::Heatmap(mat, 
#                               name = "beta",
#                               col = col_fun, 
#                               right_annotation = ha, 
#                               top_annotation = HeatmapAnnotation(Condition = c("Normal","Cancer"),
#                                                                  col = list(Condition=c("Normal"="#828282","Cancer"="#a84343")),
#                                                                  which = "column",
#                                                                  annotation_name_side = "left",
#                                                                  annotation_name_gp = grid::gpar(fontsize = 8)),
#                               cluster_columns = FALSE,
#                               heatmap_legend_param = list(legend_direction = "vertical"),
#                               column_title_gp = grid::gpar(fontsize = 8),
#                               column_names_gp = grid::gpar(fontsize = 8),
#                               column_title = paste0( round(encodePct*100,2),"% of ",
#                                                      length(unique(p2g.up.gene$peakName)),
#                                                      " unique peaks and\n",
#                                                      as.numeric(round((table(p2g.up.gene$ENCODE_Peak_Overlap)[2]/nrow(mat))*100,2)),
#                                                      "% of ", nrow(mat), " total P2Gs\noverlap with ENCODE annotations |\n",
#                                                      round(brcaPct*100,2),"% of ",
#                                                      length(unique(p2g.up.gene$peakName)),
#                                                      " unique peaks and\n",
#                                                      as.numeric(round((table(p2g.up.gene$TCGA_BRCA_Peak_Overlap)[2]/nrow(mat))*100,2)),
#                                                      "% of ", nrow(mat), " total P2Gs\noverlap with TCGA BRCA bulk ATAC peaks"))
# draw(ht,heatmap_legend_side = "right",annotation_legend_side = "right")
# 
# dev.off()
# 
# # Scatter plots of subsetted P2Gs
# p1 <- ggplot(p2g.up.gene,aes(LME_cell_type_1_model_peak_effect_size,
#                              LME_cell_type_2_model_peak_effect_size))+geom_hex(bins = 200) +
#   scale_fill_continuous(type = "viridis") +
#   theme_bw()+
#   xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
#   ggtitle(paste(nrow(p2g.up.gene)," significant differential P2Gs"))
# 
# 
# p2 <- ggplot(p2g.up.gene,aes(LME_cell_type_1_model_peak_effect_size,
#                              LME_cell_type_2_model_peak_effect_size,
#                              color=Both))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw()+
#   xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   ggtitle(paste(nrow(p2g.up.gene)," significant differential P2Gs"))
# 
# p3 <- ggplot(p2g.up.gene,aes(LME_cell_type_1_model_peak_effect_size,
#                              LME_cell_type_2_model_peak_effect_size,
#                              color=DiffClass))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw()+
#   xlab(expression("P2G" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("P2G" ~ beta ~ "effect size in normal univariate model"))+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   ggtitle(paste(nrow(p2g.up.gene)," significant differential P2Gs"))
# 
# 
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs_0_+_upgenes.pdf",width=18,height=4)
# p1+p2+p3
# dev.off()
# 
# # Gene set ORA 
# gset = msigdbr(species = "Homo sapiens", category = "H") %>% 
#   dplyr::distinct(gs_name, gene_symbol)
# # Do 0/+ first
# em.p2g.up.gene <- clusterProfiler::enricher(gene=unique(p2g.up.gene$LME_cell_type_1_model_geneName), 
#                                             TERM2GENE=gset,
#                                             pvalueCutoff=1,
#                                             qvalueCutoff=1)
# 
# saveRDS(em.p2g.up.gene,"./Basal_Cohort_Results-updates-Revised-updated/em_p2g_0_+_up_gene.rds")
# 
# enrichplot::dotplot(em.p2g.up.gene,showCategory=3)+
#   ggtitle(paste0(length(unique(p2g.up.gene$LME_cell_type_1_model_geneName))," unique genes"))
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_0_+_upgenes-DotPlot.pdf", 
#        width = 8,
#        height = 6)
# 
# barplot(em.p2g.up.gene,showCategory=3)+
#   ggtitle(paste0(length(unique(p2g.up.gene$LME_cell_type_1_model_geneName))," unique genes"))
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_0_+_upgenes-BarPlot_Count.pdf", 
#        width = 8,
#        height = 6)
# 
# mutate(em.p2g.up.gene, qscore = -log(p.adjust, base=10)) %>% 
#   barplot(x="qscore",showCategory=3)
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_0_+_upgenes-BarPlot_QScore.pdf", 
#        width = 8,
#        height = 6)
# 
# # Screen for up genes and peaks
# 
# # Get up genes in normal v. cancer
# daps <- readRDS(paste0("./Basal_Cohort_Results-updates-Revised-updated/DAPs_padj_",padj_cutoff,"_log2FC_",log2fc_cutoff,".rds"))
# up.peaks <- dplyr::filter(daps, log2FoldChange <= -log2fc_cutoff)
# 
# p2g.up.gene.up.peak <- p2g[p2g$LME_cell_type_1_model_geneName %in% up.genes$gene &
#                              p2g$LME_cell_type_1_model_peakName %in% up.peaks$gene,]
# 
# saveRDS(p2g.up.gene.up.peak,"./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction-0_+-upgenes_uppeaks.rds")
# 
# # Annotate ENCODE and TCGA regions
# 
# # Downloaded TCGA BRCA peaks from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG 
# brca <- read.delim("./miscellaneous/BRCA_peakCalls.txt",sep = "\t")
# brca.gr <- GRanges(paste0(brca$seqnames,":",brca$start,"-",brca$end))
# 
# # Downloaded ENCODE cCREs from https://screen.encodeproject.org 
# encode <- read.delim("./miscellaneous/GRCh38-cCREs.bed",header = F)
# encode.gr <- GRanges(paste0(encode$V1,":",encode$V2,"-",encode$V3))
# 
# p2g.up.gene.up.peak.gr <- GRanges(unique(p2g.up.gene.up.peak$peakName))
# 
# # Overlap with references
# overlaps <- subsetByOverlaps(x = p2g.up.gene.up.peak.gr,
#                              ranges = brca.gr)
# 
# brcaPct <- length(overlaps) / length(p2g.up.gene.up.peak.gr)
# 
# p2g.up.gene.up.peak$TCGA_BRCA_Peak_Overlap <- ifelse(p2g.up.gene.up.peak$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
# 
# overlaps <- subsetByOverlaps(x = p2g.up.gene.up.peak.gr,
#                              ranges = encode.gr)
# 
# encodePct <- length(overlaps) / length(p2g.up.gene.up.peak.gr)
# 
# p2g.up.gene.up.peak$ENCODE_Peak_Overlap <- ifelse(p2g.up.gene.up.peak$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
# 
# # Plot Heatmap of effect sizes
# mat <- as.matrix(data.frame(normal=p2g.up.gene.up.peak$LME_cell_type_2_model_peak_effect_size,
#                             cancer=p2g.up.gene.up.peak$LME_cell_type_1_model_peak_effect_size))
# 
# ha = HeatmapAnnotation(
#   ENCODE = p2g.up.gene.up.peak$ENCODE_Peak_Overlap, 
#   TCGA_BRCA = p2g.up.gene.up.peak$TCGA_BRCA_Peak_Overlap,
#   col = list(
#     ENCODE=c("TRUE"="grey70","FALSE"="gray10"),
#     TCGA_BRCA=c("TRUE"="gray70","FALSE"="gray10")
#   ),
#   which="row",
#   annotation_name_side = "top",
#   annotation_name_gp = grid::gpar(fontsize = 8)
# )
# 
# library(circlize)
# col_fun = colorRamp2(c(quantile(mat, probs = seq(.01, .99, by = .01))[1],
#                        0,
#                        quantile(mat, probs = seq(.01, .99, by = .01))[99]),
#                      c("#CCCCFF", "white", "red"))
# 
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs-Heatmap_0_+_upgenes_uppeaks.pdf",width=4.15,height=5.3)
# ht <- ComplexHeatmap::Heatmap(mat, 
#                               name = "beta",
#                               col = col_fun, 
#                               right_annotation = ha, 
#                               top_annotation = HeatmapAnnotation(Condition = c("Normal","Cancer"),
#                                                                  col = list(Condition=c("Normal"="#828282","Cancer"="#a84343")),
#                                                                  which = "column",
#                                                                  annotation_name_side = "left",
#                                                                  annotation_name_gp = grid::gpar(fontsize = 8)),
#                               cluster_columns = FALSE,
#                               heatmap_legend_param = list(legend_direction = "vertical"),
#                               column_title_gp = grid::gpar(fontsize = 8),
#                               column_names_gp = grid::gpar(fontsize = 8),
#                               column_title = paste0( round(encodePct*100,2),"% of ",
#                                                      length(unique(p2g.up.gene.up.peak$peakName)),
#                                                      " unique peaks and\n",
#                                                      as.numeric(round((table(p2g.up.gene.up.peak$ENCODE_Peak_Overlap)[2]/nrow(mat))*100,2)),
#                                                      "% of ", nrow(mat), " total P2Gs\noverlap with ENCODE annotations |\n",
#                                                      round(brcaPct*100,2),"% of ",
#                                                      length(unique(p2g.up.gene.up.peak$peakName)),
#                                                      " unique peaks and\n",
#                                                      as.numeric(round((table(p2g.up.gene.up.peak$TCGA_BRCA_Peak_Overlap)[2]/nrow(mat))*100,2)),
#                                                      "% of ", nrow(mat), " total P2Gs\noverlap with TCGA BRCA bulk ATAC peaks"))
# draw(ht,heatmap_legend_side = "right",annotation_legend_side = "right")
# 
# dev.off()
# 
# # Scatter plots of subsetted P2Gs
# p1 <- ggplot(p2g.up.gene.up.peak,aes(LME_cell_type_1_model_peak_effect_size,
#                                      LME_cell_type_2_model_peak_effect_size))+geom_hex(bins = 200) +
#   scale_fill_continuous(type = "viridis") +
#   theme_bw()+
#   xlab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in normal univariate model"))+
#   ggtitle(paste(nrow(p2g.up.gene.up.peak)," significant differential P2Gs"))
# 
# 
# p2 <- ggplot(p2g.up.gene.up.peak,aes(LME_cell_type_1_model_peak_effect_size,
#                                      LME_cell_type_2_model_peak_effect_size,
#                                      color=Both))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw()+
#   xlab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in normal univariate model"))+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   ggtitle(paste(nrow(p2g.up.gene.up.peak)," significant differential P2Gs"))
# 
# p3 <- ggplot(p2g.up.gene.up.peak,aes(LME_cell_type_1_model_peak_effect_size,
#                                      LME_cell_type_2_model_peak_effect_size,
#                                      color=DiffClass))+geom_point(size=0.5,alpha=0.5)+
#   theme_bw()+
#   xlab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in cancer univariate model"))+
#   ylab(expression("p2g.up.gene.up.peak" ~ beta ~ "effect size in normal univariate model"))+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   ggtitle(paste(nrow(p2g.up.gene.up.peak)," significant differential P2Gs"))
# 
# 
# pdf("Basal_Cohort_Results-updates-Revised-updated/Significant_Differential_P2Gs_0_+_upgenes_uppeaks.pdf",width=18,height=6)
# p1+p2+p3
# dev.off()
# 
# # Gene set ORA 
# gset = msigdbr(species = "Homo sapiens", category = "H") %>% 
#   dplyr::distinct(gs_name, gene_symbol)
# # Do 0/+ first
# em.p2g.up.gene.up.peak <- clusterProfiler::enricher(gene=unique(p2g.up.gene.up.peak$LME_cell_type_1_model_geneName), 
#                                                     TERM2GENE=gset,
#                                                     pvalueCutoff=1,
#                                                     qvalueCutoff=1)
# 
# saveRDS(em.p2g.up.gene.up.peak,"./Basal_Cohort_Results-updates-Revised-updated/em_p2g_0_+_up_gene_up_peak.rds")
# 
# enrichplot::dotplot(em.p2g.up.gene.up.peak,showCategory=3)+
#   ggtitle(paste0(length(unique(p2g.up.gene.up.peak$LME_cell_type_1_model_geneName))," unique genes"))
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_0_+_upgenes_uppeaks-DotPlot.pdf")
# 
# barplot(em.p2g.up.gene.up.peak,showCategory=3)+
#   ggtitle(paste0(length(unique(p2g.up.gene.up.peak$LME_cell_type_1_model_geneName))," unique genes"))
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_0_+_upgenes_uppeaks-BarPlot_Count.pdf")
# 
# mutate(em.p2g.up.gene.up.peak, qscore = -log(p.adjust, base=10)) %>% 
#   barplot(x="qscore",showCategory=3)
# ggsave("./Basal_Cohort_Results-updates-Revised-updated/GeneSet_ORA_Hallmark_MSigDB_0_+_upgenes_uppeaks-BarPlot_QScore.pdf")

################################################################################
# Plot Browser Tracks, Expression Boxplots, and P2G scatter plots 
################################################################################

# Function to plot browser track, expression boxplot and P2G scatter plot
plotP2Gs <- function(ArchRProj, 
                     p2g,
                     meta_cell_data,
                     encode.gr,
                     brca.gr,
                     cancerCol,
                     normalCol,
                     patientCols,
                     patientScatterCols,
                     allLevels,
                     cancerLevels,
                     normalLevels,
                     cancerSpecificEnh,
                     cancerEnh,
                     normalEnh,
                     dds,
                     valueMin,
                     valueMax,
                     outpath,
                     identifier,
                     seed,
                     highlightDistance
                     ){
  
  # Get loops 
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
  o$P2G <- paste0(o$peakName,"|",o$geneName)
  o <- o[o$P2G %in% p2g$P2G,]
  
  o <- o[order(match(o$P2G,p2g$P2G)),]
  
  p2g$idxATAC <- o$idxATAC
  p2g$idxRNA <- o$idxRNA
  
  sub <- p2g
  
  # Make loop track
  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  #Gene Info
  geneSet <- .getFeatureDF(getArrowFiles(ArchRProj),"GeneIntegrationMatrix", threads = max(floor(getArchRThreads() / 2), 1))
  
  peakSummits <- resize(peakSet, 1, "center")
  geneStarts <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
  
  resolution <- 1
  summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
  geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
  
  # Cancer loops
  cancerLoops <- .constructGR(
    seqnames = seqnames(peakSummits[ sub$idxATAC ]),
    start = summitTiles[ sub$idxATAC ],
    end = geneTiles[ sub$idxRNA ]
  )
  mcols(cancerLoops)$value <- sub$LME_cell_type_1_model_peak_effect_size
  mcols(cancerLoops)$FDR <- sub$Cancer_FDR
  
  cancerLoops <- cancerLoops[order(mcols(cancerLoops)$value, decreasing=TRUE)]
  cancerLoops <- unique(cancerLoops)
  cancerLoops <- cancerLoops[width(cancerLoops) > 0]
  cancerLoops <- sort(sortSeqlevels(cancerLoops))
  
  # Normal loops
  normalLoops <- .constructGR(
    seqnames = seqnames(peakSummits[ sub$idxATAC ]),
    start = summitTiles[ sub$idxATAC ],
    end = geneTiles[ sub$idxRNA ]
  )
  mcols(normalLoops)$value <- sub$LME_cell_type_2_model_peak_effect_size
  mcols(normalLoops)$FDR <- sub$Normal_FDR
  
  normalLoops <- normalLoops[order(mcols(normalLoops)$value, decreasing=TRUE)]
  normalLoops <- unique(normalLoops)
  normalLoops <- normalLoops[width(normalLoops) > 0]
  normalLoops <- sort(sortSeqlevels(normalLoops))
  
  returnList <- list(cancerLoops,normalLoops)
  
  set.seed(seed)
  plotbrowserBoth <- plotBrowserTrack.mod.3(ArchRProj,
                                            geneSymbol = unique(sub$geneName),
                                            upstream = ifelse(sub[grep(min(sub$peakName),sub$peakName),]$downstreamGene == TRUE,sub[grep(min(sub$peakName),sub$peakName),]$distToLinkedGeneStart+2e+04,2e+04),
                                            downstream = ifelse(sub[grep(max(sub$peakName),sub$peakName),]$downstreamGene == TRUE,2e+04,sub[grep(max(sub$peakName),sub$peakName),]$distToLinkedGeneStart+2e+04),
                                            groupBy = "Patient",
                                            groupByTop = "Condition",
                                            valueMin=valueMin,
                                            valueMax=valueMax,
                                            features = GRangesList(ENCODE_cCREs=encode.gr,
                                                                   TCGA_BRCA=brca.gr,
                                                                   MACS2_Peaks=getPeakSet(ArchRProj),
                                                                   # Diff_Peaks=GRanges(unique(sub$peakName)),
                                                                   cancer_specific_Enh=GRanges(sort(unique(cancerSpecificEnh$peakName))),
                                                                   all_cancer_Enh=GRanges(sort(unique(cancerEnh$peakName))),
                                                                   all_normal_Enh=GRanges(sort(unique(normalEnh$peakName)))
                                                                   
                                            ),
                                            loops.1 = returnList[[1]],
                                            loops.2 = returnList[[2]],
                                            # sizes = c(14, 10, 5, 5, 5, 16),
                                            sizes = c(8, 12, 4, 3, 3, 12),
                                            plotSummary = c("bulktracktop", "bulktrack", "featureTrack", "looptrack1",
                                                            "looptrack2",  "geneTrack"),
                                            palTop = c(cancerCol,normalCol),
                                            pal = patientCols,
                                            useMatrix = "GeneIntegrationMatrix",
                                            dds = dds,
                                            seed = seed,
                                            highlightDistance = highlightDistance,
                                            peakCoords = GRanges(sort(unique(sub$peakName)))
  )
  
  set.seed(seed)
  plotbrowserPatient <- plotBrowserTrack.mod.3(ArchRProj,
                                               geneSymbol = unique(sub$geneName),
                                               upstream = ifelse(sub[grep(min(sub$peakName),sub$peakName),]$downstreamGene == TRUE,sub[grep(min(sub$peakName),sub$peakName),]$distToLinkedGeneStart+2e+04,2e+04),
                                               downstream = ifelse(sub[grep(max(sub$peakName),sub$peakName),]$downstreamGene == TRUE,2e+04,sub[grep(max(sub$peakName),sub$peakName),]$distToLinkedGeneStart+2e+04),
                                               groupBy = "Patient",
                                               groupByTop = NULL,
                                               valueMin=valueMin,
                                               valueMax=valueMax,
                                               features = GRangesList(ENCODE_cCREs=encode.gr,
                                                                      TCGA_BRCA=brca.gr,
                                                                      MACS2_Peaks=getPeakSet(ArchRProj),
                                                                      # Diff_Peaks=GRanges(unique(sub$peakName)),
                                                                      cancer_specific_Enh=GRanges(sort(unique(cancerSpecificEnh$peakName))),
                                                                      all_cancer_Enh=GRanges(sort(unique(cancerEnh$peakName))),
                                                                      all_normal_Enh=GRanges(sort(unique(normalEnh$peakName)))
                                                                      
                                               ),
                                               loops.1 = returnList[[1]],
                                               loops.2 = returnList[[2]],
                                               # sizes = c(12, 6, 4, 4, 12),
                                               sizes = c(12, 4, 3, 3, 12),
                                               plotSummary = c("BulkTrack", "featureTrack", "looptrack1",
                                                               "looptrack2",  "geneTrack"),
                                               pal = patientCols,
                                               palTop = NULL,
                                               useMatrix = "GeneIntegrationMatrix",
                                               dds = dds,
                                               seed = seed,
                                               highlightDistance = highlightDistance,
                                               peakCoords = GRanges(sort(unique(sub$peakName)))
  )
  
  set.seed(seed)
  plotbrowserCondition <- plotBrowserTrack.mod.3(ArchRProj,
                                                 geneSymbol = unique(sub$geneName),
                                                 upstream = ifelse(sub[grep(min(sub$peakName),sub$peakName),]$downstreamGene == TRUE,sub[grep(min(sub$peakName),sub$peakName),]$distToLinkedGeneStart+2e+04,2e+04),
                                                 downstream = ifelse(sub[grep(max(sub$peakName),sub$peakName),]$downstreamGene == TRUE,2e+04,sub[grep(max(sub$peakName),sub$peakName),]$distToLinkedGeneStart+2e+04),
                                                 groupBy = "Condition",
                                                 groupByTop = NULL,
                                                 valueMin=valueMin,
                                                 valueMax=valueMax,
                                                 features = GRangesList(ENCODE_cCREs=encode.gr,
                                                                        TCGA_BRCA=brca.gr,
                                                                        MACS2_Peaks=getPeakSet(ArchRProj),
                                                                        # Diff_Peaks=GRanges(unique(sub$peakName)),
                                                                        cancer_specific_Enh=GRanges(sort(unique(cancerSpecificEnh$peakName))),
                                                                        all_cancer_Enh=GRanges(sort(unique(cancerEnh$peakName))),
                                                                        all_normal_Enh=GRanges(sort(unique(normalEnh$peakName)))
                                                                        
                                                 ),
                                                 loops.1 = returnList[[1]],
                                                 loops.2 = returnList[[2]],
                                                 # sizes = c(12, 6, 4, 4, 12),
                                                 sizes = c(12, 4, 3, 3, 12),
                                                 plotSummary = c("BulkTrack", "featureTrack", "looptrack1",
                                                                 "looptrack2",  "geneTrack"),
                                                 palTop = NULL,
                                                 pal = c(cancerCol,normalCol),
                                                 useMatrix = "GeneIntegrationMatrix",
                                                 dds = dds,
                                                 seed = seed,
                                                 highlightDistance = highlightDistance,
                                                 peakCoords = GRanges(sort(unique(sub$peakName)))
  )
  
  # Scatter plots
  # Create RNA and ATAC matrices
  ATAC <- assay(meta_cell_data[[2]])
  rownames(ATAC) <- paste0(meta_cell_data[[2]]@rowRanges@seqnames,":",meta_cell_data[[2]]@rowRanges@ranges)
  colnames(ATAC) <- paste0("metacell_",1:ncol(ATAC))
  
  RNA <- assay(meta_cell_data[[1]])
  rownames(RNA) <- meta_cell_data[[1]]@rowRanges$name
  colnames(RNA) <- paste0("metacell_",1:ncol(RNA))
  
  # # Make metadata
  meta <- data.frame(cell_type = meta_cell_data[[3]],
                     patient = meta_cell_data[[5]])
  rownames(meta) <- paste0("metacell_",1:ncol(ATAC))# identical metacells between RNA and ATAC
  
  peak_mat = ATAC
  gene_mat = RNA
  meta_data = meta
  covariates = NULL
  random_effect = "patient"
  cell_type = "cell_type"
  cell_type_1 = "Cancer"
  cell_type_2 = "Normal"
  meta_cells = TRUE
  id = "metacells"
  
  gene <- unique(sub$geneName)
  PlotsList <- vector(mode="list",length=length(unique(sub$peakName)))
  for ( i in 1:length(sort(unique(sub$peakName))) ){
    peak <- sort(unique(sub$peakName))[i]
    # DF
    df <- data.frame(gene = gene_mat[grep(paste0("^",gene,"$"),rownames(gene_mat)),],
                     peak = peak_mat[grep(paste0("^",peak,"$"),rownames(peak_mat)),])
    df <- cbind(df,meta_data[,colnames(meta_data) %in% c(cell_type,random_effect)])
    
    df1 <- df[df$cell_type == cell_type_1,]
    df1$gene <- base::scale(df1$gene)
    if(meta_cells){
      df1$peak <- base::scale(df1$peak)
    }
    df2 <- df[df$cell_type == cell_type_2,]
    df2$gene <- base::scale(df2$gene)
    if(meta_cells){
      df2$peak <- base::scale(df2$peak)
    }
    
    df[[cell_type]] <- ifelse(df[[cell_type]] == cell_type_1,1,0)
    df$gene <- base::scale(df$gene)
    if(meta_cells){
      df$peak <- base::scale(df$peak)
    }
    
    interaction_formula_h1 <- as.formula(paste0("gene ~ peak + ",cell_type,paste0(" + (1|",random_effect,")"," + ","peak:",cell_type)))
    interaction_formula_h1
    mod <- lmerTest::lmer(interaction_formula_h1,data = df,REML = FALSE)
    df$fit <- predict(mod)
    
    cols <- patientScatterCols
    
    df$patient <- factor(df$patient,levels=allLevels)
    plotscatterInteraction <- ggplot(df,aes(peak, gene, color = factor(patient))) +
      facet_grid(~factor(patient)) +
      geom_line(aes(y=fit), size=0.5) +
      scale_color_manual(values=cols)+
      geom_point(alpha = 0.4,size=0.3)+
      ylab(paste0(gene))+
      xlab(paste0(peak))+theme_bw()
    
    cell_type_1_formula_h1 <- as.formula(paste0("gene ~ peak + ",paste0("(1|",random_effect,")")))
    cell_type_1_formula_h1
    mod <- lmerTest::lmer(cell_type_1_formula_h1,data = df1,REML = FALSE)
    df1$fit <- predict(mod)
    
    df1$patient <- factor(df1$patient,levels=cancerLevels)
    plotscatterCellType1 <- ggplot(df1,aes(peak, gene, color = factor(patient))) +
      facet_grid(~factor(patient)) +
      geom_line(aes(y=fit), size=0.5) +
      scale_color_manual(values=rep(cancerCol,length(unique(df1$patient))))+
      geom_point(alpha = 0.4,size=0.3)+
      ylab(paste0(gene))+
      xlab(paste0(peak))+theme_bw()
    
    cell_type_2_formula_h1 <- as.formula(paste0("gene ~ peak + ",paste0("(1|",random_effect,")")))
    cell_type_2_formula_h1
    mod <- lmerTest::lmer(cell_type_2_formula_h1,data = df2,REML = FALSE)
    df2$fit <- predict(mod)
    
    df2$patient <- factor(df2$patient,levels=normalLevels)
    plotscatterCellType2 <- ggplot(df2,aes(peak, gene, color = factor(patient))) +
      facet_grid(~factor(patient)) +
      geom_line(aes(y=fit), size=0.5) +
      scale_color_manual(values=rep(normalCol,length(unique(df2$patient))))+
      geom_point(alpha = 0.4,size=0.3)+
      ylab(paste0(gene))+
      xlab(paste0(peak))+theme_bw()
    
    PlotsList[[i]] <- list(plotscatterInteraction,
                           plotscatterCellType1, 
                           plotscatterCellType2
    )
  }
  
  pdf(paste0(outpath,"P2G_Condition_and_Patient_BrowserTracks_ExpressionPlots_ScatterPlots-",
             unique(sub$geneName),"-", identifier,".pdf"),width=8,height=11)
  
  grid::grid.draw(plotbrowserBoth[[unique(sub$geneName)]])
  for(i in 1:length(PlotsList)){
    
    for( j in 1:length(PlotsList[[i]])){
      print(PlotsList[[i]][[j]])
    }
    
  }
  dev.off()
  
  pdf(paste0(outpath,"P2G_Condition_BrowserTracks_ExpressionPlots_ScatterPlots-",
             unique(sub$geneName),"-", identifier,".pdf"),width=8,height=5.5)
  
  grid::grid.draw(plotbrowserCondition[[unique(sub$geneName)]])
  for(i in 1:length(PlotsList)){
    
    for( j in 1:length(PlotsList[[i]])){
      print(PlotsList[[i]][[j]])
    }
    
  }
  dev.off()
  
  pdf(paste0(outpath,"P2G_Patient_BrowserTracks_ExpressionPlots_ScatterPlots-",
             unique(sub$geneName),"-", identifier,".pdf"),width=8,height=7)
  
  grid::grid.draw(plotbrowserPatient[[unique(sub$geneName)]])
  for(i in 1:length(PlotsList)){
    
    for( j in 1:length(PlotsList[[i]])){
      print(PlotsList[[i]][[j]])
    }
    
  }
  dev.off()
  
  for( i in 1:length(sort(unique(sub$peakName))) ){
    
    pdf(paste0(outpath,"P2G_ScatterPlots-", 
               sort(unique(sub$peakName))[i],"-",
               unique(sub$geneName),"-", identifier,"-Interaction.pdf"),width=6.5,height=2.5)
    print(PlotsList[[i]][[1]])
    dev.off()
    
    pdf(paste0(outpath,"P2G_ScatterPlots-", 
               sort(unique(sub$peakName))[i],"-",
               unique(sub$geneName),"-", identifier,"-Cancer.pdf"),width=6.5,height=2.5)
    print(PlotsList[[i]][[2]])
    dev.off()
    
    pdf(paste0(outpath,"P2G_ScatterPlots-", 
               sort(unique(sub$peakName))[i],"-",
               unique(sub$geneName),"-", identifier,"-Normal.pdf"),width=6.5,height=2.5)
    print(PlotsList[[i]][[3]])
    dev.off()
    
  }
  
}

# Read in ArchR project
proj <- loadArchRProject(path =  "./Basal_TN_Samples_scATAC-TESTING3")

idxSample <- BiocGenerics::which(proj$predictedGroup %in% c("1","2","3"))
cellsToUse <- proj$cellNames[idxSample]
proj <- proj[cellsToUse, ]

proj$Patient <- proj$Sample
proj$Patient <- plyr::mapvalues(x = proj$Patient, 
                                from = levels(factor(proj$Patient)),
                                to = c("1-35A4AL",
                                       "3-49758L",
                                       "4-49CFCL",
                                       "5-4AF75L",
                                       "6-4B146L",
                                       "2-4C2E5L")
)

proj$Condition <- ifelse(proj$Sample %in% c("49758L","49CFCL","4AF75L","4B146L"),"2-Normal","1-Cancer")

# Read in dds 
dds <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/pseudobulk_expression_dds.rds")

colData(dds)$Patient <- plyr::mapvalues(x =colData(dds)$sample_id, 
                                        from = levels(factor(colData(dds)$sample_id)),
                                        to = c("1-35A4AL",
                                               "3-49758L",
                                               "4-49CFCL",
                                               "5-4AF75L",
                                               "6-4B146L",
                                               "2-4C2E5L")
)

colData(dds)$Patient <- factor(colData(dds)$Patient, levels = c("1-35A4AL",
                                                                "2-4C2E5L",
                                                                "3-49758L",
                                                                "4-49CFCL",
                                                                "5-4AF75L",
                                                                "6-4B146L")
                               )

colData(dds)$Condition <- ifelse(colData(dds)$sample_id %in% c("49758L","49CFCL","4AF75L","4B146L"),"Normal","Cancer")

# Read in diff P2Gs with up genes
p2g <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction-+_0-upgenes.rds")
p2g$geneName <- p2g$LME_cell_type_1_model_geneName

# Add gene info

p2g$nearest <- ifelse(p2g$LME_cell_type_1_model_geneName == p2g$nearestGene,TRUE,FALSE)

genes <- getGeneAnnotation(proj)

genes <- as.data.frame(GenomicRanges::resize(genes$genes, 1, "start"))

genes$geneStart <- paste0(genes$seqnames,":",genes$start,"-",genes$end)

genes$geneName <- genes$symbol

p2g$geneName <- p2g$LME_cell_type_1_model_geneName

p2g <- merge(p2g,genes,by="geneName")

p2g$downstreamGene <- ifelse(GRanges(p2g$geneStart) > GenomicRanges::resize(GRanges(p2g$peakName),1,"center"),TRUE,FALSE)

p2g$distToLinkedGeneStart <- GenomicRanges::distance( GenomicRanges::resize(GRanges(p2g$peakName),1,"center"),
                                                      GRanges(p2g$geneStart))

p2g.HEY1 <- p2g[p2g$geneName == "HEY1",]

# Read in cancerSpecificEnh 
cancerSpecificEnh <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction-+_0-upgenes.rds")
cancerSpecificEnh <- cancerSpecificEnh[cancerSpecificEnh$LME_cell_type_1_model_geneName == "HEY1",]

# Read in cancerEnh and normalEnh 
allP2Gs <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/p2gs_univariate-postfilter.rds")
cancerEnh <- allP2Gs[allP2Gs$DiffClass %in% c("+/+","+/0","+/-") & allP2Gs$peakType %in% c("Intronic","Distal") & allP2Gs$LME_cell_type_1_model_geneName == "HEY1",]
normalEnh <- allP2Gs[allP2Gs$DiffClass %in% c("+/+","0/+","-/+") & allP2Gs$peakType %in% c("Intronic","Distal") & allP2Gs$LME_cell_type_1_model_geneName == "HEY1",]

# Establish colors for cancer and normal conditions, and patient colors based on factor leveling 
cancerCol<- "#a84343"
normalCol<- "#828282"
patientCols <- c(rep(cancerCol,2),rep(normalCol,4))
patientScatterCols <- c(rep(normalCol,4),rep(cancerCol,2))

# Read in Basal meta_cell_data 
meta_cell_data <- readRDS("./LME_Basal_out-SingFits_OLS/meta_cell_data.rds")

# Downloaded TCGA BRCA peaks from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG 
brca <- read.delim("./miscellaneous/BRCA_peakCalls.txt",sep = "\t")
brca.gr <- GRanges(paste0(brca$seqnames,":",brca$start,"-",brca$end))

# Downloaded ENCODE cCREs from https://screen.encodeproject.org 
encode <- read.delim("./miscellaneous/GRCh38-cCREs.bed",header = F)
encode.gr <- GRanges(paste0(encode$V1,":",encode$V2,"-",encode$V3))

# Set up output folders for different sets of P2Gs
outpath = "./Basal_Cohort_Results-updates-Revised-updated/"

# Plots for HEY1
plotP2Gs(ArchRProj = proj, 
         p2g = p2g.HEY1,
         meta_cell_data = meta_cell_data,
         encode.gr = encode.gr,
         brca.gr = brca.gr,
         cancerCol = cancerCol,
         normalCol = normalCol,
         patientCols = patientCols,
         patientScatterCols = patientScatterCols,
         allLevels = c("49758L","49CFCL","4AF75L","4B146L","35A4AL","4C2E5L"),
         cancerLevels = c("35A4AL","4C2E5L"),
         normalLevels = c("49758L","49CFCL","4AF75L","4B146L"),
         cancerSpecificEnh = cancerSpecificEnh,
         cancerEnh = cancerEnh,
         normalEnh = normalEnh,
         dds = dds,
         valueMin = -0.45,
         valueMax = 0.45,
         outpath = outpath,
         seed = 4321,
         highlightDistance = 2000,
         identifier = "All_Cancer_Specific_P2Gs")

# Subset to highest effect size P2G
p2g.HEY1.high <- p2g.HEY1[p2g.HEY1$LME_cell_type_1_model_peak_effect_size == max(p2g.HEY1$LME_cell_type_1_model_peak_effect_size),]

# Plots for HEY1 (highest effect size)
plotP2Gs(ArchRProj = proj, 
         p2g = p2g.HEY1.high,
         meta_cell_data = meta_cell_data,
         encode.gr = encode.gr,
         brca.gr = brca.gr,
         cancerCol = cancerCol,
         normalCol = normalCol,
         patientCols = patientCols,
         patientScatterCols = patientScatterCols,
         allLevels = c("49758L","49CFCL","4AF75L","4B146L","35A4AL","4C2E5L"),
         cancerLevels = c("35A4AL","4C2E5L"),
         normalLevels = c("49758L","49CFCL","4AF75L","4B146L"),
         cancerSpecificEnh = cancerSpecificEnh,
         cancerEnh = cancerEnh,
         normalEnh = normalEnh,
         dds = dds,
         valueMin = -0.45,
         valueMax = 0.45,
         outpath = outpath,
         seed = 4321,
         highlightDistance = 300,
         identifier = "Highest_Effect_Size_P2G")

# Subset to highest effect size P2G and neighboring
p2g.HEY1.high.neighbor <- p2g.HEY1[p2g.HEY1$peakName %in% c("chr8:79793567-79794067","chr8:79794946-79795446", "chr8:79799786-79800286"), ]

# Plots for HEY1 (highest effect size and neighboring)
plotP2Gs(ArchRProj = proj, 
         p2g = p2g.HEY1.high.neighbor,
         meta_cell_data = meta_cell_data,
         encode.gr = encode.gr,
         brca.gr = brca.gr,
         cancerCol = cancerCol,
         normalCol = normalCol,
         patientCols = patientCols,
         patientScatterCols = patientScatterCols,
         allLevels = c("49758L","49CFCL","4AF75L","4B146L","35A4AL","4C2E5L"),
         cancerLevels = c("35A4AL","4C2E5L"),
         normalLevels = c("49758L","49CFCL","4AF75L","4B146L"),
         cancerSpecificEnh = cancerSpecificEnh,
         cancerEnh = cancerEnh,
         normalEnh = normalEnh,
         dds = dds,
         valueMin = -0.45,
         valueMax = 0.45,
         outpath = outpath,
         seed = 4321,
         highlightDistance = 300,
         identifier = "Highest_Effect_Size_P2G_and_Neighbors")

