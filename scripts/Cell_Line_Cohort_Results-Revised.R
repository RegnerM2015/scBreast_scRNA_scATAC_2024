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
library(eulerr)
library(patchwork)
source("scripts/plotBrowserTrack-modified-2_LoopTracks-ExprBoxPlots-TwoBulkTracks.R")
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

# Make poster figure output folder
if(dir.exists("Cell_Line_Cohort_Results-Revised")==TRUE){
  print("Directory Cell_Line_Cohort_Results-Revised already exists!")
}else{
  dir.create("Cell_Line_Cohort_Results-Revised")
}

# UMAPs

# Plot seurat UMAP
seurat <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_CellLine_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

# Plot UMAP colored by cell type cluster
cols <- c("#662377","#EF5064","#A82973","#FC867D")

df <- as.data.frame(seurat@reductions$umap@cell.embeddings)
all.equal(rownames(df),rownames(seurat@meta.data))
df$cellType <- as.character(seurat$orig.ident)

p1 <- ggplot(df,aes(x =UMAP_1,y=UMAP_2,color = cellType))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values=cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle(paste0("scRNA-seq by cell type cluster\nn=",nrow(df)," cells"))

# Plot UMAP colored by cell type cluster
cols <- c("#662377","#EF5064","#A82973","#FC867D")

proj <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")
df <- plotEmbedding(proj,colorBy = "cellColData",name = "predictedGroup",embedding = "UMAP")
df <- as.data.frame(df$data)
df$color <- sub(".*?-","",df$color)
df$cellType <- as.character(df$color)

p2 <- ggplot(df,aes(x =x,y=y,color = cellType))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values=cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle(paste0("scATAC-seq by cell type cluster\nn=",nrow(df)," cells"))

p1 + p2

ggsave(paste0("./Cell_Line_Cohort_Results-Revised/UMAPs_",
              "CellLines",".pdf"),
       width = 14,
       height = 6 )



# Plot QC histograms

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_CellLine_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

df <- as.data.frame(rna@meta.data)
qcRNAnCount <- ggplot(df, aes(x=log10(nCount_RNA))) + 
  geom_histogram(color="gray20", fill="gray75",bins=60,boundary=log10(min(df$nCount_RNA)))+
  theme_bw()+
  geom_vline(xintercept = log10(min(df$nCount_RNA)),linetype="dashed",col="darkred")+
  geom_vline(xintercept = log10(median(df$nCount_RNA)),linetype="dashed",col="black")+
  ggtitle(paste0( paste0("median # RNA counts: ",round( median(df$nCount_RNA), 2) ), "\n",
                  paste0( "min # RNA counts: ",min(df$nCount_RNA)), "\n",
                  paste0("n=",nrow(df)," cells")))

qcRNAnFeature <- ggplot(df, aes(x=nFeature_RNA)) + 
  geom_histogram(color="gray20", fill="gray75",bins=60,boundary=min(df$nFeature_RNA))+
  theme_bw()+
  geom_vline(xintercept = min(df$nFeature_RNA),linetype="dashed",col="darkred")+
  geom_vline(xintercept = median(df$nFeature_RNA),linetype="dashed",col="black")+
  ggtitle(paste0( paste0("median # RNA features: ",round( median(df$nFeature_RNA), 2)) , "\n",
                  paste0( "min # RNA features: ",min(df$nFeature_RNA)), "\n",
                  paste0("n=",nrow(df)," cells")))

# Read in scATAC
atac <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")

df <- as.data.frame(atac@cellColData)
qcATACTSS <- ggplot(df, aes(x=TSSEnrichment)) + 
  geom_histogram(color="gray20", fill="gray75",bins=60,boundary=min(df$TSSEnrichment))+
  theme_bw()+
  geom_vline(xintercept = min(df$TSSEnrichment),linetype="dashed",col="darkred")+
  geom_vline(xintercept = median(df$TSSEnrichment),linetype="dashed",col="black")+
  ggtitle(paste0( paste0("median TSSEnrichment: ",round( median(df$TSSEnrichment), 2)), "\n",
                  paste0( "min TSSEnrichment: ",min(df$TSSEnrichment)), "\n",
                  paste0("n=",nrow(df)," cells")))

qcATACnFrags <- ggplot(df, aes(x=log10(nFrags))) + 
  geom_histogram(color="gray20", fill="gray75",bins=60,boundary=log10(min(df$nFrags)))+
  theme_bw()+
  geom_vline(xintercept = log10(min(df$nFrags)),linetype="dashed",col="darkred")+
  geom_vline(xintercept = log10(median(df$nFrags)),linetype="dashed",col="black")+
  ggtitle(paste0( paste0("median # unique fragments: ",round( median(df$nFrags), 2)), "\n",
                  paste0( "min # unique fragments: ",min(df$nFrags)), "\n",
                  paste0("n=",nrow(df)," cells")))

qcRNAnCount+qcRNAnFeature+qcATACnFrags+qcATACTSS+plot_layout(ncol=2)
ggsave("./Cell_Line_Cohort_Results-Revised/QC_Histograms_RNA_ATAC.pdf",width = 9,height = 8)

# Plot QC metrics stratified by sample
# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_CellLine_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

df <- as.data.frame(rna@meta.data)

cols <- rev(c("#662377","#A82973","#EF5064","#FC867D"))

df$Sample <- factor(df$orig.ident,levels=rev(c("HCC1143",
                                               "SUM149PT",
                                               "MCF7",
                                               "T47D")))

qcRNAnCountBySample <- ggplot(df, aes(y=Sample,
                                      x=log10(nCount_RNA),
                                      fill=Sample)) + 
  scale_fill_manual(values = cols) + 
  geom_boxplot()+
  theme_bw()+
  NoLegend()+
  ggtitle(paste("POST-QC log10(RNA counts) by sample","\n","n=",nrow(df)," cells"))

qcRNAnFeatureBySample <- ggplot(df, aes(y=Sample,
                                        x=nFeature_RNA,
                                        fill=Sample)) + 
  scale_fill_manual(values = cols) + 
  geom_boxplot()+
  theme_bw()+
  ggtitle(paste("POST-QC RNA features by sample","\n","n=",nrow(df)," cells"))

# Read in scATAC
atac <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")

df <- as.data.frame(atac@cellColData)

cols <- rev(c("#662377","#A82973","#EF5064","#FC867D"))

df$Sample <- factor(df$Sample,levels=rev(c("HCC1143",
                                               "SUM149PT",
                                               "MCF7",
                                               "T47D")))

qcATACTSSBySample <- ggplot(df, aes(y=Sample,
                                    x=TSSEnrichment,
                                    fill=Sample)) + 
  scale_fill_manual(values = cols) + 
  geom_boxplot()+
  theme_bw()+
  ggtitle(paste("POST-QC TSS enrichment by sample","\n","n=",nrow(df)," cells"))

qcATACnFragsBySample <- ggplot(df, aes(y=Sample,
                                       x=log10(nFrags),
                                       fill=Sample)) + 
  scale_fill_manual(values = cols) + 
  geom_boxplot()+
  theme_bw()+
  NoLegend()+
  ggtitle(paste("POST-QC log10(nFrags) by sample","\n","n=",nrow(df)," cells"))

qcRNAnCountBySample+qcRNAnFeatureBySample+qcATACnFragsBySample+qcATACTSSBySample+plot_layout(ncol=2)
ggsave("./Cell_Line_Cohort_Results-Revised/QC_Boxplots_RNA_ATAC.pdf",width = 9,height = 8)

# Find P2Gs 
proj <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")

peakSet <- getPeakSet(proj)
names(peakSet) <- NULL
peakSet <- as.data.frame(peakSet)
peakSet$peakName <- paste0(peakSet$seqnames,":",peakSet$start,"-",peakSet$end)

# Unbalanced
meta <- fread("LME_CellLines_out-SingFits_OLS/interactionLMM_univariateLMM_and_OLS_results-metacells.csv")
meta <- as.data.frame(meta)
colnames(meta) <- meta[2,]
meta <- meta[-2,]
names <- grep("Name",colnames(meta))
singFits <- grep("singular_fit",colnames(meta))
names <- unique(c(names,singFits))
meta[,-names] <- sapply(meta[,-names], as.numeric)
meta <- meta[complete.cases(meta),]
meta$P2G <- paste0(meta$LME_cell_type_1_model_peakName,"|",meta$LME_cell_type_1_model_geneName)

meta.int <- meta
meta.int$peakName <- meta.int$LME_cell_type_2_model_peakName
meta.int <- merge(meta.int,peakSet,by="peakName")

meta.int$Basal_FDR <- p.adjust(meta.int$LME_cell_type_1_model_peak_pval,method = "fdr")
meta.int$Luminal_FDR <- p.adjust(meta.int$LME_cell_type_2_model_peak_pval,method = "fdr")

p2g <- meta.int

classes = rep(0, nrow(p2g))

pvA = (p2g$Basal_FDR < 1e-04)
pvB = (p2g$Luminal_FDR < 1e-04)
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
                                        "-/0"),"Basal",p2g$Both)
p2g$Both <- ifelse(p2g$DiffClass %in% c("0/+",
                                        "0/-"),"Luminal",p2g$Both)
p2g$Both <- ifelse(p2g$DiffClass %in% c("0/0"),"Neither",p2g$Both)

table(p2g$DiffClass)
table(p2g$Both)

saveRDS(p2g,"./Cell_Line_Cohort_Results-Revised/p2gs_univariate-prefilter.rds")

p2g <- p2g[p2g$DiffClass != "0/0",]

saveRDS(p2g,"./Cell_Line_Cohort_Results-Revised/p2gs_univariate-postfilter.rds")

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

# Set factor level order for Both
p2g$Both <- factor(p2g$Both,levels=c("Basal","Luminal","Shared"))

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
  ggtitle(paste0("Number of significant P2Gs in each condition\nBasal: ",table(p2g$Both)[1],"  |  Luminal: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

ggsave("./Cell_Line_Cohort_Results-Revised/peakType_proportionBarchart-univariate.pdf",width = 4, height = 4)

peakTypePropUni <- df %>%
  ggplot(aes(fill=peakType, y=Pct, x= conditionSignif, label = Pct))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=rev(RColorBrewer::brewer.pal(5,"Set1")[-1]))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nBasal: ",table(p2g$Both)[1],"  |  Luminal: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
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
  ggtitle(paste0("Number of significant P2Gs in each condition\nBasal: ",table(p2g$Both)[1],"  |  Luminal: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

ggsave("./Cell_Line_Cohort_Results-Revised/encode_proportionBarchart-univariate.pdf",width = 4,height = 4)

encodePropUni <- df %>%
  ggplot(aes(fill=ENCODE, y=Pct, x= conditionSignif))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
            position=position_stack(vjust=0.5)) +
  ylab("Fraction of P2Gs")+
  xlab("Condition in which P2G is significant")+
  scale_fill_manual(values=c("gray10","gray65"))+
  ggtitle(paste0("Number of significant P2Gs in each condition\nBasal: ",table(p2g$Both)[1],"  |  Luminal: ",table(p2g$Both)[2],"  |  Shared: ",table(p2g$Both)[3]))+
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"))

pdf("Cell_Line_Cohort_Results-Revised/peakType_and_ENCODE_overlap_proportion_bar_charts.pdf",width=7,height=4)
cowplot::plot_grid(peakTypePropUni,encodePropUni, ncol = 2, nrow = 1)
dev.off()

pdf("Cell_Line_Cohort_Results-Revised/peakType_and_ENCODE_overlap_proportion_bar_charts-NoLegend.pdf",width=5,height=4)
cowplot::plot_grid(peakTypePropUni+NoLegend(),
                   encodePropUni+NoLegend(), ncol = 2, nrow = 1)
dev.off()

# Plot Euler diagrams, histograms, and proportion bar charts to summarize the P2Gs 

makePlots <- function(group1ArchRProj,
                       group2ArchRProj,
                       group1P2Gs,
                       group2P2Gs,
                       group1Color,
                       group2Color,
                       group1ColorAlt,
                       group2ColorAlt,
                       histGenesPerPeakylim,
                       histPeaksPerGeneylim,
                       outPath1,
                       outPath2,
                       outPath3,
                       outPath4,
                       outPath5,
                       outPath6,
                       metaPath1,
                       metaPath2){
  
  # Find Overlap of P2Gs
  g1P2Gs <- group1P2Gs
  g2P2Gs <- group2P2Gs
  
  group1Peaks.gr <- GRanges(g1P2Gs$peakName)
  group2Peaks.gr <- GRanges(g2P2Gs$peakName)
  
  overlaps <- findOverlaps(query = group1Peaks.gr,
                           subject = group2Peaks.gr)
  overlaps.df <- as.data.frame(overlaps)
  
  g1P2Gs$queryHits <- 1:nrow(g1P2Gs)
  g1P2Gs.annotated <- merge(g1P2Gs,overlaps.df,by = "queryHits")
  
  g2P2Gs$subjectHits <- 1:nrow(g2P2Gs)
  
  annotated <- merge(g1P2Gs.annotated,g2P2Gs,by = "subjectHits")
  
  annotated$Match <- ifelse(annotated$LME_cell_type_1_model_geneName.x == annotated$LME_cell_type_1_model_geneName.y,TRUE,FALSE)
  
  annotated.match <- annotated[annotated$Match == TRUE,]
  
  g1P2Gs$P2G_In_Group2 <- ifelse(g1P2Gs$P2G %in% annotated.match$P2G.x,TRUE,FALSE)
  
  g2P2Gs$P2G_In_Group1 <- ifelse(g2P2Gs$P2G %in% annotated.match$P2G.y,TRUE,FALSE)
  
  g1ProportionP2G <- round(table(g1P2Gs$P2G_In_Group2)[["TRUE"]]/nrow(g1P2Gs),2)
  
  g2ProportionP2G <- round(table(g2P2Gs$P2G_In_Group1)[["TRUE"]]/nrow(g2P2Gs),2)
  
  fit3 <- euler(c("group1" = table(g1P2Gs$P2G_In_Group2)[["FALSE"]], 
                  "group2" = table(g2P2Gs$P2G_In_Group1)[["FALSE"]],
                  "group1&group2" = nrow(annotated.match)))
  
  # Find overlap of enhancers
  g1P2Gs <- group1P2Gs
  g2P2Gs <- group2P2Gs
  
  group1Peaks.gr <- GRanges(unique(g1P2Gs$peakName))
  group2Peaks.gr <- GRanges(unique(g2P2Gs$peakName))
  
  g1ProportionPeak <- round(length(subsetByOverlaps(group1Peaks.gr,
                                            group2Peaks.gr))/length(group1Peaks.gr),2)
  
  g2ProportionPeak <- round(length(subsetByOverlaps(group2Peaks.gr,
                                                    group1Peaks.gr))/length(group2Peaks.gr),2)
  
  fit2 <- euler(c("group1" = length(subsetByOverlaps(group1Peaks.gr,
                                                         group2Peaks.gr,invert = TRUE)),
                  "group2" = length(subsetByOverlaps(group2Peaks.gr,
                                                       group1Peaks.gr,invert = TRUE)),
                  "group1&group2" = length(findOverlaps(query = group1Peaks.gr,
                                                             subject = group2Peaks.gr))))
  
  # Find overlap of genes
  g1P2Gs <- group1P2Gs
  g2P2Gs <- group2P2Gs
  
  group1Genes <- unique(g1P2Gs$LME_cell_type_1_model_geneName)
  group2Genes <- unique(g2P2Gs$LME_cell_type_1_model_geneName)
  
  g1ProportionGene <- round(length(intersect(group1Genes,group2Genes))/length(group1Genes),2)
  
  g2ProportionGene <- round(length(intersect(group2Genes,group1Genes))/length(group2Genes),2)
  
  fit1 <- euler(c("group1" = length(setdiff(group1Genes,
                                            group2Genes)),
                  "group2" = length(setdiff(group2Genes,
                                            group1Genes)),
                  "group1&group2" = length(intersect(group1Genes,
                                                     group2Genes))
                  ))
  
  # Write output
  pdf(outPath1,width=24,height = 10)
  gridExtra::grid.arrange(plot(fit1,
                               fills = c(group1Color,group2Color),
                               edges = FALSE,
                               quantities=TRUE,
                               main=paste0("Overlap of putative enhancer-regulated genes\n(",g1ProportionGene*100,"% of group1 | ",g2ProportionGene*100,"% of group2)")),
                          plot(fit2,
                               fills = c(group1Color,group2Color),
                               edges = FALSE,
                               quantities=TRUE,
                               main=paste0("Overlap of putative enhancers\n(",g1ProportionPeak*100,"% of group1 | ",g2ProportionPeak*100,"% of group2)")),
                          plot(fit3,
                               fills = c(group1Color,group2Color),
                               edges = FALSE,
                               quantities=TRUE,
                               main=paste0("Overlap of putative enhancer-target gene pairs\n(",g1ProportionP2G*100,"% of group1 | ",g2ProportionP2G*100,"% of group2)")),
                          nrow=1)
  dev.off()
  
  # Find overlap of raw peaksets
  
  group1Peaks.gr <- getPeakSet(group1ArchRProj)
  group2Peaks.gr <- getPeakSet(group2ArchRProj)
  
  g1ProportionPeak <- round(length(subsetByOverlaps(group1Peaks.gr,
                                                    group2Peaks.gr))/length(group1Peaks.gr),2)
  
  g2ProportionPeak <- round(length(subsetByOverlaps(group2Peaks.gr,
                                                    group1Peaks.gr))/length(group2Peaks.gr),2)
  
  fit <- euler(c("group1" = length(subsetByOverlaps(group1Peaks.gr,
                                                     group2Peaks.gr,invert = TRUE)),
                  "group2" = length(subsetByOverlaps(group2Peaks.gr,
                                                     group1Peaks.gr,invert = TRUE)),
                  "group1&group2" = length(findOverlaps(query = group1Peaks.gr,
                                                        subject = group2Peaks.gr))))

  # Write output
  pdf(outPath2,width=8,height = 10)
  gridExtra::grid.arrange(plot(fit,
                          fills = c(group1Color,group2Color),
                          edges = FALSE,
                          quantities=TRUE,
                          main=paste0("Overlap of raw peaks\n(",g1ProportionPeak*100,"% of group1 | ",g2ProportionPeak*100,"% of group2)")))
  dev.off()
  
  # Perform gene set ORA with intersection of genes 
  
  # Find overlap of genes
  g1P2Gs <- group1P2Gs
  g2P2Gs <- group2P2Gs
  
  group1Genes <- unique(g1P2Gs$LME_cell_type_1_model_geneName)
  group2Genes <- unique(g2P2Gs$LME_cell_type_1_model_geneName)
  
  # Gene set ORA 
  gset = msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::distinct(gs_name, gene_symbol)
  
  em.gene <- clusterProfiler::enricher(gene=intersect(group1Genes,group2Genes), 
                                       TERM2GENE=gset,
                                       pvalueCutoff=1,
                                       qvalueCutoff=1)
  
  saveRDS(em.gene,outPath3)
  
  # Check eFDR of p.adjust < 0.01 (see alternate bootstrapping method at the
  # bottom of script that samples the same number of enhancer-regulated genes
  # with replacement from the original starting gene list, then assesses gene
  # set enrichment with overlap of sampled genes)
  
  # Find union of genes that went into the P2G analysis
  meta_cell_data <- readRDS(metaPath1)
  
  RNA <- assay(meta_cell_data[[1]])
  rownames(RNA) <- meta_cell_data[[1]]@rowRanges$name
  colnames(RNA) <- paste0("metacell_",1:ncol(RNA))
  
  RNA.1 <- RNA
  
  meta_cell_data <- readRDS(metaPath2)
  
  RNA <- assay(meta_cell_data[[1]])
  rownames(RNA) <- meta_cell_data[[1]]@rowRanges$name
  colnames(RNA) <- paste0("metacell_",1:ncol(RNA))
  
  RNA.2 <- RNA
  
  null <- c(0)
  for(i in 1:1000){
    em.null.gene <- clusterProfiler::enricher(gene=sample(x=unique(c(rownames(RNA.1),rownames(RNA.2))),
                                                          size=length(intersect(group1Genes,group2Genes)),
                                                          replace = FALSE),
                                         TERM2GENE=gset,
                                         pvalueCutoff=1,
                                         qvalueCutoff=1)

    null <- c(null,nrow(em.null.gene@result[em.null.gene@result$p.adjust < 0.01,]))

    print(paste0(i," iterations done!"))
  }
  null <- null[-1]
  
  eFDR <- median(null)/nrow(em.gene@result[em.gene@result$p.adjust < 0.01,])

  # Plot gene set ORA dot plot
  enrichplot::dotplot(em.gene,showCategory=3)+
    ggtitle(paste0(length(intersect(group1Genes,group2Genes))," shared enhancer-regulated genes\n(eFDR ",eFDR," at p.adjust < 0.01)"))
  ggsave(outPath4, 
         width = 8,
         height = 6)
  
  # Plot version of ORA dot plot without y labels
  em.gene.edit <- em.gene
  em.gene.edit@result$ID <- as.character(seq(1:length(em.gene.edit@result$ID)))
  em.gene.edit@result$Description <- as.character(seq(1:length(em.gene.edit@result$Description)))
  
  # Plot gene set ORA dot plot
  gsEnrich <- enrichplot::dotplot(em.gene.edit,showCategory=3)+
    ylab("HALLMARK Gene Sets")+
    theme(plot.title = element_text(size = 8, face = "bold"))+
    ggtitle(paste0(length(intersect(group1Genes,group2Genes))," shared enhancer-regulated genes\n(eFDR ",eFDR," at p.adjust < 0.01)"))
  
  # Generate histograms and barcharts for genes per peak
  g1P2Gs <- group1P2Gs
  g2P2Gs <- group2P2Gs
  
  df.group1 <- data.frame(num.genes = table(g1P2Gs$peakName))
  df.group1$cat <- ifelse(df.group1$num.genes.Freq < 3,"1-2",df.group1$num.genes.Freq)
  df.group1$cat <- ifelse(df.group1$num.genes.Freq >2,"3+",df.group1$cat)
  
  df.group2 <- data.frame(num.genes = table(g2P2Gs$peakName))
  df.group2$cat <- ifelse(df.group2$num.genes.Freq < 3,"1-2",df.group2$num.genes.Freq)
  df.group2$cat <- ifelse(df.group2$num.genes.Freq >2,"3+",df.group2$cat)
  
  head(df.group2)
  head(df.group1)
  
  df.group2$type <- "group2"
  df.group1$type <- "group1"
  
  df<- rbind(df.group2,df.group1)
  
  df$type <- factor(df$type,levels = c("group2","group1"))
  
  genesPerPeakHist <- ggplot(df,aes(x=num.genes.Freq,fill=type))+geom_histogram(size=0.3,boundary = 0,bins=25,alpha=0.45, position="identity")+
    theme_bw()+
    scale_fill_manual(values=c(group2ColorAlt,group1ColorAlt))+
    scale_color_manual(values=c(group2ColorAlt,group1ColorAlt))+
    scale_y_continuous(limits=c(0,histGenesPerPeakylim)) +
    scale_x_continuous(limits = c(0,(max(df$num.genes.Freq)+1)))+
    geom_vline(data=ddply(df, "type", summarise, mean=mean(num.genes.Freq)), aes(xintercept=mean, color=type),
               linetype="dashed")+
    ggtitle(paste0("Mean # of linked genes per peak: ",
                   round(mean(df.group1$num.genes.Freq),2),
                   ", ",
                   round(mean(df.group2$num.genes.Freq),2),
                   "\nMedian # of linked genes per peak: ",
                   round(median(df.group1$num.genes.Freq),2),
                   ", ",
                   round(median(df.group2$num.genes.Freq),2),
                   "\nN: ",
                   nrow(df.group1),
                   " peaks, ",
                   nrow(df.group2),
                   " peaks"))+
    theme(plot.title = element_text(size = 8, face = "bold"))
  
  # Build table for fisher exact test
  df.group1 %>% 
    count(cat) %>% 
    mutate(perc = (n / nrow(df.group1)*100)) -> group1
  df.group2 %>% 
    count(cat) %>% 
    mutate(perc = (n / nrow(df.group2)*100)) -> group2
  
  group1$type <- "group1"
  group2$type <- "group2"
  
  dat <- data.frame(
    "3plus_no" = c(group1[group1$cat != "3+",]$n,group2[group2$cat != "3+",]$n),
    "3plus_yes" = c(group1[group1$cat == "3+",]$n,group2[group2$cat == "3+",]$n),
    row.names = c("group1", "group2"),
    stringsAsFactors = FALSE
  )
  colnames(dat) <- c("not_3+", "3+")
  fisher <- fisher.test(dat)
  
  df$type <- factor(df$type,levels = c("group1","group2"))
  df$cat <- ifelse(df$num.genes.Freq >= 3,"3+","1-2")
  df.new <- df %>% dplyr::group_by_at("type") %>% dplyr::count(cat) %>% mutate(pct= prop.table(n) * 100)
  colnames(df.new) <- c("type","cat","num","Pct")
  df.new$cat <- factor(df.new$cat,levels=c("3+","1-2"))
  
  genesPerPeakBar <- df.new %>%
    ggplot(aes(fill=cat, y=Pct, x= type, label = Pct))+
    geom_bar(stat="identity")+
    geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
              position=position_stack(vjust=0.5)) +
    ylab("% of putative enhancers")+
    xlab("type")+
    scale_fill_manual(values=rev(c("#5B84C4","#FB9B50")))+
    theme_classic()+
    theme(plot.title = element_text(size = 8, face = "bold"))+
    ggtitle(paste0("Prop. of enh per # linked genes\nFisher's exact test:\nOR=",
                   round(as.numeric(fisher$estimate),2),
                   ", p-value ",ifelse(fisher$p.value < 0.01,"< 0.01","n.s.")))
  
  # Generate histograms and barcharts for peaks per gene
  g1P2Gs <- group1P2Gs
  g2P2Gs <- group2P2Gs
  
  df.group1 <- data.frame(num.peaks = table(g1P2Gs$LME_cell_type_1_model_geneName))
  df.group1$cat <- ifelse(df.group1$num.peaks.Freq < 3,"1-2",df.group1$num.peaks.Freq)
  df.group1$cat <- ifelse(df.group1$num.peaks.Freq >2,"3+",df.group1$cat)
  
  df.group2 <- data.frame(num.peaks = table(g2P2Gs$LME_cell_type_1_model_geneName))
  df.group2$cat <- ifelse(df.group2$num.peaks.Freq < 3,"1-2",df.group2$num.peaks.Freq)
  df.group2$cat <- ifelse(df.group2$num.peaks.Freq >2,"3+",df.group2$cat)
  
  head(df.group2)
  head(df.group1)
  
  df.group2$type <- "group2"
  df.group1$type <- "group1"
  
  df<- rbind(df.group2,df.group1)
  
  df$type <- factor(df$type,levels = c("group2","group1"))
  
  peaksPerGeneHist <- ggplot(df,aes(x=num.peaks.Freq,fill=type))+geom_histogram(size=0.3,boundary = 0,bins=30,alpha=0.45, position="identity")+
    theme_bw()+
    scale_fill_manual(values=c(group2ColorAlt,group1ColorAlt))+
    scale_color_manual(values=c(group2ColorAlt,group1ColorAlt))+
    scale_y_continuous(limits=c(0,histPeaksPerGeneylim)) +
    scale_x_continuous(limits = c(0,(max(df$num.peaks.Freq)+1)))+
    geom_vline(data=ddply(df, "type", summarise, mean=mean(num.peaks.Freq)), aes(xintercept=mean, color=type),
               linetype="dashed")+
    ggtitle(paste0("Mean # of linked peaks per gene: ",
                   round(mean(df.group1$num.peaks.Freq),2),
                   ", ",
                   round(mean(df.group2$num.peaks.Freq),2),
                   "\nMedian # of linked peaks per gene: ",
                   round(median(df.group1$num.peaks.Freq),2),
                   ", ",
                   round(median(df.group2$num.peaks.Freq),2),
                   "\nN: ",
                   nrow(df.group1),
                   " genes, ",
                   nrow(df.group2),
                   " genes"))+
    theme(plot.title = element_text(size = 8, face = "bold"))
  
  # Build table for fisher exact test
  df.group1 %>% 
    count(cat) %>% 
    mutate(perc = (n / nrow(df.group1)*100)) -> group1
  df.group2 %>% 
    count(cat) %>% 
    mutate(perc = (n / nrow(df.group2)*100)) -> group2
  
  group1$type <- "group1"
  group2$type <- "group2"
  
  dat <- data.frame(
    "3plus_no" = c(group1[group1$cat != "3+",]$n,group2[group2$cat != "3+",]$n),
    "3plus_yes" = c(group1[group1$cat == "3+",]$n,group2[group2$cat == "3+",]$n),
    row.names = c("group1", "group2"),
    stringsAsFactors = FALSE
  )
  colnames(dat) <- c("not_3+", "3+")
  fisher <- fisher.test(dat)
  
  df$type <- factor(df$type,levels = c("group1","group2"))
  df$cat <- ifelse(df$num.peaks.Freq >= 3,"3+","1-2")
  df.new <- df %>% dplyr::group_by_at("type") %>% dplyr::count(cat) %>% mutate(pct= prop.table(n) * 100)
  colnames(df.new) <- c("type","cat","num","Pct")
  df.new$cat <- factor(df.new$cat,levels=c("3+","1-2"))
  
  peaksPerGeneBar <- df.new %>%
    ggplot(aes(fill=cat, y=Pct, x= type, label = Pct))+
    geom_bar(stat="identity")+
    geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
              position=position_stack(vjust=0.5)) +
    ylab("% of enhancer-regulated genes")+
    xlab("type")+
    scale_fill_manual(values=rev(c("#8A6ABA","#D1B660")))+
    theme_classic()+
    theme(plot.title = element_text(size = 8, face = "bold"))+
    ggtitle(paste0("Prop. of genes per # linked peaks\nFisher's exact test:\nOR=",
                   round(as.numeric(fisher$estimate),2),
                   ", p-value ",ifelse(fisher$p.value < 0.01,"< 0.01","n.s.")))
  

  cowplot::plot_grid(gsEnrich+NoLegend(),
                     genesPerPeakHist+NoLegend(),
                     genesPerPeakBar+NoLegend(),
                     peaksPerGeneHist+NoLegend(),
                     peaksPerGeneBar+NoLegend(),
                     ncol = 5, nrow = 1,align = "h", rel_widths = c(0.24, 0.26, 0.12, 0.26, 0.12))
  ggsave(outPath5, 
         width = 24,
         height = 7)
  
  cowplot::plot_grid(gsEnrich,
                     genesPerPeakHist,
                     genesPerPeakBar,
                     peaksPerGeneHist,
                     peaksPerGeneBar,
                     ncol = 5, nrow = 1,align = "h", rel_widths = c(0.24, 0.26, 0.12, 0.26, 0.12))
  ggsave(outPath6, 
         width = 24,
         height = 7)
  
}

# Basal patients v. Basal cell lines - different peakset
patientP2Gs <- readRDS("./Basal_Cohort_Results/p2gs_univariate-postfilter.rds")
patientP2Gs <- patientP2Gs[patientP2Gs$DiffClass %in% c("+/+","+/0","+/-") & patientP2Gs$peakType %in% c("Intronic","Distal"),]

clP2Gs <- readRDS("./Cell_Line_Cohort_Results-Revised/p2gs_univariate-postfilter.rds")
clP2Gs <- clP2Gs[clP2Gs$DiffClass %in% c("+/+","+/0","+/-") & clP2Gs$peakType %in% c("Intronic","Distal"),]

proj.patient <- loadArchRProject(path =  "./Basal_TN_Samples_scATAC-TESTING3")

proj.cl <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")

makePlots(group1ArchRProj = proj.cl,
           group2ArchRProj = proj.patient,
           group1P2Gs = clP2Gs,
           group2P2Gs = patientP2Gs,
           group1Color = "gray70",
           group2Color = "gray90",
           group1ColorAlt = "gray0",
           group2ColorAlt = "gray50",
           histGenesPerPeakylim = 43000,
           histPeaksPerGeneylim = 6000,
           outPath1 = "./Cell_Line_Cohort_Results-Revised/Basal_Cancer_Patients-Basal_Cell_Lines-P2G_Overlaps.pdf",
           outPath2 = "./Cell_Line_Cohort_Results-Revised/Basal_Cancer_Patients-Basal_Cell_Lines-PeakSet_Overlaps.pdf",
           outPath3 = "./Cell_Line_Cohort_Results-Revised/Basal_Cancer_Patients-Basal_Cell_Lines-em_gene.rds",
           outPath4 = "./Cell_Line_Cohort_Results-Revised/Basal_Cancer_Patients-Basal_Cell_Lines-GeneSet_ORA_Hallmark_MSigDB-Shared_Genes-DotPlot.pdf",
           outPath5 = "./Cell_Line_Cohort_Results-Revised/Basal_Cancer_Patients-Basal_Cell_Lines-GeneSet_ORA_with_GenesPerPeak_and_PeaksPerGene_Histograms_and_Bar_Charts-NoLegend.pdf",
           outPath6 = "./Cell_Line_Cohort_Results-Revised/Basal_Cancer_Patients-Basal_Cell_Lines-GeneSet_ORA_with_GenesPerPeak_and_PeaksPerGene_Histograms_and_Bar_Charts.pdf",
           metaPath1 = "./LME_CellLines_out-SingFits_OLS/meta_cell_data.rds",
           metaPath2 = "./LME_Basal_out-SingFits_OLS/meta_cell_data.rds")

# Luminal patients v. Luminal cell lines - different peakset
patientP2Gs <- readRDS("./Luminal_Cohort_Results/p2gs_univariate-postfilter.rds")
patientP2Gs <- patientP2Gs[patientP2Gs$DiffClass %in% c("+/+","+/0","+/-") & patientP2Gs$peakType %in% c("Intronic","Distal"),]

clP2Gs <- readRDS("./Cell_Line_Cohort_Results-Revised/p2gs_univariate-postfilter.rds")
clP2Gs <- clP2Gs[clP2Gs$DiffClass %in% c("+/+","0/+","-/+") & clP2Gs$peakType %in% c("Intronic","Distal"),]

proj.patient <- loadArchRProject(path =  "./Luminal_TN_Samples_scATAC-TESTING3")

proj.cl <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")

makePlots(group1ArchRProj = proj.cl,
           group2ArchRProj = proj.patient,
           group1P2Gs = clP2Gs,
           group2P2Gs = patientP2Gs,
           group1Color = "gray70",
           group2Color = "gray90",
           group1ColorAlt = "gray0",
           group2ColorAlt = "gray50",
           histGenesPerPeakylim = 40000,
           histPeaksPerGeneylim = 10000,
           outPath1 = "./Cell_Line_Cohort_Results-Revised/Luminal_Cancer_Patients-Luminal_Cell_Lines-P2G_Overlaps.pdf",
           outPath2 = "./Cell_Line_Cohort_Results-Revised/Luminal_Cancer_Patients-Luminal_Cell_Lines-PeakSet_Overlaps.pdf",
           outPath3 = "./Cell_Line_Cohort_Results-Revised/Luminal_Cancer_Patients-Luminal_Cell_Lines-em_gene.rds",
           outPath4 = "./Cell_Line_Cohort_Results-Revised/Luminal_Cancer_Patients-Luminal_Cell_Lines-GeneSet_ORA_Hallmark_MSigDB-Shared_Genes-DotPlot.pdf",
           outPath5 = "./Cell_Line_Cohort_Results-Revised/Luminal_Cancer_Patients-Luminal_Cell_Lines-GeneSet_ORA_with_GenesPerPeak_and_PeaksPerGene_Histograms_and_Bar_Charts-NoLegend.pdf",
           outPath6 = "./Cell_Line_Cohort_Results-Revised/Luminal_Cancer_Patients-Luminal_Cell_Lines-GeneSet_ORA_with_GenesPerPeak_and_PeaksPerGene_Histograms_and_Bar_Charts.pdf",
           metaPath1 = "./LME_CellLines_out-SingFits_OLS/meta_cell_data.rds",
           metaPath2 = "./LME_Luminal_out-SingFits_OLS/meta_cell_data.rds")


# Brainstorm strategies for computing statistical significance of overlap:
# - Fisher exact test
# - Hypergenometric (phyper)
# - Frequentist/bootstrapping strategy to randomly sample two groups of genes 
# from the pool of 20k N times and compute how many times the overlap is greather 
# than or equal to the one we observe
#
# 
# null <- c(0)
# for(i in 1:100){
#   
#   sampleGene1 <- unique(sample(genePool1,length(g1P2Gs$LME_cell_type_1_model_geneName),replace = TRUE))
#   
#   sampleGene2 <- unique(sample(genePool2,length(g2P2Gs$LME_cell_type_1_model_geneName),replace = TRUE))
#   
#   
#   em.null.gene <- clusterProfiler::enricher(gene=intersect(sampleGene1,sampleGene2),
#                                             TERM2GENE=gset,
#                                             pvalueCutoff=1,
#                                             qvalueCutoff=1)
#   
#   null <- c(null,nrow(em.null.gene@result[em.null.gene@result$p.adjust < 0.01,]))
#   
#   print(paste0(i," iterations done!"))
# }
# null <- null[-1]
# 
# median(null)/nrow(em.gene@result[em.gene@result$p.adjust < 0.01,])
#
# 
# overlap_significance <- function(genes_all, gene_sets, iterations) {
#   observed <- length(reduce(gene_sets, intersect))
#   simulated <- map_dbl(seq_len(iterations), function(x) {
#     sim <- map(lengths(gene_sets), ~sample(genes_all, .x))
#     sim <- length(reduce(sim, intersect))
#     return(sim)
#   })
#   pval <- (sum(simulated >= observed) + 1) / (iterations + 1)
#   return(list(pval=pval, simulated_values=simulated, observed=observed))
# }
# 
# overlap_significance(genes_all = )
# 
# store <- c(0)
# for(i in 1:100){
#   
#   clSub <- sample(genes,length(unique(p2g.cl$LME_cell_type_2_model_geneName)))
#   patientSub <- sample(genes,length( unique(p2g.patient$LME_cell_type_1_model_geneName)))
#   
#   print(length(intersect(clSub,patientSub)))
#   if(  length(intersect(clSub,patientSub)) >= length(intersect(unique(p2g.cl$LME_cell_type_2_model_geneName),
#                                                                unique(p2g.patient$LME_cell_type_1_model_geneName))) ){
#     store <- c(store,1)
#   }else{
#     store <- c(store,0)
#   }
# 
# }
# store <- store[-1]
