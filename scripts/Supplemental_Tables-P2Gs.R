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
library(tidyr)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(ComplexHeatmap)
library(msigdbr)
library(EnhancedVolcano)
library(readxl)
library(data.table)

# Make figure output folder
if(dir.exists("Supplemental_Tables-P2Gs")==TRUE){
  print("Directory Supplemental_Tables-P2Gs already exists!")
}else{
  dir.create("Supplemental_Tables-P2Gs")
}

# Function to format P2Gs and write out to csv
writeP2Gs <- function(p2g,
                      p2gInteraction,
                      encode,
                      renameCellType1,
                      renameCellType2,
                      cellType1_FDR,
                      cellType2_FDR,
                      diffGenes,
                      diffPeaks,
                      csvPath1,
                      csvPath2,
                      csvPath3,
                      csvPath4,
                      barchartPath
                      ){
  
  # Re-annotate peaks participating in P2Gs again for ENCODE overlap as done previously
  # Downloaded ENCODE cCREs from https://screen.encodeproject.org 
  encode.gr <- GRanges(paste0(encode$V1,":",encode$V2,"-",encode$V3))
  
  p2g.gr <- GRanges(unique(p2g$peakName))
  
  overlaps <- subsetByOverlaps(x = p2g.gr,
                               ranges = encode.gr)
  
  encodePct <- length(overlaps) / length(p2g.gr)
  
  p2g$ENCODE_Peak_Overlap <- ifelse(p2g$peakName %in% paste0(overlaps@seqnames,":",overlaps@ranges),"TRUE","FALSE")
  
  # Rename columns
  p2g$geneName <- p2g$LME_cell_type_1_model_geneName
  p2g$peakName_geneName <- paste0(p2g$peakName,"_",p2g$geneName)
  
  p2g$genomic_region_of_peak <- p2g$peakType
  p2g$peak_overlap_with_ENCODE <- p2g$ENCODE_Peak_Overlap
  
  p2g$LME_cell_type_1_model_peak_pval_FDR <- p2g[[cellType1_FDR]]
  p2g$LME_cell_type_2_model_peak_pval_FDR <- p2g[[cellType2_FDR]]
  
  p2g$significant_effect_size_in_cell_type_1 <- ifelse(p2g$LME_cell_type_1_model_peak_pval_FDR < 1e-04,TRUE,FALSE)
  p2g$significant_effect_size_in_cell_type_2 <- ifelse(p2g$LME_cell_type_2_model_peak_pval_FDR < 1e-04,TRUE,FALSE)
  
  if(!is.null(p2gInteraction)){
    
    p2g$directional_change_in_effect_size_between_conditions <- p2g$DiffClass
    
    # Annotate FDR from interaction results
    # Show FDR values for all P2Gs tested in interaction model then indicate which P2Gs reached significance with TRUE/FALSE
    
    p2g.intron.distal <- p2g[p2g$genomic_region_of_peak %in% c("Intronic","Distal"),]
    p2g.intron.distal$FDR <- p.adjust(p2g.intron.distal$LME_interaction_model_peak_cell_type_interaction_pval,method = "fdr")
    
    p2g.intron.distal$LME_interaction_model_peak_cell_type_interaction_pval_FDR <- p2g.intron.distal$FDR
    p2g.intron.distal <- p2g.intron.distal[,c("P2G",
                                              "LME_interaction_model_peak_cell_type_interaction_pval_FDR")]
    
    # Quick check to confirm consistency of interaction FDR values 
    p2g.intron.distal.sig <- p2g.intron.distal[p2g.intron.distal$LME_interaction_model_peak_cell_type_interaction_pval_FDR < 1e-04,]
    
    all.equal(p2g.intron.distal.sig$P2G,p2gInteraction$P2G)
    
    all.equal(p2g.intron.distal.sig$LME_interaction_model_peak_cell_type_interaction_pval_FDR,p2gInteraction$FDR)
    
    # Quick check to confirm consistency of ENCODE annotations
    p2g.test <- p2g[p2g$P2G %in% p2gInteraction$P2G,]
    
    all.equal(p2g.test$P2G,p2gInteraction$P2G)
    
    all.equal(p2g.test$peakName,p2gInteraction$peakName)
    
    all.equal(p2g.test$peak_overlap_with_ENCODE,p2gInteraction$ENCODE_Peak_Overlap)
    
    # Merge results to annotate the interaction FDR
    # Exonic/promoter P2Gs and P2Gs that did not reach significance in the interaction
    # model phase will have NA for the interaction FDR column
    
    p2g.annotated <- merge(p2g,p2g.intron.distal,by="P2G",all.x=TRUE)
    
    p2g.annotated$significant_change_in_effect_size_between_conditions <- ifelse(p2g.annotated$LME_interaction_model_peak_cell_type_interaction_pval_FDR < 1e-04,
                                                                                 TRUE,
                                                                                 FALSE
    )
    
    # Merge diffGenes results to annotate DEGs
    colnames(diffGenes)[-grep("geneName",colnames(diffGenes))] <- paste0("differentially_expressed_gene_",colnames(diffGenes)[-grep("geneName",colnames(diffGenes))])
    p2g.annotated <- merge(p2g.annotated,diffGenes,by="geneName",all.x=TRUE)
    
    p2g.annotated$differentially_expressed_gene <- ifelse(p2g.annotated$geneName %in% diffGenes$geneName, TRUE, FALSE)
    p2g.annotated$upregulated_gene_in_cell_type_1 <- ifelse(p2g.annotated$differentially_expressed_gene_log2FoldChange >= 0.58, TRUE, FALSE)
    p2g.annotated$downregulated_gene_in_cell_type_1 <- ifelse(p2g.annotated$differentially_expressed_gene_log2FoldChange <= -0.58, TRUE, FALSE)
    
    # Merge diffPeaks results to annotate DAPs
    colnames(diffPeaks)[-grep("peakName",colnames(diffPeaks))] <- paste0("differentially_accessible_peak_",colnames(diffPeaks)[-grep("peakName",colnames(diffPeaks))])
    p2g.annotated <- merge(p2g.annotated,diffPeaks,by="peakName",all.x=TRUE)
    
    p2g.annotated$differentially_accessible_peak <- ifelse(p2g.annotated$peakName %in% diffPeaks$peakName, TRUE, FALSE)
    p2g.annotated$upregulated_peak_in_cell_type_1 <- ifelse(p2g.annotated$differentially_accessible_peak_log2FoldChange >= 0.58, TRUE, FALSE)
    p2g.annotated$downregulated_peak_in_cell_type_1 <- ifelse(p2g.annotated$differentially_accessible_peak_log2FoldChange <= -0.58, TRUE, FALSE)
    
    # Subset to relevant columns before writing out 
    p2g.annotated <- p2g.annotated[,c("peakName_geneName",
                                      "peakName",
                                      "genomic_region_of_peak",
                                      "peak_overlap_with_ENCODE",
                                      "geneName",
                                      "significant_effect_size_in_cell_type_1",
                                      "significant_effect_size_in_cell_type_2",
                                      "directional_change_in_effect_size_between_conditions",
                                      "significant_change_in_effect_size_between_conditions",
                                      # Cancer condition
                                      "LME_cell_type_1_model_peak_effect_size",
                                      "LME_cell_type_1_model_peak_se",
                                      "LME_cell_type_1_model_peak_df",
                                      "LME_cell_type_1_model_peak_tvalue",
                                      "LME_cell_type_1_model_peak_pval",
                                      "LME_cell_type_1_model_peak_pval_FDR",
                                      "LME_cell_type_1_model_peak_ANOVA_Chisq",
                                      "LME_cell_type_1_model_peak_ANOVA_pval",
                                      "LME_cell_type_1_model_peak_singular_fit",
                                      "LME_cell_type_1_model_random_effect_variance",
                                      "LME_cell_type_1_model_random_effect_sd",
                                      # Normal condition
                                      "LME_cell_type_2_model_peak_effect_size",
                                      "LME_cell_type_2_model_peak_se",
                                      "LME_cell_type_2_model_peak_df",
                                      "LME_cell_type_2_model_peak_tvalue",
                                      "LME_cell_type_2_model_peak_pval",
                                      "LME_cell_type_2_model_peak_pval_FDR",
                                      "LME_cell_type_2_model_peak_ANOVA_Chisq",
                                      "LME_cell_type_2_model_peak_ANOVA_pval",
                                      "LME_cell_type_2_model_peak_singular_fit",
                                      "LME_cell_type_2_model_random_effect_variance",
                                      "LME_cell_type_2_model_random_effect_sd",
                                      # Interaction model
                                      "LME_interaction_model_peak_effect_size",
                                      "LME_interaction_model_peak_se",
                                      "LME_interaction_model_peak_df",
                                      "LME_interaction_model_peak_tvalue",
                                      "LME_interaction_model_peak_pval",
                                      "LME_interaction_model_peak_ANOVA_Chisq",
                                      "LME_interaction_model_peak_ANOVA_pval",
                                      "LME_interaction_model_peak_cell_type_interaction_effect_size",
                                      "LME_interaction_model_peak_cell_type_interaction_se",
                                      "LME_interaction_model_peak_cell_type_interaction_df",
                                      "LME_interaction_model_peak_cell_type_interaction_tvalue",
                                      "LME_interaction_model_peak_cell_type_interaction_pval",
                                      "LME_interaction_model_peak_cell_type_interaction_pval_FDR",
                                      "LME_interaction_model_peak_cell_type_interaction_ANOVA_Chisq",
                                      "LME_interaction_model_peak_cell_type_interaction_ANOVA_pval",
                                      "LME_interaction_model_peak_cell_type_interaction_singular_fit",
                                      "LME_interaction_model_random_effect_variance",
                                      "LME_interaction_model_random_effect_sd",
                                      # Differential expression
                                      "differentially_expressed_gene",                                
                                      "upregulated_gene_in_cell_type_1",
                                      "downregulated_gene_in_cell_type_1",   
                                      "differentially_expressed_gene_baseMean",
                                      "differentially_expressed_gene_log2FoldChange",
                                      "differentially_expressed_gene_lfcSE",                          
                                      "differentially_expressed_gene_stat",
                                      "differentially_expressed_gene_pvalue",                        
                                      "differentially_expressed_gene_padj",
                                      # Differential accessibility
                                      "differentially_accessible_peak",
                                      "upregulated_peak_in_cell_type_1",                              
                                      "downregulated_peak_in_cell_type_1",
                                      "differentially_accessible_peak_baseMean",
                                      "differentially_accessible_peak_log2FoldChange",
                                      "differentially_accessible_peak_lfcSE",
                                      "differentially_accessible_peak_stat",
                                      "differentially_accessible_peak_pvalue",
                                      "differentially_accessible_peak_padj"                          
                         
                                      
    )]
  }else{
    # Subset to relevant columns before writing out 
    p2g.annotated <- p2g[,c("peakName_geneName",
                            "peakName",
                            "genomic_region_of_peak",
                            "peak_overlap_with_ENCODE",
                            "geneName",
                            "significant_effect_size_in_cell_type_1",
                            "significant_effect_size_in_cell_type_2",
                            # Cancer condition
                            "LME_cell_type_1_model_peak_effect_size",
                            "LME_cell_type_1_model_peak_se",
                            "LME_cell_type_1_model_peak_df",
                            "LME_cell_type_1_model_peak_tvalue",
                            "LME_cell_type_1_model_peak_pval",
                            "LME_cell_type_1_model_peak_pval_FDR",
                            "LME_cell_type_1_model_peak_ANOVA_Chisq",
                            "LME_cell_type_1_model_peak_ANOVA_pval",
                            "LME_cell_type_1_model_peak_singular_fit",
                            "LME_cell_type_1_model_random_effect_variance",
                            "LME_cell_type_1_model_random_effect_sd",
                            # Normal condition
                            "LME_cell_type_2_model_peak_effect_size",
                            "LME_cell_type_2_model_peak_se",
                            "LME_cell_type_2_model_peak_df",
                            "LME_cell_type_2_model_peak_tvalue",
                            "LME_cell_type_2_model_peak_pval",
                            "LME_cell_type_2_model_peak_pval_FDR",
                            "LME_cell_type_2_model_peak_ANOVA_Chisq",
                            "LME_cell_type_2_model_peak_ANOVA_pval",
                            "LME_cell_type_2_model_peak_singular_fit",
                            "LME_cell_type_2_model_random_effect_variance",
                            "LME_cell_type_2_model_random_effect_sd"
    )]
  }
  
  # Quick check to confirm consistency of genomic region proportions and encode proportions (bar charts)
  p2g.annotated.check <- p2g.annotated
  
  # genomic_region_of_peak proportion barchart
  p2g.annotated.check$Both <- ifelse(p2g.annotated.check$significant_effect_size_in_cell_type_1 == TRUE,"cellType1","cellType2")
  p2g.annotated.check$Both <- ifelse(p2g.annotated.check$significant_effect_size_in_cell_type_1 == TRUE & 
                                       p2g.annotated.check$significant_effect_size_in_cell_type_2 == TRUE,"Shared",p2g.annotated.check$Both)
  p2g.annotated.check$Both <- factor(p2g.annotated.check$Both,levels=c("cellType2","cellType1","Shared"))
  
  df <- p2g.annotated.check %>% dplyr::group_by_at("Both") %>% dplyr::count(genomic_region_of_peak) %>% mutate(pct= prop.table(n) * 100)
  colnames(df) <- c("conditionSignif","genomic_region_of_peak","Cells","Pct")
  df$genomic_region_of_peak <- factor(df$genomic_region_of_peak,levels=c("Intronic","Distal","Promoter","Exonic"))
  
  genomic_region_of_peak_PropUni <- df %>%
    ggplot(aes(fill=genomic_region_of_peak, y=Pct, x= conditionSignif, label = Pct))+
    geom_bar(stat="identity")+
    geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
              position=position_stack(vjust=0.5)) +
    ylab("Fraction of p2g.annotated.checks")+
    xlab("Condition in which p2g.annotated.check is significant")+
    scale_fill_manual(values=rev(RColorBrewer::brewer.pal(5,"Set1")[-1]))+
    ggtitle(paste0("Number of significant p2g.annotated.checks in each condition\ncellType2: ",table(p2g.annotated.check$Both)[1],"  |  cellType1: ",table(p2g.annotated.check$Both)[2],"  |  Shared: ",table(p2g.annotated.check$Both)[3]))+
    theme_classic()+
    theme(plot.title = element_text(size = 8, face = "bold"))
  
  # ENCDOE proportion barchart
  df <- p2g.annotated.check %>% dplyr::group_by_at("Both") %>% dplyr::count(peak_overlap_with_ENCODE) %>% mutate(pct= prop.table(n) * 100)
  colnames(df) <- c("conditionSignif","ENCODE","Cells","Pct")
  df$genomic_region_of_peak <- factor(df$ENCODE)
  
  encodePropUni <- df %>%
    ggplot(aes(fill=ENCODE, y=Pct, x= conditionSignif))+
    geom_bar(stat="identity")+
    geom_text(aes(label=paste0(sprintf("%1.1f", Pct),"%")),
              position=position_stack(vjust=0.5)) +
    ylab("Fraction of p2g.annotated.checks")+
    xlab("Condition in which p2g.annotated.check is significant")+
    scale_fill_manual(values=c("gray10","gray65"))+
    ggtitle(paste0("Number of significant p2g.annotated.checks in each condition\ncellType2: ",table(p2g.annotated.check$Both)[1],"  |  cellType1: ",table(p2g.annotated.check$Both)[2],"  |  Shared: ",table(p2g.annotated.check$Both)[3]))+
    theme_classic()+
    theme(plot.title = element_text(size = 8, face = "bold"))
  
  cowplot::plot_grid(genomic_region_of_peak_PropUni,encodePropUni, ncol = 2, nrow = 1)
  ggsave(barchartPath,width=7,height=4)

  # Rename columns storing statistics for the cancer (cell type 1) and normal (cell type 2) LME models 
  colnames(p2g.annotated) <- str_replace(string=colnames(p2g.annotated),
                                         pattern = "cell_type_1",
                                         replacement = renameCellType1)
  colnames(p2g.annotated) <- str_replace(string=colnames(p2g.annotated),
                                         pattern = "cell_type_2",
                                         replacement = renameCellType2)
  
  # Convert NAs to "NA"s
  p2g.annotated[is.na(p2g.annotated)] <- "NA"
  
  # Write P2Gs out to csv file
  fwrite(p2g.annotated,csvPath1,sep = ",",col.names = TRUE)
  
  if(!is.null(p2gInteraction)){
    
    # Write out P2G quantities for univariate by normal
    p2g.annotated$Both <- ifelse(p2g.annotated$significant_effect_size_in_cancer_condition == TRUE &
                                   p2g.annotated$significant_effect_size_in_normal_condition == TRUE, "Shared","fill")
    p2g.annotated$Both <- ifelse(p2g.annotated$significant_effect_size_in_cancer_condition == TRUE &
                                   p2g.annotated$significant_effect_size_in_normal_condition == FALSE, "Cancer",p2g.annotated$Both)
    p2g.annotated$Both <- ifelse(p2g.annotated$significant_effect_size_in_cancer_condition == FALSE &
                                   p2g.annotated$significant_effect_size_in_normal_condition == TRUE, "Normal",p2g.annotated$Both)
    
    fwrite(as.data.frame(table(p2g.annotated$Both)),csvPath2,sep = ",",col.names = TRUE)
    
    # Write out P2G quantities for interaction by normal, cancer, shared
    p2g.annotated.sub <- p2g.annotated[p2g.annotated$significant_change_in_effect_size_between_conditions == TRUE,]
    
    fwrite(as.data.frame(table(p2g.annotated.sub$Both)),csvPath3,sep = ",",col.names = TRUE)
    
    # Write out P2G quantities for interaction by DiffClass
    fwrite(as.data.frame(table(p2g.annotated.sub$directional_change_in_effect_size_between_conditions)),csvPath4,sep = ",",col.names = TRUE)
    
  }else{
    
    # Write out P2G quantities for univariate by luminal
    p2g.annotated$Both <- ifelse(p2g.annotated$significant_effect_size_in_basal_condition == TRUE &
                                   p2g.annotated$significant_effect_size_in_luminal_condition == TRUE, "Shared","fill")
    p2g.annotated$Both <- ifelse(p2g.annotated$significant_effect_size_in_basal_condition == TRUE &
                                   p2g.annotated$significant_effect_size_in_luminal_condition == FALSE, "Basal",p2g.annotated$Both)
    p2g.annotated$Both <- ifelse(p2g.annotated$significant_effect_size_in_basal_condition == FALSE &
                                   p2g.annotated$significant_effect_size_in_luminal_condition == TRUE, "Luminal",p2g.annotated$Both)
    
    fwrite(as.data.frame(table(p2g.annotated$Both)),csvPath2,sep = ",",col.names = TRUE)
  
  }
    
}

# Downloaded ENCODE cCREs from https://screen.encodeproject.org 
encode <- read.delim("./miscellaneous/GRCh38-cCREs.bed",header = F)

# Read in Basal P2Gs
p2g <- readRDS("./Basal_Cohort_Results/p2gs_univariate-postfilter.rds")
p2gInteraction <-  readRDS("./Basal_Cohort_Results/p2gs_interaction.rds")
diffGenes <- as.data.frame(readRDS("./Basal_Cohort_Results/resSig_genes.rds"))
diffGenes$geneName <- rownames(diffGenes)
diffPeaks <- as.data.frame(readRDS("./Basal_Cohort_Results/resSig_peaks.rds"))
diffPeaks$peakName <- rownames(diffPeaks)
writeP2Gs(p2g=p2g,
          p2gInteraction = p2gInteraction,
          encode=encode,
          renameCellType1 = "cancer_condition",
          renameCellType2 = "normal_condition",
          cellType1_FDR = "Cancer_FDR",
          cellType2_FDR = "Normal_FDR",
          diffGenes = diffGenes,
          diffPeaks = diffPeaks,
          csvPath1 = "./Supplemental_Tables-P2Gs/Basal_Cohort_P2Gs.csv",
          csvPath2 = "./Supplemental_Tables-P2Gs/Basal_Cohort_P2G_Quantities-Univariate_by_Normal_Cancer_Shared.csv",
          csvPath3 = "./Supplemental_Tables-P2Gs/Basal_Cohort_P2G_Quantities-Interaction_by_Normal_Cancer_Shared.csv",
          csvPath4 = "./Supplemental_Tables-P2Gs/Basal_Cohort_P2G_Quantities-Interaction_by_DiffClass.csv",
          barchartPath="./Supplemental_Tables-P2Gs/Basal_Cohort_Barchart_Check.pdf")

# Read in Luminal P2Gs
p2g <- readRDS("./Luminal_Cohort_Results/p2gs_univariate-postfilter.rds")
p2gInteraction <-  readRDS("./Luminal_Cohort_Results/p2gs_interaction.rds")
diffGenes <- as.data.frame(readRDS("./Luminal_Cohort_Results/resSig_genes.rds"))
diffGenes$geneName <- rownames(diffGenes)
diffPeaks <- as.data.frame(readRDS("./Luminal_Cohort_Results/resSig_peaks.rds"))
diffPeaks$peakName <- rownames(diffPeaks)
writeP2Gs(p2g=p2g,
          p2gInteraction = p2gInteraction,
          encode=encode,
          renameCellType1 = "cancer_condition",
          renameCellType2 = "normal_condition",
          cellType1_FDR = "Cancer_FDR",
          cellType2_FDR = "Normal_FDR",
          diffGenes = diffGenes,
          diffPeaks = diffPeaks,
          csvPath1 = "./Supplemental_Tables-P2Gs/Luminal_Cohort_P2Gs.csv",
          csvPath2 = "./Supplemental_Tables-P2Gs/Luminal_Cohort_P2G_Quantities-Univariate_by_Normal_Cancer_Shared.csv",
          csvPath3 = "./Supplemental_Tables-P2Gs/Luminal_Cohort_P2G_Quantities-Interaction_by_Normal_Cancer_Shared.csv",
          csvPath4 = "./Supplemental_Tables-P2Gs/Luminal_Cohort_P2G_Quantities-Interaction_by_DiffClass.csv",
          barchartPath="./Supplemental_Tables-P2Gs/Luminal_Cohort_Barchart_Check.pdf")

# Read in cell line P2Gs
p2g <- readRDS("./Cell_Line_Cohort_Results/p2gs_univariate-postfilter.rds")
writeP2Gs(p2g=p2g,
          p2gInteraction = NULL,
          encode=encode,
          renameCellType1 = "basal_condition",
          renameCellType2 = "luminal_condition",
          cellType1_FDR = "Basal_FDR",
          cellType2_FDR = "Luminal_FDR",
          diffGenes = NULL,
          diffPeaks = NULL,
          csvPath1 = "./Supplemental_Tables-P2Gs/Cell_Line_Cohort_P2Gs.csv",
          csvPath2 = "./Supplemental_Tables-P2Gs/Cell_Line_Cohort_P2G_Quantities-Univariate_by_Normal_Cancer_Shared.csv",
          csvPath3 = "./Supplemental_Tables-P2Gs/Cell_Line_Cohort_P2G_Quantities-Interaction_by_Normal_Cancer_Shared.csv",
          csvPath4 = "./Supplemental_Tables-P2Gs/Cell_Line_Cohort_P2G_Quantities-Interaction_by_DiffClass.csv",
          barchartPath="./Supplemental_Tables-P2Gs/Cell_Line_Cohort_Barchart_Check.pdf")


print("csv files and barchart checks successfully generated!")
