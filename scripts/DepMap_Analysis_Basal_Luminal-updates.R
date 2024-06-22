################################################################################
# Matt Regner
# Franco Lab
# Reorder boxplots
# Fill color of boxplots (use same default colors in corresponding cell line plots)
# Remove _BREAST from cell line names
# Rename axis labels
# Do global DepMap analysis
# Do stats for boxplots (kruskal wallis or ANOVA) 
################################################################################
library(depmap)
library(ExperimentHub)
library(dplyr)
library(cowplot)
library(Seurat)
library(patchwork)
library(BSDA)
library(data.table)
library(ggplot2)
library(ArchR)
# Make output folder
if(dir.exists("DepMap_Results-updates")){
  print("Directory DepMap_Results-updates already exists!")
}else{
  dir.create("DepMap_Results-updates")
}

# Read in subtype info for cell lines
subtypes <- read.delim("./miscellaneous/Subtype_Call_Cell_Lines_for_Matt_Regner_03152024.txt")
subtypes$cell_line <- paste0(subtypes$cell_line,"_BREAST")

## create ExperimentHub query object
eh <- ExperimentHub()

rnai <- eh[["EH3080"]]
# mutationCalls <- eh[["EH3085"]]
# metadata <- eh[["EH3086"]]
tpm<- eh[["EH3084"]]
# copyNumber <- eh[["EH3082"]]
crispr <- eh[["EH3081"]]
# drug_sensitivity <- eh[["EH3087"]]

# Function to plot DepMap and expr for gene of interest
plotGeneDepMap <- function(eh,
                           subtypes,
                           gene,
                           rnai,
                           crispr,
                           tpm,
                           rnaiOutPath,
                           crisprOutPath,
                           crisprOutPath2,
                           crisprOutPath3){
  
  # CCLE expression by subtype and by cell line (ranked)
  # tpm.sub <- tpm[tpm$cell_line %in% subtypes$cell_line,]
  tpm.sub <- dplyr::filter(tpm,cell_line %in% subtypes$cell_line)
  subtyped.tpm <- merge(tpm.sub,subtypes,by="cell_line")
  subtyped.tpm$cell_line <- gsub("\\_.*","",subtyped.tpm$cell_line)
  
  kw <- kruskal.test(expression ~ SUBTYPE.SUGGESTED.OVERVIEW,
                     data=subtyped.tpm[subtyped.tpm$gene_name == gene,] )
  
  exprSubtype <- ggplot(subtyped.tpm[subtyped.tpm$gene_name == gene,],
                        aes(reorder(SUBTYPE.SUGGESTED.OVERVIEW,expression),expression,fill=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_boxplot()+
    theme_bw()+
    ggtitle(paste0("Expression of ",gene," by BC subtype\nKW test p-value=",kw$p.val))+
    ylab("TPM")+
    xlab("Subtype")
  
  exprCellLine <- ggplot( subtyped.tpm[subtyped.tpm$gene_name == gene,],
                          aes(expression,
                              reorder(cell_line,expression),
                              color=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_point()+
    theme_bw()+
    ggtitle(paste0("Expression of ",gene," in BC cell lines"))+
    ylab("Cell Line")
  
  print("Made expression plots")
  
  # RNAi by subtype and by cell line (ranked)
  rnai$geneName <- rnai$gene_name
  # rnai.sub <- rnai[rnai$cell_line %in% subtypes$cell_line,]
  rnai.sub <- dplyr::filter(rnai,cell_line %in% subtypes$cell_line)
  subtyped.rnai <- merge(rnai.sub,subtypes,by="cell_line")
  subtyped.rnai$cell_line <- gsub("\\_.*","",subtyped.rnai$cell_line)
  
  kw <- kruskal.test(dependency ~ SUBTYPE.SUGGESTED.OVERVIEW,
                     data=subtyped.rnai[subtyped.rnai$gene_name == gene,] )
  
  rnaiSubtype <- ggplot(subtyped.rnai[subtyped.rnai$gene_name == gene,],
                        aes(reorder(SUBTYPE.SUGGESTED.OVERVIEW,dependency),dependency,fill=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_boxplot()+
    theme_bw()+
    ggtitle(paste0("RNAi of ",gene," by BC subtype\nKW test p-value=",kw$p.val))+
    ylab("Dependency")+
    xlab("Subtype")
  
  rnaiCellLine <- ggplot( subtyped.rnai[subtyped.rnai$gene_name == gene,],
                          aes(x=dependency,
                              y=reorder(cell_line,
                                        dependency),
                              color=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_point()+
    theme_bw()+
    ggtitle(paste0("RNAi of ",gene," in BC cell lines"))+
    ylab("Cell Line")
  
  print("Made RNAi plots")
  
  # CRISPR by subtype and by cell line
  crispr$geneName <- crispr$gene_name
  # crispr.sub <- crispr[crispr$cell_line %in% subtypes$cell_line,]
  crispr.sub <- dplyr::filter(crispr,cell_line %in% subtypes$cell_line)
  subtyped.crispr <- merge(crispr.sub,subtypes,by="cell_line")
  subtyped.crispr$cell_line <- gsub("\\_.*","",subtyped.crispr$cell_line)
  
  kw <- kruskal.test(dependency ~ SUBTYPE.SUGGESTED.OVERVIEW,
                     data=subtyped.crispr[subtyped.crispr$gene_name == gene,] )
  
  crisprSubtype <- ggplot(subtyped.crispr[subtyped.crispr$gene_name == gene,],
                          aes(reorder(SUBTYPE.SUGGESTED.OVERVIEW,dependency),dependency,fill=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_boxplot()+
    theme_bw()+
    ggtitle(paste0("CRISPR of ",gene," by BC subtype\nKW test p-value=",kw$p.val))+
    ylab("Dependency")+
    xlab("Subtype")
  
  crisprCellLine <- ggplot( subtyped.crispr[subtyped.crispr$gene_name == gene,],
                            aes(x=dependency,
                                y=reorder(cell_line,
                                          dependency),
                                color=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_point()+
    theme_bw()+
    ggtitle(paste0("CRISPR of ",gene," in BC cell lines"))+
    ylab("Cell Line")
  
  print("Made CRISPR plots")
  
  rnaiSubtype+rnaiCellLine+exprSubtype+exprCellLine+plot_layout(ncol=2)
  ggsave(rnaiOutPath,width=12,height = 14)
  
  crisprSubtype+crisprCellLine+exprSubtype+exprCellLine+plot_layout(ncol=2)
  ggsave(crisprOutPath,width=12,height = 14)
  
  crisprCellLine+exprSubtype+plot_layout(ncol=1)
  ggsave(crisprOutPath2,width=10,height = 12)
  
  crisprCellLine+NoLegend()+exprSubtype+NoLegend()+plot_layout(ncol=1)
  ggsave(crisprOutPath3,width=8,height = 12)
  
}

# Plot DepMap and Expr for HEY1
plotGeneDepMap(eh=eh,
               subtypes=subtypes,
               gene="HEY1",
               rnai=as.data.frame(rnai),
               crispr=as.data.frame(crispr),
               tpm=as.data.frame(tpm),
               rnaiOutPath="./DepMap_Results-updates/HEY1_rnai_and_tpm.pdf",
               crisprOutPath="./DepMap_Results-updates/HEY1_crispr_and_tpm.pdf",
               crisprOutPath2="./DepMap_Results-updates/HEY1_crispr_and_tpm-vertical.pdf",
               crisprOutPath3="./DepMap_Results-updates/HEY1_crispr_and_tpm-vertical-NoLegend.pdf")

# Plot DepMap and Expr for CRABP2
plotGeneDepMap(eh=eh,
               subtypes=subtypes,
               gene="CRABP2",
               rnai=as.data.frame(rnai),
               crispr=as.data.frame(crispr),
               tpm=as.data.frame(tpm),
               rnaiOutPath="./DepMap_Results-updates/CRABP2_rnai_and_tpm.pdf",
               crisprOutPath="./DepMap_Results-updates/CRABP2_crispr_and_tpm.pdf",
               crisprOutPath2="./DepMap_Results-updates/CRABP2_crispr_and_tpm-vertical.pdf",
               crisprOutPath3="./DepMap_Results-updates/CRABP2_crispr_and_tpm-vertical-NoLegend.pdf")


# Function to DepMap and Expr for HEY1/CRABP2 side by side
plot2GeneDepMap <- function(eh,
                            subtypes,
                            gene1,
                            gene2,
                            crispr,
                            tpm,
                            outPath1,
                            outPath2){
  
  # CRISPR by subtype and by cell line
  crispr$geneName <- crispr$gene_name
  # crispr.sub <- crispr[crispr$cell_line %in% subtypes$cell_line,]
  crispr.sub <- dplyr::filter(crispr,cell_line %in% subtypes$cell_line)
  subtyped.crispr <- merge(crispr.sub,subtypes,by="cell_line")
  subtyped.crispr$cell_line <- gsub("\\_.*","",subtyped.crispr$cell_line)
  
  crisprG1CellLine <- ggplot( subtyped.crispr[subtyped.crispr$gene_name == gene1,],
                              aes(x=dependency,
                                  y=reorder(cell_line,
                                            dependency),
                                  color=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_point(size=4)+
    theme_bw()+
    ggtitle(paste0("CRISPR of ",gene1," in BC cell lines","\nN=",
                   length(unique(subtyped.crispr[subtyped.crispr$gene_name == gene1,]$cell_line))," cell lines"))+
    ylab("Cell Line")
  
  print("Made CRISPR plots gene1")
  
  # CCLE expression by subtype and by cell line (ranked)
  # tpm.sub <- tpm[tpm$cell_line %in% subtypes$cell_line,]
  tpm.sub <- dplyr::filter(tpm,cell_line %in% subtypes$cell_line)
  subtyped.tpm <- merge(tpm.sub,subtypes,by="cell_line")
  subtyped.tpm$cell_line <- gsub("\\_.*","",subtyped.tpm$cell_line)
  subtyped.tpm <- subtyped.tpm[ subtyped.tpm$cell_line %in% unique(subtyped.crispr[subtyped.crispr$gene_name == gene1,]$cell_line),]
  
  kw <- kruskal.test(expression ~ SUBTYPE.SUGGESTED.OVERVIEW,
                     data=subtyped.tpm[subtyped.tpm$gene_name == gene1,] )
  
  exprG1Subtype <- ggplot(subtyped.tpm[subtyped.tpm$gene_name == gene1,],
                          aes(reorder(SUBTYPE.SUGGESTED.OVERVIEW,expression),expression,fill=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_boxplot()+
    theme_bw()+
    ggtitle(paste0("Expression of ",gene1," by BC subtype\nKW test p-value=",
                   round(kw$p.val,3),"\nN=",
                   length(unique(subtyped.tpm[subtyped.tpm$gene_name == gene1,]$cell_line))," cell lines"))+
    ylab("TPM")+
    xlab("Subtype")
  
  print("Made expression plots gene1")
  
  
  # CRISPR by subtype and by cell line
  crispr$geneName <- crispr$gene_name
  # crispr.sub <- crispr[crispr$cell_line %in% subtypes$cell_line,]
  crispr.sub <- dplyr::filter(crispr,cell_line %in% subtypes$cell_line)
  subtyped.crispr <- merge(crispr.sub,subtypes,by="cell_line")
  subtyped.crispr$cell_line <- gsub("\\_.*","",subtyped.crispr$cell_line)
  
  crisprG2CellLine <- ggplot( subtyped.crispr[subtyped.crispr$gene_name == gene2,],
                              aes(x=dependency,
                                  y=reorder(cell_line,
                                            dependency),
                                  color=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_point(size=4)+
    theme_bw()+
    ggtitle(paste0("CRISPR of ",gene2," in BC cell lines","\nN=",
                   length(unique(subtyped.crispr[subtyped.crispr$gene_name == gene2,]$cell_line))," cell lines"))+
    ylab("Cell Line")
  
  print("Made CRISPR plots gene2")
  
  # CCLE expression by subtype and by cell line (ranked)
  # tpm.sub <- tpm[tpm$cell_line %in% subtypes$cell_line,]
  tpm.sub <- dplyr::filter(tpm,cell_line %in% subtypes$cell_line)
  subtyped.tpm <- merge(tpm.sub,subtypes,by="cell_line")
  subtyped.tpm$cell_line <- gsub("\\_.*","",subtyped.tpm$cell_line)
  subtyped.tpm <- subtyped.tpm[ subtyped.tpm$cell_line %in% unique(subtyped.crispr[subtyped.crispr$gene_name == gene2,]$cell_line),]
  
  kw <- kruskal.test(expression ~ SUBTYPE.SUGGESTED.OVERVIEW,
                     data=subtyped.tpm[subtyped.tpm$gene_name == gene2,] )
  
  exprG2Subtype <- ggplot(subtyped.tpm[subtyped.tpm$gene_name == gene2,],
                          aes(reorder(SUBTYPE.SUGGESTED.OVERVIEW,expression),expression,fill=SUBTYPE.SUGGESTED.OVERVIEW))+
    geom_boxplot()+
    theme_bw()+
    ggtitle(paste0("Expression of ",gene2," by BC subtype\nKW test p-value=",
                   round(kw$p.val,3),"\nN=",
                   length(unique(subtyped.tpm[subtyped.tpm$gene_name == gene2,]$cell_line))," cell lines"))+
    ylab("TPM")+
    xlab("Subtype")
  
  print("Made expression plots gene2")
  
  crisprG1CellLine+crisprG2CellLine+exprG1Subtype+exprG2Subtype+plot_layout(ncol=2)
  ggsave(outPath1,width=20,height = 18)
  
  crisprG1CellLine+NoLegend()+crisprG2CellLine+exprG1Subtype+NoLegend()+exprG2Subtype+plot_layout(ncol=2)
  ggsave(outPath2,width=16,height = 18)
  
}

# Plot DepMap and Expr for HEY1/CRABP2 side by side
plot2GeneDepMap(eh=eh,
                subtypes=subtypes,
                gene1="HEY1",
                gene2="CRABP2",
                crispr=as.data.frame(crispr),
                tpm=as.data.frame(tpm),
                outPath1="./DepMap_Results-updates/HEY1_CRABP2_crispr_and_tpm.pdf",
                outPath2="./DepMap_Results-updates/HEY1_CRABP2_crispr_and_tpm-NoLegend.pdf")


# # Functions to plot series of visualizations for patients
# 
# plotFreqDistb <- function(genes,
#                           crispr,
#                           proj,
#                           gene,
#                           subtype,
#                           subtypes,
#                           hjust,
#                           outPath){
#   
#   # CRISPR by subtype and by cell line
#   crispr$geneName <- crispr$gene_name
#   crispr.sub <- dplyr::filter(crispr,cell_line %in% subtypes$cell_line)
#   # crispr.sub <- crispr[crispr$cell_line %in% subtypes$cell_line,]
#   subtyped.crispr <- merge(crispr.sub,subtypes,by="cell_line")
#   # Subset to genes cell lines and calculate mean dependency for each gene across cell lines
#   subtyped.crispr.genes <- subtyped.crispr[subtyped.crispr$SUBTYPE.SUGGESTED.OVERVIEW == subtype,]
#   genes.mean <- subtyped.crispr.genes %>% group_by(geneName) %>%
#     summarize(mean_dependency=mean(dependency))
#   genes.mean <- genes.mean[complete.cases(genes.mean),]
#   
#   genes.mean.sorted <- dplyr::arrange(genes.mean[genes.mean$geneName %in% genes$geneName,],mean_dependency)
#   genes.mean.sorted$rank <- 1:nrow(genes.mean.sorted)
#   
#   ggplot(genes.mean.sorted,aes(x=rank,y=mean_dependency,label=geneName))+
#     geom_point()+
#     geom_text_repel(aes(label=ifelse(rank<6 | geneName == gene,as.character(geneName),'')),hjust=hjust,vjust=0)+
#     theme_bw()+
#     ylab(paste0("Mean dependency score across ", subtype ," cell lines"))+
#     ggtitle(paste0("Freq. distribution plot of ",nrow(genes.mean.sorted), " enh-regulated genes ranked by dependency score"))
#   ggsave(outPath,width=7,height = 7)
#   
#   print("Made freq. distrb. rank plot")
#   
#   
# }
# 
# 
# plotRandBgd <- function(genes,
#                         crispr,
#                         proj,
#                         gene,
#                         subtype,
#                         subtypes,
#                         hjust,
#                         outPath){
#   # CRISPR by subtype and by cell line
#   crispr$geneName <- crispr$gene_name
#   crispr.sub <- dplyr::filter(crispr,cell_line %in% subtypes$cell_line)
#   # crispr.sub <- crispr[crispr$cell_line %in% subtypes$cell_line,]
#   subtyped.crispr <- merge(crispr.sub,subtypes,by="cell_line")
#   # Subset to genes cell lines and calculate mean dependency for each gene across cell lines
#   subtyped.crispr.genes <- subtyped.crispr[subtyped.crispr$SUBTYPE.SUGGESTED.OVERVIEW == subtype,]
#   genes.mean <- subtyped.crispr.genes %>% group_by(geneName) %>%
#     summarize(mean_dependency=mean(dependency))
#   genes.mean <- genes.mean[complete.cases(genes.mean),]
#   
#   # Get ArchR proj genes
#   geneInitial <- getGenes(proj)
#   
#   genes.mean.geneSet <- genes.mean[genes.mean$geneName %in% geneInitial$symbol,]
#   
#   # Generate background distribution of dependency scores for random gene sets
#   store <-c(0)
#   emp <- c()
#   for(i in 1:1000){
#     # Cancer enhancer-regualted genes
#     cancer <- genes.mean.geneSet[genes.mean.geneSet$geneName %in% genes$geneName,]
#     cancer$group <- "cancer_enh_genes"
#     
#     # Random genes
#     random <- genes.mean.geneSet[genes.mean.geneSet$geneName %in% sample(genes.mean.geneSet$geneName,nrow(cancer)),]
#     random$group <- "random_genes"
#     
#     # store <- c(store,res1$p.value)
#     store <- c(store,mean(random$mean_dependency))
#     emp <- c(emp,mean(random$mean_dependency) <= mean(cancer$mean_dependency))
#     
#     print(i)
#   }
#   store <- store[-1]
#   emp <- emp[-1]
#   empPval <- sum(emp)/length(emp)
#   
#   # Do stats
#   # Two sample
#   zres <- z.test(genes.mean.geneSet[genes.mean.geneSet$geneName %in% genes$geneName,]$mean_dependency,
#                  store, mu=0, alternative = "two.sided",
#                  sigma.x = sd(genes.mean.geneSet[genes.mean.geneSet$geneName %in% genes$geneName,]$mean_dependency),
#                  sigma.y = sd(store)
#   )
#   # One sample
#   # zres <- z.test(genes.mean.geneSet[genes.mean.geneSet$geneName %in% genes$geneName,]$mean_dependency,
#   #                mu=mean(store),
#   #                sigma.x=sd(store))
#   
#   store.df <- data.frame(mean_dependency_score=store,fill="fill")
#   
#   ggplot(store.df,aes(mean_dependency_score))+
#     geom_histogram(color="gray20", fill="gray75",bins=60)+
#     theme_bw()+
#     geom_vline(xintercept = mean(genes.mean.geneSet[genes.mean.geneSet$geneName %in% genes$geneName,]$mean_dependency), color="red", linetype = "dashed")+
#     geom_vline(xintercept = mean(store), color="blue", linetype = "dashed")+
#     ggtitle(paste0("Distribution of mean dependency scores for 1,000 randomly sampled gene sets of size ",length(unique(genes.mean.geneSet[genes.mean.geneSet$geneName %in% genes$geneName,]$mean_dependency)),"\nRed: mean dependency score of ",length(unique(genes.mean.geneSet[genes.mean.geneSet$geneName %in% genes$geneName,]$mean_dependency))," Enh-regulated genes (",
#                    round(mean(genes.mean.geneSet[genes.mean.geneSet$geneName %in% genes$geneName,]$mean_dependency),3),")\nBlue: mean dependency score of null distribution (",round(mean(store),3),")",
#                    "\nZ test p-value=",round(zres$p.value,3),
#                    "\nEmpirical p-value=",empPval)
#     )
#   
#   ggsave(outPath,width=8,height = 8)
#   
#   print("Made bootstraping histogram")
#   
# }
# 
# plotCompareBoxPlots <- function(group1,
#                                 group2,
#                                 crispr,
#                                 proj,
#                                 subtype,
#                                 subtypes,
#                                 outPath){
#   
#   # CRISPR by subtype and by cell line
#   crispr$geneName <- crispr$gene_name
#   crispr.sub <- dplyr::filter(crispr,cell_line %in% subtypes$cell_line)
#   # crispr.sub <- crispr[crispr$cell_line %in% subtypes$cell_line,]
#   subtyped.crispr <- merge(crispr.sub,subtypes,by="cell_line")
#   # Subset to genes cell lines and calculate mean dependency for each gene across cell lines
#   subtyped.crispr.genes <- subtyped.crispr[subtyped.crispr$SUBTYPE.SUGGESTED.OVERVIEW == subtype,]
#   genes.mean <- subtyped.crispr.genes %>% group_by(geneName) %>%
#     summarize(mean_dependency=mean(dependency))
#   genes.mean <- genes.mean[complete.cases(genes.mean),]
#   
#   # Get ArchR proj genes
#   geneInitial <- getGenes(proj)
#   
#   genes.mean.geneSet <- genes.mean[genes.mean$geneName %in% geneInitial$symbol,]
#   
#   # Cancer
#   cancer <- genes.mean.geneSet[genes.mean.geneSet$geneName %in% group1$geneName,]
#   cancer$group <- "cancer_enh_genes"
#   
#   # Normal
#   normal <- genes.mean.geneSet[genes.mean.geneSet$geneName %in% group2$geneName,]
#   normal$group <- "normal_enh_genes"
#   
#   # Combine
#   comb <- rbind(cancer,normal)
#   
#   res1 <- wilcox.test(mean_dependency ~ group,data=comb,correct=FALSE)
#   
#   res2 <- wilcox.test(cancer$mean_dependency,normal$mean_dependency,correct=FALSE)
#   
#   if(res1$statistic == res2$statistic){
#     print("wilcox test statistics are equal")
#   }else{
#     print("wilcox test statistics are not equal")
#   }
#   
#   means <- aggregate(mean_dependency ~  group, comb, mean)
#   p <- ggplot(data=comb, aes(x=group, y=mean_dependency, fill=group)) + geom_boxplot() +
#     geom_text(data = means, aes(label = round(mean_dependency,3), y = mean_dependency + 0.5))+theme_bw()+
#     ggtitle(paste0(subtype," subtype patient analysis\nWilcox test p-value: ",round(res1$p.value,3),"\n",nrow(cancer)," enh-regulated genes in cancer\n",nrow(normal)," enh-regulated genes in normal"))
#   
#   ggsave(outPath,width=8,height = 8)
#   
#   print("Made boxplots")
# }
# 
# 
# # Plots for Basal patients
# basal <- readRDS("./Basal_Cohort_Results-updates/p2gs_interaction-+_0-upgenes.rds")
# basal$geneName <- basal$LME_cell_type_1_model_geneName
# 
# proj <- loadArchRProject(path =  "./Basal_TN_Samples_scATAC-TESTING3")
# 
# cancer <- readRDS("./Basal_Cohort_Results-updates/p2gs_univariate-postfilter.rds")
# cancer <- cancer[cancer$DiffClass %in% c("+/+","+/0","+/-") & cancer$peakType %in% c("Intronic","Distal"),]
# cancer$geneName <- cancer$LME_cell_type_1_model_geneName
# 
# normal <- readRDS("./Basal_Cohort_Results-updates/p2gs_univariate-postfilter.rds")
# normal <- normal[normal$DiffClass %in% c("+/+","0/+","-/+") & normal$peakType %in% c("Intronic","Distal"),]
# normal$geneName <- normal$LME_cell_type_2_model_geneName
# 
# plotFreqDistb(genes=basal,
#               crispr=crispr,
#               proj=proj,
#               subtype="Basal",
#               subtypes=subtypes,
#               gene="HEY1",
#               hjust=-0.1,
#               outPath="./DepMap_Results-updates/Freq_Dist_Rank-Basal_Cancer_Specific_Enh_Upregulated_Genes.pdf")
# 
# plotRandBgd(genes=basal,
#             crispr=crispr,
#             proj=proj,
#             subtype="Basal",
#             subtypes=subtypes,
#             outPath="./DepMap_Results-updates/Histogram_of_Rand_Background_Dependency_Scores-Basal_Cancer_Specific_Enh_Upregulated_Genes.pdf")
# 
# plotRandBgd(genes=cancer,
#             crispr=crispr,
#             proj=proj,
#             subtype="Basal",
#             subtypes=subtypes,
#             outPath="./DepMap_Results-updates/Histogram_of_Rand_Background_Dependency_Scores-Basal_Cancer_Enh_Regulated_Genes.pdf")
# 
# plotCompareBoxPlots(group1=cancer,
#                     group2=normal,
#                     crispr=crispr,
#                     proj=proj,
#                     subtype="Basal",
#                     subtypes=subtypes,
#                     outPath="./DepMap_Results-updates/Boxplot_Compare-Basal_Cancer_Enh_Regulated_Genes_Versus_LP_Normal_Enh_Regulated_Genes.pdf")
# 
# 
# # Boxplots for condition-specific and upregulated genes
# degs <- readRDS("Basal_Cohort_Results-updates/DEGs_padj_0.05_log2FC_0.58.rds")
# up.genes <- degs[degs$log2FoldChange > 0.58,]
# down.genes <- degs[degs$log2FoldChange < -0.58,]
# 
# cancer <- readRDS("./Basal_Cohort_Results-updates/p2gs_interaction.rds")
# cancer$geneName <- cancer$LME_cell_type_1_model_geneName
# cancer <- cancer[cancer$geneName %in% up.genes$gene & cancer$DiffClass == "+/0",]
# 
# normal <- readRDS("./Basal_Cohort_Results-updates/p2gs_interaction.rds")
# normal$geneName <- normal$LME_cell_type_1_model_geneName
# normal <- normal[normal$geneName %in% down.genes$gene & normal$DiffClass == "0/+",]
# 
# plotCompareBoxPlots(group1=cancer,
#                     group2=normal,
#                     crispr=crispr,
#                     proj=proj,
#                     subtype="Basal",
#                     subtypes=subtypes,
#                     outPath="./DepMap_Results-updates/Boxplot_Compare-Basal_Cancer_Specific_Enh_Upregulated_Genes_Versus_LP_Normal_Specific_Enh_Upregulated_Genes.pdf")
# 
# 
# # Plots for Luminal patients
# Luminal <- readRDS("./Luminal_Cohort_Results-updates/p2gs_interaction-+_0-upgenes.rds")
# Luminal$geneName <- Luminal$LME_cell_type_1_model_geneName
# 
# proj <- loadArchRProject(path =  "./Luminal_TN_Samples_scATAC-TESTING3")
# 
# cancer <- readRDS("./Luminal_Cohort_Results-updates/p2gs_univariate-postfilter.rds")
# cancer <- cancer[cancer$DiffClass %in% c("+/+","+/0","+/-") & cancer$peakType %in% c("Intronic","Distal"),]
# cancer$geneName <- cancer$LME_cell_type_1_model_geneName
# 
# normal <- readRDS("./Luminal_Cohort_Results-updates/p2gs_univariate-postfilter.rds")
# normal <- normal[normal$DiffClass %in% c("+/+","0/+","-/+") & normal$peakType %in% c("Intronic","Distal"),]
# normal$geneName <- normal$LME_cell_type_2_model_geneName
# 
# plotFreqDistb(genes=Luminal,
#               crispr=crispr,
#               proj=proj,
#               subtype="Luminal",
#               subtypes=subtypes,
#               gene="CRABP2",
#               hjust=-0.1,
#               outPath="./DepMap_Results-updates/Freq_Dist_Rank-Luminal_Cancer_Specific_Enh_Upregulated_Genes.pdf")
# 
# plotRandBgd(genes=Luminal,
#             crispr=crispr,
#             proj=proj,
#             subtype="Luminal",
#             subtypes=subtypes,
#             outPath="./DepMap_Results-updates/Histogram_of_Rand_Background_Dependency_Scores-Luminal_Cancer_Specific_Enh_Upregulated_Genes.pdf")
# 
# plotRandBgd(genes=cancer,
#             crispr=crispr,
#             proj=proj,
#             subtype="Luminal",
#             subtypes=subtypes,
#             outPath="./DepMap_Results-updates/Histogram_of_Rand_Background_Dependency_Scores-Luminal_Cancer_Enh_Regulated_Genes.pdf")
# 
# plotCompareBoxPlots(group1=cancer,
#                     group2=normal,
#                     crispr=crispr,
#                     proj=proj,
#                     subtype="Luminal",
#                     subtypes=subtypes,
#                     outPath="./DepMap_Results-updates/Boxplot_Compare-Luminal_Cancer_Enh_Regulated_Genes_Versus_LP_Normal_Enh_Regulated_Genes.pdf")
# 
# 
# # Boxplots for condition-specific and upregulated genes
# degs <- readRDS("Luminal_Cohort_Results-updates/DEGs_padj_0.05_log2FC_0.58.rds")
# up.genes <- degs[degs$log2FoldChange > 0.58,]
# down.genes <- degs[degs$log2FoldChange < -0.58,]
# 
# cancer <- readRDS("./Luminal_Cohort_Results-updates/p2gs_interaction.rds")
# cancer$geneName <- cancer$LME_cell_type_1_model_geneName
# cancer <- cancer[cancer$geneName %in% up.genes$gene & cancer$DiffClass == "+/0",]
# 
# normal <- readRDS("./Luminal_Cohort_Results-updates/p2gs_interaction.rds")
# normal$geneName <- normal$LME_cell_type_1_model_geneName
# normal <- normal[normal$geneName %in% down.genes$gene & normal$DiffClass == "0/+",]
# 
# plotCompareBoxPlots(group1=cancer,
#                     group2=normal,
#                     crispr=crispr,
#                     proj=proj,
#                     subtype="Luminal",
#                     subtypes=subtypes,
#                     outPath="./DepMap_Results-updates/Boxplot_Compare-Luminal_Cancer_Specific_Enh_Upregulated_Genes_Versus_LP_Normal_Specific_Enh_Upregulated_Genes.pdf")
# 
# 
# # Plots for Basal cell lines
# basal <- fread("./Supplemental_Tables-P2Gs/Cell_Line_Cohort_P2Gs.csv")
# basal <- basal[basal$significant_effect_size_in_basal_condition == TRUE &
#                  basal$LME_basal_condition_model_peak_effect_size > 0 &
#                  basal$genomic_region_of_peak %in% c("Intronic","Distal"),]
# 
# proj <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")
# 
# plotFreqDistb(genes=basal,
#               crispr=crispr,
#               proj=proj,
#               subtype="Basal",
#               subtypes=subtypes,
#               gene="HEY1",
#               hjust=-1,
#               outPath="./DepMap_Results-updates/Freq_Dist_Rank-Basal_Cell_Line_Enh_Regulated_Genes.pdf")
# 
# plotRandBgd(genes=basal,
#             crispr=crispr,
#             proj=proj,
#             subtype="Basal",
#             subtypes=subtypes,
#             outPath="./DepMap_Results-updates/Histogram_of_Rand_Background_Dependency_Scores-Basal_Cell_Line_Enh_Regulated_Genes.pdf")
# 
# 
# # Plots for Luminal cell lines
# luminal <- fread("./Supplemental_Tables-P2Gs/Cell_Line_Cohort_P2Gs.csv")
# luminal <- luminal[luminal$significant_effect_size_in_luminal_condition == TRUE &
#                      luminal$LME_luminal_condition_model_peak_effect_size > 0 &
#                      luminal$genomic_region_of_peak %in% c("Intronic","Distal"),]
# 
# proj <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")
# 
# plotFreqDistb(genes=luminal,
#               crispr=crispr,
#               proj=proj,
#               subtype="Luminal",
#               subtypes=subtypes,
#               gene="CRABP2",
#               hjust=-2,
#               outPath="./DepMap_Results-updates/Freq_Dist_Rank-Luminal_Cell_Line_Enh_Regulated_Genes.pdf")
# 
# plotRandBgd(genes=luminal,
#             crispr=crispr,
#             proj=proj,
#             subtype="Luminal",
#             subtypes=subtypes,
#             outPath="./DepMap_Results-updates/Histogram_of_Rand_Background_Dependency_Scores-Luminal_Cell_Line_Enh_Regulated_Genes.pdf")
