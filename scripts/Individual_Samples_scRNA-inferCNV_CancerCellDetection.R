################################################################################
# Matt Regner and Aatish Thennavan
# Franco Lab
# Description: This script performs the following tasks  
#         1) Read in Seurat object with cluster cell type annotations
#         2) Subset to epithelial and immune/endothelial fraction
#         3) Run inferCNV
#         4) Run correlation scatter plot method to infer cancer cells
################################################################################
library(Seurat)
library(tidyverse)
library(dplyr)
library(scales)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(stringi)
library(infercnv)
library(EnsDb.Hsapiens.v86)
library(parallel)
library(future)
library(cluster)
library(fpc)
options(future.globals.maxSize = 6000 * 1024^2)
plan("multicore")
set.seed(1)

args=(commandArgs(TRUE))
ID <- args[[1]]
print(ID)

rna <- readRDS(paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering-CellTypeAnnotations.rds"))

Idents(rna) <- "Cluster_Majority_No_CCA_Reference_Cell_Type_Minor"

if(dir.exists("inferCNV_scRNA")==TRUE){
  print("Directory inferCNV_scRNA already exists!")
}else{
  dir.create("inferCNV_scRNA")
}

# Subset Seurat object to epithelial fraction and immune/endothelial
idx.1 <- grep("T cells",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.2 <- grep("T-cells",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.3 <- grep("B cells",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.4 <- grep("NK cells",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.5 <- grep("NKT cells",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.6 <- grep("Macrophage",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.7 <- grep("Monocyte",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.8 <- grep("Cycling_Myeloid",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.9 <- grep("DCs",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.10 <- grep("Plasmablasts",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.11 <- grep("Mast",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.12 <- grep("Endothelial",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.13 <- grep("ithelial",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.14 <- grep("Luminal",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.15 <- grep("Basal",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.16 <- grep("Cancer",unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx <- unique(c(idx.1,idx.2,idx.3,idx.4,idx.5,idx.6,
                idx.7,idx.8,idx.9,idx.10,idx.11,idx.12,
                idx.13,idx.14,idx.15,idx.16))
clusters.to.use <- unique(rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor)[idx]
rna.subset <- subset(x = rna, idents = clusters.to.use)

# Set reference group names
idx.1 <- grep("T cells",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.2 <- grep("T-cells",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.3 <- grep("B cells",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.4 <- grep("NK cells",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.5 <- grep("NKT cells",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.6 <- grep("Macrophage",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.7 <- grep("Monocyte",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.8 <- grep("Cycling_Myeloid",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.9 <- grep("DCs",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.10 <- grep("Plasmablasts",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.11 <- grep("Mast",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx.12 <- grep("Endothelial",unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor))
idx <- unique(c(idx.1,idx.2,idx.3,idx.4,idx.5,idx.6,
                idx.7,idx.8,idx.9,idx.10,idx.11,idx.12))
ref_group_names <- unique(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor)[idx]

rna.subset$cell.type <- ifelse(rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor %in% ref_group_names,
                               rna.subset$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor,ID)

# Set annotations file
annotations_df <- data.frame(cell.name=rownames(rna.subset@meta.data),
                             cell.type=rna.subset@meta.data$cell.type)
write.table(annotations_df,
            paste0("./inferCNV_scRNA/",ID,"_annotations_file.txt"),
            sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)

# Create inferCNV object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=rna.subset@assays$RNA@counts[,colnames(rna.subset)],
                                    annotations_file=paste0("./inferCNV_scRNA/",ID,"_annotations_file.txt"),
                                    delim="\t",
                                    gene_order_file="./miscellaneous/refdata-cellranger-GRCh38-3_0_0-gene-ordering-file.txt",
                                    ref_group_names=ref_group_names) 
# Run inferCNV
infercnv_obj = infercnv::run(infercnv_obj,
                             out_dir= paste0("./inferCNV_scRNA/",ID,"_inferCNV_output"),
                             cutoff=0.1,
                             window_length=101,
                             max_centered_threshold=3,
                             cluster_by_groups=F,
                             plot_steps=F,
                             denoise=T,
                             sd_amplifier=1.3,
                             analysis_mode = "samples",
                             HMM=F)

# Read in denoised values and perform scatter plot method for cancer cell detection
##Correlation to itself
#1) Read in your CNA values file for the TUMOR
DCIS<-read.table(paste0("./inferCNV_scRNA/",ID,"_inferCNV_output/infercnv.observations.txt"), sep="", header = TRUE)
infercnv_output <- as.data.frame(t(DCIS))

#2) Rescale all the CNA values for the TUMOR
scaled_df <- as.data.frame(rescale(as.matrix(infercnv_output), c(-1,1)))

#3) Correlate values to the top 5% of cells with high CNA values in the TUMOR
CNA_values <- apply(scaled_df, 1, function(y) {
  #y[is.na(y)] <- 0
  #scaled_y <- rescale(y, c(-1, 1))
  return(mean(y^2))
})
CNA_value_df <- data.frame(
  row.names = names(CNA_values),
  CNA_value = CNA_values
)
CNA_order <- order(CNA_value_df$CNA_value, decreasing=T)
ordered_CNA_values  <- data.frame(
  row.names = rownames(CNA_value_df)[CNA_order],
  CNA_value = CNA_value_df[CNA_order,]
)
top_cancer <- head(ordered_CNA_values, nrow(ordered_CNA_values)*0.05)
top_cancer_CNV_average <- apply(infercnv_output[rownames(top_cancer),], 2, mean)
cancer_correlations <- apply(infercnv_output, 1, function(x) {
  if (length(unique(as.numeric(x))) == 1) {
    cor_result <- data.frame(cor.estimate="no_CNVs_recorded",
                             cor.p.value="no_CNVs_recorded")
  } else {
    cor <- cor.test(as.numeric(x), top_cancer_CNV_average, method = "kendall")
    cor_result <- data.frame(cor$estimate, cor$p.value)
  }
  return(cor_result)
})
correlation_df <- do.call("rbind", cancer_correlations)
epithelial_metadata<-cbind(CNA_value_df, correlation_df)

print(paste0("Rownames check CNA_value_df with correlation_df: ",identical(rownames(CNA_value_df), rownames(correlation_df))))

print(paste0("Rownames check infercnv_output with epithelial_metadata: ",identical(rownames(infercnv_output), rownames(epithelial_metadata))))

if(dir.exists("inferCNV_scRNA")==TRUE){
  print("Directory inferCNV_scRNA already exists!")
}else{
  dir.create("inferCNV_scRNA")
}

# Plot scatter plot
ggplot(epithelial_metadata,aes(x=CNA_value,y=cor.estimate))+
  geom_point()+
  geom_hline(yintercept = 0.4,linetype="dashed")+
  geom_vline(xintercept = 0.02,linetype="dashed")+
  ggtitle(paste0("inferCNV cancer cell detection: ",ID))+
  theme_bw()
ggsave(paste0("./inferCNV_scRNA/ScatterPlot_CancerCellDetection_",ID,".pdf"),
       width = 6,height = 5)

# Label cancer cells
quad_df <- data.frame(
  row.names = rownames(epithelial_metadata),
  CNA_value = epithelial_metadata$CNA_value, 
  cor.estimate = epithelial_metadata$cor.estimate
)
# scale data:
scaled_quad_df <- scale(quad_df) %>% as.data.frame()
# run silhouette cluster analysis to determine clusters and thresholds:
pamk_result <- pamk(scaled_quad_df, krange=1:4)
print(pamk_result$nc)
silhouette_result <- pam(scaled_quad_df, pamk_result$nc)
#saveRDS(silhouette_result, paste0(Robject_dir, "silhouette_result.Rdata"))



# if no. clusters estimated to be > 1, use cluster information to set 
# normal vs cancer thresholds:
if (pamk_result$nc > 1) {
  
  cancer_x_threshold_sd_multiplier <- 2
  cancer_y_threshold_sd_multiplier <- 1.5
  
  sil_values <- as.data.frame(silhouette_result$silinfo$widths)
  sil_result <- data.frame(row.names=names(silhouette_result$clustering),
                           cluster=silhouette_result$clustering,
                           sil_width=sil_values$sil_width)
  
  # add sil_result to epithelial_metadata:
  epithelial_metadata <- cbind(epithelial_metadata, sil_result)
  
  # determine normal and cancer clusters by determining the max CNA values and
  # correlation with top 5% cancer:
  cluster_split <- split(epithelial_metadata, epithelial_metadata$cluster)
  names(cluster_split) <- paste0("cluster_", names(cluster_split))
  # determine order of clusters by adding mean CNA and correlation values:
  cluster_means <- sort(
    unlist(
      lapply(cluster_split, function(x) {
        return(mean(x$CNA_value) + mean(x$cor.estimate))
      })
    )
  )
  # determine second cluster from axes as cancer cluster closest to axes:
  first_cancer_cluster <- names(cluster_means[2])
  first_cancer_df <- eval(parse(text=paste0("cluster_split$", first_cancer_cluster)))

  # define x-axis as 2 std devs left of mean:
  CNA_mean <- mean(first_cancer_df$CNA_value)
  CNA_std_dev <- sd(first_cancer_df$CNA_value)
  x_int <- round(CNA_mean - (cancer_x_threshold_sd_multiplier*CNA_std_dev), 3)
  # define y-axis as 1.5 std devs below mean:
  cor_mean <- mean(first_cancer_df$cor.estimate)
  cor_std_dev <- sd(first_cancer_df$cor.estimate)
  y_int <- round(cor_mean - (cancer_y_threshold_sd_multiplier*cor_std_dev), 3)
  
  # if (adjust_normal_thresholds) {
  #   if (sample_name %in% adjustment_df$sample_id) {
  #     x_int <- adjustment_df$x_int[adjustment_df$sample_id == sample_name]
  #     y_int <- adjustment_df$y_int[adjustment_df$sample_id == sample_name]
  #   }
  # }
  
  # define normal and cancer cells:
  epithelial_metadata$normal_cell_call <- "cancer"
  
  epithelial_metadata$normal_cell_call[
    epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
  ] <- "normal"
  epithelial_metadata$normal_cell_call[
    epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate > y_int
  ] <- "unassigned"
  epithelial_metadata$normal_cell_call[
    epithelial_metadata$CNA_value > x_int & epithelial_metadata$cor.estimate < y_int
  ] <- "unassigned"
  
  # create quad plot:
  ggplot(epithelial_metadata, 
         aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))+
    geom_point()+
    scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                       labels=c("Cancer", "Normal", "Unassigned"))+
    xlab("Infercnv level")+
    ylab("Corr. with top 5% cancer")+
    theme(legend.title = element_blank())+
    geom_vline(xintercept = x_int)+
    geom_hline(yintercept = y_int)
  ggsave(paste0("./inferCNV_scRNA/ScatterPlot_CancerCellDetection_",ID,".pdf"),
         width = 6,height = 5)
  
  saveRDS(epithelial_metadata,paste0("./inferCNV_scRNA/Epithelial_Metadata_",ID,".rds"))
} else {
  cancer_x_threshold_sd_multiplier <- 1
  cancer_y_threshold_sd_multiplier <- 1.25
  
  # define x-axis as 1 std devs left of mean:
  CNA_mean <- mean(epithelial_metadata$CNA_value)
  CNA_std_dev <- sd(epithelial_metadata$CNA_value)
  x_int <- round(CNA_mean - (cancer_x_threshold_sd_multiplier*CNA_std_dev), 3)
  # define y-axis as 1.25 std devs below mean:
  cor_mean <- mean(epithelial_metadata$cor.estimate)
  cor_std_dev <- sd(epithelial_metadata$cor.estimate)
  y_int <- round(cor_mean - (cancer_y_threshold_sd_multiplier*cor_std_dev), 3)
  
  # if (adjust_normal_thresholds) {
  #   if (sample_name %in% adjustment_df$sample_id) {
  #     x_int <- adjustment_df$x_int[adjustment_df$sample_id == sample_name]
  #     y_int <- adjustment_df$y_int[adjustment_df$sample_id == sample_name]
  #   }
  # }
  
  # define normal and cancer cells:
  epithelial_metadata$normal_cell_call <- "cancer"
  
  epithelial_metadata$normal_cell_call[
    epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
  ] <- "normal"
  epithelial_metadata$normal_cell_call[
    epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate > y_int
  ] <- "unassigned"
  epithelial_metadata$normal_cell_call[
    epithelial_metadata$CNA_value > x_int & epithelial_metadata$cor.estimate < y_int
  ] <- "unassigned"
  
  # create quad plot:
  ggplot(epithelial_metadata, 
         aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))+
    geom_point()+
    scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                       labels=c("Cancer", "Normal", "Unassigned"))+
    xlab("Infercnv level")+
    ylab("Corr. with top 5% cancer")+
    theme(legend.title = element_blank())+
    geom_vline(xintercept = x_int)+
    geom_hline(yintercept = y_int)
  ggsave(paste0("./inferCNV_scRNA/ScatterPlot_CancerCellDetection_",ID,".pdf"),
         width = 6,height = 5)
  
  saveRDS(epithelial_metadata,paste0("./inferCNV_scRNA/Epithelial_Metadata_",ID,".rds"))
}

# Append cancer/normal predicitons to original Seurat object
rna$barcode <- rownames(rna@meta.data)
epithelial_metadata$barcode <- gsub("\\.","-",rownames(epithelial_metadata))
rna@meta.data <- merge(rna@meta.data,epithelial_metadata,by="barcode",all.x=TRUE)
rownames(rna@meta.data) <- rna$barcode

if(dir.exists("Processed_Seurat_Objects_scRNA")==TRUE){
  print("Directory Processed_Seurat_Objects_scRNA already exists!")
}else{
  dir.create("Processed_Seurat_Objects_scRNA")
}
saveRDS(rna,paste0("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_",ID,"-MultiKClustering-CellTypeAnnotations-inferCNV.rds"))
