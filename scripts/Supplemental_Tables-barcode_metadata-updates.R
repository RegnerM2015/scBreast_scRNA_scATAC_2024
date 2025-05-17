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
library(writexl)
library(data.table)

# Make output folder
if(dir.exists("Supplemental_Tables-barcode_metadata-updates")==TRUE){
  print("Directory Supplemental_Tables-barcode_metadata-updates already exists!")
}else{
  dir.create("Supplemental_Tables-barcode_metadata-updates")
}

################################################################################
# Write out scRNA barcode metadata from:
# - full cohort
# - basal cohort
# - luminal cohort
# - cell line cohort
################################################################################

# Full cohort

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
rna$barcode <- rownames(rna@meta.data)

# Make inferCNV column
rna$normal_cell_call[is.na(rna$normal_cell_call)] <- "not_tested"
rna$inferCNV <- plyr::mapvalues(rna$normal_cell_call,
                                from=c("cancer","not_tested","normal","unassigned"),
                                to=c("high","not_tested","low","ambiguous"))
# Make SCSubtype column
rna$SCSubtype[is.na(rna$SCSubtype)] <- "not_tested"

# Make patient sample columns
rna$sample <- rna$orig.ident
rna$patient <- plyr::mapvalues(rna$sample,
                               from=c("35A4AL",
                                      "4C2E5L",
                                      "35EE8L",
                                      "3821AL",
                                      "3B3E9L",
                                      "3C7D1L",
                                      "3D388L",
                                      "3FCDEL",
                                      "43E7BL",
                                      "43E7CL",
                                      "44F0AL",
                                      "45CB0L",
                                      "49758L",
                                      "49CFCL",
                                      "4AF75L",
                                      "4B146L"),
                               to=c("Patient_5",
                                    "Patient_6",
                                    "Patient_7",
                                    "Patient_8",
                                    "Patient_9",
                                    "Patient_10",
                                    "Patient_13",
                                    "Patient_11",
                                    "Patient_14-1",
                                    "Patient_14-2",
                                    "Patient_12",
                                    "Patient_15",
                                    "Patient_1",
                                    "Patient_2",
                                    "Patient_3",
                                    "Patient_4"))

# RNA
levels(factor(rna$RNA_snn_res.0.4))
rna$cluster <- rna$RNA_snn_res.0.4
rna$cell.type <- gsub('[[:digit:]]+', '', rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor)
rna$cell.type <- gsub("^\\-","",rna$cell.type)

Idents(object = rna) <- "RNA_snn_res.0.4"
# Assign celltype label to each cluster based on majority label
cells <- rna@meta.data %>% dplyr::group_by(RNA_snn_res.0.4) %>% dplyr::count(cell.type) 

cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in levels(factor(cells$RNA_snn_res.0.4))){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.4 ==i) %>% dplyr::arrange(desc(n))
  print(cells.sub)
  cluster.ids[[i]] <- cells.sub$cell.type[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$cell.type <- Idents(rna)

rna$cellTypeCluster <- paste0(rna$cluster,"-",rna$cell.type)

# Rename majority cell type cluster labels if they have inferCNV high cells
# (for example, some Luminal BC cells received a label of Mature Luminal, but were inferCNV high)
table(rna$inferCNV,rna$cellTypeCluster)

rna$cellTypeCluster <- gsub(" .*","",rna$cellTypeCluster)
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "11-Cancer",replacement = "11-Basal-like_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "5-PVL",replacement = "5-Perivascular-like_cells")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "1-T",replacement = "1-T_cells")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "18-Cancer",replacement = "18-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "0-CAFs",replacement = "0-Fibroblasts")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "10-CAFs",replacement = "10-Fibroblasts")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "16-B",replacement = "16-B_cells")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "13-CAFs",replacement = "13-Fibroblasts")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "12-Mature",replacement = "12-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "14-Cancer",replacement = "14-Basal-like_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "7-Luminal",replacement = "7-Luminal_Progenitors")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "3-CAFs",replacement = "3-Fibroblasts")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "15-Mast",replacement = "15-Mast_cells")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "21-Luminal",replacement = "21-Luminal_Progenitors")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "9-Mature",replacement = "9-Mature_Luminal")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "4-Myoepithelial",replacement = "4-Basal_Epithelial_cells")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "22-Myoepithelial",replacement = "22-Basal_Epithelial_cells")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "6-Cancer",replacement = "6-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "23-Cancer",replacement = "23-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "20-Mature",replacement = "20-Luminal_BC")

all.equal(as.character(rna$RNA_snn_res.0.4),as.character(gsub("\\-.*","",rna$cellTypeCluster)))

rna$cell_type_label <- as.character(sub(".*?-", "", rna$cellTypeCluster))
rna$cluster_number_RNA_snn_res.0.4 <- rna$RNA_snn_res.0.4
rna$cluster_number_RNA_snn_res.0.4_cell_type_label  <- rna$cellTypeCluster

# Subset to relevant columns
rna.meta <- rna@meta.data[,c("barcode",
                             "sample",
                             "patient",
                             "nCount_RNA",
                             "nFeature_RNA",
                             "percent.mt",
                             "cluster_number_RNA_snn_res.0.4",
                             "cell_type_label",
                             "cluster_number_RNA_snn_res.0.4_cell_type_label",
                             "inferCNV",
                             "SCSubtype"
                             )]


# Append UMAP coordinates
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna.meta)))
all.equal(rownames(rna.df),rownames(rna.meta))
all.equal(rownames(rna.df),colnames(rna))

rna.meta$UMAP_1 <- rna.df$UMAP_1
rna.meta$UMAP_2 <- rna.df$UMAP_2

full_cohort_rna_metadata <- rna.meta


# Basal cohort 

rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Basal_TN_Subset-TESTING.rds")
rna$barcode <- rownames(rna@meta.data)

# Make inferCNV column
rna$normal_cell_call[is.na(rna$normal_cell_call)] <- "not_tested"
rna$inferCNV <- plyr::mapvalues(rna$normal_cell_call,
                                from=c("cancer","not_tested","normal","unassigned"),
                                to=c("high","not_tested","low","ambiguous"))
# Make SCSubtype column
rna$SCSubtype[is.na(rna$SCSubtype)] <- "not_tested"

# Make patient sample columns
rna$sample <- rna$orig.ident
rna$patient <- plyr::mapvalues(rna$sample,
                               from=c("35A4AL",
                                      "4C2E5L",
                                      "49758L",
                                      "49CFCL",
                                      "4AF75L",
                                      "4B146L"),
                               to=c("Patient_5",
                                    "Patient_6",
                                    "Patient_1",
                                    "Patient_2",
                                    "Patient_3",
                                    "Patient_4"))

# RNA
levels(factor(rna$RNA_snn_res.0.015))
rna$cluster <- rna$RNA_snn_res.0.015
rna$cell.type <- gsub('[[:digit:]]+', '', rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor)
rna$cell.type <- gsub("^\\-","",rna$cell.type)

Idents(object = rna) <- "RNA_snn_res.0.015"
# Assign celltype label to each cluster based on majority label
cells <- rna@meta.data %>% dplyr::group_by(RNA_snn_res.0.015) %>% dplyr::count(cell.type) 

cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in levels(factor(cells$RNA_snn_res.0.015))){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.015 ==i) %>% dplyr::arrange(desc(n))
  print(cells.sub)
  cluster.ids[[i]] <- cells.sub$cell.type[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$cell.type <- Idents(rna)

rna$cellTypeCluster <- paste0(rna$cluster,"-",rna$cell.type)

# Rename majority cell type cluster labels if they have inferCNV high cells
# (for example, some Luminal BC cells received a label of Mature Luminal, but were inferCNV high)
table(rna$inferCNV,rna$cellTypeCluster)

rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "2-Cancer Basal SC",replacement = "2-Basal-like_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "3-Cancer Basal SC",replacement = "3-Basal-like_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "1-Luminal Progenitors",replacement = "1-Normal_Luminal_Progenitor")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "0-Myoepithelial",replacement = "0-Normal_Basal_Epithelial")

all.equal(as.character(rna$RNA_snn_res.0.015),as.character(gsub("\\-.*","",rna$cellTypeCluster)))

rna$cell_type_label <- as.character(sub(".*?-", "", rna$cellTypeCluster))
rna$cluster_number_RNA_snn_res.0.015 <- rna$RNA_snn_res.0.015
rna$cluster_number_RNA_snn_res.0.015_cell_type_label  <- rna$cellTypeCluster

# Subset to relevant columns
rna.meta <- rna@meta.data[,c("barcode",
                             "sample",
                             "patient",
                             "nCount_RNA",
                             "nFeature_RNA",
                             "percent.mt",
                             "cluster_number_RNA_snn_res.0.015",
                             "cell_type_label",
                             "cluster_number_RNA_snn_res.0.015_cell_type_label",
                             "inferCNV",
                             "SCSubtype"
)]

# Append UMAP coordinates
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna.meta)))
all.equal(rownames(rna.df),rownames(rna.meta))
all.equal(rownames(rna.df),colnames(rna))

rna.meta$UMAP_1 <- rna.df$UMAP_1
rna.meta$UMAP_2 <- rna.df$UMAP_2

basal_cohort_rna_metadata <- rna.meta


# Luminal cohort 

rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Luminal_TN_Subset-TESTING.rds")
rna$barcode <- rownames(rna@meta.data)

# Make inferCNV column
rna$normal_cell_call[is.na(rna$normal_cell_call)] <- "not_tested"
rna$inferCNV <- plyr::mapvalues(rna$normal_cell_call,
                                from=c("cancer","not_tested","normal","unassigned"),
                                to=c("high","not_tested","low","ambiguous"))
# Make SCSubtype column
rna$SCSubtype[is.na(rna$SCSubtype)] <- "not_tested"

# Make patient sample columns
rna$sample <- rna$orig.ident
rna$patient <- plyr::mapvalues(rna$sample,
                               from=c("35EE8L",
                                      "3821AL",
                                      "3B3E9L",
                                      "3C7D1L",
                                      "3D388L",
                                      "3FCDEL",
                                      "43E7BL",
                                      "43E7CL",
                                      "44F0AL",
                                      "45CB0L",
                                      "49758L",
                                      "49CFCL",
                                      "4AF75L",
                                      "4B146L"),
                               to=c("Patient_7",
                                    "Patient_8",
                                    "Patient_9",
                                    "Patient_10",
                                    "Patient_13",
                                    "Patient_11",
                                    "Patient_14-1",
                                    "Patient_14-2",
                                    "Patient_12",
                                    "Patient_15",
                                    "Patient_1",
                                    "Patient_2",
                                    "Patient_3",
                                    "Patient_4"))

# RNA
levels(factor(rna$RNA_snn_res.0.015))
rna$cluster <- rna$RNA_snn_res.0.015
rna$cell.type <- gsub('[[:digit:]]+', '', rna$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor)
rna$cell.type <- gsub("^\\-","",rna$cell.type)

Idents(object = rna) <- "RNA_snn_res.0.015"
# Assign celltype label to each cluster based on majority label
cells <- rna@meta.data %>% dplyr::group_by(RNA_snn_res.0.015) %>% dplyr::count(cell.type) 

cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in levels(factor(cells$RNA_snn_res.0.015))){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.015 ==i) %>% dplyr::arrange(desc(n))
  print(cells.sub)
  cluster.ids[[i]] <- cells.sub$cell.type[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$cell.type <- Idents(rna)

rna$cellTypeCluster <- paste0(rna$cluster,"-",rna$cell.type)

# Rename majority cell type cluster labels if they have inferCNV high cells
# (for example, some Luminal BC cells received a label of Mature Luminal, but were inferCNV high)
table(rna$inferCNV,rna$cellTypeCluster)

rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "3-Cancer LumA SC",replacement = "3-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "0-Mature Luminal",replacement = "0-Normal_Mature_Luminal")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "1-Cancer LumA SC",replacement = "1-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "5-Cancer LumA SC",replacement = "5-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "2-Mature Luminal",replacement = "2-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "6-Cancer LumA SC",replacement = "6-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "4-Mature Luminal",replacement = "4-Luminal_BC")
rna$cellTypeCluster <- str_replace(rna$cellTypeCluster,pattern = "7-Cancer LumA SC",replacement = "7-Luminal_BC")

all.equal(as.character(rna$RNA_snn_res.0.015),as.character(gsub("\\-.*","",rna$cellTypeCluster)))

rna$cell_type_label <- as.character(sub(".*?-", "", rna$cellTypeCluster))
rna$cluster_number_RNA_snn_res.0.015 <- rna$RNA_snn_res.0.015
rna$cluster_number_RNA_snn_res.0.015_cell_type_label  <- rna$cellTypeCluster

# Subset to relevant columns
rna.meta <- rna@meta.data[,c("barcode",
                             "sample",
                             "patient",
                             "nCount_RNA",
                             "nFeature_RNA",
                             "percent.mt",
                             "cluster_number_RNA_snn_res.0.015",
                             "cell_type_label",
                             "cluster_number_RNA_snn_res.0.015_cell_type_label",
                             "inferCNV",
                             "SCSubtype"
)]

# Append UMAP coordinates
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna.meta)))
all.equal(rownames(rna.df),rownames(rna.meta))
all.equal(rownames(rna.df),colnames(rna))

rna.meta$UMAP_1 <- rna.df$UMAP_1
rna.meta$UMAP_2 <- rna.df$UMAP_2

luminal_cohort_rna_metadata <- rna.meta


# Cell Line cohort 
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_CellLine_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
rna$barcode <- rownames(rna@meta.data)

# Make patient sample columns
rna$cell_line <- rna$orig.ident

# Subset to relevant columns
rna.meta <- rna@meta.data[,c("barcode",
                             "cell_line",
                             "nCount_RNA",
                             "nFeature_RNA",
                             "percent.mt"
)]

# Append UMAP coordinates
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna.meta)))
all.equal(rownames(rna.df),rownames(rna.meta))
all.equal(rownames(rna.df),colnames(rna))

rna.meta$UMAP_1 <- rna.df$UMAP_1
rna.meta$UMAP_2 <- rna.df$UMAP_2

cl_cohort_rna_metadata <- rna.meta



# Quick check on UMAPs

# Full Cohort
p1 <- ggplot(full_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=patient))+geom_point(size=0.1)+theme_classic()
p2 <- ggplot(full_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=sample))+geom_point(size=0.1)+theme_classic()
p3 <- ggplot(full_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=cluster_number_RNA_snn_res.0.4_cell_type_label))+geom_point(size=0.1)+theme_classic()
p4 <- ggplot(full_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=cluster_number_RNA_snn_res.0.4))+geom_point(size=0.1)+theme_classic()
p5 <- ggplot(full_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=inferCNV))+geom_point(size=0.1)+theme_classic()
p6 <- ggplot(full_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=SCSubtype))+geom_point(size=0.1)+theme_classic()
cowplot::plot_grid(p1,
                   p2,
                   p3,
                   p4,
                   p5,
                   p6,
                   ncol = 3, nrow = 2)
ggsave("./Supplemental_Tables-barcode_metadata-updates/Full_Cohort_scRNA_UMAP_check.pdf",width=16,height=9)

# Basal Cohort
p1 <- ggplot(basal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=patient))+geom_point(size=0.1)+theme_classic()
p2 <- ggplot(basal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=sample))+geom_point(size=0.1)+theme_classic()
p3 <- ggplot(basal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=cluster_number_RNA_snn_res.0.015_cell_type_label))+geom_point(size=0.1)+theme_classic()
p4 <- ggplot(basal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=cluster_number_RNA_snn_res.0.015))+geom_point(size=0.1)+theme_classic()
p5 <- ggplot(basal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=inferCNV))+geom_point(size=0.1)+theme_classic()
p6 <- ggplot(basal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=SCSubtype))+geom_point(size=0.1)+theme_classic()
cowplot::plot_grid(p1,
                   p2,
                   p3,
                   p4,
                   p5,
                   p6,
                   ncol = 3, nrow = 2)
ggsave("./Supplemental_Tables-barcode_metadata-updates/Basal_Cohort_scRNA_UMAP_check.pdf",width=16,height=9)


# Luminal Cohort
p1 <- ggplot(luminal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=patient))+geom_point(size=0.1)+theme_classic()
p2 <- ggplot(luminal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=sample))+geom_point(size=0.1)+theme_classic()
p3 <- ggplot(luminal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=cluster_number_RNA_snn_res.0.015_cell_type_label))+geom_point(size=0.1)+theme_classic()
p4 <- ggplot(luminal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=cluster_number_RNA_snn_res.0.015))+geom_point(size=0.1)+theme_classic()
p5 <- ggplot(luminal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=inferCNV))+geom_point(size=0.1)+theme_classic()
p6 <- ggplot(luminal_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=SCSubtype))+geom_point(size=0.1)+theme_classic()
cowplot::plot_grid(p1,
                   p2,
                   p3,
                   p4,
                   p5,
                   p6,
                   ncol = 3, nrow = 2)
ggsave("./Supplemental_Tables-barcode_metadata-updates/Luminal_Cohort_scRNA_UMAP_check.pdf",width=16,height=9)

# Cell Line Cohort
p1 <- ggplot(cl_cohort_rna_metadata,aes(UMAP_1,UMAP_2,color=cell_line))+geom_point(size=0.1)+theme_classic()
cowplot::plot_grid(p1,
                   ncol = 1, nrow = 1)
ggsave("./Supplemental_Tables-barcode_metadata-updates/Cell_Line_Cohort_scRNA_UMAP_check.pdf",width=6,height=6)


# Write metadata from each cohort into sheets
excel <- writexl::write_xlsx(x = list("Full_Cohort_barcode_meta" = full_cohort_rna_metadata,
                                  "Basal_Cohort_barcode_meta" = basal_cohort_rna_metadata,
                                  "Luminal_Cohort_barcode_meta" = luminal_cohort_rna_metadata,
                                  "Cell_Line_Cohort_barcode_meta" = cl_cohort_rna_metadata),
                             path = "./Supplemental_Tables-barcode_metadata-updates/scRNA_barcode_metadata.xlsx")


################################################################################
# Write out scATAC barcode metadata from:
# - full cohort
# - basal cohort
# - luminal cohort
# - cell line cohort
################################################################################

# Full cohort

# ATAC
atac <- loadArchRProject(path =  "./Patient_Samples_scATAC")
atac$barcode <- rownames(atac@cellColData)

# Grab UMAP coords
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "cluster"
atac.df$cluster <- sub(".*?-","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac.df$UMAP_1 <- atac.df$x
atac.df$UMAP_2 <- atac.df$y

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

all.equal(atac.df$cluster,atac$predictedGroup)

atac$UMAP_1 <- atac.df$UMAP_1
atac$UMAP_2 <- atac.df$UMAP_2

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
rna$barcode <- rownames(rna@meta.data)

# Make inferCNV column in RNA
rna$normal_cell_call[is.na(rna$normal_cell_call)] <- "not_tested"
rna$inferCNV <- plyr::mapvalues(rna$normal_cell_call,
                                from=c("cancer","not_tested","normal","unassigned"),
                                to=c("high","not_tested","low","ambiguous"))
# Make SCSubtype column in RNA
rna$SCSubtype[is.na(rna$SCSubtype)] <- "not_tested"

atac.meta <- as.data.frame(atac@cellColData)

all.equal(rownames(atac.meta),atac$cellNames)
atac.meta$cellNames <- rownames(atac.meta)
all.equal(atac.meta$cellNames,atac$cellNames)
length(intersect(atac.meta$predictedCell,colnames(rna)))
length(intersect(atac.meta$predictedCell,rownames(rna@meta.data)))
length(intersect(atac.meta$predictedCell,rna$barcode))
atac.meta$barcode <- atac.meta$predictedCell
atac.meta$id  <- 1:nrow(atac.meta)

atac.meta <- merge(atac.meta,rna@meta.data,by="barcode")
atac.meta <- atac.meta[order(atac.meta$id), ]

all.equal(atac.meta$cellNames,atac$cellNames)

# Make patient sample columns
atac.meta$sample <- atac.meta$Sample
atac.meta$patient <- plyr::mapvalues(atac.meta$sample,
                               from=c("35A4AL",
                                      "4C2E5L",
                                      "35EE8L",
                                      "3821AL",
                                      "3B3E9L",
                                      "3C7D1L",
                                      "3D388L",
                                      "3FCDEL",
                                      "43E7BL",
                                      "43E7CL",
                                      "44F0AL",
                                      "45CB0L",
                                      "49758L",
                                      "49CFCL",
                                      "4AF75L",
                                      "4B146L"),
                               to=c("Patient_5",
                                    "Patient_6",
                                    "Patient_7",
                                    "Patient_8",
                                    "Patient_9",
                                    "Patient_10",
                                    "Patient_13",
                                    "Patient_11",
                                    "Patient_14-1",
                                    "Patient_14-2",
                                    "Patient_12",
                                    "Patient_15",
                                    "Patient_1",
                                    "Patient_2",
                                    "Patient_3",
                                    "Patient_4"))

# Rename predictedGroup
atac.meta$predicted_cluster_number_RNA_snn_res.0.4 <- atac.meta$predictedGroup

# Rename cellNames back to barcode
atac.meta$barcode <- atac.meta$cellNames

# Subset to relevant columns
atac.meta <- atac.meta[,c("barcode",
                          "sample",
                          "patient",
                          "TSSEnrichment",
                          "nFrags",
                          "predicted_cluster_number_RNA_snn_res.0.4", 
                          "inferCNV",
                          "SCSubtype",
                          "UMAP_1",
                          "UMAP_2")]

full_cohort_atac_metadata <- atac.meta


# Basal cohort

# ATAC
atac <- loadArchRProject(path =  "./Basal_TN_Samples_scATAC-TESTING3")
atac$barcode <- rownames(atac@cellColData)

# Grab UMAP coords
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "cluster"
atac.df$cluster <- sub(".*?-","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac.df$UMAP_1 <- atac.df$x
atac.df$UMAP_2 <- atac.df$y

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

all.equal(atac.df$cluster,atac$predictedGroup)

atac$UMAP_1 <- atac.df$UMAP_1
atac$UMAP_2 <- atac.df$UMAP_2

# Read in scRNA
rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Basal_TN_Subset-TESTING.rds")
rna$barcode <- rownames(rna@meta.data)

# Make inferCNV column in RNA
rna$normal_cell_call[is.na(rna$normal_cell_call)] <- "not_tested"
rna$inferCNV <- plyr::mapvalues(rna$normal_cell_call,
                                from=c("cancer","not_tested","normal","unassigned"),
                                to=c("high","not_tested","low","ambiguous"))
# Make SCSubtype column in RNA
rna$SCSubtype[is.na(rna$SCSubtype)] <- "not_tested"

atac.meta <- as.data.frame(atac@cellColData)

all.equal(rownames(atac.meta),atac$cellNames)
atac.meta$cellNames <- rownames(atac.meta)
all.equal(atac.meta$cellNames,atac$cellNames)
length(intersect(atac.meta$predictedCell,colnames(rna)))
length(intersect(atac.meta$predictedCell,rownames(rna@meta.data)))
length(intersect(atac.meta$predictedCell,rna$barcode))
atac.meta$barcode <- atac.meta$predictedCell
atac.meta$id  <- 1:nrow(atac.meta)

atac.meta <- merge(atac.meta,rna@meta.data,by="barcode")
atac.meta <- atac.meta[order(atac.meta$id), ]

all.equal(atac.meta$cellNames,atac$cellNames)

# Make patient sample columns
atac.meta$sample <- atac.meta$Sample
atac.meta$patient <- plyr::mapvalues(atac.meta$sample,
                                     from=c("35A4AL",
                                            "4C2E5L",
                                            "49758L",
                                            "49CFCL",
                                            "4AF75L",
                                            "4B146L"),
                                     to=c("Patient_5",
                                          "Patient_6",
                                          "Patient_1",
                                          "Patient_2",
                                          "Patient_3",
                                          "Patient_4"))

# Rename predictedGroup
atac.meta$predicted_cluster_number_RNA_snn_res.0.015 <- atac.meta$predictedGroup

# Rename cellNames back to barcode
atac.meta$barcode <- atac.meta$cellNames

# Subset to relevant columns
atac.meta <- atac.meta[,c("barcode",
                          "sample",
                          "patient",
                          "TSSEnrichment",
                          "nFrags",
                          "predicted_cluster_number_RNA_snn_res.0.015", 
                          "inferCNV",
                          "SCSubtype",
                          "UMAP_1",
                          "UMAP_2")]

basal_cohort_atac_metadata <- atac.meta


# Luminal Cohort

# ATAC
atac <- loadArchRProject(path =  "./Luminal_TN_Samples_scATAC-TESTING3")
atac$barcode <- rownames(atac@cellColData)

# Grab UMAP coords
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "cluster"
atac.df$cluster <- sub(".*?-","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac.df$UMAP_1 <- atac.df$x
atac.df$UMAP_2 <- atac.df$y

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

all.equal(atac.df$cluster,atac$predictedGroup)

atac$UMAP_1 <- atac.df$UMAP_1
atac$UMAP_2 <- atac.df$UMAP_2

# Read in scRNA
rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Luminal_TN_Subset-TESTING.rds")
rna$barcode <- rownames(rna@meta.data)

# Make inferCNV column in RNA
rna$normal_cell_call[is.na(rna$normal_cell_call)] <- "not_tested"
rna$inferCNV <- plyr::mapvalues(rna$normal_cell_call,
                                from=c("cancer","not_tested","normal","unassigned"),
                                to=c("high","not_tested","low","ambiguous"))
# Make SCSubtype column in RNA
rna$SCSubtype[is.na(rna$SCSubtype)] <- "not_tested"

atac.meta <- as.data.frame(atac@cellColData)

all.equal(rownames(atac.meta),atac$cellNames)
atac.meta$cellNames <- rownames(atac.meta)
all.equal(atac.meta$cellNames,atac$cellNames)
length(intersect(atac.meta$predictedCell,colnames(rna)))
length(intersect(atac.meta$predictedCell,rownames(rna@meta.data)))
length(intersect(atac.meta$predictedCell,rna$barcode))
atac.meta$barcode <- atac.meta$predictedCell
atac.meta$id  <- 1:nrow(atac.meta)

atac.meta <- merge(atac.meta,rna@meta.data,by="barcode")
atac.meta <- atac.meta[order(atac.meta$id), ]

all.equal(atac.meta$cellNames,atac$cellNames)

# Make patient sample columns
atac.meta$sample <- atac.meta$Sample
atac.meta$patient <- plyr::mapvalues(atac.meta$sample,
                                     from=c("35EE8L",
                                            "3821AL",
                                            "3B3E9L",
                                            "3C7D1L",
                                            "3D388L",
                                            "3FCDEL",
                                            "43E7BL",
                                            "43E7CL",
                                            "44F0AL",
                                            "45CB0L",
                                            "49758L",
                                            "49CFCL",
                                            "4AF75L",
                                            "4B146L"),
                                     to=c("Patient_7",
                                          "Patient_8",
                                          "Patient_9",
                                          "Patient_10",
                                          "Patient_13",
                                          "Patient_11",
                                          "Patient_14-1",
                                          "Patient_14-2",
                                          "Patient_12",
                                          "Patient_15",
                                          "Patient_1",
                                          "Patient_2",
                                          "Patient_3",
                                          "Patient_4"))


# Rename predictedGroup
atac.meta$predicted_cluster_number_RNA_snn_res.0.015 <- atac.meta$predictedGroup

# Rename cellNames back to barcode
atac.meta$barcode <- atac.meta$cellNames

# Subset to relevant columns
atac.meta <- atac.meta[,c("barcode",
                          "sample",
                          "patient",
                          "TSSEnrichment",
                          "nFrags",
                          "predicted_cluster_number_RNA_snn_res.0.015", 
                          "inferCNV",
                          "SCSubtype",
                          "UMAP_1",
                          "UMAP_2")]

luminal_cohort_atac_metadata <- atac.meta


# Cell Line cohort 


# ATAC
atac <- loadArchRProject(path =  "./CellLine_Samples_scATAC-update")
atac$barcode <- rownames(atac@cellColData)

# Grab UMAP coords
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "cluster"
atac.df$cluster <- sub(".*?-","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac.df$UMAP_1 <- atac.df$x
atac.df$UMAP_2 <- atac.df$y

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

all.equal(atac.df$cluster,atac$predictedGroup)

atac$UMAP_1 <- atac.df$UMAP_1
atac$UMAP_2 <- atac.df$UMAP_2

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_CellLine_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
rna$barcode <- rownames(rna@meta.data)

atac.meta <- as.data.frame(atac@cellColData)

all.equal(rownames(atac.meta),atac$cellNames)
atac.meta$cellNames <- rownames(atac.meta)
all.equal(atac.meta$cellNames,atac$cellNames)
length(intersect(atac.meta$predictedCell,colnames(rna)))
length(intersect(atac.meta$predictedCell,rownames(rna@meta.data)))
length(intersect(atac.meta$predictedCell,rna$barcode))
atac.meta$barcode <- atac.meta$predictedCell
atac.meta$id  <- 1:nrow(atac.meta)

atac.meta <- merge(atac.meta,rna@meta.data,by="barcode")
atac.meta <- atac.meta[order(atac.meta$id), ]

all.equal(atac.meta$cellNames,atac$cellNames)

# Make patient sample columns
atac.meta$cell_line <- atac.meta$Sample

# Rename cellNames back to barcode
atac.meta$barcode <- atac.meta$cellNames

# Subset to relevant columns
atac.meta <- atac.meta[,c("barcode",
                          "cell_line",
                          "TSSEnrichment",
                          "nFrags",
                          "UMAP_1",
                          "UMAP_2")]

cl_cohort_atac_metadata <- atac.meta


# Quick check on UMAPs

# Full Cohort
p1 <- ggplot(full_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=patient))+geom_point(size=0.1)+theme_classic()
p2 <- ggplot(full_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=sample))+geom_point(size=0.1)+theme_classic()
p3 <- ggplot(full_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=predicted_cluster_number_RNA_snn_res.0.4))+geom_point(size=0.1)+theme_classic()
p4 <- ggplot(full_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=inferCNV))+geom_point(size=0.1)+theme_classic()
p5 <- ggplot(full_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=SCSubtype))+geom_point(size=0.1)+theme_classic()
cowplot::plot_grid(p1,
                   p2,
                   p3,
                   p4,
                   p5,
                   ncol = 3, nrow = 2)
ggsave("./Supplemental_Tables-barcode_metadata-updates/Full_Cohort_scATAC_UMAP_check.pdf",width=16,height=9)

# Basal Cohort
p1 <- ggplot(basal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=patient))+geom_point(size=0.1)+theme_classic()
p2 <- ggplot(basal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=sample))+geom_point(size=0.1)+theme_classic()
p3 <- ggplot(basal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=predicted_cluster_number_RNA_snn_res.0.015))+geom_point(size=0.1)+theme_classic()
p4 <- ggplot(basal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=inferCNV))+geom_point(size=0.1)+theme_classic()
p5 <- ggplot(basal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=SCSubtype))+geom_point(size=0.1)+theme_classic()
cowplot::plot_grid(p1,
                   p2,
                   p3,
                   p4,
                   p5,
                   ncol = 3, nrow = 2)
ggsave("./Supplemental_Tables-barcode_metadata-updates/Basal_Cohort_scATAC_UMAP_check.pdf",width=16,height=9)


# Luminal Cohort
p1 <- ggplot(luminal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=patient))+geom_point(size=0.1)+theme_classic()
p2 <- ggplot(luminal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=sample))+geom_point(size=0.1)+theme_classic()
p3 <- ggplot(luminal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=predicted_cluster_number_RNA_snn_res.0.015))+geom_point(size=0.1)+theme_classic()
p4 <- ggplot(luminal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=inferCNV))+geom_point(size=0.1)+theme_classic()
p5 <- ggplot(luminal_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=SCSubtype))+geom_point(size=0.1)+theme_classic()
cowplot::plot_grid(p1,
                   p2,
                   p3,
                   p4,
                   p5,
                   ncol = 3, nrow = 2)
ggsave("./Supplemental_Tables-barcode_metadata-updates/Luminal_Cohort_scATAC_UMAP_check.pdf",width=16,height=9)

# Cell Line Cohort
p1 <- ggplot(cl_cohort_atac_metadata,aes(UMAP_1,UMAP_2,color=cell_line))+geom_point(size=0.1)+theme_classic()
cowplot::plot_grid(p1,
                   ncol = 1, nrow = 1)
ggsave("./Supplemental_Tables-barcode_metadata-updates/Cell_Line_Cohort_scATAC_UMAP_check.pdf",width=6,height=6)


# Write metadata from each cohort into sheets
excel <- writexl::write_xlsx(x = list("Full_Cohort_barcode_meta" = full_cohort_atac_metadata,
                                      "Basal_Cohort_barcode_meta" = basal_cohort_atac_metadata,
                                      "Luminal_Cohort_barcode_meta" = luminal_cohort_atac_metadata,
                                      "Cell_Line_Cohort_barcode_meta" = cl_cohort_atac_metadata),
                             path = "./Supplemental_Tables-barcode_metadata-updates/scATAC_barcode_metadata.xlsx")

print("xlsx files and UMAP checks successfully generated!")
