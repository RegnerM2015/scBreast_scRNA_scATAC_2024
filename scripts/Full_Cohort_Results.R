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
library(patchwork)
source("scripts/Matrix.utils.R")

set.seed(1234)
addArchRThreads(threads = 32)
addArchRGenome("hg38")

# Make output folder
if(dir.exists("Full_Cohort_Results")){
  print("Directory Full_Cohort_Results already exists!")
}else{
  dir.create("Full_Cohort_Results")
}

# Plot RNA and ATAC UMAP by cell type cluster 

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

# Plot cell type UMAPs for RNA/ATAC
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
all.equal(rownames(rna.df),rownames(rna@meta.data))
all.equal(rownames(rna.df),colnames(rna))
rna.df$Sample <- rna$orig.ident

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

# Cancer cell check
# Annotate cancer normal cells
rna$SCSubtype[is.na(rna$SCSubtype)] <- "notTested_or_predictedNormal_or_unassigned"
rna$normal_cell_call[is.na(rna$normal_cell_call)] <- "not_tested"

rna$SCSubtype_Cancer_Normal <- ifelse(rna$SCSubtype ==  "notTested_or_predictedNormal_or_unassigned",
                                      rna$normal_cell_call,rna$SCSubtype)

confusionMatrix(rna$SCSubtype_Cancer_Normal,rna$cellTypeCluster)


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

# Help sort the cluster numbers:
###############################
basalCancer <- grep("Basal-like_BC",levels(factor(rna$cellTypeCluster)))
luminalCancer <- grep("Luminal_BC",levels(factor(rna$cellTypeCluster)))
matureLuminal <- grep("\\bMature_Luminal\\b",levels(factor(rna$cellTypeCluster)))
normalBasal <- grep("Basal_Epithelial",levels(factor(rna$cellTypeCluster)))
normalLP <- grep("Luminal_Progenitors",levels(factor(rna$cellTypeCluster)))
fibroblasts <- grep("Fibroblasts",levels(factor(rna$cellTypeCluster)))
tCells <- grep("T_cells",levels(factor(rna$cellTypeCluster)))
bCells <- grep("B_cells",levels(factor(rna$cellTypeCluster)))
macrophage <- grep("Macrophage",levels(factor(rna$cellTypeCluster)))
plasmablasts <- grep("Plasmablasts",levels(factor(rna$cellTypeCluster)))
mastCells <- grep("Mast_cells",levels(factor(rna$cellTypeCluster)))
pvl <- grep("Perivascular-like_cells",levels(factor(rna$cellTypeCluster)))
endothelial <- grep("Endothelial",levels(factor(rna$cellTypeCluster)))


cell.types.idx <- c(basalCancer,normalLP,normalBasal,
                    luminalCancer,matureLuminal,
                    fibroblasts,pvl,endothelial,
                    tCells,bCells,plasmablasts,macrophage,mastCells)

store <- numeric(0)
for(i in 1:length(cell.types.idx)){
  name <- levels(factor(rna$cellTypeCluster))[cell.types.idx[i]]
  print(gsub("-.*","",name))
  new.name <- gsub("-.*","",name)
  new.num <- as.numeric(new.name)
  store[i] <- new.num
}
print(store)
#####################################################

my_levels <- store


# Relevel object@ident
rna$clusterNumber <- factor(x = rna$RNA_snn_res.0.4, levels = my_levels)

basal.cols <- colorRampPalette(c("#A50505", "#FF9B9B"))
basal.cols <- basal.cols(6)

luminal.cols <- colorRampPalette(c("#096B02", "#A4DCA0"))
luminal.cols <- luminal.cols(6)

stromal.cols <- colorRampPalette(c("#8B059F", "#EFBEF6"))
stromal.cols <- stromal.cols(7)

t.cols <- "#FF8000"
b.cols <- c("#404040","#C0C0C0")
myeloid <- c("#244BCB","#AEBFF6")

cols <- c(basal.cols,luminal.cols,stromal.cols,t.cols,b.cols,myeloid)

all.equal(rownames(rna.df),rownames(rna@meta.data))
rna.df$cellTypeCluster <- factor(rna$cellTypeCluster,levels=c("11-Basal-like_BC",
                                                              "14-Basal-like_BC",
                                                              "21-Luminal_Progenitors",
                                                              "7-Luminal_Progenitors",
                                                              "22-Basal_Epithelial_cells",
                                                              "4-Basal_Epithelial_cells",
                                                              "12-Luminal_BC",
                                                              "18-Luminal_BC",
                                                              "20-Luminal_BC",
                                                              "23-Luminal_BC",
                                                              "6-Luminal_BC",
                                                              "9-Mature_Luminal",
                                                              "0-Fibroblasts",
                                                              "10-Fibroblasts",
                                                              "13-Fibroblasts",
                                                              "3-Fibroblasts",
                                                              "5-Perivascular-like_cells",
                                                              "17-Endothelial",
                                                              "2-Endothelial",
                                                              "1-T_cells",
                                                              "16-B_cells",
                                                              "19-Plasmablasts",
                                                              "8-Macrophage",
                                                              "15-Mast_cells"
                                                              ))
rna.df$clusterNumber <- rna$clusterNumber

ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = clusterNumber))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle(paste0("scRNA-seq colored by cluster (n=",nrow(rna.df)," cells)"))
ggsave("./Full_Cohort_Results/UMAP_RNA_Cluster_Number.pdf",width = 8,height = 7)

p1 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = clusterNumber))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()+
  ggtitle(paste0("scRNA-seq colored by cluster (n=",nrow(rna.df)," cells)"))
LabelClusters(p1,id="clusterNumber",color="black",repel = T,size=6)
ggsave("./Full_Cohort_Results/UMAP_RNA_Cluster_Number-Labels.pdf",width = 8,height = 7)

ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cellTypeCluster))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle(paste0("scRNA-seq colored by cell type cluster (n=",nrow(rna.df)," cells)"))
ggsave("./Full_Cohort_Results/UMAP_RNA_Cell_Type_Cluster.pdf",width = 8,height = 7)

p1 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cellTypeCluster))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()+
  ggtitle(paste0("scRNA-seq colored by cell type cluster (n=",nrow(rna.df)," cells)"))
LabelClusters(p1,id="cellTypeCluster",color="black",repel = T,size=6)
ggsave("./Full_Cohort_Results/UMAP_RNA_Cell_Type_Cluster-Labels.pdf",width = 8,height = 7)

# ATAC
atac <- loadArchRProject(path =  "./Patient_Samples_scATAC")

atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$cell.type <- sub(".*?-","",atac.df$color)
atac.df$cluster <- as.factor(atac.df$cell.type)

# Relevel object@ident
atac.df$cluster.new <- factor(x = atac.df$cluster, levels = my_levels)# use my_levels from RNA above

basal.cols <- colorRampPalette(c("#A50505", "#FF9B9B"))
basal.cols <- basal.cols(6)

luminal.cols <- colorRampPalette(c("#096B02", "#A4DCA0"))
luminal.cols <- luminal.cols(6)

stromal.cols <- colorRampPalette(c("#8B059F", "#EFBEF6"))
stromal.cols <- stromal.cols(7)

t.cols <- "#FF8000"
b.cols <- c("#404040","#C0C0C0")
myeloid <- c("#244BCB","#AEBFF6")

cols <- c(basal.cols,luminal.cols,stromal.cols,t.cols,b.cols,myeloid)

p1 <- ggplot(atac.df,aes(x =x,y=y,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle(paste0("scATAC-seq colored by inferred cluster (n=",nrow(atac.df)," cells)"))
ggsave("./Full_Cohort_Results/UMAP_ATAC_Cluster_Number.pdf",width = 8,height = 7)

p1 <- ggplot(atac.df,aes(x = x,y=y,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()+
  ggtitle(paste0("scATAC-seq colored by inferred cluster (n=",nrow(atac.df)," cells)"))
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)
ggsave("./Full_Cohort_Results/UMAP_ATAC_Cluster_Number-Labels.pdf",width = 8,height = 7)

# Plot RNA and ATAC UMAP by sample

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
all.equal(rownames(rna.df),rownames(rna@meta.data))
all.equal(rownames(rna.df),colnames(rna))
rna.df$Sample <- factor(rna$orig.ident,levels=c("35A4AL",
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
                                                "4B146L"))

library(RColorBrewer)
cols <- RColorBrewer::brewer.pal(12,"Paired")
cols <- c(cols[1:8],"#a97ac2","#ccb4d9",cols[c(10,12)],"#B4B4B4","#828282","#505050","#323232")

set.seed(2)
shuffle <- sample(1:nrow(rna.df),nrow(rna.df)) # shuffle df rows to avoid overplotting
ggplot(rna.df[shuffle,],aes(x = UMAP_1,y=UMAP_2,color = Sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle(paste0("scRNA-seq colored by patient sample (n=",nrow(rna.df)," cells)"))
ggsave("./Full_Cohort_Results/UMAP_RNA_Patient_Sample.pdf",width = 8,height = 7)

p1 <- ggplot(rna.df[shuffle,],aes(x = UMAP_1,y=UMAP_2,color = Sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()+
  ggtitle(paste0("scRNA-seq colored by patient sample (n=",nrow(rna.df)," cells)"))
LabelClusters(p1,id="Sample",color="black",repel = T,size=6)
ggsave("./Full_Cohort_Results/UMAP_RNA_Patient_Sample-Labels.pdf",width = 8,height = 7)

# ATAC
atac <- loadArchRProject(path =  "./Patient_Samples_scATAC")

atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "Sample",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$Sample <- sub(".*?-","",atac.df$color)
atac.df$Sample <- factor(atac.df$Sample,levels=c("35A4AL",
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
                                                "4B146L"))

library(RColorBrewer)
cols <- RColorBrewer::brewer.pal(12,"Paired")
cols <- c(cols[1:8],"#a97ac2","#ccb4d9",cols[c(10,12)],"#B4B4B4","#828282","#505050","#323232")

set.seed(2)
shuffle <- sample(1:nrow(atac.df),nrow(atac.df)) # shuffle df rows to avoid overplotting
ggplot(atac.df[shuffle,],aes(x =x,y=y,color = Sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggtitle(paste0("scATAC-seq colored by patient sample (n=",nrow(atac.df)," cells)"))
ggsave("./Full_Cohort_Results/UMAP_ATAC_Patient_Sample.pdf",width = 8,height = 7)

p1 <- ggplot(atac.df[shuffle,],aes(x = x,y=y,color = Sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()+
  ggtitle(paste0("scATAC-seq colored by patient sample (n=",nrow(atac.df)," cells)"))
LabelClusters(p1,id="Sample",color="black",repel = T,size=8)
ggsave("./Full_Cohort_Results/UMAP_ATAC_Patient_Sample-Labels.pdf",width = 8,height = 7)

# Proportion barcharts

rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

meta <- rna@meta.data

meta$Cluster <- factor(meta$RNA_snn_res.0.4,levels=rev(c("11","14",
                                                     "6","12","18","20","23",
                                                     "7","21",
                                                     "4","22",
                                                     "9",
                                                     "0","3","10","13",
                                                     "5",
                                                     "2","17",
                                                     "1",
                                                     "16","19",
                                                     "8","15")))

meta$Sample <- factor(meta$orig.ident,levels=c("35A4AL",
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
                                                "4B146L"))

library(RColorBrewer)
cols <- RColorBrewer::brewer.pal(12,"Paired")
cols <- c(cols[1:8],"#a97ac2","#ccb4d9",cols[c(10,12)],"#B4B4B4","#828282","#505050","#323232")

# Begin plotting
df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(SCSubtype)
colnames(df) <- c("Cluster","SCSubtype","Cells")
p3 <- df %>%
  ggplot(aes(fill=SCSubtype, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()+
  ggtitle("scRNA-seq")

df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(normal_cell_call)
colnames(df) <- c("Cluster","CNV_status","Cells")
df$CNV_status[is.na(df$CNV_status)] <- "NA"
df$CNV_status <- as.factor(df$CNV_status)
df$CNV_status <- plyr::mapvalues(df$CNV_status,
                                 from=c("cancer","NA","normal","unassigned"),
                                 to=c("high","NA","low","ambiguous"))
df$CNV_status <- factor(df$CNV_status,
                        levels=rev(c("NA","low","ambiguous","high") ))

p2 <- df %>%
  ggplot(aes(fill=CNV_status, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  scale_fill_manual(values=rev(c("#807a7a","#faa0a0","#de4b4b","#b80404")))+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()+
  ggtitle("scRNA-seq")

df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")
p1 <- df %>%
  ggplot(aes(fill=Sample, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  scale_fill_manual(values=cols)+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()+
  ggtitle("scRNA-seq")

p1+p2+p3
ggsave(paste0("./Full_Cohort_Results/Proportion_Bar_Charts_Sample_inferCNV_SCSubtype_RNA",".pdf"),
       width = 14,
       height = 7)

meta$cell.type <- gsub('[[:digit:]]+', '', meta$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor)
meta$cell.type <- gsub("^\\-","",meta$cell.type)
df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(cell.type)
colnames(df) <- c("Cluster","refCellType","Cells")
p4 <- df %>%
  ggplot(aes(fill=refCellType, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()+
  ggtitle("scRNA-seq")
p4
ggsave(paste0("./Full_Cohort_Results/Proportion_Bar_Charts_Reference_Cell_Type_Annotations_RNA",".pdf"),
       width = 10,
       height = 7)

# ATAC
atac <- loadArchRProject(path =  "./Patient_Samples_scATAC")

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
rna$barcode <- rownames(rna@meta.data)

meta <- as.data.frame(atac@cellColData)

all.equal(rownames(meta),atac$cellNames)
meta$cellNames <- rownames(meta)
all.equal(meta$cellNames,atac$cellNames)
length(intersect(meta$predictedCell,colnames(rna)))
length(intersect(meta$predictedCell,rownames(rna@meta.data)))
length(intersect(meta$predictedCell,rna$barcode))
meta$barcode <- meta$predictedCell
meta$id  <- 1:nrow(meta)

meta <- merge(meta,rna@meta.data,by="barcode")
meta <- meta[order(meta$id), ]

all.equal(meta$cellNames,atac$cellNames)

meta$Cluster <- factor(meta$predictedGroup,levels=rev(c("11","14",
                                                     "6","12","18","20","23",
                                                     "7","21",
                                                     "4","22",
                                                     "9",
                                                     "0","3","10","13",
                                                     "5",
                                                     "2","17",
                                                     "1",
                                                     "16","19",
                                                     "8","15")))

meta$Sample <- factor(meta$Sample,levels=c("35A4AL",
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
                                               "4B146L"))


library(RColorBrewer)
cols <- RColorBrewer::brewer.pal(12,"Paired")
cols <- c(cols[1:8],"#a97ac2","#ccb4d9",cols[c(10,12)],"#B4B4B4","#828282","#505050","#323232")

# Begin plotting
df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(SCSubtype)
colnames(df) <- c("Cluster","SCSubtype","Cells")
p3 <- df %>%
  ggplot(aes(fill=SCSubtype, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()+
  ggtitle("scATAC-seq")

df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(normal_cell_call)
colnames(df) <- c("Cluster","CNV_status","Cells")
df$CNV_status[is.na(df$CNV_status)] <- "NA"
df$CNV_status <- as.factor(df$CNV_status)
df$CNV_status <- plyr::mapvalues(df$CNV_status,
                                 from=c("cancer","NA","normal","unassigned"),
                                 to=c("high","NA","low","ambiguous"))
df$CNV_status <- factor(df$CNV_status,
                        levels=rev(c("NA","low","ambiguous","high") ))

p2 <- df %>%
  ggplot(aes(fill=CNV_status, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  scale_fill_manual(values=rev(c("#807a7a","#faa0a0","#de4b4b","#b80404")))+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()+
  ggtitle("scATAC-seq")

df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")
p1 <- df %>%
  ggplot(aes(fill=Sample, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  scale_fill_manual(values=cols)+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  theme_classic()+
  ggtitle("scATAC-seq")

p1+p2+p3
ggsave(paste0("./Full_Cohort_Results/Proportion_Bar_Charts_Sample_inferCNV_SCSubtype_ATAC",".pdf"),
       width = 14,
       height = 7)


# Plot QC histograms

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

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
atac <- loadArchRProject(path =  "./Patient_Samples_scATAC")

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
ggsave("./Full_Cohort_Results/QC_Histograms_RNA_ATAC.pdf",width = 9,height = 8)

# Plot prediction score ATAC UMAP
atac <- loadArchRProject(path =  "./Patient_Samples_scATAC")
atac$barcode <- rownames(atac@cellColData)
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "cluster"
atac.df$cluster <- sub(".*?-","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac$idx <- 1:nrow(atac@cellColData)

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

all.equal(atac.df$cluster,atac$predictedGroup)
atac.df$predictedScore <- as.numeric(atac$predictedScore)

umapATACpredScores <- ggplot(atac.df,aes(x = x, y= y,color = predictedScore))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_gradient2(midpoint=0.6,mid="white",low="blue",high="red3")+
  ggtitle(paste0("scATAC-seq colored by prediction score (n=",nrow(atac.df)," cells)"))
ggsave("./Full_Cohort_Results/UMAP_ATAC_Prediction_Scores.pdf",width = 8,height = 7)

predScoreHist <- ggplot(atac.df,aes(x=predictedScore))+
  geom_histogram(binwidth = 0.025,fill="gray", color="black", alpha=0.9)+
  theme_classic()+
  ggtitle(paste0("Prediction score distribution across scATAC-seq cells (n=",nrow(atac.df)," cells)"))
ggsave("./Full_Cohort_Results/Histogram_ATAC_Prediction_Scores.pdf",width = 8,height = 12)

# Proportion barcharts

rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

meta <- rna@meta.data

meta$Cluster <- factor(meta$RNA_snn_res.0.4,levels=rev(c("11","14",
                                                         "6","12","18","20","23",
                                                         "7","21",
                                                         "4","22",
                                                         "9",
                                                         "0","3","10","13",
                                                         "5",
                                                         "2","17",
                                                         "1",
                                                         "16","19",
                                                         "8","15")))

meta$Sample <- factor(meta$orig.ident,levels=c("35A4AL",
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
                                               "4B146L"))

library(RColorBrewer)
cols <- RColorBrewer::brewer.pal(12,"Paired")
cols <- c(cols[1:8],"#a97ac2","#ccb4d9",cols[c(10,12)],"#B4B4B4","#828282","#505050","#323232")

# Begin plotting
df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(SCSubtype)
colnames(df) <- c("Cluster","SCSubtype","Cells")
pSubtype.RNA <- df %>%
  ggplot(aes(fill=SCSubtype, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  ggtitle("scRNA-seq")

df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(normal_cell_call)
colnames(df) <- c("Cluster","CNV_status","Cells")
df$CNV_status[is.na(df$CNV_status)] <- "NA"
df$CNV_status <- as.factor(df$CNV_status)
df$CNV_status <- plyr::mapvalues(df$CNV_status,
                                 from=c("cancer","NA","normal","unassigned"),
                                 to=c("high","NA","low","ambiguous"))
df$CNV_status <- factor(df$CNV_status,
                        levels=rev(c("NA","low","ambiguous","high") ))

pCNV.RNA <- df %>%
  ggplot(aes(fill=CNV_status, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  scale_fill_manual(values=rev(c("#807a7a","#faa0a0","#de4b4b","#b80404")))+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  ggtitle("scRNA-seq")

df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")
pSample.RNA <- df %>%
  ggplot(aes(fill=Sample, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  scale_fill_manual(values=cols)+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  ggtitle("scRNA-seq")

meta$cell.type <- gsub('[[:digit:]]+', '', meta$Cluster_Majority_No_CCA_Reference_Cell_Type_Minor)
meta$cell.type <- gsub("^\\-","",meta$cell.type)
df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(cell.type)
colnames(df) <- c("Cluster","refCellType","Cells")
df$refCellType <- factor(df$refCellType,levels=c("Cancer Basal SC",
                                                 "Cancer Cycling",
                                                 "Cancer LumA SC",
                                                 "Luminal Progenitors",
                                                 "Myoepithelial",
                                                 "Mature Luminal",
                                                 "CAFs myCAF-like",
                                                 "CAFs MSC iCAF-like",
                                                 "PVL Immature",
                                                 "PVL Differentiated",
                                                 "Endothelial ACKR",
                                                 "Endothelial CXCL",
                                                 "Endothelial Lymphatic LYVE",
                                                 "Endothelial RGS",
                                                 "T cells CD+",
                                                 "Cycling T-cells",
                                                 "B cells Memory",
                                                 "Plasmablasts",
                                                 "Macrophage",
                                                 "Monocyte",
                                                 "DCs",
                                                 "Mast cell"
))
pCellType.RNA <- df %>%
  ggplot(aes(fill=refCellType, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  ggtitle("scRNA-seq")

# ATAC
atac <- loadArchRProject(path =  "./Patient_Samples_scATAC")

# Read in scRNA
rna <- readRDS("./Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")
rna$barcode <- rownames(rna@meta.data)

meta <- as.data.frame(atac@cellColData)

all.equal(rownames(meta),atac$cellNames)
meta$cellNames <- rownames(meta)
all.equal(meta$cellNames,atac$cellNames)
length(intersect(meta$predictedCell,colnames(rna)))
length(intersect(meta$predictedCell,rownames(rna@meta.data)))
length(intersect(meta$predictedCell,rna$barcode))
meta$barcode <- meta$predictedCell
meta$id  <- 1:nrow(meta)

meta <- merge(meta,rna@meta.data,by="barcode")
meta <- meta[order(meta$id), ]

all.equal(meta$cellNames,atac$cellNames)

meta$Cluster <- factor(meta$predictedGroup,levels=rev(c("11","14",
                                                        "6","12","18","20","23",
                                                        "7","21",
                                                        "4","22",
                                                        "9",
                                                        "0","3","10","13",
                                                        "5",
                                                        "2","17",
                                                        "1",
                                                        "16","19",
                                                        "8","15")))

meta$Sample <- factor(meta$Sample,levels=c("35A4AL",
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
                                           "4B146L"))

library(RColorBrewer)
cols <- RColorBrewer::brewer.pal(12,"Paired")
cols <- c(cols[1:8],"#a97ac2","#ccb4d9",cols[c(10,12)],"#B4B4B4","#828282","#505050","#323232")

# Begin plotting
df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(SCSubtype)
colnames(df) <- c("Cluster","SCSubtype","Cells")
pSubtype.ATAC <- df %>%
  ggplot(aes(fill=SCSubtype, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  ggtitle("scATAC-seq")

df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(normal_cell_call)
colnames(df) <- c("Cluster","CNV_status","Cells")
df$CNV_status[is.na(df$CNV_status)] <- "NA"
df$CNV_status <- as.factor(df$CNV_status)
df$CNV_status <- plyr::mapvalues(df$CNV_status,
                                 from=c("cancer","NA","normal","unassigned"),
                                 to=c("high","NA","low","ambiguous"))
df$CNV_status <- factor(df$CNV_status,
                        levels=rev(c("NA","low","ambiguous","high") ))

pCNV.ATAC <- df %>%
  ggplot(aes(fill=CNV_status, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  scale_fill_manual(values=rev(c("#807a7a","#faa0a0","#de4b4b","#b80404")))+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  ggtitle("scATAC-seq")

df <- meta %>% dplyr::group_by_at("Cluster") %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")
pSample.ATAC <- df %>%
  ggplot(aes(fill=Sample, y=Cells, x= Cluster))+
  geom_bar(position="fill", stat="identity")+NoLegend()+
  theme(axis.title.x = element_text(size = 9))+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.text.y = element_text(size = 8))+
  scale_fill_manual(values=cols)+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Fraction of cells")+
  ggtitle("scATAC-seq")

pCellType.RNA+pSample.RNA+pCNV.RNA+pSubtype.RNA+pSample.ATAC+pCNV.ATAC+pSubtype.ATAC+plot_layout(ncol = 7)
ggsave(paste0("./Full_Cohort_Results/Proportion_Bar_Charts_Cell_Type_inferCNV_Subtype",".pdf"),
       width = 28,
       height = 10)

# Make pseudobulk hierarchical clustering heatmap and PCA of subtype-specific cancer and normal profiles

# Subset to Basal SC, normal basal, and normal LP

# Read in scRNA
rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

rna <- rna[,rna$RNA_snn_res.0.4 %in% c("11",
                                       "14",
                                       "4",
                                       "7")]
normals <- c("49758L","4AF75L","4B146L","49CFCL")
rna$From_Normal_Sample <- ifelse(rna$orig.ident %in% normals,"Yes","No")

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 %in% c("4","7") & rna$From_Normal_Sample == "No","Drop","Keep" )

rna <- rna[,rna$Keep.Cells != "Drop"]

basals <- c("35A4AL","4C2E5L")
rna$From_Basal_Sample <- ifelse(rna$orig.ident %in% basals,"Yes","No")

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 %in% c("11","14") & rna$From_Basal_Sample == "No","Drop","Keep" )

rna <- rna[,rna$Keep.Cells != "Drop"]

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 %in% c("11","14") & rna$SCSubtype != "Basal_SC","Drop","Keep" )
basal <- rna[,rna$Keep.Cells != "Drop"]

# Subset to Luminal SC and normal ML

# Read in scRNA
rna <- readRDS("Processed_Seurat_Objects_scRNA/Processed_Seurat_Object_Patient_Cohort-MultiKClustering-CellTypeAnnotations-inferCNV-SCSubtype.rds")

# Subset
rna <- rna[,rna$RNA_snn_res.0.4 %in% c("20",
                                       "6",
                                       "23",
                                       "18",
                                       "12",
                                       "9")]
normals <- c("49758L","4AF75L","4B146L","49CFCL")
rna$From_Normal_Sample <- ifelse(rna$orig.ident %in% normals,"Yes","No")

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 == "9" & rna$From_Normal_Sample == "No","Drop","Keep" )

rna <- rna[,rna$Keep.Cells != "Drop"]

luminals <- c("35EE8L", "3821AL", "3B3E9L", "3C7D1L", "3D388L", "3FCDEL", "43E7BL", "43E7CL", "44F0AL", "45CB0L")
rna$From_Luminal_Sample <- ifelse(rna$orig.ident %in% luminals,"Yes","No")

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 %in% c("20",
                                                    "6",
                                                    "23",
                                                    "18",
                                                    "12") & rna$From_Luminal_Sample == "No","Drop","Keep" )

rna <- rna[,rna$Keep.Cells != "Drop"]

rna$Keep.Cells <- ifelse(rna$RNA_snn_res.0.4 %in% c("20",
                                                    "6",
                                                    "23",
                                                    "18",
                                                    "12") & rna$SCSubtype %ni% c("LumA_SC","LumB_SC"),"Drop","Keep" )
lum <- rna[,rna$Keep.Cells != "Drop"]

rna <- merge(lum,basal)

rna$SCSubtype <- rna$SCSubtype %>% tidyr::replace_na('NA')

rna$SCSubtype <- ifelse(rna$SCSubtype %in% c("LumA_SC","LumB_SC"),"Luminal_SC",rna$SCSubtype)
rna$Pseudobulk <- paste0(rna$SCSubtype,"_",rna$orig.ident)
rna$Pseudobulk <- str_replace_all(string = rna$Pseudobulk,pattern = "NA",replacement = "Normal")

rna$Pseudobulk <- ifelse(rna$RNA_snn_res.0.4 == "9",paste0("Normal_ML_",rna$orig.ident),rna$Pseudobulk)
rna$Pseudobulk <- ifelse(rna$RNA_snn_res.0.4 == "4",paste0("Normal_Basal_",rna$orig.ident),rna$Pseudobulk)
rna$Pseudobulk <- ifelse(rna$RNA_snn_res.0.4 == "7",paste0("Normal_LP_",rna$orig.ident),rna$Pseudobulk)

table(rna$Pseudobulk)

# This custom function generates a heatmap with hierarchical clustering
pseudobulk_heatmap <- function(seurat,
                               groupBy,
                               clusters,
                               topProp,
                               sigclust_params,
                               resList){
  
  Idents(seurat) <- groupBy
  seurat <- subset(x = seurat, idents = clusters)
  
  # Extract raw counts and metadata to create SingleCellExperiment object
  counts <- seurat@assays$RNA@counts
  
  metadata <- seurat@meta.data
  
  # Set up metadata as desired for aggregation 
  
  # Create single cell experiment object
  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = metadata)
  
  groups <- colData(sce)[, groupBy ]
  unique(groups)
  
  # Aggregate across cluster-sample groups
  # transposing row/columns to have cell_ids as row names matching those of groups
  aggr_counts <- aggregate.Matrix(t(counts(sce)),
                                  groupings = groups, fun = "sum")
  
  # Explore output matrix
  class(aggr_counts)
  dim(aggr_counts)
  aggr_counts <- t(aggr_counts)
  
  cluster_metadata <- data.frame(groupBy = colnames(aggr_counts),
                                 cellType = sub("^([^_]*_[^_]*).*", "\\1", colnames(aggr_counts)),
                                 sample = sub(".*\\_", "", colnames(aggr_counts)),
                                 cluster = c(rep("2",2),rep("1",10),rep("2",8),rep("1",4)),# sigclust2 cluster annotations added retrospectively
                                 status = ifelse(sub(".*\\_", "", colnames(aggr_counts)) %in% c("49758L", "49CFCL", "4AF75L", "4B146L"),"normal","cancer"))
  colnames(cluster_metadata)[1] <- groupBy
  
  all.equal(cluster_metadata[[groupBy]],names(table(colData(sce)[[groupBy]])))
  
  cluster_metadata$cell_count <- as.numeric(table(colData(sce)[[groupBy]]))
  
  # Read in clincal info
  
  library(readxl)
  # Read in clinical data
  clinical <- read_excel("miscellaneous/Franco-Perou-SingleCellBreastCancerDataset-July2022_MR-MJR_KW.xlsx")
  clinical$sample <- clinical$Sample_ID
  clinical$Patient_Age_at_Surgery <- as.numeric(clinical$Patient_Age_at_Surgery)
  clinical$BMI <- as.numeric(clinical$BMI)
  colnames(clinical)[16] <- "Menopause_Status"
  clinical$Menopause_Status <- gsub("\\ .*","",clinical$Menopause_Status)
  clinical <- clinical[!duplicated(clinical),]
  
  cluster_metadata <- merge(cluster_metadata,clinical,by="sample")
  
  cluster_metadata <- cluster_metadata[order(match(cluster_metadata[[groupBy]],colnames(aggr_counts))),]
  
  all.equal(cluster_metadata[[groupBy]],colnames(aggr_counts))
  
  row_sub = apply(aggr_counts, 1, function(row) all(row !=0 ))
  
  dds <- DESeqDataSetFromMatrix(aggr_counts[row_sub,],
                                colData = cluster_metadata,
                                design = ~ 1)
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  
  #Data transformation & Scaling
  rlog.rna <- rlog(dds, blind=T)
  
  # Sigclust denodgram
  rv <- rowVars(assay(rlog.rna))
  select <- order(rv, decreasing = TRUE)[ seq_len( round(topProp*length(rv)) ) ]
  
  shc_result <- sigclust2::shc(as.matrix(t(assay(rlog.rna)[select, ])),
                               metric = sigclust_params[[1]],
                               linkage = sigclust_params[[2]],
                               alpha = sigclust_params[[3]],
                               n_min = sigclust_params[[4]])
  
  resList[[1]] <- plot(shc_result)
  
  hclust <- shc_result$hc_dat
  
  ha = HeatmapAnnotation(
    
    Status = rlog.rna$status,
    Type = rlog.rna$cellType,
    Sample = rlog.rna$sample,
    col = list(Status = c("cancer" = "#b50404", "normal" = "#706f6f" ),
               Type = c("Basal_SC" = "#A50505", "Luminal_SC" = "#479841",
                        "Normal_Basal" = "#FF9B9B", "Normal_LP" = "#DB5F5F",
                        "Normal_ML" = "#A4DCA0"),
               Sample = c("35A4AL"= "#A6CEE3",
                          "4C2E5L" = "#1F78B4",
                          "35EE8L" = "#B2DF8A",
                          "3821AL" = "#33A02C",
                          "3B3E9L" = "#FB9A99",
                          "3C7D1L" = "#E31A1C",
                          "3D388L" = "#FDBF6F",
                          "3FCDEL" = "#FF7F00",
                          "43E7BL" = "#a97ac2",
                          "43E7CL" = "#ccb4d9",
                          "44F0AL" = "#6A3D9A",
                          "45CB0L" = "#B15928",
                          "49758L" = "#B4B4B4",
                          "49CFCL" = "#828282",
                          "4AF75L" = "#505050",
                          "4B146L" = "#323232" ))
    
  )
  
  ht <- Heatmap(as.matrix(t(scale(t(assay(rlog.rna)[select, ])))),cluster_columns = hclust,column_split = 2,
                top_annotation = ha,show_row_names = TRUE,show_column_names = TRUE)
  
  resList[[2]] <- ht
  
  resList[[3]] <- rlog.rna
  
  resList[[4]] <- assay(rlog.rna)[select, ]
  
  return(resList)
}

res <- pseudobulk_heatmap(seurat = rna,
                          groupBy = "Pseudobulk",
                          clusters = unique(rna$Pseudobulk),
                          topProp = 0.10,
                          sigclust_params = list("euclidean",
                                                 "ward.D2",
                                                 0.05,
                                                 24),
                          resList = vector(mode='list', length=4))

pdf("./Full_Cohort_Results/Pseudobulk_Hierarchical_Clustering_Dendrogram.pdf",width=8,height = 9)
res[[1]]
dev.off()

pdf("./Full_Cohort_Results/Pseudobulk_Hierarchical_Clustering_Heatmap.pdf",width=10,height = 8)
res[[2]]
dev.off()

saveRDS(res[[3]],"./Full_Cohort_Results/rlog_object.rds")

saveRDS(res[[4]],"./Full_Cohort_Results/rlog_object-top_variable_features.rds")

# Plot PCA of discrete variables
topProp <- 0.10
n_features = round(topProp*length(res[[3]]))

p1 <- plotPCA(res[[3]],intgroup = "cluster",ntop = n_features)+
  ggtitle(paste0("PCA colored by SigClust2 cluster\nTop 10% variable genes: ",n_features))+
  theme_bw()+
  scale_color_manual(values=c("#2081f9","#f99820"))

p2 <- plotPCA(res[[3]],intgroup = "status",ntop = n_features)+
  ggtitle(paste0("PCA colored by status\nTop 10% variable genes: ",n_features))+
  theme_bw()+
  scale_color_manual(values=c("#b50404","#706f6f"))

p3 <- plotPCA(res[[3]],intgroup = "cellType",ntop = n_features)+
  ggtitle(paste0("PCA colored by cell type\nTop 10% variable genes: ",n_features))+
  theme_bw()+
  scale_color_manual(values=c("#A50505","#479841","#FF9B9B","#DB5F5F","#A4DCA0"))

res[[3]]$sample <- factor(res[[3]]$sample,levels=c("35A4AL",
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
                                                   "4B146L"))


p4 <- plotPCA(res[[3]],intgroup = "sample",ntop = n_features)+
  ggtitle(paste0("PCA colored by sample\nTop 10% variable genes: ",n_features))+
  theme_bw()+
  scale_color_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                              "#FDBF6F", "#FF7F00", "#a97ac2", "#ccb4d9", "#6A3D9A",
                              "#B15928", "#B4B4B4", "#828282", "#505050", "#323232"))

p1+p2+p3+p4+plot_layout(ncol = 2)
ggsave("./Full_Cohort_Results/PCA_cellType_status_cluster_sample.pdf",
       width = 11,
       height = 11)


# Plot PCA of continuous variables 
topProp <- 0.10
n_features = round(topProp*length(res[[3]]))

p1 <- plotPCA(res[[3]],intgroup = "cell_count",ntop = n_features)+
  scale_color_gradient(low="grey",high="blue")+
  ggtitle(paste0("PCA colored by cell count\nTop 10% variable genes: ",n_features))+
  theme_bw()

p2 <- plotPCA(res[[3]],intgroup = "sizeFactor",ntop = n_features)+
  scale_color_gradient(low="grey",high="purple")+
  ggtitle(paste0("PCA colored by sizeFactor\nTop 10% variable genes: ",n_features))+
  theme_bw()

p3 <- plotPCA(res[[3]],intgroup = "BMI",ntop = n_features)+
  scale_color_gradient(low="grey",high="red")+
  ggtitle(paste0("PCA colored by BMI\nTop 10% variable genes: ",n_features))+
  theme_bw()

p4 <- plotPCA(res[[3]],intgroup = "Patient_Age_at_Surgery",ntop = n_features)+
  scale_color_gradient(low="grey",high="orange")+
  ggtitle(paste0("PCA colored by age\nTop 10% variable genes: ",n_features))+
  theme_bw()

p1+p2+p3+p4+plot_layout(ncol = 2)
ggsave("./Full_Cohort_Results/PCA_cellCount_sizeFactor_BMI_Age.pdf",
       width = 11,
       height = 11)
