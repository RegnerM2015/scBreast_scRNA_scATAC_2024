# miscellaneous

refdata-cellranger-GRCh38-3_0_0-gene-ordering-file.txt was made using the following command:

```
python ../scripts/gtf_to_position_file.py --attribute_name gene_name ./genes.gtf ./refdata-cellranger-GRCh38-3_0_0-gene-ordering-file.txt
```
Output:
```
Number of lines read: 2565061
Number of comments: 5
Number of entries: 33514
Number of duplicate entries: 2531537
Number of entries written: 33514
```

genes.gtf is listed in the .gitignore as this file is too large to track in a git repo. It was originally retrieved from the cellranger reference located under the following path: 
`/datastore/nextgenout5/share/labs/francolab/Data/refdata-cellranger-GRCh38-3.0.0/genes`

gtf_to_position_file.py was provided by the inferCNV team/developers:
https://github.com/broadinstitute/infercnv/blob/master/scripts/gtf_to_position_file.py

crispr_mean_median_depmap_scores.csv and rnai_mean_median_depmap_scores.csv were created using the following R script below: 

```
################################################################################
# Matt Regner
# Franco Lab
################################################################################
library(depmap)
library(ExperimentHub)
library(dplyr)
set.seed(1)

subtypes <- read.delim("./miscellaneous/PAM50_call_comparison-CMP.txt")
subtypes$cell_line <- paste0(subtypes$X,"_BREAST")

## create ExperimentHub query object
eh <- ExperimentHub()

rnai <- eh[["EH2260"]]
rnai$geneName <- rnai$gene_name
rnai.sub <- rnai[rnai$cell_line %in% subtypes$cell_line,]
subtyped.rnai <- merge(rnai.sub,subtypes,by="cell_line")

# Subset to basal cell lines and calculate mean dependency for each gene across cell lines
subtyped.rnai.basal <- subtyped.rnai[subtyped.rnai$SUBTYPE.SUGGESTED.OVERVIEW == "Basal",]
basal.mean <- subtyped.rnai.basal %>% group_by(geneName) %>% 
  summarize(mean_dependency=mean(dependency))

# Subset to luminal cell lines and calculate mean dependency for each gene across cell lines
subtyped.rnai.luminal <- subtyped.rnai[subtyped.rnai$SUBTYPE.SUGGESTED.OVERVIEW == "Luminal",]
luminal.mean <- subtyped.rnai.luminal %>% group_by(geneName) %>% 
  summarize(mean_dependency=mean(dependency))

# Subset to basal cell lines and calculate median dependency for each gene across cell lines
subtyped.rnai.basal <- subtyped.rnai[subtyped.rnai$SUBTYPE.SUGGESTED.OVERVIEW == "Basal",]
basal.median <- subtyped.rnai.basal %>% group_by(geneName) %>% 
  summarize(median_dependency=median(dependency))

# Subset to luminal cell lines and calculate median dependency for each gene across cell lines
subtyped.rnai.luminal <- subtyped.rnai[subtyped.rnai$SUBTYPE.SUGGESTED.OVERVIEW == "Luminal",]
luminal.median <- subtyped.rnai.luminal %>% group_by(geneName) %>% 
  summarize(median_dependency=median(dependency))

all.equal(basal.mean$geneName,
          luminal.mean$geneName)
all.equal(basal.mean$geneName,
          luminal.median$geneName)
all.equal(basal.median$geneName,
          luminal.median$geneName)
all.equal(basal.median$geneName,
          luminal.mean$geneName)

rnai_gene_effect_scores <- data.frame(geneName = basal.mean$geneName,
                                      mean_Basal = basal.mean$mean_dependency,
                                      mean_Luminal = luminal.mean$mean_dependency,
                                      median_Basal = basal.median$median_dependency,
                                      median_Luminal = luminal.median$median_dependency,
                                      cell_lines = paste0("Basal: ",
                                                          toString(unique(subtyped.rnai.basal$cell_line)),
                                                          " | Luminal: ",toString(unique(subtyped.rnai.luminal$cell_line))),
                                      R_package_versions = paste0(R.Version()$version.string,
                                                                   " | depmap: ",
                                                                   packageVersion("depmap"),
                                                                   " | ExperimentHub: ",
                                                                   packageVersion("ExperimentHub")),
                                      ExperimentHub_accession_number = "EH2260")

write.table(rnai_gene_effect_scores,"./miscellaneous/rnai_mean_median_depmap_scores.csv",
            row.names = F,col.names = T,sep = ",")

################################################################################
crispr <- eh[["EH2261"]]

crispr$geneName <- crispr$gene_name
crispr.sub <- crispr[crispr$cell_line %in% subtypes$cell_line,]
subtyped.crispr <- merge(crispr.sub,subtypes,by="cell_line")

# Subset to basal cell lines and calculate mean dependency for each gene across cell lines
subtyped.crispr.basal <- subtyped.crispr[subtyped.crispr$SUBTYPE.SUGGESTED.OVERVIEW == "Basal",]
basal.mean <- subtyped.crispr.basal %>% group_by(geneName) %>% 
  summarize(mean_dependency=mean(dependency))

# Subset to luminal cell lines and calculate mean dependency for each gene across cell lines
subtyped.crispr.luminal <- subtyped.crispr[subtyped.crispr$SUBTYPE.SUGGESTED.OVERVIEW == "Luminal",]
luminal.mean <- subtyped.crispr.luminal %>% group_by(geneName) %>% 
  summarize(mean_dependency=mean(dependency))

# Subset to basal cell lines and calculate median dependency for each gene across cell lines
subtyped.crispr.basal <- subtyped.crispr[subtyped.crispr$SUBTYPE.SUGGESTED.OVERVIEW == "Basal",]
basal.median <- subtyped.crispr.basal %>% group_by(geneName) %>% 
  summarize(median_dependency=median(dependency))

# Subset to luminal cell lines and calculate median dependency for each gene across cell lines
subtyped.crispr.luminal <- subtyped.crispr[subtyped.crispr$SUBTYPE.SUGGESTED.OVERVIEW == "Luminal",]
luminal.median <- subtyped.crispr.luminal %>% group_by(geneName) %>% 
  summarize(median_dependency=median(dependency))

all.equal(basal.mean$geneName,
          luminal.mean$geneName)
all.equal(basal.mean$geneName,
          luminal.median$geneName)
all.equal(basal.median$geneName,
          luminal.median$geneName)
all.equal(basal.median$geneName,
          luminal.mean$geneName)

crispr_gene_effect_scores <- data.frame(geneName = basal.mean$geneName,
                                      mean_Basal = basal.mean$mean_dependency,
                                      mean_Luminal = luminal.mean$mean_dependency,
                                      median_Basal = basal.median$median_dependency,
                                      median_Luminal = luminal.median$median_dependency,
                                      cell_lines = paste0("Basal: ",
                                                          toString(unique(subtyped.crispr.basal$cell_line)),
                                                          " | Luminal: ",toString(unique(subtyped.crispr.luminal$cell_line))),
                                      R_package_versions = paste0(R.Version()$version.string,
                                                                   " | depmap: ",
                                                                   packageVersion("depmap"),
                                                                   " | ExperimentHub: ",
                                                                   packageVersion("ExperimentHub")),
                                      ExperimentHub_accession_number = "EH2261")

write.table(crispr_gene_effect_scores,"./miscellaneous/crispr_mean_median_depmap_scores.csv",
            row.names = F,col.names = T,sep = ",")

```
