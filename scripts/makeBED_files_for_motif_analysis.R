################################################################################
# Matt Regner
# Franco Lab

################################################################################
library(ExperimentHub)
library(dplyr)
library(cowplot)
library(patchwork)
library(BSDA)
library(data.table)
library(ggplot2)
library(stringr)

# Make output folder
if(dir.exists("Basal_P2G_BED_files")){
  print("Directory Basal_P2G_BED_files already exists!")
}else{
  dir.create("Basal_P2G_BED_files")
}

# Read in all differential P2Gs
p2g <- readRDS("./Basal_Cohort_Results-updates-Revised/p2gs_interaction.rds")

for(i in unique(p2g$DiffClass)){
  
  sub <- p2g[p2g$DiffClass == i,]
  bed <- data.frame(chrom=gsub("\\:.*","",sub$peakName),
                    chromStart=gsub(".*:(.+)-.*", "\\1", sub$peakName),
                    chromEnd=sub('.*-', '', sub$peakName),
                    name=sub$peakName,
                    column5="NA",
                    strand="both")
  bed <- bed[!duplicated(bed),]
  
  j <- str_replace(i,"/","_")
  fwrite(bed,file=paste0("Basal_P2G_BED_files/unique_peaks_in_",j,"_P2Gs.bed"),
         sep = "\t")
  
}

p2g.cancer <- readRDS("./Basal_Cohort_Results-updates-Revised/p2gs_interaction-+_0-upgenes.rds")
bed <- data.frame(chrom=gsub("\\:.*","",p2g.cancer$peakName),
                  chromStart=gsub(".*:(.+)-.*", "\\1", p2g.cancer$peakName),
                  chromEnd=sub('.*-', '', p2g.cancer$peakName),
                  name=p2g.cancer$peakName,
                  column5="NA",
                  strand="both")
bed <- bed[!duplicated(bed),]
fwrite(bed,file=paste0("Basal_P2G_BED_files/unique_peaks_in_+_0_P2Gs-upgenes.bed"),
       sep = "\t")

p2g.normal <- readRDS("./Basal_Cohort_Results-updates-Revised/p2gs_interaction-0_+-upgenes.rds")
bed <- data.frame(chrom=gsub("\\:.*","",p2g.normal$peakName),
                  chromStart=gsub(".*:(.+)-.*", "\\1", p2g.normal$peakName),
                  chromEnd=sub('.*-', '', p2g.normal$peakName),
                  name=p2g.normal$peakName,
                  column5="NA",
                  strand="both")
bed <- bed[!duplicated(bed),]
fwrite(bed,file=paste0("Basal_P2G_BED_files/unique_peaks_in_0_+_P2Gs-upgenes.bed"),
       sep = "\t")

# Make output folder
if(dir.exists("Luminal_P2G_BED_files")){
  print("Directory Luminal_P2G_BED_files already exists!")
}else{
  dir.create("Luminal_P2G_BED_files")
}

# Read in all differential P2Gs
p2g <- readRDS("./Luminal_Cohort_Results-updates-Revised/p2gs_interaction.rds")

for(i in unique(p2g$DiffClass)){
  
  sub <- p2g[p2g$DiffClass == i,]
  bed <- data.frame(chrom=gsub("\\:.*","",sub$peakName),
                    chromStart=gsub(".*:(.+)-.*", "\\1", sub$peakName),
                    chromEnd=sub('.*-', '', sub$peakName),
                    name=sub$peakName,
                    column5="NA",
                    strand="both")
  bed <- bed[!duplicated(bed),]
  
  j <- str_replace(i,"/","_")
  fwrite(bed,file=paste0("Luminal_P2G_BED_files/unique_peaks_in_",j,"_P2Gs.bed"),
         sep = "\t")
  
}

p2g.cancer <- readRDS("./Luminal_Cohort_Results-updates-Revised/p2gs_interaction-+_0-upgenes.rds")
bed <- data.frame(chrom=gsub("\\:.*","",p2g.cancer$peakName),
                  chromStart=gsub(".*:(.+)-.*", "\\1", p2g.cancer$peakName),
                  chromEnd=sub('.*-', '', p2g.cancer$peakName),
                  name=p2g.cancer$peakName,
                  column5="NA",
                  strand="both")
bed <- bed[!duplicated(bed),]
fwrite(bed,file=paste0("Luminal_P2G_BED_files/unique_peaks_in_+_0_P2Gs-upgenes.bed"),
       sep = "\t")

p2g.normal <- readRDS("./Luminal_Cohort_Results-updates-Revised/p2gs_interaction-0_+-upgenes.rds")
bed <- data.frame(chrom=gsub("\\:.*","",p2g.normal$peakName),
                  chromStart=gsub(".*:(.+)-.*", "\\1", p2g.normal$peakName),
                  chromEnd=sub('.*-', '', p2g.normal$peakName),
                  name=p2g.normal$peakName,
                  column5="NA",
                  strand="both")
bed <- bed[!duplicated(bed),]
fwrite(bed,file=paste0("Luminal_P2G_BED_files/unique_peaks_in_0_+_P2Gs-upgenes.bed"),
       sep = "\t")

